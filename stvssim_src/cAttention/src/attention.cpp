//#include <afxwin.h>
#include <stdio.h>
#include <stack>

#include "cxcore.h"
#include "highgui.h"
#include "cv.h"
#include "attention.h"
//#include "cvgabor.h"

using namespace std;

#ifndef max 
#define max(a,b)            (((a) > (b)) ? (a) : (b)) 
#endif 

#ifndef max3
#define max3(a, b, c)		(a > max(b, c)) ? (a) : (max(b,c))
#endif

#ifndef min 
#define min(a,b)            (((a) < (b)) ? (a) : (b)) 
#endif 

#ifndef min3
#define min3(a, b, c)       (a < min(b, c)) ? (a) : (min(b,c))
#endif 

// Inital Cluster Center using histogram. Preprocess for color quantization
// img:			the input image
// matCenter:	the result
// maxCluster:	the maximum of number of cluster
// CoverRate:	
int InitialClusterCenter(const IplImage* img, CvMat* matCenter, double CoverRate,
						 int R_bin, int G_bin, int B_bin)
{
	cvZero(matCenter);
	// Calcualte image histogram
	int hist_size[] = {R_bin, G_bin, B_bin};
	float r_ranges[] = { 0, 255 };
	float g_ranges[] = { 0, 255 };
	float b_ranges[] = { 0, 255 };
	float* ranges[] = { r_ranges, g_ranges, b_ranges };
	CvHistogram	*h	= cvCreateHist( 3, hist_size, CV_HIST_ARRAY, ranges, 1 );
	cvClearHist(h);

	IplImage	*plane_r, *plane_g, *plane_b;
	plane_r		= cvCreateImage(cvGetSize(img), 8, 1);
	plane_g		= cvCreateImage(cvGetSize(img), 8, 1);
	plane_b		= cvCreateImage(cvGetSize(img), 8, 1);

	IplImage* planes[] = { plane_r, plane_g, plane_b };

	cvCvtPixToPlane( img, plane_b, plane_g, plane_r, 0 );
	cvCalcHist(planes, h, 0, NULL);
	cvNormalizeHist(h, 1);

	// Choose the bins of higher values as cluster center
	int NumCluster = 0, maxCluster	= matCenter->rows;
	double Rate	= 0;
	int maxIndex[3], r_Index, g_Index, b_Index;
	float	minValue = 0, maxValue = 0;
	int	k = 0, nDim = 3;
	while (Rate < CoverRate && NumCluster < maxCluster)
	{
		cvGetMinMaxHistValue(h, &minValue, &maxValue, 0, maxIndex);
		r_Index = maxIndex[0];
		g_Index = maxIndex[1];
		b_Index = maxIndex[2];

		// The order is B-G-R
		CV_MAT_ELEM(*matCenter, double, NumCluster, 0) = (double)(2 * b_Index + 1) * 256 / (2 * B_bin);
		CV_MAT_ELEM(*matCenter, double, NumCluster, 1) = (double)(2 * g_Index + 1) * 256 / (2 * G_bin);
		CV_MAT_ELEM(*matCenter, double, NumCluster, 2) = (double)(2 * r_Index + 1) * 256 / (2 * R_bin);

		cvSetReal3D(h->bins, r_Index, g_Index, b_Index, 0);
		NumCluster ++;
		Rate += maxValue;
	}

	// Release memory
	cvReleaseImage(&plane_r);
	cvReleaseImage(&plane_g);
	cvReleaseImage(&plane_b);
	cvReleaseHist(&h);

	return NumCluster;
}

// img:			Input image
// imgResult:	color quantization result
// ClusterCount: The number of cluster in previous result, used only when "UsePrevious", 
// labels:		an image, same size with img, 32U
// center_arr:	the center of each cluster, center_arr.rows = maxCluster, center_arr.cols = img.dims
// the color in center_arr is consistent with img, i.e.:BGR
// ElemCount:	The number of elements of each cluster
// Useprevious: if 1, the current cluster center will be used; if 2, the current label will be used; if 0, cluster from scratch
// minDis:		the clusters of distance smaller than minDis will be merged.
// minCluster:	The min number of cluster element, if the number of element of a cluster is smaller than minCluster, it will be merged to others.
int GLA( const IplImage* img, IplImage* imgResult, 
		int ClusterCount, int maxCluster,
		IplImage* labels, CvMat* centers_arr, CvMat* ElemCount,
		CvTermCriteria termcrit, int UsePrevious, 
		double minDis, int minCluster )
{
	int cluster_count = maxCluster;
	double minDistance	= minDis * minDis;
	CvMat* centers = 0;
	CvMat* old_centers = 0;
	CvMat* counters = 0;

	CvRNG rng = CvRNG(-1);
	int i, j, x, y, k;
	int iter;
	double max_dist;

	int dims = 3;

	// parameter check

	if( maxCluster < 1 )
	{
		//		AfxMessageBox("Number of clusters should be positive");
		return 0;
	}

	//    termcrit = cvCheckTermCriteria( termcrit, 1e-6, 100 );

	termcrit.epsilon *= termcrit.epsilon;
	int sample_count = img->height * img->width;

	if( cluster_count > sample_count )
		cluster_count = sample_count;

	if ( UsePrevious == 1 && centers_arr == 0 )
	{
		//		AfxMessageBox("if use previous cluster centers, centers should be a matrix!");
		return 0;
	}

	if (centers_arr != 0) 
	{
		centers = (CvMat*)centers_arr;

		if (centers->rows != maxCluster || centers->cols != dims) 
		{
			//			AfxMessageBox("centers must be a matrix of size maxCluster * dims!");
			return 0;
		}

		if (CV_MAT_DEPTH(centers->type) != CV_64F)
		{
			//			AfxMessageBox("centers must be a matrix of 64F!");
			return 0;
		}

		old_centers = cvCreateMat( maxCluster, dims, CV_64FC1 );
	}
	else
	{
		centers = cvCreateMat( cluster_count, dims, CV_64FC1 );
		old_centers = cvCreateMat( cluster_count, dims, CV_64FC1 );
	}

	if (ElemCount != 0)
	{
		if (ElemCount->rows != maxCluster)
		{
			//			AfxMessageBox("ElemCount must be a matrix of maxCluster columns!");
			return 0;
		}
		counters = ElemCount;
	}
	else
	{
		counters = cvCreateMat( cluster_count, 1, CV_32SC1 );
	}

	// initialize cluster centers
	if (UsePrevious != 1)
	{
		if (UsePrevious != 2)	// cluster from scratch
		{
			for (y = 0; y < img->height; y ++)
			{
				for (x = 0; x < img->width; x ++)
				{
					CV_IMAGE_ELEM(labels, int, y, x) = cvRandInt(&rng) % cluster_count;
				}
			}
		}
		// compute centers
		cvZero( centers );
		cvZero( counters );

		for (y = 0; y < img->height; y ++)
		{
			for (x = 0; x < img->width; x ++)
			{
				unsigned char * s = &CV_IMAGE_ELEM(img, unsigned char, y, x * 3);
				k	= CV_IMAGE_ELEM(labels, int, y, x);
				double*	c = (double*)(centers->data.ptr + k*centers->step);
				for( j = 0; j < dims; j++ )
				{
					c[j] += (double)(s[j]);
				}
				CV_MAT_ELEM(*counters, int, k, 0) ++;
			}
		}

		for( k = 0; k < cluster_count; k++ )
		{
			double* c = (double*)(centers->data.ptr + k*centers->step);
			if( counters->data.i[k] != 0 )
			{
				double scale = 1.0 / CV_MAT_ELEM(*counters, int, k, 0);
				for( j = 0; j < dims; j++ )
				{
					c[j] *= scale;
				}
			}
			else
			{
				x = cvRandInt( &rng ) % img->width;
				y = cvRandInt( &rng ) % img->height;
				unsigned char* s = &CV_IMAGE_ELEM(img, unsigned char, y, x * 3);
				for( j = 0; j < dims; j++ )
				{
					c[j] = (double)(s[j]);
				}
			}
		}
	}
	else	// use the previous cluster center
	{
		for (i = ClusterCount; i < cluster_count; i ++)
		{
			x = cvRandInt(&rng) % img->width;
			y = cvRandInt(&rng) % img->height;
			unsigned char* s = &CV_IMAGE_ELEM(img, unsigned char, y, x * 3);
			for (j = 0; j < dims; j ++)
			{
				cvmSet(centers, i, j, (double)(s[j]));
			}
		}
	}

	max_dist = termcrit.epsilon * 2;

	// Start to cluster
	for( iter = 0; iter < termcrit.max_iter; iter++ )
	{
		cvCopy(centers, old_centers);
		cvZero( centers );
		cvZero( counters );

		// assign labels and update the centers at the same time
		for(  y = 0; y < img->height; y ++ )
		{
			for (x = 0; x < img->width; x ++)
			{
				unsigned char* s = &CV_IMAGE_ELEM(img, unsigned char, y, x * 3);
				int k_best = 0;
				double min_dist = DBL_MAX;

				// clustering using the old_centers
				for( k = 0; k < cluster_count; k++ )
				{
					double* c = (double*)(old_centers->data.ptr + k * old_centers->step);
					double dist = 0;

					for( j = 0; j < dims; j++ )
					{
						double t = c[j] - (double)(s[j]);
						dist += t*t;
					}

					if( min_dist > dist )
					{
						min_dist = dist;
						k_best = k;
					}
				}

				CV_IMAGE_ELEM(labels, int, y, x) = k_best;

				// update centers
				double* c = (double*)(centers->data.ptr + k_best * centers->step);
				for( j = 0; j < dims; j++ )
				{
					c[j] += s[j];
				}
				counters->data.i[k_best]++;
			}
		}

		max_dist = 0;
		for( k = 0; k < cluster_count; k++ )
		{
			double* c = (double*)(centers->data.ptr + k * centers->step);
			if( counters->data.i[k] != 0 )
			{
				double scale = 1./counters->data.i[k];
				for( j = 0; j < dims; j++ )
				{
					c[j] *= scale;
				}
			}
			else
			{
				x = cvRandInt( &rng ) % img->width;
				y = cvRandInt( &rng ) % img->height;
				unsigned char* s = &CV_IMAGE_ELEM(img, unsigned char, y, x * 3);
				for( j = 0; j < dims; j++ )
				{
					c[j] = (double)(s[j]);
				}
			}

			double dist = 0;
			double* c_o = (double*)(old_centers->data.ptr + k*old_centers->step);
			for( j = 0; j < dims; j++ )
			{
				double t = c[j] - c_o[j];
				dist += t*t;
			}
			if( k == 0 || max_dist < dist )
				max_dist = dist;
		}

		if( max_dist < termcrit.epsilon )
			break;
	}

	// ensure that we do not have empty clusters
	// Merge
	bool flag = true;
	int k1, k2;
	while (flag && cluster_count > 1)
	{
		flag = false;
		//		k = 0;
		k1 = 0;
		while (k1 < cluster_count)
		{
			k2 = k1 + 1;
			while (k2 < cluster_count)
			{
				if ( counters->data.i[k2] < minCluster )
				{
					double* c = (double*)(centers->data.ptr + k2 * centers->step);
					double* c_o = (double*)(centers->data.ptr + (cluster_count - 1) * centers->step);
					for (i = 0; i < dims; i ++)
					{
						c[i] = c_o[i];
					}
					counters->data.i[k2] = counters->data.i[cluster_count - 1];
					cluster_count --;
					flag = true;
				}
				else
				{
					double* c = (double*)(centers->data.ptr + k1 * centers->step);
					double* c_o = (double*)(centers->data.ptr + k2 * centers->step);
					double d = 0;
					for (i = 0; i < dims; i ++)
					{
						d += (c[i] - c_o[i]) * (c[i] - c_o[i]);
					}
					if (d < minDistance)
					{
						double* c = (double*)(centers->data.ptr + k2 * centers->step);
						double* c_o = (double*)(centers->data.ptr + (cluster_count - 1) * centers->step);
						for (i = 0; i < dims; i ++)
						{
							c[i] = c_o[i];
						}
						counters->data.i[k2] = counters->data.i[cluster_count - 1];
						cluster_count --;
						flag = true;
					}
					else
					{
						k2 ++;
					}
				}
			}
			k1 ++;
		}

		if (!flag) break;

		cvCopy(centers, old_centers);
		cvZero( centers );
		cvZero( counters );
		// assign labels
		for (y = 0; y < img->height; y ++)
		{
			for (x = 0; x < img->width; x ++)
			{
				unsigned char* s = &CV_IMAGE_ELEM(img, unsigned char, y, x * 3);
				int k_best = 0;
				double min_dist = DBL_MAX;

				for( k = 0; k < cluster_count; k++ )
				{
					double* c = (double*)(old_centers->data.ptr + k * old_centers->step);
					double dist = 0;

					for( j = 0; j < dims; j++ )
					{
						double t = c[j] - (double)(s[j]);
						dist += t*t;
					}

					if( min_dist > dist )
					{
						min_dist = dist;
						k_best = k;
					}
				}

				CV_IMAGE_ELEM(labels, int, y, x) = k_best;

				// update cluster centers
				double* c = (double*)(centers->data.ptr + k_best * centers->step);
				for( j = 0; j < dims; j++ )
				{
					c[j] += (double)(s[j]);
				}
				counters->data.i[k_best]++;
			}
		}

		for( k = 0; k < cluster_count; k++ )
		{
			double* c = (double*)(centers->data.ptr + k*centers->step);
			if( counters->data.i[k] != 0 )
			{
				double scale = 1./counters->data.i[k];
				for( j = 0; j < dims; j++ )
				{
					c[j] *= scale;
				}
			}
		}
	}

	// Release memory and exit
	if (centers_arr == 0)
		cvReleaseMat( &centers );
	if (ElemCount == 0)
		cvReleaseMat( &counters);

	cvReleaseMat( &old_centers );

	if (0 != imgResult)
	{
		int pos = 0;
		int x = 0, y = 0, i = 0, j = 0, k = 0;
		for (y = 0; y < imgResult->height; y ++)
		{
			for (x = 0; x < imgResult->width; x ++)
			{
				k = CV_IMAGE_ELEM(labels, int, y, x);
				double* dp = &CV_MAT_ELEM(*centers_arr, double, k, 0);
				unsigned char* uc = &CV_IMAGE_ELEM(imgResult, unsigned char, y, x * 3);
				uc[0] = (int)(dp[0]);
				uc[1] = (int)(dp[1]);
				uc[2] = (int)(dp[2]);
			}
		}
	}

	return cluster_count;
}

// label the regions from color quantization result
int GetRegion(IplImage* imgSeg, IplImage* RegionLabel, int Noise)
{
	int pos = 0;
	int x = 0, y = 0, i = 0, j = 0, k = 0;

	// mark region
	int	nNumRegion = 0;
	unsigned char* pDst;

	for (y = 0; y < RegionLabel->height; y ++)
	{
		for (x = 0; x < RegionLabel->width; x ++)
		{
			CV_IMAGE_ELEM(RegionLabel, int, y, x) = -1;
		}
	}

	// Set the mask
	IplImage* mask	= cvCreateImage(cvSize(imgSeg->width + 2, imgSeg->height + 2), 8 ,1);
	cvSetZero(mask);
	pDst = &CV_IMAGE_ELEM(mask, unsigned char, 0, 0);
	memset(pDst, 255, mask->width * sizeof(unsigned char));
	pDst = &CV_IMAGE_ELEM(mask, unsigned char, mask->height - 1, 0);
	memset(pDst, 255, mask->width * sizeof(unsigned char));
	for ( y = 1; y < mask->height - 1; y ++)
	{
		CV_IMAGE_ELEM(mask, unsigned char, y, 0) = 255;
		CV_IMAGE_ELEM(mask, unsigned char, y, mask->width - 1) = 255;
	}

	CvConnectedComp		comp;

	// Floodfill
	CvPoint ptSeed = cvPoint(-1, -1);

	// Mark the regions

	int n = 0;
	for (y = 0; y < RegionLabel->height; y ++)
	{
		for (x = 0; x < RegionLabel->width; x ++)
		{
			if ( CV_IMAGE_ELEM(RegionLabel, int, y, x) == -1)
			{
				ptSeed = cvPoint(x, y);
				cvFloodFill(imgSeg, ptSeed, cvScalarAll(255), cvScalarAll(0), cvScalarAll(0), &comp,
					CV_FLOODFILL_MASK_ONLY , mask);

				if (comp.area < Noise)
				{
					n = 0;
				}
				else
				{
					n = ++ nNumRegion;
				}
				for ( i = comp.rect.y; i < comp.rect.y + comp.rect.height; i ++)
				{
					for (j = comp.rect.x; j < comp.rect.x + comp.rect.width; j ++)
					{
						if (CV_IMAGE_ELEM(RegionLabel, int, i, j) == -1 &&
							CV_IMAGE_ELEM(mask, unsigned char, i+1, j+1) != 0 )
						{
							CV_IMAGE_ELEM(RegionLabel, int, i, j) = n; // regionlable is "1" based, 0 means "noise"
						}
					}
				}
			}
		}
	}

	cvWatershed(imgSeg, RegionLabel);
	// mark the edge points
	double dis = 0, minDis = -1;
	unsigned char *up1, *up2;
	bool flag = 1;
	while (flag)
	{
		flag = 0;
		for (y = 0; y < RegionLabel->height; y ++)
		{
			for (x = 0; x < RegionLabel->width; x ++)
			{
				if (CV_IMAGE_ELEM(RegionLabel, int, y, x) <= 0)
				{
					minDis	= 3 * 255 * 255;
					up1 = &CV_IMAGE_ELEM(imgSeg, unsigned char, y, 3 * x);
					n = -1;
					// left neighbor
					if ((x - 1) >= 0 && (CV_IMAGE_ELEM(RegionLabel, int, y, x - 1) > 0))
					{
						up2 = &CV_IMAGE_ELEM(imgSeg, unsigned char, y, 3 * (x - 1));
						dis = 0;
						dis += (double)(up1[0] - up2[0]) * (up1[0] - up2[0]);
						dis += (double)(up1[1] - up2[1]) * (up1[1] - up2[1]);
						dis += (double)(up1[2] - up2[2]) * (up1[2] - up2[2]);
						if (dis < minDis)
						{
							minDis	= dis;
							n		= CV_IMAGE_ELEM(RegionLabel, int, y, x - 1);
						}
					}
					// right neighbor
					if ((x + 1) < imgSeg->width && (CV_IMAGE_ELEM(RegionLabel, int, y, x + 1) > 0))
					{
						up2 = &CV_IMAGE_ELEM(imgSeg, unsigned char, y, 3 * (x + 1));
						dis = 0;
						dis += (double)(up1[0] - up2[0]) * (up1[0] - up2[0]);
						dis += (double)(up1[1] - up2[1]) * (up1[1] - up2[1]);
						dis += (double)(up1[2] - up2[2]) * (up1[2] - up2[2]);
						if (dis < minDis)
						{
							minDis	= dis;
							n		= CV_IMAGE_ELEM(RegionLabel, int, y, x + 1);
						}
					}
					// up neighbor
					if ((y - 1) >= 0 && (CV_IMAGE_ELEM(RegionLabel, int, y - 1, x) > 0))
					{
						up2 = &CV_IMAGE_ELEM(imgSeg, unsigned char, y - 1, 3 * x);
						dis = 0;
						dis += (double)(up1[0] - up2[0]) * (up1[0] - up2[0]);
						dis += (double)(up1[1] - up2[1]) * (up1[1] - up2[1]);
						dis += (double)(up1[2] - up2[2]) * (up1[2] - up2[2]);
						if (dis < minDis)
						{
							minDis	= dis;
							n		= CV_IMAGE_ELEM(RegionLabel, int, y - 1, x);
						}
					}
					// down neighbor
					if ((y + 1) < imgSeg->height && (CV_IMAGE_ELEM(RegionLabel, int, y + 1, x) > 0))
					{
						up2 = &CV_IMAGE_ELEM(imgSeg, unsigned char, y + 1, 3 * x);
						dis = 0;
						dis += (double)(up1[0] - up2[0]) * (up1[0] - up2[0]);
						dis += (double)(up1[1] - up2[1]) * (up1[1] - up2[1]);
						dis += (double)(up1[2] - up2[2]) * (up1[2] - up2[2]);
						if (dis < minDis)
						{
							minDis	= dis;
							n		= CV_IMAGE_ELEM(RegionLabel, int, y + 1, x);
						}
					}

					if (n == -1)
					{
						flag = 1;
					}
					else
					{
						CV_IMAGE_ELEM(RegionLabel, int, y, x) = n;
					}
				}
			}
		}
	}

	nNumRegion	= 0;
	for (y = 0; y < RegionLabel->height; y ++)
	{
		for (x = 0; x < RegionLabel->width; x ++)
		{
			nNumRegion	= max(nNumRegion, CV_IMAGE_ELEM(RegionLabel, int, y, x));
			CV_IMAGE_ELEM(RegionLabel, int, y, x) = CV_IMAGE_ELEM(RegionLabel, int, y, x) - 1;
		}
	}

	cvReleaseImage(&mask);
	return nNumRegion;
}

bool CalcRegionContrast(const IplImage* img, IplImage* RegionLabel, IplImage* Contrast, 
						int nNumRegion, double alpha, double lambda)
{
	cvZero(Contrast);

	int RegionCount = nNumRegion;
	int NumOfPixel	= RegionLabel->height * RegionLabel->width;

	CvMat* matRadius	= cvCreateMat(RegionCount, 1, CV_64FC1);
	CvMat* matDistance	= cvCreateMat(RegionCount, RegionCount, CV_64FC1);
	CvMat* matFeatureDis= cvCreateMat(RegionCount, RegionCount, CV_64FC1);
	CvMat* matCenter	= cvCreateMat(RegionCount, 2, CV_32SC1);
	CvMat* matArea		= cvCreateMat(RegionCount, 1, CV_64FC1);
	CvMat* matContrast	= cvCreateMat(RegionCount, 1, CV_64FC1);
	CvMat* matRegionFea	= cvCreateMat(RegionCount, 3, CV_64FC1);
	cvZero(matRadius);
	cvZero(matDistance);
	cvZero(matFeatureDis);
	cvZero(matCenter);
	cvZero(matArea);
	cvZero(matRegionFea);
	cvZero(matContrast);

	int		x = 0, y = 0, i = 0, j = 0, k1 = 0, k2 = 0;
	double	dRadius_2	= 0;
	double	dDis_2		= 0;
	double	dFeaDis_2	= 0;
	CvPoint pt1, pt2;
	double dArea = 0;
	double *dpFea1 = 0, *dpFea2 = 0;

	unsigned char *ucp;
	for(y = 0; y < RegionLabel->height; y ++)
	{
		for(x = 0; x < RegionLabel->width; x ++)
		{
			k1 = CV_IMAGE_ELEM(RegionLabel, int, y, x);
			if (k1 < 0)
			{
				cvReleaseMat(&matRadius);
				cvReleaseMat(&matDistance);
				cvReleaseMat(&matFeatureDis);
				cvReleaseMat(&matCenter);
				cvReleaseMat(&matArea);
				cvReleaseMat(&matContrast);
				cvReleaseMat(&matRegionFea);
				return 0;
			}

			CV_MAT_ELEM(*matArea, double, k1, 0)	= CV_MAT_ELEM(*matArea, double, k1, 0) + 1;
			CV_MAT_ELEM(*matCenter, int, k1, 0)		= CV_MAT_ELEM(*matCenter, int, k1, 0) + x;
			CV_MAT_ELEM(*matCenter, int, k1, 1)		= CV_MAT_ELEM(*matCenter, int, k1, 1) + y;

			ucp			= &CV_IMAGE_ELEM(img, unsigned char, y, x * 3);
			dpFea1		= &CV_MAT_ELEM(*matRegionFea, double, k1, 0);
			dpFea1[0]	+= ucp[0];
			dpFea1[1]	+= ucp[1];
			dpFea1[2]	+= ucp[2];
		}
	}

	// Calculate the area and center of each region
	// Calculate the spatial distance and feature distance of each region
	for(i = 0; i < RegionCount; i ++)
	{
		dArea	= CV_MAT_ELEM(*matArea, double, i, 0);
		dRadius_2	= dArea / pi;
		CV_MAT_ELEM(*matRadius, double, i, 0)	= dRadius_2;

		if (0 == dArea)
		{
			//			AfxMessageBox("dArea == 0");
			continue;
		}
		CV_MAT_ELEM(*matCenter, int, i, 0) = (int)(CV_MAT_ELEM(*matCenter, int, i, 0) / dArea);
		CV_MAT_ELEM(*matCenter, int, i, 1) = (int)(CV_MAT_ELEM(*matCenter, int, i, 1) / dArea);

		dpFea1	= &CV_MAT_ELEM(*matRegionFea, double, i, 0);
		dpFea1[0] = dpFea1[0] / dArea;
		dpFea1[1] = dpFea1[1] / dArea;
		dpFea1[2] = dpFea1[2] / dArea;

		pt1 = cvPoint(CV_MAT_ELEM(*matCenter, int, i, 0),  CV_MAT_ELEM(*matCenter, int, i, 1));

		for(j = 0; j < i; j ++)
		{
			pt2			= cvPoint(CV_MAT_ELEM(*matCenter, int, j, 0),  CV_MAT_ELEM(*matCenter, int, j, 1));
			dpFea2		= &CV_MAT_ELEM(*matRegionFea, double, j, 0);

			dDis_2		= pow((double)(pt1.x - pt2.x), 2) + pow((double)(pt1.y - pt2.y), 2);
			CV_MAT_ELEM(*matDistance, double, i, j)		= CV_MAT_ELEM(*matDistance, double, j, i)	= dDis_2;

			dFeaDis_2	= pow(dpFea1[0] - dpFea2[0], 2) + pow(dpFea1[1] - dpFea2[1], 2) + pow(dpFea1[2] - dpFea2[2], 2);
			CV_MAT_ELEM(*matFeatureDis, double, i, j)	= CV_MAT_ELEM(*matFeatureDis, double, j, i)	=dFeaDis_2;
		}
	}

	const double lambda_2	= lambda * lambda;	// the ratio between center and surround
	//	const double alpha_2	= alpha * alpha;	// the ratio between region radius and sigma
	const double alpha_2	= (1 - 1 / lambda_2) / (2 * log(lambda_2));

	cvScale(matArea, matArea, 1 / (double)NumOfPixel, 0);

	// Calculate the contrast of each region
	double sigma_2, dContrast	= 0, minDis;
	for(i = 0; i < RegionCount; i ++)
	{
		dRadius_2	= CV_MAT_ELEM(*matRadius, double, i, 0);
		if (dRadius_2 == 0)
		{
			continue;
		}

		sigma_2 = dRadius_2 * alpha_2;
		dContrast	= 0;
		minDis = 4 * log(lambda_2) * sigma_2 / ( 1 - 1 / lambda_2);

		for(j = 0; j < RegionCount; j ++)
		{
			dFeaDis_2			= CV_MAT_ELEM(*matFeatureDis, double, i, j);
			dArea				= CV_MAT_ELEM(*matArea, double, j, 0);
			if (dArea == 0)
			{
				continue;
			}
			double dRadius_j_2	= CV_MAT_ELEM(*matRadius, double, j, 0);

			dDis_2		= CV_MAT_ELEM(*matDistance, double, i, j);
			dDis_2		= pow(sqrt(dDis_2) - sqrt(dRadius_j_2), 2);
			if (dDis_2 < minDis) 
			{
				dDis_2 = minDis;
			}
			dDis_2		= exp(-dDis_2 / (2 * sigma_2 * lambda_2)) / lambda_2 - exp(-dDis_2 / (2 * sigma_2));
			//			CV_MAT_ELEM(*matDistance, double, j, i)	= dDis_2;
			dContrast	+= dFeaDis_2 * dDis_2 * dArea;
		}
		dContrast	/= sigma_2;
		CV_MAT_ELEM(*matContrast, double, i, 0)		= dContrast;
	}

	double minValue = 0, maxValue = 0, scale, shift;
	cvMinMaxLoc(matContrast, &minValue, &maxValue);

	if (minValue == maxValue)
	{
		cvReleaseMat(&matRadius);
		cvReleaseMat(&matDistance);
		cvReleaseMat(&matFeatureDis);
		cvReleaseMat(&matCenter);
		cvReleaseMat(&matArea);
		cvReleaseMat(&matContrast);
		cvReleaseMat(&matRegionFea);
		return 0;
	}

	scale	= 1 / (maxValue - minValue);
	shift	= - minValue * scale;
	cvScale(matContrast, matContrast, scale, shift);
	double MeanContrast = 0, TotalArea = 0;
	for (int k = 0; k < RegionCount; k ++)
	{
		dArea		= CV_MAT_ELEM(*matArea, double, k, 0);
		MeanContrast+=CV_MAT_ELEM(*matContrast, double, k, 0) * dArea;
		TotalArea	+= dArea;
	}
	MeanContrast	/= TotalArea;
	scale	= 255 * 0.5 / MeanContrast;
	cvScale(matContrast, matContrast, scale, 0);

	for(y = 0; y < RegionLabel->height; y ++)
	{
		for(x = 0; x < RegionLabel->width; x ++)
		{
			k1	= CV_IMAGE_ELEM(RegionLabel, int, y, x);

			CV_IMAGE_ELEM(Contrast, unsigned char, y, x)	= max(0, min(255, (int)(CV_MAT_ELEM(*matContrast, double, k1, 0))));
		}
	}

	// Release memory
	cvReleaseMat(&matRadius);
	cvReleaseMat(&matDistance);
	cvReleaseMat(&matFeatureDis);
	cvReleaseMat(&matCenter);
	cvReleaseMat(&matArea);
	cvReleaseMat(&matContrast);
	cvReleaseMat(&matRegionFea);

	return 1;

}

bool CalcRegionContrast_float(const IplImage* img, IplImage* RegionLabel, IplImage* Contrast, 
							  int nNumRegion, double alpha, double lambda)
{
	cvZero(Contrast);

	int RegionCount = nNumRegion;
	int NumOfPixel	= RegionLabel->height * RegionLabel->width;

	CvMat* matRadius	= cvCreateMat(RegionCount, 1, CV_64FC1);
	CvMat* matDistance	= cvCreateMat(RegionCount, RegionCount, CV_64FC1);
	CvMat* matFeatureDis= cvCreateMat(RegionCount, RegionCount, CV_64FC1);
	CvMat* matCenter	= cvCreateMat(RegionCount, 2, CV_32SC1);
	CvMat* matArea		= cvCreateMat(RegionCount, 1, CV_64FC1);
	CvMat* matContrast	= cvCreateMat(RegionCount, 1, CV_64FC1);
	CvMat* matRegionFea	= cvCreateMat(RegionCount, 3, CV_64FC1);
	cvZero(matRadius);
	cvZero(matDistance);
	cvZero(matFeatureDis);
	cvZero(matCenter);
	cvZero(matArea);
	cvZero(matRegionFea);
	cvZero(matContrast);

	int		x = 0, y = 0, i = 0, j = 0, k1 = 0, k2 = 0;
	double	dRadius_2	= 0;
	double	dDis_2		= 0;
	double	dFeaDis_2	= 0;
	CvPoint pt1, pt2;
	double dArea = 0;
	double *dpFea1 = 0, *dpFea2 = 0;

	unsigned char *ucp;
	for(y = 0; y < RegionLabel->height; y ++)
	{
		for(x = 0; x < RegionLabel->width; x ++)
		{
			k1 = CV_IMAGE_ELEM(RegionLabel, int, y, x);
			if (k1 < 0)
			{
				cvReleaseMat(&matRadius);
				cvReleaseMat(&matDistance);
				cvReleaseMat(&matFeatureDis);
				cvReleaseMat(&matCenter);
				cvReleaseMat(&matArea);
				cvReleaseMat(&matContrast);
				cvReleaseMat(&matRegionFea);
				return 0;
			}

			CV_MAT_ELEM(*matArea, double, k1, 0)	= CV_MAT_ELEM(*matArea, double, k1, 0) + 1;
			CV_MAT_ELEM(*matCenter, int, k1, 0)		= CV_MAT_ELEM(*matCenter, int, k1, 0) + x;
			CV_MAT_ELEM(*matCenter, int, k1, 1)		= CV_MAT_ELEM(*matCenter, int, k1, 1) + y;

			ucp			= &CV_IMAGE_ELEM(img, unsigned char, y, x * 3);
			dpFea1		= &CV_MAT_ELEM(*matRegionFea, double, k1, 0);
			dpFea1[0]	+= ucp[0];
			dpFea1[1]	+= ucp[1];
			dpFea1[2]	+= ucp[2];
		}
	}

	// Calculate the area and center of each region
	// Calculate the spatial distance and feature distance of each region
	for(i = 0; i < RegionCount; i ++)
	{
		dArea	= CV_MAT_ELEM(*matArea, double, i, 0);
		dRadius_2	= dArea / pi;
		CV_MAT_ELEM(*matRadius, double, i, 0)	= dRadius_2;

		if (0 == dArea)
		{
			//			AfxMessageBox("dArea == 0");
			continue;
		}
		CV_MAT_ELEM(*matCenter, int, i, 0) = (int)(CV_MAT_ELEM(*matCenter, int, i, 0) / dArea);
		CV_MAT_ELEM(*matCenter, int, i, 1) = (int)(CV_MAT_ELEM(*matCenter, int, i, 1) / dArea);

		dpFea1	= &CV_MAT_ELEM(*matRegionFea, double, i, 0);
		dpFea1[0] = dpFea1[0] / dArea;
		dpFea1[1] = dpFea1[1] / dArea;
		dpFea1[2] = dpFea1[2] / dArea;

		pt1 = cvPoint(CV_MAT_ELEM(*matCenter, int, i, 0),  CV_MAT_ELEM(*matCenter, int, i, 1));

		for(j = 0; j < i; j ++)
		{
			pt2			= cvPoint(CV_MAT_ELEM(*matCenter, int, j, 0),  CV_MAT_ELEM(*matCenter, int, j, 1));
			dpFea2		= &CV_MAT_ELEM(*matRegionFea, double, j, 0);

			dDis_2		= pow((double)(pt1.x - pt2.x), 2) + pow((double)(pt1.y - pt2.y), 2);
			CV_MAT_ELEM(*matDistance, double, i, j)		= CV_MAT_ELEM(*matDistance, double, j, i)	= dDis_2;

			dFeaDis_2	= pow(dpFea1[0] - dpFea2[0], 2) + pow(dpFea1[1] - dpFea2[1], 2) + pow(dpFea1[2] - dpFea2[2], 2);
			CV_MAT_ELEM(*matFeatureDis, double, i, j)	= CV_MAT_ELEM(*matFeatureDis, double, j, i)	=dFeaDis_2;
		}
	}

	const double lambda_2	= lambda * lambda;	// the ratio between center and surround
	//	const double alpha_2	= alpha * alpha;	// the ratio between region radius and sigma
	const double alpha_2	= (1 - 1 / lambda_2) / (2 * log(lambda_2));

	cvScale(matArea, matArea, 1 / (double)NumOfPixel, 0);

	// Calculat the contrast of each region
	double sigma_2, dContrast	= 0, minDis;
	for(i = 0; i < RegionCount; i ++)
	{
		dRadius_2	= CV_MAT_ELEM(*matRadius, double, i, 0);
		if (dRadius_2 == 0)
		{
			continue;
		}

		sigma_2 = dRadius_2 * alpha_2;
		dContrast	= 0;
		minDis = 4 * log(lambda_2) * sigma_2 / ( 1 - 1 / lambda_2);

		for(j = 0; j < RegionCount; j ++)
		{
			dFeaDis_2			= CV_MAT_ELEM(*matFeatureDis, double, i, j);
			dArea				= CV_MAT_ELEM(*matArea, double, j, 0);
			if (dArea == 0)
			{
				continue;
			}
			double dRadius_j_2	= CV_MAT_ELEM(*matRadius, double, j, 0);

			dDis_2		= CV_MAT_ELEM(*matDistance, double, i, j);
			dDis_2		= pow(sqrt(dDis_2) - sqrt(dRadius_j_2), 2);
			if (dDis_2 < minDis) 
			{
				dDis_2 = minDis;
			}
			dDis_2		= exp(-dDis_2 / (2 * sigma_2 * lambda_2)) / lambda_2 - exp(-dDis_2 / (2 * sigma_2));
			//			CV_MAT_ELEM(*matDistance, double, j, i)	= dDis_2;
			dContrast	+= dFeaDis_2 * dDis_2 * dArea;
		}
		dContrast	/= sigma_2;
		CV_MAT_ELEM(*matContrast, double, i, 0)		= dContrast;
	}

	double minValue = 0, maxValue = 0;
	cvMinMaxLoc(matContrast, &minValue, &maxValue);

	if (minValue == maxValue)
	{
		cvReleaseMat(&matRadius);
		cvReleaseMat(&matDistance);
		cvReleaseMat(&matFeatureDis);
		cvReleaseMat(&matCenter);
		cvReleaseMat(&matArea);
		cvReleaseMat(&matContrast);
		cvReleaseMat(&matRegionFea);
		return 0;
	}

	for(y = 0; y < RegionLabel->height; y ++)
	{
		for(x = 0; x < RegionLabel->width; x ++)
		{
			k1	= CV_IMAGE_ELEM(RegionLabel, int, y, x);

			CV_IMAGE_ELEM(Contrast, float, y, x)	= CV_MAT_ELEM(*matContrast, double, k1, 0);
		}
	}

	// Release memory
	cvReleaseMat(&matRadius);
	cvReleaseMat(&matDistance);
	cvReleaseMat(&matFeatureDis);
	cvReleaseMat(&matCenter);
	cvReleaseMat(&matArea);
	cvReleaseMat(&matContrast);
	cvReleaseMat(&matRegionFea);

	return 1;

}

bool CalcInformation(IplImage* FeatureLabel, CvMat* ElemCount, int cluster_count,
					 IplImage* InformationMap)
{
	CvMat*	matInformation = cvCreateMat(cluster_count, 1, CV_64FC1);
	cvZero(matInformation);
	int x = 0, y = 0, k = 0;
	double dTotalNum = (double)(FeatureLabel->height * FeatureLabel->width);
	double p = 0;
	double info = 0;

	for(k = 0; k < cluster_count; k ++)
	{
		if (CV_MAT_ELEM(*ElemCount, int, k, 0) < 0)	return false;

		p		= CV_MAT_ELEM(*ElemCount, int, k, 0) / dTotalNum;
		info = 0;
		if (!(p == 0))
		{
			info	= - log(p);
		}
		CV_MAT_ELEM(*matInformation, double, k, 0) = info;
	}
	double minValue = 0, maxValue = 0, scale, shift;

	cvMinMaxLoc(matInformation, &minValue, &maxValue);
	if (minValue == maxValue)
	{
		cvZero(matInformation);
		cvReleaseMat(&matInformation);
		return false;
	}
	scale	= 255 / (maxValue - minValue);
	shift	= - minValue * scale;
	cvScale(matInformation, matInformation, scale, shift);

	for(y = 0; y < FeatureLabel->height; y ++)
	{
		for(x = 0; x < FeatureLabel->width; x ++)
		{
			k = CV_IMAGE_ELEM(FeatureLabel, int, y, x);
			CV_IMAGE_ELEM(InformationMap, unsigned char, y, x) = (int)(CV_MAT_ELEM(*matInformation, double, k, 0));
		}
	}

	cvReleaseMat(&matInformation);
	return true;
}

//////////////////////////////////////////////////////////////////////////

int	MinVar(int hist[256])
{
	int threshold = 128;
	int i=0, j=0;
	int num = 0, num1=0, num2=0;
	int sum = 0, sum1=0, sum2=0;
	double mean = 0, mean1=0, mean2=0;
	double var=0;
	double MinVariance = 0;
	for (i=0; i<256; i++)
	{
		num  += hist[i];
		sum  += hist[i] * i;
	}
	if (num == 0)
	{
		return -1;
	}
	mean = sum/num;
	for (i=0; i<256; i++)
	{
		MinVariance += hist[i] * (i-mean) * (i-mean);
		threshold = -1;
	}
	for (i=0; i<255; i++)
	{
		if (hist[i] == 0) continue;
		sum1 += i*hist[i];
		sum2 = sum - sum1;
		num1 += hist[i];
		num2 = num - num1;
		if (num1 == 0 || num2 == 0) continue;
		mean1=sum1/num1;
		mean2=sum2/num2;
		var=0;
		for(j=0; j<=i; j++)
		{
			if(hist[j] == 0) continue;
			var += hist[j] * (j-mean1) * (j-mean1);
		}
		for(j=i+1; j<256; j++)
		{
			if(hist[j] == 0) continue;
			var += hist[j] * (j-mean2) * (j-mean2);
		}
		if (var <= MinVariance )
		{
			MinVariance = var;
			threshold = i;
		}
	}
	return threshold;
}
int		MinVar(int* hist, int bin, double &Sw)
{
	int threshold = 128;
	int i=0, j=0;
	int num = 0, num1=0, num2=0;
	int sum = 0, sum1=0, sum2=0;
	double mean = 0, mean1=0, mean2=0;
	double var=0;
	double MinVariance = 0;
	for (i=0; i<bin; i++)
	{
		num  += hist[i];
		sum  += hist[i] * i;
	}
	if (num == 0)
	{
		Sw = -1;
		return -1;
	}
	mean = sum/num;
	for (i=0; i<bin; i++)
	{
		MinVariance += hist[i] * (i-mean) * (i-mean);
		threshold = -1;
	}
	for (i=0; i < bin - 1; i++)
	{
		if (hist[i] == 0) continue;
		sum1 += i*hist[i];
		sum2 = sum - sum1;
		num1 += hist[i];
		num2 = num - num1;
		if (num1 == 0 || num2 == 0) continue;
		mean1=sum1/num1;
		mean2=sum2/num2;
		var=0;
		for(j=0; j<=i; j++)
		{
			if(hist[j] == 0) continue;
			var += hist[j] * (j-mean1) * (j-mean1);
		}
		for(j=i+1; j<bin; j++)
		{
			if(hist[j] == 0) continue;
			var += hist[j] * (j-mean2) * (j-mean2);
		}
		if (var <= MinVariance )
		{
			MinVariance = var;
			threshold = i;
		}
	}
	Sw = MinVariance	/ num;
	return threshold;

}
int		MinVar(CvHistogram *hist, int bin, double &Sw)
{
	int threshold = 128;
	int i=0, j=0;
	double num = 0, num1=0, num2=0;
	double sum = 0, sum1=0, sum2=0;
	double mean = 0, mean1=0, mean2=0;
	double var=0;
	double MinVariance = 0;
	for (i=0; i<bin; i++)
	{
		num  += cvQueryHistValue_1D(hist, i);
		sum  += cvQueryHistValue_1D(hist, i) * i;
	}
	if (num == 0)
	{
		Sw = -1;
		return -1;
	}
	mean = sum/num;
	for (i=0; i < bin; i++)
	{
		MinVariance += cvQueryHistValue_1D(hist, i) * (i-mean) * (i-mean);
		threshold = -1;
	}
	for (i=0; i < bin - 1; i++)
	{
		if (cvQueryHistValue_1D(hist, i) == 0) continue;
		sum1 += i * cvQueryHistValue_1D(hist, i);
		sum2 = sum - sum1;
		num1 += cvQueryHistValue_1D(hist, i);
		num2 = num - num1;
		if (num1 == 0 || num2 == 0) continue;
		mean1=sum1/num1;
		mean2=sum2/num2;
		var=0;
		for(j=0; j<=i; j++)
		{
			if(cvQueryHistValue_1D(hist, j) == 0) continue;
			var += cvQueryHistValue_1D(hist, j) * (j-mean1) * (j-mean1);
		}
		for(j=i+1; j<256; j++)
		{
			if(cvQueryHistValue_1D(hist, j) == 0) continue;
			var += cvQueryHistValue_1D(hist, j) * (j-mean2) * (j-mean2);
		}
		if (var <= MinVariance )
		{
			MinVariance = var;
			threshold = i;
		}
	}

	Sw = MinVariance / num;
	return threshold;
}

// Calculate the maximum variance between clusters
// The return value is the threshold
int		MaxSb(int* hist, int bin, double &Sb)
{
	int k = 0, threshold = 0;
	int total_num	= 0, total_num_1	= 0, total_num_2 = 0;;
	double miu = 0, miu_1 = 0, miu_2 = 0;
	double max_sdv = 0, sdv	= 0;

	for (k = 0; k < bin; k ++)
	{
		miu			+= k * hist[k];
		total_num	+= hist[k];
	}
	miu_2 = miu;
	miu	/= total_num;
	for (k = 0; k < bin; k ++)
	{
		miu_1		+= k * hist[k];
		total_num_1	+= hist[k];
		miu_2		-=  k * hist[k];
		total_num_2	=  total_num - total_num_1;

		sdv			= total_num_1 * pow(miu - miu_1 / total_num_1, 2.0) + total_num_2 * pow(miu - miu_2 / total_num_2, 2.0);
		if (sdv	> max_sdv)
		{
			max_sdv		= sdv;
			threshold	= k;
		}
	}

	Sb	= max_sdv / total_num;
	if (0 == Sb)
	{
		return -1;
	}
	return threshold;
}
// <= threshold
int		MaxSb(CvHistogram* hist, int bin, double &Sb)
{
	int k = 0, threshold = 0;
	double total_num	= 0, total_num_1	= 0, total_num_2 = 0;;
	double miu = 0, miu_1 = 0, miu_2 = 0;
	double max_sdv = 0, sdv	= 0;

	for (k = 0; k < bin; k ++)
	{
		miu			+= k * cvQueryHistValue_1D(hist, k);
		total_num	+= cvQueryHistValue_1D(hist, k);
	}
	miu_2	= miu;
	miu	/= total_num;
	for (k = 0; k < bin; k ++)
	{
		miu_1		+= k * cvQueryHistValue_1D(hist, k);
		total_num_1	+= cvQueryHistValue_1D(hist, k);
		miu_2		-=  k * cvQueryHistValue_1D(hist, k);
		total_num_2	=  total_num - total_num_1;

		sdv			= total_num_1 * pow(miu - miu_1 / total_num_1, 2.0) + total_num_2 * pow(miu - miu_2 / total_num_2, 2.0);
		if (sdv	> max_sdv)
		{
			max_sdv		= sdv;
			threshold	= k;
		}
	}

	Sb	= max_sdv / total_num;
	if (0 == Sb)
	{
		return -1;
	}
	return threshold;
}

// This funciton detect the ROI from the saliency map, by claculating the original moment and central moment
// reference: Contrast based ******,MSRA, Ma Yufei
CvRect	GetROIMoment(IplImage* img, double alpha)
{
	CvRect rect = cvRect(0, 0, 0, 0);
	if (1 != img->nChannels)
	{
		return rect;
	}
	int x = 0, y = 0;
	double x0 = 0, y0 = 0;
	double total = 0;
	double w = 0, h = 0;
	double saliency = 0;
	for (y = 0; y < img->height; y ++)
	{
		for (x = 0; x < img->width; x ++)
		{
			saliency	= (double)(CV_IMAGE_ELEM(img, unsigned char, y, x));
			total		+= saliency;
			x0			+= saliency * x;
			y0			+= saliency * y;
		}
	}
	if (total == 0)
	{
		return rect;
	}
	x0 /= total;
	y0 /= total;

	for (y = 0; y < img->height; y ++)
	{
		for (x = 0; x < img->width; x ++)
		{
			saliency	= (double)(CV_IMAGE_ELEM(img, unsigned char, y, x));
			w			+= saliency * fabs((double)x - x0);
			h			+= saliency * fabs((double)y - y0);
		}
	}
	w	/= total;
	h	/= total;

	int left	= (int)(max(x0 - w * alpha, 0));
	int right	= (int)(min(x0 + w * alpha, img->width - 1));
	int top		= (int)(max(y0 - h * alpha, 0));
	int bottom	= (int)(min(y0 + h * alpha, img->height - 1));
	rect.x = left;
	rect.y = top;
	rect.width	= right - left + 1;
	rect.height = bottom- top + 1;

	return rect;
}

bool	ExpandRect(IplImage* img, CvRect &rect, double thres_1, double thres_2, int step)
{
	//	IplImage* imgColor = cvCreateImage(cvGetSize(img), 8, 3);
	//	cvMerge(img, img, img, 0, imgColor);

	if (1 != img->nChannels)
	{
		return 0;
	}

	// Initialize the rect using moment method
	if (rect.width == 0 || rect.height == 0)
	{
		return 0;
	}
	int left	= rect.x;
	int right	= left + rect.width - 1;
	int	top		= rect.y;
	int bottom	= top + rect.height -1;
	//	cvRectangle(imgColor, cvPoint(left, top), cvPoint(right, bottom), CV_RGB(255, 255, 255));
	//	cvNamedWindow("image", 1);
	//	cvShowImage("image", imgColor);
	//	cvWaitKey(0);

	//FILE *fp;
	//	fp = fopen("E:\\rect.txt", "w");
	//	fprintf(fp, "%d\t%d%d\t%d", left, top, right, bottom);
	//	fclose(fp);

	// Expand
	//	double MeanSaliency = 0;
	int x = 0, y = 0;
	//	for (y = top; y <= bottom; y ++)
	//	{
	//		for (x = left; x <= right; x ++)
	//		{
	//			MeanSaliency += (double)(CV_IMAGE_ELEM(img, unsigned char, y, x));
	//		}
	//	}
	//	MeanSaliency /= (double)(rect.width * rect.height);
	//	double thres = (MeanSaliency * thres_1);
	//	thres = thres * thres;
	double thres	= thres_1 * thres_1;

	double E_In = 0, E_Out = 0, EP = 0;

	// expand left
	//fp	= fopen("E:\\left.txt", "a");
	for (y = top; y <= bottom; y ++)
	{
		for (x = left; x < min(left + step, img->width); x ++)
		{
			E_In	+= (double)(CV_IMAGE_ELEM(img, unsigned char, y, x)) / 255;
		}
	}
	E_In /= ((double)step * (bottom - top + 1));
	while (left > 0)
	{
		E_Out	= 0;
		for (y = top; y <= bottom; y ++)
		{
			for (x = left - 1; x >= max(left - step, 0); x --)
			{
				E_Out	+= (double)(CV_IMAGE_ELEM(img, unsigned char, y, x)) / 255;
			}
		}
		E_Out /= ((double)step * (bottom - top + 1));
		EP = E_In * E_Out;
		if (EP < thres || (E_Out / E_In < thres_2)) 
		{
			//fprintf(fp, "%f\t%f\n", E_In * E_Out, E_Out / E_In);
			break;
		}
		left -= step;
		E_In = E_Out;
		//		cvLine(imgColor, cvPoint(left, top), cvPoint(left, bottom), CV_RGB(255, 0, 0));
		//		cvNamedWindow("image", 1);
		//		cvShowImage("image", imgColor);
		//		cvWaitKey(0);
	}
	//	fclose(fp);

	// expand right
	//	fp = fopen("E:\\right", "w");
	E_In	= 0;
	for (y = top; y <= bottom; y ++)
	{
		for (x = right; x >= max(right - step + 1, 0); x --)
		{
			E_In	+= (double)(CV_IMAGE_ELEM(img, unsigned char, y, x)) / 255;
		}
	}
	E_In /= ((double)step * (bottom - top + 1));
	while (right < img->width)
	{
		E_Out	= 0;
		for (y = top; y <= bottom; y ++)
		{
			for (x = right + 1; x < min(right + step + 1, img->width); x ++)
			{
				E_Out	+=(double)(CV_IMAGE_ELEM(img, unsigned char, y, x)) / 255;
			}
		}
		E_Out /= ((double)step * (bottom - top + 1));
		EP = E_In * E_Out;
		if (EP < thres || (E_Out / E_In < thres_2))
		{
			//fprintf(fp, "%f\t%f\n", E_In * E_Out, E_Out / E_In);
			break;
		}
		right += step;
		E_In = E_Out;
		//		fprintf(fp, "%f\t%f\t%f\t%f\n", E_In, E_Out, E_In * E_Out, E_Out / E_In);
		//		cvLine(imgColor, cvPoint(right, top), cvPoint(right, bottom), CV_RGB(255, 0, 0));
		//		cvNamedWindow("image", 1);
		//		cvShowImage("image", imgColor);
		//		cvWaitKey(0);
	}
	//	fclose(fp);

	// expand up
	//	fp = fopen("E:\\up.txt", "w");
	E_In	= 0;
	for (y = top; y < min(top + step, img->height); y ++)
	{
		for (x = left; x <= right; x ++)
		{
			E_In	+= (double)(CV_IMAGE_ELEM(img, unsigned char, y, x)) / 255;
		}
	}
	E_In /= ((double)step * (right - left + 1));
	while (top > 0)
	{
		E_Out	= 0;
		for (y = top - 1; y >= max(0, top - step); y --)
		{
			for (x = left; x <= right; x ++)
			{
				E_Out	+= (double)(CV_IMAGE_ELEM(img, unsigned char, y, x)) / 255;
			}
		}
		E_Out /= ((double)step * (right - left + 1));
		EP = E_In * E_Out;

		if (EP < thres || (E_Out / E_In < thres_2))
		{
			//fprintf(fp, "%f\t%f\n", E_In * E_Out, E_Out / E_In);
			break;
		}
		top -= step;
		E_In = E_Out;
		//		fprintf(fp, "%f\t%f\t%f\t%f\n", E_In, E_Out, E_In * E_Out, E_Out / E_In);
		//		cvLine(imgColor, cvPoint(left, top), cvPoint(right, top), CV_RGB(255, 0, 0));
		//		cvNamedWindow("image", 1);
		//		cvShowImage("image", imgColor);
		//		cvWaitKey(0);
	}
	//	fclose(fp);

	// expand down
	//	fp = fopen("E:\\down.txt", "w");
	E_In	= 0;
	for (y = bottom; y >= max(bottom - step + 1, 0); y --)
	{
		for (x = left; x <= right; x ++)
		{
			E_In	+= (double)(CV_IMAGE_ELEM(img, unsigned char, y, x)) / 255;
		}
	}
	E_In /= ((double)step * (right - left + 1));
	while (bottom < img->height)
	{
		E_Out	= 0;
		for (y = bottom + 1; y < min(img->height, bottom + step + 1); y ++)
		{
			for (x = left; x <= right; x ++)
			{
				E_Out	+= (double)(CV_IMAGE_ELEM(img, unsigned char, y, x)) / 255;
			}
		}
		E_Out /= ((double)step * (right - left + 1));
		EP = E_In * E_Out;
		if (EP < thres || (E_Out / E_In < thres_2))
		{
			//fprintf(fp, "%f\t%f\n", E_In * E_Out, E_Out / E_In);
			break;
		}
		bottom += step;
		E_In = E_Out;
		//		fprintf(fp, "%f\t%f\t%f\t%f\n", E_In, E_Out, E_In * E_Out, E_Out / E_In);
		//		cvLine(imgColor, cvPoint(left, bottom), cvPoint(right, bottom), CV_RGB(255, 0, 0));
		//		cvNamedWindow("image", 1);
		//		cvShowImage("image", imgColor);
	}
	//fclose(fp);

	rect.x = left;
	rect.y = top;
	rect.width = right - left + 1;
	rect.height= bottom - top + 1;

	//	cvWaitKey(0);
	//	cvReleaseImage(&imgColor);
	return 1;
}

CvRect	GetROICover(IplImage* img, double alpha)
{
	CvRect rect = cvRect(0, 0, 0, 0);
	if (1 != img->nChannels)
	{
		return rect;
	}
	int x = 0, y = 0;
	double x0 = 0, y0 = 0;
	double total = 0;
	double w = 0, h = 0;
	double saliency = 0;
	for (y = 0; y < img->height; y ++)
	{
		for (x = 0; x < img->width; x ++)
		{
			saliency	= (double)(CV_IMAGE_ELEM(img, unsigned char, y, x));
			total		+= saliency;
			x0			+= saliency * x;
			y0			+= saliency * y;
		}
	}
	if (total == 0)
	{
		return rect;
	}
	x0 /= total;
	y0 /= total;

	for (y = 0; y < img->height; y ++)
	{
		for (x = 0; x < img->width; x ++)
		{
			saliency	= (double)(CV_IMAGE_ELEM(img, unsigned char, y, x));
			w			+= saliency * fabs((double)x - x0);
			h			+= saliency * fabs((double)y - y0);
		}
	}
	w	/= total;
	h	/= total;

	int left	= (int)(max(x0 - w , 0));
	int right	= (int)(min(x0 + w , img->width - 1));
	int top		= (int)(max(y0 - h , 0));
	int bottom	= (int)(min(y0 + h , img->height - 1));

	//
	saliency = 0;
	for (y = top; y <= bottom; y ++)
	{
		for (x = left; x <= right; x ++)
		{
			saliency	+= (double)(CV_IMAGE_ELEM(img, unsigned char, y, x));
		}
	}

	double dl = 0, dr = 0, dt = 0, db = 0;
	while (saliency / total < alpha)
	{
		dl = 0;
		dr = 0;
		dt = 0;
		db = 0;
		if (left > 0)
		{
			for (y = top; y <= bottom; y ++)
			{
				dl	+= (double)(CV_IMAGE_ELEM(img, unsigned char, y, left - 1));
			}
			dl /= (double)(bottom - top + 1);
		}

		if (right < img->width - 1)
		{
			for (y = top; y <= bottom; y ++)
			{
				dr	+= (double)(CV_IMAGE_ELEM(img, unsigned char, y, right + 1));
			}
			dr /= (double)(bottom - top + 1);
		}

		if (top > 0)
		{
			for (x = left; x <= right; x ++)
			{
				dt	+= (double)(CV_IMAGE_ELEM(img, unsigned char, top - 1, x));
			}
			dt /= (double)(right - left + 1);
		}

		if (bottom < img->height - 1)
		{
			for (x = left; x <= right; x ++)
			{
				db	+= (double)(CV_IMAGE_ELEM(img, unsigned char, bottom + 1, x));
			}
			db /= (double)(right - left + 1);
		}

		double dmax = max(dl, max(dr, max(dt, db)));
		if (dmax == 0)
		{
			break;
		}
		if (dl == dmax)
		{
			left --;
			saliency += dl * (double)(bottom - top + 1);
		}
		else if (dr == dmax)
		{
			right ++;
			saliency += dr * (double)(bottom - top + 1);
		}
		else if (dt == dmax)
		{
			top --;
			saliency += dt * (double)(right - left + 1);
		}
		else
		{
			bottom ++;
			saliency += db * (double)(right - left + 1);
		}
	}

	rect.x = left;
	rect.y = top;
	rect.width	= right - left + 1;
	rect.height = bottom- top + 1;

	return rect;
}

CvRect	GetROI(IplImage	*img)
{
	CvRect rect = cvRect(0, 0, 0, 0);
	if (1 != img->nChannels)
	{
		return rect;
	}

	int hist_size[] = {256};
	float r_ranges[] = { 0, 255 };
	float* ranges[] = { r_ranges };
	CvHistogram	*h	= cvCreateHist( 1, hist_size, CV_HIST_ARRAY, ranges, 1 );
	cvClearHist(h);

	IplImage* planes[] = { img};

	cvCalcHist(planes, h, 0, NULL);

	CvMat* mat = (CvMat*)(h->bins);
	int	thres = MinVar((int*)(mat->data.ptr));

	IplImage* mask = cvCreateImage(cvGetSize(img), 8 ,1);
	cvThreshold(img, mask, (double)thres, 255, CV_THRESH_BINARY);

	rect	= cvBoundingRect(mask);

	cvReleaseHist(&h);
	cvReleaseImage(&mask);
	return	rect;
}

// Choose the two thresholds for fuzzy growing by minimize the difference between
// the prior entropy of the two fuzzy sets
// 
void DoubleThreshold(CvHistogram* hist, int bin, int& low, int& high, bool bNorm)
{
	CvHistogram* h;

	if (!bNorm)
	{
		int hist_size[] = {bin};
		float rangesGray[] = { 0, 255 };
		float* ranges[] = { rangesGray };
		CvHistogram* h = cvCreateHist( 1, hist_size, CV_HIST_ARRAY, ranges, 1 );
		cvCopyHist(hist, &h);
		cvNormalizeHist(h, 1);
	}
	else
	{
		h = hist;
	}
	//////////////////////////////////////////////////////////////////////////

	int i, j, k;
	int t1 = 0, t2 = bin - 1, d;
	double p1,  p2, p;
	double p10 =0, p20= 0;
	double pk1, pk2;
	double H1 = 0, H2 = 0;
	double H10= 0, H20= 0;
	double MinErr = 0, Err;

	d = t2 - t1;
	p1 = 0;
	p2 = 0;
	for (k = 0; k < bin; k ++)
	{
		p = cvQueryHistValue_1D(h, k);
		if (p == 0) continue;

		pk1	 = p * (t2 - k) / d;
		p1  += pk1;
	}
	p2 = 1 - p1;
	for (k = 0; k < bin; k ++)
	{
		p = cvQueryHistValue_1D(h, k);
		if (p == 0) continue;

		pk1	 = p * (t2 - k) / d;
		pk2	 = p - pk1;
		if (pk1 != 0)	
		{
			H1  -= (pk1 / p1) * log(pk1 / p1);
		}
		if (pk2 != 0)
		{
			H2  -= (pk2 / p2) * log(pk2 / p2);
		}
	}
	MinErr	= (H1 - H2) * (H1 - H2);

	for (i = 0; i < bin; i ++)
	{
		p = cvQueryHistValue_1D(h, i);
		if (p == 0) continue;

		p10 += p;
		for (j = bin - 1; j > i; j --)
		{
			d	= j - i;
			p1	= p10;
			for (k = i + 1; k < j; k ++)
			{
				p = cvQueryHistValue_1D(h, k);
				if (p == 0) continue;

				pk1	 = p * (j - k) / d;
				p1  += pk1;
			}
			p2 = 1 - p1;

			H1 = 0;
			H2 = 0;
			for (k = 0; k < bin; k ++)
			{
				p = cvQueryHistValue_1D(h, k);
				if (p == 0) continue;

				pk1	= p * min(d, max(0, j - k)) / d;
				pk2 = p - pk1;
				if (pk1 != 0)	
				{
					H1  -= (pk1 / p1) * log(pk1 / p1);
				}
				if (pk2 != 0)
				{
					H2  -= (pk2 / p2) * log(pk2 / p2);
				}
			}

			Err	= (H1 - H2) * (H1 - H2);

			if (Err < MinErr)
			{
				MinErr	= Err;
				t1 = i;
				t2 = j;
			}
		}
	}

	low = t1;
	high= t2;
	//////////////////////////////////////////////////////////////////////////
	if (!bNorm)
	{
		cvReleaseHist(&h);
	}
}

CvRect	GetROIFuzzyGrowing(IplImage* img)
{
	CvRect rect = cvRect(0, 0, 0, 0);
	if (1 != img->nChannels)
	{
		return rect;
	}

	int hist_size[] = {256};
	float r_ranges[] = { 0, 255 };
	float* ranges[] = { r_ranges };
	CvHistogram	*h	= cvCreateHist( 1, hist_size, CV_HIST_ARRAY, ranges, 1 );
	cvClearHist(h);

	IplImage* planes[] = { img};

	cvCalcHist(planes, h, 0, NULL);
	cvNormalizeHist(h, 1);
	int low = 0, high = 0;
	DoubleThreshold(h, 256, low, high, true);

	int x = 0, y = 0;

	IplImage *mask	= cvCloneImage(img);
	int thres_low = max(0, low - 1), thres_high = min(254, high - 1);
	int CurrValue = 0;
	for (y = 0; y < mask->height; y ++)
	{
		for (x = 0; x < mask->width; x ++)
		{
			CurrValue = CV_IMAGE_ELEM(mask, unsigned char, y, x);
			if (CurrValue <= low)
			{
				CV_IMAGE_ELEM(mask, unsigned char, y, x) = 0;
			}
			else if (CurrValue < high)
			{
				CV_IMAGE_ELEM(mask, unsigned char, y, x) = 64;
			}
			else
			{
				CV_IMAGE_ELEM(mask, unsigned char, y, x) = 128;
			}
		}
	}

	// Floodfill

	//CvPoint ptSeed = cvPoint(-1, -1);
	//CvConnectedComp		comp;
	//double area = 0;
	for (y = 0; y < mask->height; y ++)
	{
		for (x = 0; x < mask->width; x ++)
		{
			CurrValue	= CV_IMAGE_ELEM(mask, unsigned char, y, x);
			if (CurrValue == 128 )
			{
				cvFloodFill(mask, cvPoint(x, y), cvScalarAll(255), cvScalarAll(70), cvScalarAll(255), 
					0, CV_FLOODFILL_FIXED_RANGE);

				//ptSeed = cvPoint(x, y);
				//cvFloodFill(mask, ptSeed, cvScalarAll(255), cvScalarAll(70), cvScalarAll(255), 
				//			&comp, CV_FLOODFILL_FIXED_RANGE);
				//if (comp.area > area)
				//{
				//	area = comp.area;
				//	rect = comp.rect;
				//}
			}
		}
	}

	for (y = 0; y < mask->height; y ++)
	{
		for (x = 0; x < mask->width; x ++)
		{
			if (CV_IMAGE_ELEM(mask, unsigned char, y, x) < 255)
			{
				CV_IMAGE_ELEM(mask, unsigned char, y, x) = 0;
			}
		}
	}
	rect = cvBoundingRect(mask);


	cvReleaseImage(&mask);
	cvReleaseHist(&h);
	return rect;
}

CvRect GetOptimalROI(IplImage* img, double alpha1, double alpha2)
{
	CvRect rect = cvRect(0, 0, 0, 0);
	if (1 != img->nChannels)
	{
		return rect;
	}

	int x = 0, y = 0;
	double CoverRate = 1;
	double density = 0;
	double total = 0, total_rect = 0;
	for (y = 0; y < img->height; y ++)
	{
		for (x = 0; x < img->width; x ++)
		{
			total += CV_IMAGE_ELEM(img, unsigned char, y, x);
		}
	}
	if (total == 0)
	{
		return rect;
	}
	double area = img->width * img->height;
	density		= (double)total / (255 * area);
	total_rect	= total;

	int left = 0, top = 0, right = img->width - 1, bottom = img->height - 1;
	double dis = (1 - CoverRate) * (1 - CoverRate) * alpha1 + (1 - density) * (1 - density) * alpha2;
	double dis_temp, total_temp, area_temp;

	bool flag = 1;
	// for debug //////////////////////////////////////////////////////////////////////////
	//	IplImage *ImgDebug	= cvCreateImage(cvGetSize(img), 8, 3);
	//	cvMerge(img, img, img, 0, ImgDebug);
	//	cvNamedWindow("Debug", 1);
	//	cvShowImage("Debug", ImgDebug);
	//	cvWaitKey(0);
	//////////////////////////////////////////////////////////////////////////

	while (flag && left < right && top < bottom)
	{
		flag	= 0;
		// left
		total_temp	= total_rect;
		area_temp	= area;
		for (y = top; y <= bottom; y ++)
		{
			total_temp	-= CV_IMAGE_ELEM(img, unsigned char, y, left);
			area_temp	--;
		}
		CoverRate	= total_temp / total;
		density		= total_temp / (area_temp * 255);
		dis_temp	= (1 - CoverRate) * (1 - CoverRate) * alpha1 + (1 - density) * (1 - density) * alpha2;
		if (dis_temp <= dis)
		{
			left ++;
			total_rect	= total_temp;
			area		= area_temp;
			dis			= dis_temp;
			flag		= 1;
			// for debug //////////////////////////////////////////////////////////////////////////
			//			cvLine(ImgDebug, cvPoint(left, top), cvPoint(left, bottom), CV_RGB(255, 0, 0));
			//			cvShowImage("Debug", ImgDebug);
			//			cvWaitKey(0);
		}

		// top
		total_temp	= total_rect;
		area_temp	= area;
		for	(x = left; x <= right; x ++)
		{
			total_temp	-= CV_IMAGE_ELEM(img, unsigned char, top, x);
			area_temp	--;
		}
		CoverRate	= total_temp / total;
		density		= total_temp / (area_temp * 255);
		dis_temp	= (1 - CoverRate) * (1 - CoverRate) * alpha1 + (1 - density) * (1 - density) * alpha2;
		if (dis_temp <= dis)
		{
			top ++;
			total_rect	= total_temp;
			area		= area_temp;
			dis			= dis_temp;
			flag		= 1;
			// for debug //////////////////////////////////////////////////////////////////////////
			//			cvLine(ImgDebug, cvPoint(left, top), cvPoint(right, top), CV_RGB(255, 0, 0));
			//			cvShowImage("Debug", ImgDebug);
			//			cvWaitKey(0);
		}

		// right
		total_temp	= total_rect;
		area_temp	= area;
		for (y = top; y <= bottom; y ++)
		{
			total_temp	-= CV_IMAGE_ELEM(img, unsigned char, y, right);
			area_temp	--;
		}
		CoverRate	= total_temp / total;
		density		= total_temp / (area_temp * 255);
		dis_temp	= (1 - CoverRate) * (1 - CoverRate) * alpha1 + (1 - density) * (1 - density) * alpha2;
		if (dis_temp <= dis)
		{
			right --;
			total_rect	= total_temp;
			area		= area_temp;
			dis			= dis_temp;
			flag		= 1;
			// for debug //////////////////////////////////////////////////////////////////////////
			//			cvLine(ImgDebug, cvPoint(right, top), cvPoint(right, bottom), CV_RGB(255, 0, 0));
			//			cvShowImage("Debug", ImgDebug);
			//			cvWaitKey(0);
		}

		//bottom
		total_temp	= total_rect;
		area_temp	= area;
		for	(x = left; x <= right; x ++)
		{
			total_temp	-= CV_IMAGE_ELEM(img, unsigned char, bottom, x);
			area_temp	--;
		}
		CoverRate	= total_temp / total;
		density		= total_temp / (area_temp * 255);
		dis_temp	= (1 - CoverRate) * (1 - CoverRate) * alpha1 + (1 - density) * (1 - density) * alpha2;
		if (dis_temp <= dis)
		{
			bottom --;
			total_rect	= total_temp;
			area		= area_temp;
			dis			= dis_temp;
			flag		= 1;
			// for debug //////////////////////////////////////////////////////////////////////////
			//			cvLine(ImgDebug, cvPoint(left, bottom), cvPoint(right, bottom), CV_RGB(255, 0, 0));
			//			cvShowImage("Debug", ImgDebug);
			//			cvWaitKey(0);
		}
	}
	rect.x	= left;
	rect.y	= top;
	rect.width	= right - left + 1;
	rect.height = bottom- top + 1;

	// for debug //////////////////////////////////////////////////////////////////////////
	// 	cvReleaseImage(&ImgDebug);

	return rect;
}

// Expand the rect by optimize cover rate and density
bool	OptimizeRect(IplImage* img, CvRect &rect, double alpha1, double alpha2)
{
	if (1 != img->nChannels)
	{
		return 0;
	}

	int left	= rect.x;
	int right	= left + rect.width - 1;
	int	top		= rect.y;
	int bottom	= top + rect.height -1;

	int x = 0, y = 0;
	double CoverRate = 1;
	double density = 0;
	double total = 0, total_rect = 0;
	for (y = 0; y < img->height; y ++)
	{
		for (x = 0; x < img->width; x ++)
		{
			total += CV_IMAGE_ELEM(img, unsigned char, y, x);
			if (x >= left && x <= right && y >= top && y <= bottom)
			{
				total_rect += CV_IMAGE_ELEM(img, unsigned char, y, x);
			}
		}
	}
	if (total == 0)
	{
		return 0;
	}
	double area = rect.width * rect.height;
	density		= (double)total_rect / (255 * area);
	CoverRate	= total_rect / total;

	double dis = (1 - CoverRate) * (1 - CoverRate) * alpha1 + (1 - density) * (1 - density) * alpha2;
	double dis_temp, total_temp, area_temp;

	bool flag = 1;
	// for debug //////////////////////////////////////////////////////////////////////////
	//	IplImage *ImgDebug	= cvCreateImage(cvGetSize(img), 8, 3);
	//	cvMerge(img, img, img, 0, ImgDebug);
	//	cvNamedWindow("Debug", 1);
	//	cvShowImage("Debug", ImgDebug);
	//	cvWaitKey(0);
	//////////////////////////////////////////////////////////////////////////

	while (flag)
	{
		flag	= 0;
		// left
		if (left > 0)
		{
			total_temp	= total_rect;
			area_temp	= area;
			for (y = top; y <= bottom; y ++)
			{
				total_temp	+= CV_IMAGE_ELEM(img, unsigned char, y, left - 1);
				area_temp	++;
			}
			CoverRate	= total_temp / total;
			density		= total_temp / (area_temp * 255);
			dis_temp	= (1 - CoverRate) * (1 - CoverRate) * alpha1 + (1 - density) * (1 - density) * alpha2;
			if (dis_temp <= dis)
			{
				left --;
				total_rect	= total_temp;
				area		= area_temp;
				dis			= dis_temp;
				flag		= 1;
				// for debug //////////////////////////////////////////////////////////////////////////
				//				cvLine(ImgDebug, cvPoint(left, top), cvPoint(left, bottom), CV_RGB(255, 0, 0));
				//				cvShowImage("Debug", ImgDebug);
				//				cvWaitKey(0);
			}
		}

		// top
		if (top > 0)
		{
			total_temp	= total_rect;
			area_temp	= area;
			for	(x = left; x <= right; x ++)
			{
				total_temp	+= CV_IMAGE_ELEM(img, unsigned char, top - 1, x);
				area_temp	++;
			}
			CoverRate	= total_temp / total;
			density		= total_temp / (area_temp * 255);
			dis_temp	= (1 - CoverRate) * (1 - CoverRate) * alpha1 + (1 - density) * (1 - density) * alpha2;
			if (dis_temp <= dis)
			{
				top --;
				total_rect	= total_temp;
				area		= area_temp;
				dis			= dis_temp;
				flag		= 1;
				// for debug //////////////////////////////////////////////////////////////////////////
				//				cvLine(ImgDebug, cvPoint(left, top), cvPoint(right, top), CV_RGB(255, 0, 0));
				//				cvShowImage("Debug", ImgDebug);
				//				cvWaitKey(0);
			}
		}

		// right
		if (right < img->width - 1)
		{
			total_temp	= total_rect;
			area_temp	= area;
			for (y = top; y <= bottom; y ++)
			{
				total_temp	+= CV_IMAGE_ELEM(img, unsigned char, y, right + 1);
				area_temp	++;
			}
			CoverRate	= total_temp / total;
			density		= total_temp / (area_temp * 255);
			dis_temp	= (1 - CoverRate) * (1 - CoverRate) * alpha1 + (1 - density) * (1 - density) * alpha2;
			if (dis_temp <= dis)
			{
				right ++;
				total_rect	= total_temp;
				area		= area_temp;
				dis			= dis_temp;
				flag		= 1;
				// for debug //////////////////////////////////////////////////////////////////////////
				//				cvLine(ImgDebug, cvPoint(right, top), cvPoint(right, bottom), CV_RGB(255, 0, 0));
				//				cvShowImage("Debug", ImgDebug);
				//				cvWaitKey(0);
			}
		}

		//bottom
		if (bottom < img->height - 1)
		{
			total_temp	= total_rect;
			area_temp	= area;
			for	(x = left; x <= right; x ++)
			{
				total_temp	+= CV_IMAGE_ELEM(img, unsigned char, bottom + 1, x);
				area_temp	++;
			}
			CoverRate	= total_temp / total;
			density		= total_temp / (area_temp * 255);
			dis_temp	= (1 - CoverRate) * (1 - CoverRate) * alpha1 + (1 - density) * (1 - density) * alpha2;
			if (dis_temp <= dis)
			{
				bottom ++;
				total_rect	= total_temp;
				area		= area_temp;
				dis			= dis_temp;
				flag		= 1;
				// for debug //////////////////////////////////////////////////////////////////////////
				//				cvLine(ImgDebug, cvPoint(left, bottom), cvPoint(right, bottom), CV_RGB(255, 0, 0));
				//				cvShowImage("Debug", ImgDebug);
				//				cvWaitKey(0);
			}
		}
	}
	rect.x	= left;
	rect.y	= top;
	rect.width	= right - left + 1;
	rect.height = bottom- top + 1;

	// for debug //////////////////////////////////////////////////////////////////////////
	// 	cvReleaseImage(&ImgDebug);

	return 1;
}

// Compare two rects
// The first one is the "Ground Truth"
// Precision = (rect1 & rect2) / rect2
// Recall = (rect1 & rect2) / rect1
bool	Cmp2Rects(CvRect rectGT, CvRect rect, double &Precision, double &Recall)
{
	int	L	= max(rectGT.x, rect.x);
	int	T	= max(rectGT.y, rect.y);
	int	R	= min(rectGT.x + rectGT.width - 1, rect.x + rect.width - 1);
	int	B	= min(rectGT.y + rectGT.height - 1, rect.y + rect.height - 1);
	if (L > R || T > B)
	{
		Precision	= 0;
		Recall		= 0;
	}
	else
	{
		double dArea		= (double)(rect.width) * (rect.height);
		double dAreaGT		= (double)(rectGT.width) * (rectGT.height);
		if (dArea == 0 || dAreaGT == 0)
		{
			return 0;
		}
		double	dAreaUnion	= (R - L + 1) * (B - T + 1);
		Precision	= dAreaUnion / dArea;
		Recall		= dAreaUnion / dAreaGT;
	}
	return 1;
}

// Calcualte the BDE(Boundary Displacement Error), CVPR 2007, MSRA, Tie Liu
// There are bugs not solved yet.
bool	BDERect(CvRect rectGT, CvRect rect, double &BDE_Mean, double &BDE_Sdv)
{
	if (0 == rect.height || 0 == rect.width || 0 == rectGT.width || 0 == rectGT.height)
	{
		return false;
	}
	int n = (rect.width + rect.height) * 2 - 4;
	if (n < 1)
	{
		return false;
	}
	CvMat* dis = cvCreateMat(n, 1, CV_64FC1); // To store the distance
	cvZero(dis);
	int k = 0, x = 0, y = 0;
	int G_x = 0, G_y = 0;	// the ground truth point
	int left = rect.x, right = left + rect.width - 1;
	int top  = rect.y, bottom = top + rect.height - 1;
	double minDistance, distance = 0;

	// The top boundary
	y = top;
	for(x = rect.x; x < rect.x + rect.width; x ++)
	{
		// initialize the distance as the its distance to the left-top corner
		minDistance	= pow(rectGT.x - x, 2.0) + pow(rectGT.y - y, 2.0);
		if (0 != minDistance)
		{
			G_y	= rectGT.y;
			for (G_x = rectGT.x; G_x < rectGT.x + rectGT.width; G_x ++)
			{
				distance	= pow(G_x - x, 2.0) + pow(G_y - y, 2.0);
				if (distance < minDistance)
				{
					minDistance	= distance;
				}
			}
		}
		CV_MAT_ELEM(*dis, double, k, 0)	= sqrt(minDistance);
		k ++;
	}
	// the bottom boundary
	y = bottom;
	for(x = rect.x; x < rect.x + rect.width; x ++)
	{
		// initialize the distance as the its distance to the left-top corner
		minDistance	= pow(rectGT.x - x, 2.0) + pow(rectGT.y - y, 2.0);
		if (0 != minDistance)
		{
			G_y	= rectGT.y + rectGT.height - 1;
			for (G_x = rectGT.x; G_x < rectGT.x + rectGT.width; G_x ++)
			{
				distance	= pow(G_x - x, 2.0) + pow(G_y - y, 2.0);
				if (distance < minDistance)
				{
					minDistance	= distance;
				}
			}
		}
		CV_MAT_ELEM(*dis, double, k, 0)	= sqrt(minDistance);
		k ++;
	}
	// the left boundary
	// the four corner points will not be calculated again
	x = left;
	for (y = rect.y + 1; y < rect.y + rect.height - 1; y ++)
	{
		// initialize the distance as the its distance to the left-top corner
		minDistance	= pow(rectGT.x - x, 2.0) + pow(rectGT.y - y, 2.0);
		if (0 != minDistance)
		{
			G_x	= rectGT.x;
			for (G_y = rectGT.y; G_y < rectGT.y + rectGT.height; G_y ++)
			{
				distance	= pow(G_x - x, 2.0) + pow(G_y - y, 2.0);
				if (distance < minDistance)
				{
					minDistance	= distance;
				}
			}
		}
		CV_MAT_ELEM(*dis, double, k, 0)	= sqrt(minDistance);
		k ++;
	}
	// the left boundary
	// the four corner points will not be calculated again
	x = right;
	for (y = rect.y + 1; y < rect.y + rect.height - 1; y ++)
	{
		// initialize the distance as the its distance to the left-top corner
		minDistance	= pow(rectGT.x - x, 2.0) + pow(rectGT.y - y, 2.0);
		if (0 != minDistance)
		{
			G_x	= rectGT.x + rectGT.width - 1;
			for (G_y = rectGT.y; G_y < rectGT.y + rectGT.height; G_y ++)
			{
				distance	= pow(G_x - x, 2.0) + pow(G_y - y, 2.0);
				if (distance < minDistance)
				{
					minDistance	= distance;
				}
			}
		}
		CV_MAT_ELEM(*dis, double, k, 0)	= sqrt(minDistance);
		k ++;
	}

	CvScalar m, v;
	cvAvgSdv(dis, &m, &v);
	BDE_Mean	= m.val[0];
	BDE_Sdv		= v.val[0];

	cvReleaseMat(&dis);
	return 1;
}


bool	CalcRegionContrastMul(CvMat* feature, IplImage* RegionLabel, IplImage* Contrast, 
							  int nNumRegion, double alpha, double lambda)
{
	cvZero(Contrast);

	int RegionCount = nNumRegion;
	int NumOfPixel	= RegionLabel->height * RegionLabel->width;

	CvMat* matRadius	= cvCreateMat(RegionCount, 1, CV_64FC1);
	CvMat* matDistance	= cvCreateMat(RegionCount, RegionCount, CV_64FC1);
	CvMat* matFeatureDis= cvCreateMat(RegionCount, RegionCount, CV_64FC1);
	CvMat* matCenter	= cvCreateMat(RegionCount, 2, CV_32SC1);
	CvMat* matArea		= cvCreateMat(RegionCount, 1, CV_64FC1);
	CvMat* matContrast	= cvCreateMat(RegionCount, 1, CV_64FC1);
	CvMat* matRegionFea	= cvCreateMat(RegionCount, feature->cols, CV_64FC1);
	cvZero(matRadius);
	cvZero(matDistance);
	cvZero(matFeatureDis);
	cvZero(matCenter);
	cvZero(matArea);
	cvZero(matRegionFea);
	cvZero(matContrast);

	int		x = 0, y = 0, i = 0, j = 0, k1 = 0, k2 = 0;
	double	dRadius_2	= 0;
	double	dDis_2		= 0;
	double	dFeaDis_2	= 0;
	CvPoint pt1, pt2;
	double dArea = 0;
	double *dpFea1 = 0, *dpFea2 = 0;

	//unsigned char *ucp;
	float *ucp;
	for(y = 0; y < RegionLabel->height; y ++)
	{
		for(x = 0; x < RegionLabel->width; x ++)
		{
			k1 = CV_IMAGE_ELEM(RegionLabel, int, y, x);
			if (k1 < 0)
			{
				cvReleaseMat(&matRadius);
				cvReleaseMat(&matDistance);
				cvReleaseMat(&matFeatureDis);
				cvReleaseMat(&matCenter);
				cvReleaseMat(&matArea);
				cvReleaseMat(&matContrast);
				cvReleaseMat(&matRegionFea);
				return 0;
			}

			CV_MAT_ELEM(*matArea, double, k1, 0)	= CV_MAT_ELEM(*matArea, double, k1, 0) + 1;
			CV_MAT_ELEM(*matCenter, int, k1, 0)		= CV_MAT_ELEM(*matCenter, int, k1, 0) + x;
			CV_MAT_ELEM(*matCenter, int, k1, 1)		= CV_MAT_ELEM(*matCenter, int, k1, 1) + y;

			ucp			= &CV_MAT_ELEM(*feature, float, y * RegionLabel->width + x, 0);
			dpFea1		= &CV_MAT_ELEM(*matRegionFea, double, k1, 0);
			for (int dim = 0; dim < feature->cols; dim ++)
			{
				dpFea1[dim]	+= ucp[dim];
			}
		}
	}

	// Calculate the area and center of each region
	// Calculate the spatial distance and feature distance of each region
	for(i = 0; i < RegionCount; i ++)
	{
		dArea	= CV_MAT_ELEM(*matArea, double, i, 0);
		dRadius_2	= dArea / pi;
		CV_MAT_ELEM(*matRadius, double, i, 0)	= dRadius_2;

		if (0 == dArea)
		{
			//			AfxMessageBox("dArea == 0");
			continue;
		}
		CV_MAT_ELEM(*matCenter, int, i, 0) = (int)(CV_MAT_ELEM(*matCenter, int, i, 0) / dArea);
		CV_MAT_ELEM(*matCenter, int, i, 1) = (int)(CV_MAT_ELEM(*matCenter, int, i, 1) / dArea);

		dpFea1	= &CV_MAT_ELEM(*matRegionFea, double, i, 0);
		for (int dim = 0; dim < feature->cols; dim ++)
		{
			dpFea1[dim] = dpFea1[dim] / dArea;
		}

		pt1 = cvPoint(CV_MAT_ELEM(*matCenter, int, i, 0),  CV_MAT_ELEM(*matCenter, int, i, 1));

		for(j = 0; j < i; j ++)
		{
			pt2			= cvPoint(CV_MAT_ELEM(*matCenter, int, j, 0),  CV_MAT_ELEM(*matCenter, int, j, 1));
			dpFea2		= &CV_MAT_ELEM(*matRegionFea, double, j, 0);

			dDis_2		= pow((double)(pt1.x - pt2.x), 2) + pow((double)(pt1.y - pt2.y), 2);
			CV_MAT_ELEM(*matDistance, double, i, j)		= CV_MAT_ELEM(*matDistance, double, j, i)	= dDis_2;

			dFeaDis_2	= 0;
			for (int dim = 0; dim < feature->cols; dim ++)
			{
				dFeaDis_2 += pow(dpFea1[dim] - dpFea2[dim], 2);
			}
			CV_MAT_ELEM(*matFeatureDis, double, i, j)	= CV_MAT_ELEM(*matFeatureDis, double, j, i)	=dFeaDis_2;
		}
	}

	const double lambda_2	= lambda * lambda;	// the ratio between center and surround
	//	const double alpha_2	= alpha * alpha;	// the ratio between region radius and sigma
	const double alpha_2	= (1 - 1 / lambda_2) / (2 * log(lambda_2));

	cvScale(matArea, matArea, 1 / (double)NumOfPixel, 0);

	// Calculat the contrast of each region
	double sigma_2, dContrast	= 0, minDis;
	for(i = 0; i < RegionCount; i ++)
	{
		dRadius_2	= CV_MAT_ELEM(*matRadius, double, i, 0);
		if (dRadius_2 == 0)
		{
			continue;
		}

		sigma_2 = dRadius_2 * alpha_2;
		dContrast	= 0;
		minDis = 4 * log(lambda_2) * sigma_2 / ( 1 - 1 / lambda_2);

		for(j = 0; j < RegionCount; j ++)
		{
			dFeaDis_2			= CV_MAT_ELEM(*matFeatureDis, double, i, j);
			dArea				= CV_MAT_ELEM(*matArea, double, j, 0);
			if (dArea == 0)
			{
				continue;
			}
			double dRadius_j_2	= CV_MAT_ELEM(*matRadius, double, j, 0);

			dDis_2		= CV_MAT_ELEM(*matDistance, double, i, j);
			dDis_2		= pow(sqrt(dDis_2) - sqrt(dRadius_j_2), 2);
			if (dDis_2 < minDis) 
			{
				dDis_2 = minDis;
			}
			dDis_2		= exp(-dDis_2 / (2 * sigma_2 * lambda_2)) / lambda_2 - exp(-dDis_2 / (2 * sigma_2));
			//			CV_MAT_ELEM(*matDistance, double, j, i)	= dDis_2;
			dContrast	+= dFeaDis_2 * dDis_2 * dArea;
		}
		dContrast	/= sigma_2;
		CV_MAT_ELEM(*matContrast, double, i, 0)		= dContrast;
	}

	double minValue = 0, maxValue = 0, scale, shift;
	cvMinMaxLoc(matContrast, &minValue, &maxValue);

	if (minValue == maxValue)
	{
		cvReleaseMat(&matRadius);
		cvReleaseMat(&matDistance);
		cvReleaseMat(&matFeatureDis);
		cvReleaseMat(&matCenter);
		cvReleaseMat(&matArea);
		cvReleaseMat(&matContrast);
		cvReleaseMat(&matRegionFea);
		return 0;
	}

	scale	= 1 / (maxValue - minValue);
	shift	= - minValue * scale;
	cvScale(matContrast, matContrast, scale, shift);
	double MeanContrast = 0, TotalArea = 0;
	for (int k = 0; k < RegionCount; k ++)
	{
		dArea		= CV_MAT_ELEM(*matArea, double, k, 0);
		MeanContrast+=CV_MAT_ELEM(*matContrast, double, k, 0) * dArea;
		TotalArea	+= dArea;
	}
	MeanContrast	/= TotalArea;
	scale	= 255 * 0.5 / MeanContrast;
	cvScale(matContrast, matContrast, scale, 0);

	for(y = 0; y < RegionLabel->height; y ++)
	{
		for(x = 0; x < RegionLabel->width; x ++)
		{
			k1	= CV_IMAGE_ELEM(RegionLabel, int, y, x);

			CV_IMAGE_ELEM(Contrast, unsigned char, y, x)	= min(255, (int)(CV_MAT_ELEM(*matContrast, double, k1, 0)));
		}
	}

	// Release memory
	cvReleaseMat(&matRadius);
	cvReleaseMat(&matDistance);
	cvReleaseMat(&matFeatureDis);
	cvReleaseMat(&matCenter);
	cvReleaseMat(&matArea);
	cvReleaseMat(&matContrast);
	cvReleaseMat(&matRegionFea);

	return 1;
}


//////////////////////////////////////////////////////////////////////////
// RGB information density
// RGB color contrast
// Map fusion: within-variance
//////////////////////////////////////////////////////////////////////////
bool	ImgAttention(const IplImage* image, IplImage* map)
{
	// Check input parameters
	if (!(3 == image->nChannels) || !(1 == map->nChannels))
	{
		return false;
	}
	if (!(image->width == map->width) || !(image->height == map->height))
	{
		return false;
	}

	// Begin
	int MaxClusterNum		= 8;
	int ClusterNum			= 0;
	int RegionNum			= 0;
	int x = 0, y = 0;
	CvMat* ClusterCenter	= cvCreateMat(MaxClusterNum, 3, CV_64FC1);
	CvMat* ElemCount		= cvCreateMat(MaxClusterNum, 1, CV_32SC1);
	IplImage* FeaLabel		= cvCreateImage(cvGetSize(image), IPL_DEPTH_32S, 1);
	IplImage* RegionLabel	= cvCreateImage(cvGetSize(image), IPL_DEPTH_32S, 1);
	IplImage* imgSeg		= cvCreateImage(cvGetSize(image), 8, 3);
	IplImage* Contrast_map			= cvCreateImage(cvGetSize(image), 8, 1);
	IplImage* Information_map		= cvCreateImage(cvGetSize(image), 8, 1);

	ClusterNum				= InitialClusterCenter(image, ClusterCenter, 0.95);
	ClusterNum				= GLA( image, imgSeg, ClusterNum, MaxClusterNum, FeaLabel, ClusterCenter, ElemCount, 
		cvTermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 100, sqrt(3.0) ), 1, 32, 0);
	RegionNum				= GetRegion(imgSeg, RegionLabel, 100);

	CalcRegionContrast(image, RegionLabel, Contrast_map, RegionNum);

	CalcInformation(FeaLabel, ElemCount, ClusterNum, Information_map);

	// fuse information map and contrast map
	double dtemp1, dtemp2;
	double w_information, w_contrast;
	int	thres = 0;

	int hist_size[] = {256};
	float r_ranges[] = { 0, 255 };
	float* ranges[] = { r_ranges };
	CvHistogram	*h	= cvCreateHist( 1, hist_size, CV_HIST_ARRAY, ranges, 1 );

	// calculate the within-cluster variance of information map
	cvClearHist(h);
	IplImage* planes[] = { Information_map};
	cvCalcHist(planes, h, 0, NULL);
	thres	= MinVar(h, 256, w_information);

	// calculate the within-cluster variance of contrast map
	cvClearHist(h);
	planes[0] = Contrast_map;
	cvCalcHist(planes, h, 0, NULL);
	thres	= MinVar(h, 256, w_contrast);

	w_information	= w_contrast / (w_contrast + w_information);
	w_contrast		= 1 - w_information;

	for (y = 0; y < image->height; y ++)
	{
		for (x = 0; x < image->width; x ++)
		{
			dtemp1 = (double)(CV_IMAGE_ELEM(Contrast_map, unsigned char, y, x));
			dtemp2 = (double)(CV_IMAGE_ELEM(Information_map, unsigned char, y, x));
			CV_IMAGE_ELEM(map, unsigned char, y, x) = (int)(dtemp1 * w_contrast + dtemp2 * w_information);
		}
	}

	// Release memory
	cvReleaseImage(&FeaLabel);
	cvReleaseImage(&RegionLabel);
	cvReleaseImage(&imgSeg);
	cvReleaseImage(&Contrast_map);
	cvReleaseImage(&Information_map);
	cvReleaseHist(&h);
	cvReleaseMat(&ClusterCenter);
	cvReleaseMat(&ElemCount);

	return true;
}


bool	Binarize(IplImage	*img)
{
	if (1 != img->nChannels)
	{
		return false;
	}

	int hist_size[] = {256};
	float r_ranges[] = { 0, 255 };
	float* ranges[] = { r_ranges };
	CvHistogram	*h	= cvCreateHist( 1, hist_size, CV_HIST_ARRAY, ranges, 1 );
	cvClearHist(h);

	IplImage* planes[] = { img};

	cvCalcHist(planes, h, 0, NULL);

	double Sb = 0;
	int	thres = MaxSb(h, 256, Sb);

	if (-1 == thres)
	{
		cvReleaseHist(&h);
		return	false;
	}
	cvThreshold(img, img, (double)thres, 255, CV_THRESH_BINARY);

	cvReleaseHist(&h);

	return	true;
}