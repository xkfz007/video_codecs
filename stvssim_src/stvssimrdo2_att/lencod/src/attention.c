//#include <afxwin.h>
#include <stdio.h>
//#include <stack>

#include "cxcore.h"
#include "highgui.h"
#include "cv.h"
#include "attention.h"
//#include "cvgabor.h"

// Inital Cluster Center using histogram. Preprocess for color quantization
// img:			the input image
// matCenter:	the result
// maxCluster:	the maximum of number of cluster
// CoverRate:	
int InitialClusterCenter(const IplImage* img, CvMat* matCenter, double CoverRate)
{
int R_bin = 8;
int G_bin = 8;
int B_bin = 8;

	IplImage	*plane_r		= cvCreateImage(cvGetSize(img), 8, 1);
	IplImage	*plane_g		= cvCreateImage(cvGetSize(img), 8, 1);
	IplImage	*plane_b		= cvCreateImage(cvGetSize(img), 8, 1);
	IplImage* planes[] = { plane_r, plane_g, plane_b };

	int NumCluster, maxCluster;
	double Rate;
	int maxIndex[3], r_Index, g_Index, b_Index;
	float	minValue , maxValue;
	int	k , nDim ;

	// Calcualte image histogram
	int hist_size[] = {R_bin, G_bin, B_bin};
	float r_ranges[] = { 0, 255 };
	float g_ranges[] = { 0, 255 };
	float b_ranges[] = { 0, 255 };
	float* ranges[] = { r_ranges, g_ranges, b_ranges };
	CvHistogram	*h	= cvCreateHist( 3, hist_size, CV_HIST_ARRAY, ranges, 1 );
	cvClearHist(h);
	cvZero(matCenter);

	cvCvtPixToPlane( img, plane_b, plane_g, plane_r, 0 );
	cvCalcHist(planes, h, 0, NULL);
	cvNormalizeHist(h, 1);

	// Choose the bins of higher values as cluster center
	NumCluster = 0;
	maxCluster	= matCenter->rows;
	Rate	= 0;
	minValue = 0;
	maxValue = 0;

	k = 0;
	nDim = 3;
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

	CvRNG rng = (CvRNG)-1;
	int i, j, x, y, k;
	int iter;
	double max_dist;

	int dims = 3;
	int sample_count;

	int flag;
	int k1, k2;
	// parameter check

	if( maxCluster < 1 )
	{
		//		AfxMessageBox("Number of clusters should be positive");
		return 0;
	}

	termcrit.epsilon *= termcrit.epsilon;
	sample_count = img->height * img->width;

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
				double *c;
				k	= CV_IMAGE_ELEM(labels, int, y, x);
				c = (double*)(centers->data.ptr + k*centers->step);
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
				unsigned char*s;
				x = cvRandInt( &rng ) % img->width;
				y = cvRandInt( &rng ) % img->height;
				 s = &CV_IMAGE_ELEM(img, unsigned char, y, x * 3);
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
			unsigned char*s;

			x = cvRandInt(&rng) % img->width;
			y = cvRandInt(&rng) % img->height;
			s = &CV_IMAGE_ELEM(img, unsigned char, y, x * 3);
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
		cvCopy(centers, old_centers,0);
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
				double* c;

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
				c = (double*)(centers->data.ptr + k_best * centers->step);
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
			double dist;
			double* c_o;

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
			{unsigned char*s;
				x = cvRandInt( &rng ) % img->width;
				y = cvRandInt( &rng ) % img->height;
				s = &CV_IMAGE_ELEM(img, unsigned char, y, x * 3);
				for( j = 0; j < dims; j++ )
				{
					c[j] = (double)(s[j]);
				}
			}


			dist = 0;
			c_o = (double*)(old_centers->data.ptr + k*old_centers->step);
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
	flag = 1;
	while (flag && cluster_count > 1)
	{
		flag = 0;
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
					flag = 1;
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
						flag =1;
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

		cvCopy(centers, old_centers,0);
		cvZero( centers );
		cvZero( counters );
		// assign labels
		for (y = 0; y < img->height; y ++)
		{
			for (x = 0; x < img->width; x ++)
			{
				unsigned char* s = &CV_IMAGE_ELEM(img, unsigned char, y, x * 3);
				double *c;
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
				c = (double*)(centers->data.ptr + k_best * centers->step);
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
				double *dp;
				unsigned char* uc;
				k = CV_IMAGE_ELEM(labels, int, y, x);
				dp = &CV_MAT_ELEM(*centers_arr, double, k, 0);
				uc = &CV_IMAGE_ELEM(imgResult, unsigned char, y, x * 3);
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
	IplImage* mask;
	CvConnectedComp		comp;
	CvPoint ptSeed ;
	int n;

	double dis,minDis;
	unsigned char*up1,*up2;
	int flag;

	for (y = 0; y < RegionLabel->height; y ++)
	{
		for (x = 0; x < RegionLabel->width; x ++)
		{
			CV_IMAGE_ELEM(RegionLabel, int, y, x) = -1;
		}
	}

	// Set the mask
	mask	= cvCreateImage(cvSize(imgSeg->width + 2, imgSeg->height + 2), 8 ,1);
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

	// Floodfill
	ptSeed = cvPoint(-1, -1);

	// Mark the regions

	n = 0;
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
	dis = 0;
	minDis = -1;
	flag = 1;
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

int CalcRegionContrast(const IplImage* img, IplImage* RegionLabel, IplImage* Contrast, 
						int nNumRegion)
{
	double alpha = 0.33;
	double lambda = 9.2410;
    int		x = 0, y = 0, i = 0, j = 0, k1 = 0, k2 = 0;
	double	dRadius_2	= 0;
	double	dDis_2		= 0;
	double	dFeaDis_2	= 0;
	CvPoint pt1, pt2;
	double dArea = 0;
	double *dpFea1 = 0, *dpFea2 = 0;

	unsigned char *ucp;

	double lambda_2,alpha_2;
	double sigma_2,dContrast,minDis;
	double minValue, maxValue, scale, shift;
	double MeanContrast, TotalArea;
	int k;

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

	cvZero(Contrast);
	
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

	lambda_2	= lambda * lambda;	// the ratio between center and surround
	//	const double alpha_2	= alpha * alpha;	// the ratio between region radius and sigma
	alpha_2	= (1 - 1 / lambda_2) / (2 * log(lambda_2));

	cvScale(matArea, matArea, 1 / (double)NumOfPixel, 0);

	// Calculate the contrast of each region
	dContrast	= 0;
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
			double dRadius_j_2;
			dFeaDis_2			= CV_MAT_ELEM(*matFeatureDis, double, i, j);
			dArea				= CV_MAT_ELEM(*matArea, double, j, 0);
			if (dArea == 0)
			{
				continue;
			}
			dRadius_j_2	= CV_MAT_ELEM(*matRadius, double, j, 0);

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

	minValue = 0;
	maxValue = 0;
	cvMinMaxLoc(matContrast, &minValue, &maxValue,0,0,0);

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

	MeanContrast = 0;
	TotalArea = 0;
	for (k = 0; k < RegionCount; k ++)
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

int CalcInformation(IplImage* FeatureLabel, CvMat* ElemCount, int cluster_count,
					 IplImage* InformationMap)
{
	CvMat*	matInformation = cvCreateMat(cluster_count, 1, CV_64FC1);
	int x = 0, y = 0, k = 0;
	double dTotalNum = (double)(FeatureLabel->height * FeatureLabel->width);
	double p = 0;
	double info = 0;
	double minValue, maxValue, scale, shift;

	cvZero(matInformation);

	for(k = 0; k < cluster_count; k ++)
	{
		if (CV_MAT_ELEM(*ElemCount, int, k, 0) < 0)
			return 0;

		p		= CV_MAT_ELEM(*ElemCount, int, k, 0) / dTotalNum;
		info = 0;
		if (!(p == 0))
		{
			info	= - log(p);
		}
		CV_MAT_ELEM(*matInformation, double, k, 0) = info;
	}
	minValue = 0;
	maxValue = 0;

	cvMinMaxLoc(matInformation, &minValue, &maxValue,0,0,0);
	if (minValue == maxValue)
	{
		cvZero(matInformation);
		cvReleaseMat(&matInformation);
		return 0;
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
	return 1;
}

//////////////////////////////////////////////////////////////////////////

int		MinVar(CvHistogram *hist, int bin, double *Sw)
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
		*Sw = -1;
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

	*Sw = MinVariance / num;
	return threshold;
}

// <= threshold
int		MaxSb(CvHistogram* hist, int bin, double *Sb)
{
	int k = 0, threshold = 0;
	double total_num	= 0, total_num_1	= 0, total_num_2 = 0;
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

	*Sb	= max_sdv / total_num;
	if (0 == *Sb)
	{
		return -1;
	}
	return threshold;
}

// 
int	CalcRegionContrastMul(CvMat* feature, IplImage* RegionLabel, IplImage* Contrast, 
							  int nNumRegion)
{
 double alpha = 0.33;
 double lambda = 9.2410;
	int		x = 0, y = 0, i = 0, j = 0, k1 = 0, k2 = 0;
	double	dRadius_2	= 0;
	double	dDis_2		= 0;
	double	dFeaDis_2	= 0;
	CvPoint pt1, pt2;
	double dArea = 0;
	double *dpFea1 = 0, *dpFea2 = 0;
	float *ucp;
	int dim;

	double lambda_2,alpha_2;
	double sigma_2,dContrast,minDis;
	double minValue, maxValue, scale, shift;
	double MeanContrast, TotalArea;
	int k;

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

	cvZero(Contrast);

	//unsigned char *ucp;
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
			for (dim = 0; dim < feature->cols; dim ++)
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
		for (dim = 0; dim < feature->cols; dim ++)
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
			for (dim = 0; dim < feature->cols; dim ++)
			{
				dFeaDis_2 += pow(dpFea1[dim] - dpFea2[dim], 2);
			}
			CV_MAT_ELEM(*matFeatureDis, double, i, j)	= CV_MAT_ELEM(*matFeatureDis, double, j, i)	=dFeaDis_2;
		}
	}

	lambda_2	= lambda * lambda;	// the ratio between center and surround
	//	const double alpha_2	= alpha * alpha;	// the ratio between region radius and sigma
	alpha_2	= (1 - 1 / lambda_2) / (2 * log(lambda_2));

	cvScale(matArea, matArea, 1 / (double)NumOfPixel, 0);

	// Calculat the contrast of each region
	dContrast	= 0;
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
			double dRadius_j_2;
			dFeaDis_2			= CV_MAT_ELEM(*matFeatureDis, double, i, j);
			dArea				= CV_MAT_ELEM(*matArea, double, j, 0);
			if (dArea == 0)
			{
				continue;
			}
			dRadius_j_2	= CV_MAT_ELEM(*matRadius, double, j, 0);

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

	minValue=0;
	maxValue=0;
	cvMinMaxLoc(matContrast, &minValue, &maxValue,0,0,0);

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
	MeanContrast = 0;
	TotalArea = 0;
	for (k = 0; k < RegionCount; k ++)
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
