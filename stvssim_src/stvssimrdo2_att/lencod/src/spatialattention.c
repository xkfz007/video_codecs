#include <stdio.h>
#include <stdlib.h>

#include "saptialattention.h"
#include "attention.h"
#include "gabor.h"
#include "utility.h"
#include "global.h"
#include "memalloc.h"

extern InputParameters *params;
extern ImageParameters *img;

byte **attention_data;
//float **att_wgt;
float **att_mbWgt;


static IplImage* imgYUV;
static IplImage* imgRGB;

static double fps;
static int MaxClusterNum;

static IplImage*Saliency_Static;
static IplImage*Saliency_Motion; 
static IplImage*Novelty_RGB;
static IplImage*AttentionMap;

static int ys;
static int xs;

static IplImage*FeaLabel;
static CvMat *ElemCount;
static IplImage*RegionLabel;

static IplImage*imgHSV;
static IplImage*imgCurrent;
static IplImage*imgPrevious;

static CvMat*matColorDistance;
static CvMat *matMean;
static CvHistogram	*h;

static int hist_size[] = {256};
static float r_ranges[] = { 0, 255 };
static float* ranges[] = { r_ranges };
static int FeatureDim ;
static IplImage*imgSeg;
static CvMat*ClusterCenter;
static IplImage* Information_RGB;
static IplImage* Contrast_HSV;

static CvMat* Feature;
static IplImage* imgGabor;
static IplImage* Contrast_Gabor;

static CvMat *Real;
static CvMat *Imag;

static CvMat *matVel_X;
static CvMat *matVel_Y;
static CvMat *mat_Motion_Feature;

void IntialAttentionVar()
{
	ys=params->source.height;
	xs=params->source.width;
	
	get_mem2D(&attention_data,ys,xs);
//	get_mem2Dfloat(&att_wgt,ys,xs);
	get_mem2Dfloat(&att_mbWgt,ys/16,xs/16);

	imgYUV = cvCreateImage(cvSize(xs,ys),8,3);
	imgRGB = cvCreateImage(cvSize(xs,ys),8,3);
	fps=30;
	MaxClusterNum=4;

	Saliency_Static	= cvCreateImage(cvSize(xs, ys), 8, 1);
	Saliency_Motion	= cvCreateImage(cvSize(xs, ys), 8, 1);
	Novelty_RGB	= cvCreateImage(cvSize(xs, ys), 8, 1);
	AttentionMap = cvCreateImage(cvSize(xs, ys), 8, 1);
	
	FeaLabel		= cvCreateImage(cvSize(xs, ys), IPL_DEPTH_32S, 1);
	cvZero(FeaLabel);

	ElemCount		= cvCreateMat(MaxClusterNum, 1, CV_32SC1);
	RegionLabel		= cvCreateImage(cvSize(xs, ys), IPL_DEPTH_32S, 1);

	imgHSV	= cvCreateImage(cvSize(xs, ys), 8, 3);
	imgCurrent		= cvCreateImage(cvSize(xs, ys), 8, 1);
	imgPrevious		= cvCreateImage(cvSize(xs, ys), 8, 1);
	cvZero(imgPrevious);

	matColorDistance	= cvCreateMat(ys, xs, CV_64FC1);
	cvZero(matColorDistance);
	matMean		= cvCreateMat(ys, xs, CV_64FC3);
	cvZero(matMean);

	h	= cvCreateHist( 1, hist_size, CV_HIST_ARRAY, ranges, 1 );

	imgSeg			= cvCreateImage(cvSize(xs, ys), 8, 3);
	ClusterCenter	= cvCreateMat(MaxClusterNum, 3, CV_64FC1);

	Information_RGB	= cvCreateImage(cvSize(xs, ys), 8, 1);
	Contrast_HSV	= cvCreateImage(cvSize(xs, ys), 8, 1);

    FeatureDim = 16;

	Feature	= cvCreateMat(ys*xs, FeatureDim, CV_32F);
	imgGabor	= cvCreateImage(cvSize(xs, ys), IPL_DEPTH_32F, 1);
	Contrast_Gabor	= cvCreateImage(cvSize(xs, ys), 8, 1);
 	Real = cvCreateMat( xs, xs, CV_32FC1);
 	Imag = cvCreateMat( xs, xs, CV_32FC1);
	matVel_X		= cvCreateMat(ys, xs, CV_32FC1);
	matVel_Y		= cvCreateMat(ys, xs, CV_32FC1);
	mat_Motion_Feature	= cvCreateMat(ys * xs, 2, CV_32FC1);
}
void ReleaseAttentionVar()
{

	free_mem2D(attention_data);
//	free_mem2Dfloat(att_wgt);
	free_mem2Dfloat(att_mbWgt);

	cvReleaseImage(&Saliency_Static);
	cvReleaseImage(&Saliency_Motion);
	cvReleaseImage(&Novelty_RGB);
	cvReleaseImage(&AttentionMap);

	cvReleaseImage(&FeaLabel);
	cvReleaseMat(&ElemCount);
	cvReleaseImage(&RegionLabel);
	cvReleaseImage(&imgHSV);
	cvReleaseImage(&imgCurrent);
	cvReleaseImage(&imgPrevious);
	cvReleaseMat(&matColorDistance);
	cvReleaseMat(&matMean);
	cvReleaseHist(&h);

	cvReleaseImage(&imgYUV);
	cvReleaseImage(&imgRGB);
	
	cvReleaseImage(&imgSeg);
	cvReleaseMat(&ClusterCenter);

	cvReleaseImage(&Information_RGB);
	cvReleaseImage(&Contrast_HSV);

	cvReleaseMat(&Feature);
	cvReleaseImage(&Contrast_Gabor);
    cvReleaseImage(&imgGabor);
	//delete gabor;
 	cvReleaseMat( &Real );
 	cvReleaseMat( &Imag );

	cvReleaseMat(&matVel_X);
	cvReleaseMat(&matVel_Y);
	cvReleaseMat(&mat_Motion_Feature);

}
void seqSpatialAttention(imgpel**imgY_att,imgpel**imgU_att,imgpel**imgV_att)//(char*yuvPath,int xs,int ys,char*savepath)
{
IplImage* imgTemp;
IplImage* attMap;
	int RegionNum=0;
	int y,x;
	int i,j;
	int flg1=0;
	int flg2=0;
	int flg3=0;
	long sum=0;
	float mbWgt_mean;

		for (y = 0; y < ys; y ++)
		{
			for (x = 0; x < xs; x ++)
			{
				((byte*)(imgYUV->imageData + y * imgYUV->widthStep))[3 * x + 0] = (byte)imgY_att[y][x];
				((byte*)(imgYUV->imageData + y * imgYUV->widthStep))[3 * x + 1] = (byte)imgU_att[y/2][x/2];
				((byte*)(imgYUV->imageData + y * imgYUV->widthStep))[3 * x + 2] = (byte)imgV_att[y/2][x/2];
			}
		}
		cvCvtColor(imgYUV,imgRGB,CV_YCrCb2RGB);

		RegionNum=ImgSegmentation(imgRGB,MaxClusterNum,FeaLabel,ElemCount,RegionLabel);

		cvCvtColor(imgRGB, imgCurrent, CV_BGR2GRAY);
		cvCvtColor(imgRGB, imgHSV, CV_BGR2HSV);

#if  _ATT_STATIC_
		StaticSaliency(Saliency_Static,imgHSV,imgCurrent,RegionNum,MaxClusterNum,FeaLabel,ElemCount,RegionLabel,h);
		flg1=1;
#endif
#if _ATT_MOTION_
		MotionSaliency(Saliency_Motion,imgCurrent,imgPrevious,RegionNum,RegionLabel);

		imgTemp	= imgCurrent;
		imgCurrent			= imgPrevious;
		imgPrevious			= imgTemp;
		flg2=1;
#endif
#if _ATT_NOVELTY_
 		StaticNovelty(Novelty_RGB,imgRGB,matColorDistance,matMean,img->number);	
		flg3=1;
#endif 
		if(flg1==1&&flg2==0&&flg3==0)
			attMap=Saliency_Static;
		else if(flg1==0&&flg2==1&&flg3==0)
			attMap=Saliency_Motion;
		else if(flg1==0&&flg2==0&&flg3==1)
			attMap=Novelty_RGB;
		else if(flg1==1&&flg2==1&&flg3==0)
		{
			GetMap2(AttentionMap,Saliency_Static,Saliency_Motion,h);
			attMap=AttentionMap;
		}
		else if(flg1==1&&flg2==0&&flg3==1)
		{
			GetMap2(AttentionMap,Saliency_Static,Novelty_RGB,h);
			attMap=AttentionMap;
		}
		else if(flg1==0&&flg2==1&&flg3==1)
		{
			GetMap2(AttentionMap,Saliency_Motion,Novelty_RGB,h);
			attMap=AttentionMap;
		}
		else if(flg1==1&&flg2==1&&flg3==1)
		{
 		GetMap3(AttentionMap,Saliency_Static,Saliency_Motion,Novelty_RGB,h);
			attMap=AttentionMap;
		}

		{
			//FILE *pf=fopen("att_orig.txt","a+");
	for (y = 0; y < ys; y ++)
	{
		for(x=0;x<xs;x++)
		{
			attention_data[y][x]=(byte)((attMap->imageData+ y * attMap->widthStep)[x]+128);
			//fprintf(pf,"%7d ",(attMap->imageData+ y * attMap->widthStep)[x]);
			sum+=attention_data[y][x];
		}
	//	fprintf(pf,"\n");
	}
//////	fclose(pf);

		}

		mbWgt_mean=sum*1.0/(xs*ys)*16*16;
//	for (y = 0; y < ys; y ++)
//		for(x=0;x<xs;x++)
//		{
//			att_wgt[y][x]=(float)attention_data[y][x];
//		}
		
	for(y=0;y<ys/16;++y)
	{
		for(x=0;x<xs/16;++x)
		{
			float mbWgt=0;
			for(j=0;j<16;++j)
				for(i=0;i<16;++i)
					mbWgt+=attention_data[y*16+j][x*16+i]*1.0;
			att_mbWgt[y][x]=mbWgt_mean/mbWgt;
		}

	}

#if _DEBUG
	{
		FILE *pstatic=fopen("attention.txt","w");
		int i,j;
		unsigned char** aa=(unsigned char**)calloc(ys,sizeof(unsigned char*));

		for(j=0;j<ys;++j)
		{
			aa[j]=(unsigned char*)calloc(xs,sizeof(unsigned char));
			for(i=0;i<xs;++i)
			{
				aa[j][i]=attention_data[j][i];
				fprintf(pstatic,"%d ",aa[j][i]);
			}
			fprintf(pstatic,"\n");
			free(aa[j]);
		}
		fclose(pstatic);
		free(aa);
	}
	{
//		FILE *patt=fopen("att.txt","w");
		FILE *pmb_att=fopen("mb_att.txt","w");
		int i,j;

//		for(j=0;j<ys;++j)
//		{
//			for(i=0;i<xs;++i)
//			{
//				fprintf(patt,"%8.7f ",att_wgt[j][i]);
//			}
//			fprintf(patt,"\n");
//		}
		for(j=0;j<ys/16;++j)
		{
			for(i=0;i<xs/16;++i)
			{
				fprintf(pmb_att,"%8.7f ",att_mbWgt[j][i]);
			}
			fprintf(pmb_att,"\n");
		}
//		fclose(patt);
		fclose(pmb_att);
	}
#endif

}
int  ImgSegmentation(IplImage*imgFrame,int MaxClusterNum,IplImage*FeaLabel,CvMat*ElemCount,IplImage*RegionLabel)
{
	int Width=imgFrame->width;
	int Height=imgFrame->height;

	int ClusterNum=0;
	int RegionNum;

	cvZero(ClusterCenter);
	// image segmentation using RGB image
	ClusterNum				= InitialClusterCenter(imgFrame, ClusterCenter, 0.95);
	ClusterNum				= GLA( imgFrame, imgSeg, ClusterNum, MaxClusterNum, FeaLabel, ClusterCenter, ElemCount, 
		cvTermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 100, sqrt(3.0) ), 1, 32, 0);
	RegionNum				= GetRegion(imgSeg, RegionLabel, 100);


	return RegionNum;
}
void StaticSaliency(IplImage* Saliency_Static,IplImage*imgHSV,IplImage*imgCurrent,int RegionNum,int MaxClusterNum,IplImage*FeaLabel,CvMat*ElemCount,IplImage*RegionLabel,CvHistogram*h)
{
	int Width=imgHSV->width;
	int Height=imgHSV->height;
	int x, y;

	// Gabor
	int scale = 4, direction = 4;
	int s = 0, d = 0;
	int pos	= 0;
	double Sigma = 2*PI;
	double F = sqrt(2.0); 

	long wd;

	CvMat*	  FeaLabel_1D, matHeader_1;

	//////////////////////////////////////////////////////////////////////////
	double	Sw_RGB = 0, Sw_HSV	= 0, Sw_Gabor = 0;
	double  weight_RGB = 0, weight_HSV = 0, weight_Gabor = 0;
	int thres;

IplImage* planes[1];
double sum_factror;

	CalcInformation(FeaLabel, ElemCount, MaxClusterNum, Information_RGB);

	for (y = 0; y < imgHSV->height; y ++)
	{
		for (x = 0; x < imgHSV->width; x ++)
		{
			double Hue			= 2 * PI * ((unsigned char*)(imgHSV->imageData + y * imgHSV->widthStep))[3 * x + 0] / 180;
			double Saturation	= ((unsigned char*)(imgHSV->imageData + y * imgHSV->widthStep))[3 * x + 1];
			double Value		= ((unsigned char*)(imgHSV->imageData + y * imgHSV->widthStep))[3 * x + 2];
			((unsigned char*)(imgHSV->imageData + y * imgHSV->widthStep))[3 * x + 0] = (int)(Saturation * 0.5 * sin(Hue) + 128);
			((unsigned char*)(imgHSV->imageData + y * imgHSV->widthStep))[3 * x + 1] = (int)(Saturation * 0.5 * cos(Hue) + 128);
			((unsigned char*)(imgHSV->imageData + y * imgHSV->widthStep))[3 * x + 2] = (int)(0.5 * Value);
		}
	}

	CalcRegionContrast(imgHSV, RegionLabel, Contrast_HSV, RegionNum);



    wd=mask_width(F,Sigma);
	Init(d*PI/direction,s,Sigma,F,wd,Real,Imag);
	for (s = 0; s < scale; s ++)
	{
		for (d = 0; d < direction; d ++)
		{
			int posTemp = 0;
			conv_img(imgCurrent, imgGabor,Real,Imag,wd,CV_GABOR_MAG);
			for (y = 0; y < imgGabor->height; y ++)
			{
				for (x = 0; x < imgGabor->width; x ++)
				{
					CV_MAT_ELEM(*Feature, float, posTemp, pos)	= CV_IMAGE_ELEM(imgGabor, float, y, x);
					posTemp ++;
				}
			}
			pos  ++;
		}
	}

	FeaLabel_1D				= cvReshape(FeaLabel, &matHeader_1, 0, Height * Width);

	CalcRegionContrastMul(Feature, RegionLabel, Contrast_Gabor, RegionNum);


	// calculate the between-cluster variance of RGB information
	cvClearHist(h);
	planes[0] = Information_RGB;
	cvCalcHist(planes, h, 0, NULL);
	thres	= MinVar(h, 256, &Sw_RGB);

	// calculate the between-cluster variance of HSV contrast
	cvClearHist(h);
	planes[0] = Contrast_HSV;
	cvCalcHist(planes, h, 0, NULL);
	thres	= MinVar(h, 256, &Sw_HSV);

	// calculate the between-cluster variance of Gabor Contrast
	cvClearHist(h);
	planes[0] = Contrast_Gabor;
	cvCalcHist(planes, h, 0, NULL);
	thres	= MinVar(h, 256, &Sw_Gabor);

	sum_factror = Sw_RGB + Sw_HSV + Sw_Gabor;
	if (!(0 == sum_factror))
	{
		weight_RGB	= 0.5 * (1 - Sw_RGB / sum_factror);
		weight_HSV	= 0.5 * (1 - Sw_HSV / sum_factror);
		weight_Gabor= 0.5 * (1 - Sw_Gabor/ sum_factror);
	}
	//IplImage* Saliency_Static	= cvCreateImage(cvSize(Width, Height), 8, 1);
	cvZero(Saliency_Static);
	cvAddWeighted(Contrast_HSV, weight_HSV, Information_RGB, weight_RGB, 0, Saliency_Static);
	cvAddWeighted(Saliency_Static, 1, Contrast_Gabor, weight_Gabor, 0, Saliency_Static);

}

void MotionSaliency(IplImage *Saliency_Motion,IplImage*imgCurrent,IplImage *imgPrevious,int RegionNum,IplImage*RegionLabel)
{

	int Width=imgCurrent->width;
	int Height=imgCurrent->height;
	int pos = 0;
	int x,y;



	cvZero(matVel_X);
	cvZero(matVel_Y);
	cvZero(mat_Motion_Feature);

	// Motion saliency
	cvCalcOpticalFlowHS( imgPrevious, imgCurrent, 1,matVel_X, matVel_Y, 1,cvTermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 100, 0.001));

	for (y = 0; y < Height; y ++)
	{
		for (x = 0; x < Width; x ++)
		{
			CV_MAT_ELEM(*mat_Motion_Feature, float, pos, 0)	= CV_MAT_ELEM(*matVel_X, float, y, x);
			CV_MAT_ELEM(*mat_Motion_Feature, float, pos, 1)	= CV_MAT_ELEM(*matVel_Y, float, y, x);
			pos ++;
		}
	}

	CalcRegionContrastMul(mat_Motion_Feature, RegionLabel, Saliency_Motion, RegionNum);


	//return Saliency_Motion;
}
void StaticNovelty(IplImage* Novelty_RGB,IplImage*imgFrame,CvMat*matColorDistance,CvMat*matMean,int n)
{

	int Width=imgFrame->width;
	int Height=imgFrame->height;

	int x = 0, y = 0, k = 0;
	double *dpM;//, *dpV;
	unsigned char* ucpC;
	double dis = 0, dColor = 0;
	double dMeanOld = 0, dVarianceOld = 0, dMeanNew = 0;

	double dMin, dMax, dscale, dshift;

	CvScalar dMean = cvScalarAll(0), dStd = cvScalarAll(0);
	CvMat	matHeader, *matNovelty;

	for (y = 0; y < imgFrame->height; y ++)
	{
		dpM = (double*)(matMean->data.ptr + y * matMean->step);
		for (x = 0; x < imgFrame->width; x ++)
		{
			dis = 0;
			ucpC= &CV_IMAGE_ELEM(imgFrame, unsigned char, y, x * 3);

			for (k = 0; k < 3; k ++)
			{
				dMeanOld	= dpM[x * 3 + k];
				dColor		= (double)ucpC[k];
				dMeanNew	= (n * dMeanOld + dColor) / (n + 1);

				dis		+= fabs(dColor - dMeanOld);
				dpM[x * 3 + k]		= dMeanNew;
			}
			dis /= 3;
			CV_MAT_ELEM(*matColorDistance, double, y, x) = dis;
		}
	}
	//int x = 0, y = 0;
	dMin = 0;
	dMax = 0;
	dscale = 1;
	dshift = 0;
	cvMinMaxLoc(matColorDistance, &dMin, &dMax, 0, 0, 0);
	cvAvgSdv(matColorDistance, &dMean, &dStd, 0);
	dMin	= max(dMin, (dMean.val[0] - 3 * dStd.val[0]) );
	dMax	= min(dMax, (dMean.val[0] + 3 * dStd.val[0]) );
	if (dMin != dMax)
	{
		dscale = 255 / (dMax - dMin);
		dshift = - dMin * dscale;
	}

	matNovelty	= cvGetMat(Novelty_RGB, &matHeader,0,0);
	cvConvertScale(matColorDistance, matNovelty, dscale, dshift);
	if (1 == imgFrame->origin)
	{
		cvFlip(Novelty_RGB,0,0);
	}
	

}
void GetMap3(IplImage *AttentionMap,IplImage*Saliency_Static,IplImage*Saliency_Motion,IplImage*Static_Novelty,CvHistogram*h)
{
	double	Sw_SS = 0, Sw_MS	= 0, Sw_SN = 0;
	double  weight_SS = 0, weight_MS = 0, weight_SN = 0;

	int thres;
	IplImage* planes[1];
	double sum_factror;
	// calculate the between-cluster variance of RGB information
	cvClearHist(h);
	//IplImage* planes[]={Saliency_Static};
	planes[0] = Saliency_Static;
	cvCalcHist(planes, h, 0, NULL);
	thres	= MinVar(h, 256, &Sw_SS);

	// calculate the between-cluster variance of HSV contrast
	cvClearHist(h);
	planes[0] = Saliency_Motion;
	cvCalcHist(planes, h, 0, NULL);
	thres	= MinVar(h, 256, &Sw_MS);

	// calculate the between-cluster variance of Gabor Contrast
	cvClearHist(h);
	planes[0] = Static_Novelty;
	cvCalcHist(planes, h, 0, NULL);
	thres	= MinVar(h, 256, &Sw_SN);

	sum_factror = Sw_SS+ Sw_MS+ Sw_SN;
	if (!(0 == sum_factror))
	{
		weight_SS	= 0.5 * (1 - Sw_SS / sum_factror);
		weight_MS	= 0.5 * (1 - Sw_MS / sum_factror);
		weight_SN	= 0.5 * (1 - Sw_SN / sum_factror);
	}

	// = cvCreateImage(cvSize(Saliency_Static), 8, 1);
	cvZero(AttentionMap);
	cvAddWeighted(Saliency_Static, weight_SS, Saliency_Motion, weight_MS, 0, AttentionMap);
	cvAddWeighted(AttentionMap, 1, Static_Novelty, weight_SN, 0, AttentionMap);

	//return AttentionMap;
}

void GetMap2(IplImage*Map,IplImage*ftA,IplImage*ftB,CvHistogram*h)
{

double	swA = 0, swB	= 0;
	double  wgtA = 0, wgtB = 0;

	int thres;
	IplImage* planes[1];
	double sum_factror;

		cvClearHist(h);
		//IplImage* planes[]={ftA};
		planes[0] = ftA;
		cvCalcHist(planes, h, 0, NULL);
		thres	= MinVar(h, 256, &swA);

		cvClearHist(h);
		planes[0] = ftB;
		cvCalcHist(planes, h, 0, NULL);
		thres	= MinVar(h, 256, &swB);

		sum_factror = swA + swB;
		if (!(0 == sum_factror))
		{
			wgtA	= (1 - swA / sum_factror);
			wgtB	= (1 - swB / sum_factror);
		}

		cvZero(Map);
		cvAddWeighted(ftA, wgtA, ftB, wgtB, 0, Map);
}