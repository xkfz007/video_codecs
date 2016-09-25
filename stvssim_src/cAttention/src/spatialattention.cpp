#include <stdio.h>
#include <stdlib.h>

#include "saptialattention.h"
#include "attention.h"
#include "gabor.h"
#include "utility.h"
char*getfilename(char*path)
{
	char *name=(char*)malloc(200*sizeof(char));
	char*p,*q;
	p=path;
	while(*p!='\0')
	{
if(*p++=='/')
			q=p;
	}
	//q++;
	int i=0;

	while(*q!='\0'&&*q!='.')
	{
		name[i++]=*q++;
	}
	name[i]='\0';
	return name;
}
void convert_data(IplImage* attMap,unsigned char**attention_data, int xs,int ys){
	for(int y=0;y<ys;y++)
		for(int x=0;x<xs;x++)
			attention_data[y][x]=(byte)((attMap->imageData+ y * attMap->widthStep)[x]+128);
}
void write_data(char *file_name,unsigned char **attention_data,int xs,int ys){
	FILE *fp=fopen(file_name,"w");
	if(fp==NULL){
		fprintf(stderr,"can't open file");
		exit(-1);
	}
	for(int y=0;y<ys;y++)
	{
		for(int x=0;x<xs;x++)
			fprintf(fp,"%4d",attention_data[y][x]);
		fprintf(fp,"\n");
	}
	fclose(fp);
}
void seqSpatialAttention(char*yuvPath,int xs,int ys,char*savepath)
{
	FILE *p_in;
	if(NULL==(p_in=fopen(yuvPath,"rb")))
	{
		printf("Can't open file:%s",yuvPath);
		exit(-1);
	}

	IplImage* imgYUV = cvCreateImage(cvSize(xs,ys),8,3);
	IplImage* imgRGB = cvCreateImage(cvSize(xs,ys),8,3);

	Sourceframe *srcframe=AllocSourceframe(xs,ys);
	byte**imgY,**imgU,**imgV;
	get_mem2D(&imgY,ys,xs);
	get_mem2D(&imgU,ys/2,xs/2);//need ensure that "ys/2" is an integer
	get_mem2D(&imgV,ys/2,xs/2);

//	char*winName="originalImg";
// 	char*winName2="Saliency_Static";
// 	char*winName3="Saliency_Motion";
// 	char*winName4="Novelty_RGB";
//	char*winName5="AttentionMap";
//	cvNamedWindow(winName,1);
// 	cvNamedWindow(winName2,0);
// 	cvNamedWindow(winName3,0);
// 	cvNamedWindow(winName4,0);
//	cvNamedWindow(winName5,1);

//	char* suffix=".avi";
	char*name=getfilename(yuvPath);

	char*staticPath=(char*)malloc(100*sizeof(char));
	char*motionPath=(char*)malloc(100*sizeof(char));
	char*noveltyPath=(char*)malloc(100*sizeof(char));
	char*attentionPath=(char*)malloc(100*sizeof(char));
	char*ssa_name=(char*)malloc(100*sizeof(char));
	char*msa_name=(char*)malloc(100*sizeof(char));
	char*na_name=(char*)malloc(100*sizeof(char));
	char*att_name=(char*)malloc(100*sizeof(char));
	strcpy(staticPath,savepath);
	strcpy(motionPath,savepath);
	strcpy(noveltyPath,savepath);
	strcpy(attentionPath,savepath);

	strcat(staticPath,name);
	strcat(motionPath,name);
	strcat(noveltyPath,name);
	strcat(attentionPath,name);

	//strcat(staticPath,"_static.avi");
	//strcat(motionPath,"_motion.avi");
	//strcat(noveltyPath,"_novelty.avi");
	//strcat(attentionPath,"_attention.avi");
	strcat(staticPath,"_SSA");
	strcat(motionPath,"_MSA");
	strcat(noveltyPath,"_NA");
	strcat(attentionPath,"_ATT");

	double fps=30;
//	CvVideoWriter *staticWriter,*motionWriter,*noveltyWriter,*attentionWriter;
//#ifdef WIN32
//#define FOURCC CV_FOURCC('I','V','5','0')
//#else
//#define FOURCC CV_FOURCC('X','2','6','4')
//#endif
//	staticWriter = cvCreateVideoWriter( staticPath, FOURCC, fps, cvSize(xs, ys), 0);
//	motionWriter = cvCreateVideoWriter( motionPath, FOURCC, fps, cvSize(xs, ys), 0);
//	noveltyWriter = cvCreateVideoWriter( noveltyPath, FOURCC, fps, cvSize(xs, ys), 0);
//	attentionWriter = cvCreateVideoWriter( attentionPath, FOURCC, fps, cvSize(xs, ys), 0);

	IplImage*Saliency_Static,*Saliency_Motion, *Novelty_RGB,*AttentionMap;
	Saliency_Static	= cvCreateImage(cvSize(xs, ys), 8, 1);
	Saliency_Motion	= cvCreateImage(cvSize(xs, ys), 8, 1);
	Novelty_RGB	= cvCreateImage(cvSize(xs, ys), 8, 1);
	AttentionMap = cvCreateImage(cvSize(xs, ys), 8, 1);
	unsigned char **SSA,**MSA,**NA,**ATT;
	get_mem2D(&SSA,ys,xs);
	get_mem2D(&MSA,ys,xs);
	get_mem2D(&NA,ys,xs);
	get_mem2D(&ATT,ys,xs);


	int MaxClusterNum=4;

	IplImage*FeaLabel		= cvCreateImage(cvSize(xs, ys), IPL_DEPTH_32S, 1);
	cvZero(FeaLabel);
	CvMat *ElemCount		= cvCreateMat(MaxClusterNum, 1, CV_32SC1);
	IplImage*RegionLabel		= cvCreateImage(cvSize(xs, ys), IPL_DEPTH_32S, 1);


	IplImage*imgHSV,*imgCurrent,*imgPrevious;
	imgHSV	= cvCreateImage(cvSize(xs, ys), 8, 3);
	imgCurrent		= cvCreateImage(cvSize(xs, ys), 8, 1);
	imgPrevious		= cvCreateImage(cvSize(xs, ys), 8, 1);
	cvZero(imgPrevious);
     
	CvMat*matColorDistance;
	matColorDistance	= cvCreateMat(ys, xs, CV_64FC1);
	cvZero(matColorDistance);
	CvMat *matMean;
	matMean		= cvCreateMat(ys, xs, CV_64FC3);
	cvZero(matMean);

	int hist_size[] = {256};
	float r_ranges[] = { 0, 255 };
	float* ranges[] = { r_ranges };
	CvHistogram	*h	= cvCreateHist( 1, hist_size, CV_HIST_ARRAY, ranges, 1 );

	int RegionNum=0;
	int i=0;
	while(ReadOneFrame(p_in,i,0,0,xs,ys,srcframe)>0)
	{
		CopyFrameToOldImg(srcframe,imgY,imgU,imgV);
		for (int y = 0; y < ys; y ++)
		{
			for (int x = 0; x < xs; x ++)
			{
				((unsigned char*)(imgYUV->imageData + y * imgYUV->widthStep))[3 * x + 0] = imgY[y][x];
				((unsigned char*)(imgYUV->imageData + y * imgYUV->widthStep))[3 * x + 1] = imgU[y/2][x/2];
				((unsigned char*)(imgYUV->imageData + y * imgYUV->widthStep))[3 * x + 2] = imgV[y/2][x/2];
			}
		}
		cvCvtColor(imgYUV,imgRGB,CV_YCrCb2RGB);

		RegionNum=ImgSegmentation(imgRGB,MaxClusterNum,FeaLabel,ElemCount,RegionLabel);

		cvCvtColor(imgRGB, imgCurrent, CV_BGR2GRAY);
		cvCvtColor(imgRGB, imgHSV, CV_BGR2HSV);

		StaticSaliency(Saliency_Static,imgHSV,imgCurrent,RegionNum,MaxClusterNum,FeaLabel,ElemCount,RegionLabel,h);
		MotionSaliency(Saliency_Motion,imgCurrent,imgPrevious,RegionNum,RegionLabel);
		IplImage* imgTemp	= imgCurrent;
		imgCurrent			= imgPrevious;
		imgPrevious			= imgTemp;


 		StaticNovelty(Novelty_RGB,imgRGB,matColorDistance,matMean,i);	
 		GetUltimateMap(AttentionMap,Saliency_Static,Saliency_Motion,Novelty_RGB,h);

		convert_data(Saliency_Static,SSA,xs,ys);
		convert_data(Saliency_Motion,MSA,xs,ys);
		convert_data(Novelty_RGB,NA,xs,ys);
		convert_data(AttentionMap,ATT,xs,ys);

		sprintf(ssa_name,"%s_%03d.txt",staticPath,i);
		sprintf(msa_name,"%s_%03d.txt",motionPath,i);
		sprintf(na_name,"%s_%03d.txt",noveltyPath,i);
		sprintf(att_name,"%s_%03d.txt",attentionPath,i);
		
		write_data(ssa_name,SSA,xs,ys);
		write_data(msa_name,MSA,xs,ys);
		write_data(na_name,NA,xs,ys);
		write_data(att_name,ATT,xs,ys);
	
// 		cvShowImage(winName,imgRGB);
// 		cvShowImage(winName2,Saliency_Static);
// 		cvShowImage(winName3,Saliency_Motion);
// 		cvShowImage(winName4,Novelty_RGB);
// 	    cvShowImage(winName5,AttentionMap);

//		cvWriteFrame(staticWriter, Saliency_Static);
//		cvWriteFrame(motionWriter, Saliency_Motion);
//		cvWriteFrame(noveltyWriter, Novelty_RGB);
//		cvWriteFrame(attentionWriter, AttentionMap);

		printf(".");
        fflush(stdout);
        if((i+1)%50==0)
        {
            printf("\n");
            fflush(stdout);
        }

		i++;
//		cvWaitKey(0);
	}
	printf("\n");
    fflush(stdout);

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

//	cvReleaseVideoWriter(&staticWriter);
//	cvReleaseVideoWriter(&motionWriter);
//	cvReleaseVideoWriter(&noveltyWriter);
//	cvReleaseVideoWriter(&attentionWriter);
	
	cvReleaseImage(&imgYUV);
	cvReleaseImage(&imgRGB);

//	cvDestroyWindow(winName);
// 	cvDestroyWindow(winName2);
// 	cvDestroyWindow(winName3);
// 	cvDestroyWindow(winName4);
// 	cvDestroyWindow(winName5);


	FreeSourceframe(srcframe);
	free_mem2D(imgY);
	free_mem2D(imgU);
	free_mem2D(imgV);
	free_mem2D(SSA);
	free_mem2D(MSA);
	free_mem2D(NA);
	free_mem2D(ATT);

	fclose(p_in);


}
// void SpatialAttention(IplImage*imgRGB,int MaxClusterNum,IplImage*FeaLabel,CvMat*ElemCount,IplImage*RegionLabel)
// {
// 	int RegionNum=0;
// 	RegionNum=ImgSegmentation(imgRGB,MaxClusterNum,FeaLabel,ElemCount,RegionLabel);
// 
// 	IplImage*imgCurrent;
// 	cvCvtColor(imgRGB, imgCurrent, CV_BGR2GRAY);
// 	cvCvtColor(imgRGB, imgHSV, CV_BGR2HSV);
// 
// 	StaticSaliency(Saliency_Static,imgHSV,imgCurrent,RegionNum,MaxClusterNum,FeaLabel,ElemCount,RegionLabel,h);
// 	MotionSaliency(Saliency_Motion,imgCurrent,imgPrevious,RegionNum,RegionLabel);
// 	IplImage* imgTemp	= imgCurrent;
// 	imgCurrent			= imgPrevious;
// 	imgPrevious			= imgTemp;
// 
// 
// 	StaticNovelty(Novelty_RGB,imgRGB,matColorDistance,matMean,i);	
// 	GetUltimateMap(AttentionMap,Saliency_Static,Saliency_Motion,Novelty_RGB,h);
// 
// }
int  ImgSegmentation(IplImage*imgFrame,int MaxClusterNum,IplImage*FeaLabel,CvMat*ElemCount,IplImage*RegionLabel)
{
	int Width=imgFrame->width;
	int Height=imgFrame->height;

	int ClusterNum=0;
	int RegionNum;

	IplImage*imgSeg			= cvCreateImage(cvSize(Width, Height), 8, 3);
	CvMat*ClusterCenter	= cvCreateMat(MaxClusterNum, 3, CV_64FC1);
	cvZero(ClusterCenter);
	// image segmentation using RGB image
	ClusterNum				= InitialClusterCenter(imgFrame, ClusterCenter, 0.95);
	ClusterNum				= GLA( imgFrame, imgSeg, ClusterNum, MaxClusterNum, FeaLabel, ClusterCenter, ElemCount, 
		cvTermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 100, sqrt(3.0) ), 1, 32, 0);
	RegionNum				= GetRegion(imgSeg, RegionLabel, 100);

	cvReleaseImage(&imgSeg);
	cvReleaseMat(&ClusterCenter);

	return RegionNum;
}
void StaticSaliency(IplImage* Saliency_Static,IplImage*imgHSV,IplImage*imgCurrent,int RegionNum,int MaxClusterNum,IplImage*FeaLabel,CvMat*ElemCount,IplImage*RegionLabel,CvHistogram*h)
{
	int Width=imgHSV->width;
	int Height=imgHSV->height;


	//	IplImage *imgHSV	= cvCreateImage(cvSize(Width, Height), 8, 3);
	//IplImage*	= cvCreateImage(cvSize(Width, Height), 8, 1);

	IplImage* Information_RGB	= cvCreateImage(cvSize(Width, Height), 8, 1);
	//	 IplImage* imgCurrent		= cvCreateImage(cvSize(Width, Height), 8, 1);
	IplImage* Contrast_HSV	= cvCreateImage(cvSize(Width, Height), 8, 1);

	// int RegionNum;
	// ImgSegmentation(imgFrame,&RegionNum,MaxClusterNum,FeaLabel,ElemCount,RegionLabel);
	// static saliency
	// calculate information density with RGB
	CalcInformation(FeaLabel, ElemCount, MaxClusterNum, Information_RGB);

	// calculate HSV contrast
	// 	cvCvtColor(imgFrame, imgCurrent, CV_BGR2GRAY);
	// 	cvCvtColor(imgFrame, imgHSV, CV_BGR2HSV);

	int x = 0, y = 0;
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

	// Gabor
	int scale = 4, direction = 4, FeatureDim = 16;
	int s = 0, d = 0;
	int pos	= 0;
	double Sigma = 2*PI;
	double F = sqrt(2.0);

	CvMat* Feature	= cvCreateMat(Height * Width, FeatureDim, CV_32F);
	IplImage* imgGabor	= cvCreateImage(cvSize(Width, Height), IPL_DEPTH_32F, 1);

	IplImage* Contrast_Gabor	= cvCreateImage(cvSize(Width, Height), 8, 1);

// 	CvGabor *gabor = new CvGabor;
// 	gabor->Init(d * PI / direction, s, Sigma, F);
    long wd=mask_width(F,Sigma);
 	CvMat *Real,*Imag;
 	Real = cvCreateMat( Width, Width, CV_32FC1);
 	Imag = cvCreateMat( Width, Width, CV_32FC1);
	Init(d*PI/direction,s,Sigma,F,wd,Real,Imag);
	for (s = 0; s < scale; s ++)
	{
		for (d = 0; d < direction; d ++)
		{
			conv_img(imgCurrent, imgGabor,Real,Imag,wd,CV_GABOR_MAG);
			int posTemp = 0;
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
	//delete gabor;
 	cvReleaseMat( &Real );
 	cvReleaseMat( &Imag );
    cvReleaseImage(&imgGabor);

	CvMat*	  FeaLabel_1D, matHeader_1;
	FeaLabel_1D				= cvReshape(FeaLabel, &matHeader_1, 0, Height * Width);

	CalcRegionContrastMul(Feature, RegionLabel, Contrast_Gabor, RegionNum);

	//////////////////////////////////////////////////////////////////////////
	double	Sw_RGB = 0, Sw_HSV	= 0, Sw_Gabor = 0;
	double  weight_RGB = 0, weight_HSV = 0, weight_Gabor = 0;
	int thres;

	// 	int hist_size[] = {256};
	// 	float r_ranges[] = { 0, 255 };
	// 	float* ranges[] = { r_ranges };
	// 	CvHistogram	*h	= cvCreateHist( 1, hist_size, CV_HIST_ARRAY, ranges, 1 );

	// calculate the between-cluster variance of RGB information
	cvClearHist(h);
	IplImage* planes[] = { Information_RGB};
	cvCalcHist(planes, h, 0, NULL);
	thres	= MinVar(h, 256, Sw_RGB);

	// calculate the between-cluster variance of HSV contrast
	cvClearHist(h);
	planes[0] = Contrast_HSV;
	cvCalcHist(planes, h, 0, NULL);
	thres	= MinVar(h, 256, Sw_HSV);

	// calculate the between-cluster variance of Gabor Contrast
	cvClearHist(h);
	planes[0] = Contrast_Gabor;
	cvCalcHist(planes, h, 0, NULL);
	thres	= MinVar(h, 256, Sw_Gabor);

	double sum_factror = Sw_RGB + Sw_HSV + Sw_Gabor;
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


	cvReleaseImage(&Information_RGB);
	cvReleaseImage(&Contrast_HSV);
	cvReleaseMat(&Feature);
	cvReleaseImage(&Contrast_Gabor);
	

//	return Saliency_Static;

}

void MotionSaliency(IplImage *Saliency_Motion,IplImage*imgCurrent,IplImage *imgPrevious,int RegionNum,IplImage*RegionLabel)
{

	int Width=imgCurrent->width;
	int Height=imgCurrent->height;


	CvMat *matVel_X, *matVel_Y, *mat_Motion_Feature;
	matVel_X		= cvCreateMat(Height, Width, CV_32FC1);
	matVel_Y		= cvCreateMat(Height, Width, CV_32FC1);
	mat_Motion_Feature	= cvCreateMat(Height * Width, 2, CV_32FC1);

	cvZero(matVel_X);
	cvZero(matVel_Y);
	cvZero(mat_Motion_Feature);

	// Motion saliency
	cvCalcOpticalFlowHS( imgPrevious, imgCurrent, 1,matVel_X, matVel_Y, 1,cvTermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 100, 0.001));

	int pos = 0;
	for (int y = 0; y < Height; y ++)
	{
		for (int x = 0; x < Width; x ++)
		{
			CV_MAT_ELEM(*mat_Motion_Feature, float, pos, 0)	= CV_MAT_ELEM(*matVel_X, float, y, x);
			CV_MAT_ELEM(*mat_Motion_Feature, float, pos, 1)	= CV_MAT_ELEM(*matVel_Y, float, y, x);
			pos ++;
		}
	}

	CalcRegionContrastMul(mat_Motion_Feature, RegionLabel, Saliency_Motion, RegionNum);

	cvReleaseMat(&matVel_X);
	cvReleaseMat(&matVel_Y);
	cvReleaseMat(&mat_Motion_Feature);

	//return Saliency_Motion;
}
void StaticNovelty(IplImage* Novelty_RGB,IplImage*imgFrame,CvMat*matColorDistance,CvMat*matMean,int n)
{

	int Width=imgFrame->width;
	int Height=imgFrame->height;

	//Novelty_RGB	= cvCreateImage(cvSize(Width, Height), 8, 1);

	
	//StaticNovelty(imgFrame, matColorDistance, matMean, FrameNo);


	// Begin to calculate novelty and update model
	int x = 0, y = 0, k = 0;
	double *dpM;//, *dpV;
	unsigned char* ucpC;
	double dis = 0, dColor = 0;
	double dMeanOld = 0, dVarianceOld = 0, dMeanNew = 0;

	// Calculate novelty
	// update parameters

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
	double dMin = 0, dMax = 0, dscale = 1, dshift = 0;
	dMin = 0;
	dMax = 0;
	dscale = 1;
	dshift = 0;
	cvMinMaxLoc(matColorDistance, &dMin, &dMax, 0, 0, 0);
	CvScalar dMean = cvScalarAll(0), dStd = cvScalarAll(0);
	cvAvgSdv(matColorDistance, &dMean, &dStd, 0);
	dMin	= max(dMin, (dMean.val[0] - 3 * dStd.val[0]) );
	dMax	= min(dMax, (dMean.val[0] + 3 * dStd.val[0]) );
	if (dMin != dMax)
	{
		dscale = 255 / (dMax - dMin);
		dshift = - dMin * dscale;
	}

	CvMat	matHeader, *matNovelty;
	matNovelty	= cvGetMat(Novelty_RGB, &matHeader);
	cvConvertScale(matColorDistance, matNovelty, dscale, dshift);
	if (1 == imgFrame->origin)
	{
		cvFlip(Novelty_RGB);
	}
	
	//cvReleaseMat(&matColorDistance);

	//return Novelty_RGB;
}
void GetUltimateMap(IplImage *AttentionMap,IplImage*Saliency_Static,IplImage*Saliency_Motion,IplImage*Novelty_RGB,CvHistogram*h)
{
	double	Sw_SS = 0, Sw_MS	= 0, Sw_SN = 0;
	double  weight_SS = 0, weight_MS = 0, weight_SN = 0;

	int thres;

	// calculate the between-cluster variance of RGB information
	cvClearHist(h);
	IplImage* planes[]={Saliency_Static};
	planes[0] = Saliency_Static;
	cvCalcHist(planes, h, 0, NULL);
	thres	= MinVar(h, 256, Sw_SS);

	// calculate the between-cluster variance of HSV contrast
	cvClearHist(h);
	planes[0] = Saliency_Motion;
	cvCalcHist(planes, h, 0, NULL);
	thres	= MinVar(h, 256, Sw_MS);

	// calculate the between-cluster variance of Gabor Contrast
	cvClearHist(h);
	planes[0] = Novelty_RGB;
	cvCalcHist(planes, h, 0, NULL);
	thres	= MinVar(h, 256, Sw_SN);

	double sum_factror;
	sum_factror = Sw_SS + Sw_MS + Sw_SN;
	if (!(0 == sum_factror))
	{
		weight_SS	= 0.5 * (1 - Sw_SS / sum_factror);
		weight_MS	= 0.5 * (1 - Sw_MS / sum_factror);
		weight_SN	= 0.5 * (1 - Sw_SN / sum_factror);
	}

	// = cvCreateImage(cvSize(Saliency_Static), 8, 1);
	cvZero(AttentionMap);
	cvAddWeighted(Saliency_Static, weight_SS, Saliency_Motion, weight_MS, 0, AttentionMap);
	cvAddWeighted(AttentionMap, 1, Novelty_RGB, weight_SN, 0, AttentionMap);

	//return AttentionMap;
}
