#ifndef SAPTATALATTENTION_H
#define SAPTATALATTENTION_H
#include "cxcore.h"
#include "highgui.h"
#include "cv.h"

void seqSpatialAttention(char*yuvPath,int xs,int ys,char*savepath);
int  ImgSegmentation(IplImage*imgFrame,int MaxClusterNum,IplImage*FeaLabel,CvMat*ElemCount,IplImage*RegionLabel);
void StaticSaliency(IplImage*,IplImage*imgHSV,IplImage*imgCurrent,int RegionNum,int MaxClusterNum,IplImage*FeaLabel,CvMat*ElemCount,IplImage*RegionLabel,CvHistogram*h);
//void MotionSaliency(IplImage*,IplImage*imgCurrent,IplImage *imgPrevious,int RegionNum,IplImage*RegionLabel);
void MotionSaliency(IplImage *Saliency_Motion,IplImage*imgCurrent,IplImage *imgPrevious,int RegionNum,IplImage*RegionLabel);
void StaticNovelty(IplImage*,IplImage*imgFrame,CvMat*,CvMat*,int n);
void GetUltimateMap(IplImage*,IplImage*Saliency_Static,IplImage*Saliency_Motion,IplImage*Novelty_RGB,CvHistogram*h);
#endif