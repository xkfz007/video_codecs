#ifndef SAPTATALATTENTION_H
#define SAPTATALATTENTION_H

#include "cxcore.h"
#include "highgui.h"
#include "cv.h"
#include "att_stv.h"
#include "global.h"


void IntialAttentionVar();
void ReleaseAttentionVar();
void seqSpatialAttention(imgpel**imgY_att,imgpel**imgU_att,imgpel**imgV_att);
int  ImgSegmentation(IplImage*imgFrame,int MaxClusterNum,IplImage*FeaLabel,CvMat*ElemCount,IplImage*RegionLabel);
void StaticSaliency(IplImage*,IplImage*imgHSV,IplImage*imgCurrent,int RegionNum,int MaxClusterNum,IplImage*FeaLabel,CvMat*ElemCount,IplImage*RegionLabel,CvHistogram*h);
void MotionSaliency(IplImage *Saliency_Motion,IplImage*imgCurrent,IplImage *imgPrevious,int RegionNum,IplImage*RegionLabel);
void StaticNovelty(IplImage*,IplImage*imgFrame,CvMat*,CvMat*,int n);
void GetMap2(IplImage*Map,IplImage*ftA,IplImage*ftB,CvHistogram*h);
void GetMap3(IplImage *AttentionMap,IplImage*Saliency_Static,IplImage*Saliency_Motion,IplImage*Static_Novelty,CvHistogram*h);
#endif