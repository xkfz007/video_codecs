#ifndef GABOR_H
#define GABOR_H
#include <cv.h>
#include <highgui.h>
#define PI 3.14159265
#define CV_GABOR_REAL 1
#define CV_GABOR_IMAG 2
#define CV_GABOR_MAG  3
#define CV_GABOR_PHASE 4
long mask_width(double Sigma,double K);
void Init(double Phi, int iNu, double Sigma, double F,long wd,CvMat*Real,CvMat*Imag);
void conv_img(IplImage *src, IplImage *dst,CvMat*Real,CvMat*Imag, long wd,int Type);
#endif