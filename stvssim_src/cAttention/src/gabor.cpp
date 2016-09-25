#include "gabor.h"
#include "utility.h"
long mask_width(double Sigma,double K)
{
	long lwd;
	double dModSigma = Sigma/K;
	double dwd = cvRound(dModSigma*6 + 1);
	//test whether dwd is an odd.
	if (fmod(dwd, 2.0)==0.0) dwd++;
	lwd = (long)dwd;

	return lwd;
}
void Init(double Phi, int iNu, double Sigma, double F,long wd,CvMat*Real,CvMat*Imag)
{

// 	bInitialised = false;
// 	bKernel = false;
	//double Sigma = dSigma;
	//double F = dF;

	double Kmax = PI/2;

	// Absolute value of K
	double K = Kmax / pow(F, (double)iNu);
	//double Phi = dPhi;
/*	bInitialised = true;*/
	//wd = mask_width();
// 	double dModSigma = Sigma/K;
// 	double dwd = cvRound(dModSigma*6 + 1);
// 	//test whether dwd is an odd.
// 	if (fmod(dwd, 2.0)==0.0) 
// 		dwd++;
// 	int wd = (long)dwd;

// 	CvMat *Real = cvCreateMat( wd, wd, CV_32FC1);
// 	CvMat *Imag = cvCreateMat( wd, wd, CV_32FC1);
	
/*	CvMat *mReal, *mImag;*/
// 	Real = cvCreateMat( wd, wd, CV_32FC1);
// 	Imag = cvCreateMat( wd, wd, CV_32FC1);

	/**************************** Gabor Function ****************************/ 
	int x, y;
	double dReal;
	double dImag;
	double dTemp1, dTemp2, dTemp3;

	for (int i = 0; i < wd; i++)
	{
		for (int j = 0; j < wd; j++)
		{
			x = i-(wd-1)/2;
			y = j-(wd-1)/2;
			dTemp1 = (pow(K,2)/pow(Sigma,2))*exp(-(pow((double)x,2)+pow((double)y,2))*pow(K,2)/(2*pow(Sigma,2)));
			dTemp2 = cos(K*cos(Phi)*x + K*sin(Phi)*y) - exp(-(pow(Sigma,2)/2));
			dTemp3 = sin(K*cos(Phi)*x + K*sin(Phi)*y);
			dReal = dTemp1*dTemp2;
			dImag = dTemp1*dTemp3; 
			//gan_mat_set_el(pmReal, i, j, dReal);
			//cvmSet( (CvMat*)mReal, i, j, dReal );
			cvSetReal2D((CvMat*)Real, i, j, dReal );
			//gan_mat_set_el(pmImag, i, j, dImag);
			//cvmSet( (CvMat*)mImag, i, j, dImag );
			cvSetReal2D((CvMat*)Imag, i, j, dImag );

		} 
	}
	/**************************** Gabor Function ****************************/
	//bKernel = true;
// 	cvCopy(mReal, Real, NULL);
// 	cvCopy(mImag, Imag, NULL);
	//printf("A %d x %d Gabor kernel with %f PI in arc is created.\n", wd, wd, Phi/PI);
// 	cvReleaseMat( &mReal );
// 	cvReleaseMat( &mImag );

}



void conv_img(IplImage *src, IplImage *dst,CvMat*Real,CvMat*Imag, long wd,int Type)
{
// 	CvMat*Real = cvCreateMat( wd, wd, CV_32FC1);
// 	CvMat*Imag = cvCreateMat( wd, wd, CV_32FC1);
	double ve;//, re,im;

	CvMat *mat = cvCreateMat(src->width, src->height, CV_32FC1);
	for (int i = 0; i < src->width; i++)
	{
		for (int j = 0; j < src->height; j++)
		{
			ve = CV_IMAGE_ELEM(src, uchar, j, i);
			CV_MAT_ELEM(*mat, float, i, j) = (float)ve;
		}
	}

	CvMat *rmat = cvCreateMat(src->width, src->height, CV_32FC1);
	CvMat *imat = cvCreateMat(src->width, src->height, CV_32FC1);

	switch (Type)
	{
	case CV_GABOR_REAL:
		cvFilter2D( (CvMat*)mat, (CvMat*)mat, (CvMat*)Real, cvPoint( (wd-1)/2, (wd-1)/2));
		break;
	case CV_GABOR_IMAG:
		cvFilter2D( (CvMat*)mat, (CvMat*)mat, (CvMat*)Imag, cvPoint( (wd-1)/2, (wd-1)/2));
		break;
	case CV_GABOR_MAG:
		cvFilter2D( (CvMat*)mat, (CvMat*)rmat, (CvMat*)Real, cvPoint( (wd-1)/2, (wd-1)/2));
		cvFilter2D( (CvMat*)mat, (CvMat*)imat, (CvMat*)Imag, cvPoint( (wd-1)/2, (wd-1)/2));

		cvPow(rmat,rmat,2); 
		cvPow(imat,imat,2);
		cvAdd(imat,rmat,mat); 
		cvPow(mat,mat,0.5); 
		break;
	case CV_GABOR_PHASE:
		break;
	}

	if (dst->depth == IPL_DEPTH_8U)
	{
		cvNormalize((CvMat*)mat, (CvMat*)mat, 0, 255, CV_MINMAX);
		for (int i = 0; i < mat->rows; i++)
		{
			for (int j = 0; j < mat->cols; j++)
			{
				ve = CV_MAT_ELEM(*mat, float, i, j);
				CV_IMAGE_ELEM(dst, uchar, j, i) = (uchar)cvRound(ve);
			}
		}
	}

	if (dst->depth == IPL_DEPTH_32F)
	{
		for (int i = 0; i < mat->rows; i++)
		{
			for (int j = 0; j < mat->cols; j++)
			{
				ve = cvGetReal2D((CvMat*)mat, i, j);
				cvSetReal2D( (IplImage*)dst, j, i, ve );
			}
		}
	}
	cvReleaseMat(&imat);
	cvReleaseMat(&rmat);
	cvReleaseMat(&mat);
// 	cvReleaseMat(&Real);
// 	cvReleaseMat(&Imag);
}
