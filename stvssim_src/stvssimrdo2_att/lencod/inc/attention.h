#ifndef ATTENTION_H
#define ATTENTION_H
#include "cxcore.h"

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

#define pi			3.14159265


int		InitialClusterCenter(const IplImage* img, CvMat* matCenter, double CoverRate);

int		GLA(	const IplImage* img, IplImage* imgResult, 
			int ClusterCount, int maxCluster,
			IplImage* labels, CvMat* centers_arr, CvMat* ElemCount,
			CvTermCriteria termcrit, int UsePrevious, 
			double minDis, int minCluster );

int		GetRegion(IplImage* imgSeg, IplImage* RegionLabel, int Noise);

int	CalcRegionContrast(const IplImage* img, IplImage* RegionLabel, IplImage* Contrast, 
						   int nNumRegion);

// calculate contrast,multi dimension feature 2009-09-29
int	CalcRegionContrastMul(CvMat* feature, IplImage* RegionLabel, IplImage* Contrast, 
							  int nNumRegion);

int	CalcInformation(IplImage* FeatureLabel, CvMat* ElemCount, int cluster_count,
						IplImage* Information);


int		MinVarArray(int hist[256]);
int		MinVar(CvHistogram *hist, int bin, double *Sw);
int		MaxSb(CvHistogram* hist, int bin, double *Sb);

#endif