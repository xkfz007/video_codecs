#include "cxcore.h"
//#include <stack>

#define pi			3.14159265

//typedef stack<CvPoint>  CPointStack;


int		InitialClusterCenter(const IplImage* img, CvMat* matCenter, double CoverRate,
							 int R_bin = 8, int G_bin = 8, int B_bin = 8);

int		GLA(	const IplImage* img, IplImage* imgResult, 
			int ClusterCount, int maxCluster,
			IplImage* labels, CvMat* centers_arr, CvMat* ElemCount,
			CvTermCriteria termcrit, int UsePrevious, 
			double minDis, int minCluster );

int		GetRegion(IplImage* imgSeg, IplImage* RegionLabel, int Noise);

bool	CalcRegionContrast(const IplImage* img, IplImage* RegionLabel, IplImage* Contrast, 
						   int nNumRegion, double alpha = 0.33, double lambda = 9.2410);

bool	CalcRegionContrast_float(const IplImage* img, IplImage* RegionLabel, IplImage* Contrast, 
								 int nNumRegion, double alpha = 0.33, double lambda = 9.2410);

// calculate contrast,multi dimension feature 2009-09-29
bool	CalcRegionContrastMul(CvMat* feature, IplImage* RegionLabel, IplImage* Contrast, 
							  int nNumRegion, double alpha = 0.33, double lambda = 9.2410);

bool	CalcInformation(IplImage* FeatureLabel, CvMat* ElemCount, int cluster_count,
						IplImage* Information);

bool	ImgAttention(const IplImage* image, IplImage* map);
/*bool	ImgAttention_HSV(const IplImage* image, IplImage* map);*/
/*bool	ImgAttention_HSV_Gabor_Late(const IplImage* image, IplImage* map);*/
/*bool	ImgAttention_HSV_Gabor_Early(const IplImage* image, IplImage* map);*/

int		MinVar(int hist[256]);
int		MinVar(int* hist, int bin, double &Sw);
int		MinVar(CvHistogram *hist, int bin, double &Sw);
int		MaxSb(int* hist, int bin, double &Sb);
int		MaxSb(CvHistogram* hist, int bin, double &Sb);

// Get the ROI according to 1-order center moment
CvRect	GetROIMoment(IplImage* img, double alpha = 1);

// Expand the rect
bool	ExpandRect(IplImage* img, CvRect &rect, double thres_1 = 0.35, double thres_2 = 0.03, int step = 1);

// Get the ROI according to cover rate
CvRect	GetROICover(IplImage* img, double alpha = 0.90);

// Get the ROI by binarize the map and get the bounding rect of the points
CvRect	GetROI(IplImage	*img);

// Get the low and high threshold using prior entropy method (MSRA, Ma Yufei)
void	DoubleThreshold(CvHistogram* hist, int bin, int& low, int& high, bool bNorm);

//	Get the ROI using fuzzy growing
CvRect	GetROIFuzzyGrowing(IplImage* img);

// Get the ROI by balance cover ratio and density
CvRect GetOptimalROI(IplImage* img, double alpha1 = 0.5, double alpha2 = 0.5);

// Expand the rect by optimize cover rate and density
bool	OptimizeRect(IplImage* img, CvRect &rect, double alpha1 = 0.5, double alpha2 = 0.5);

// Compare two rects
bool	Cmp2Rects(CvRect rectGT, CvRect rect, double &Precision, double &Recall); 

// Calculate the BDE(Boundary Displacement Error), CVPR 2007, MSRA, Tie Liu
bool	BDERect(CvRect rectGT, CvRect rect, double &BDE_Mean, double &BDE_Sdv);

// Extract the texture information of an image 2009-09-29


bool	Binarize(IplImage	*img);