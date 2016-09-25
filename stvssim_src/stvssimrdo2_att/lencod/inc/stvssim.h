#ifndef _STVSSIM_H_
#define _STVSSIM_H_
#include "global.h"
#include "att_stv.h"

#define  MODENUM 4
#define  FLTNUM 4
#define  PI 3.1415926f
#define ORIENTS 32 
#define ORIENTS2 16 

imgpel **refPicsData[REFNUM][3];
imgpel **srcPicsData[REFNUM][3];
short ****mb_mv;//8x4x4x2:7种模式，以4x4块为单位
//short ****mb_mv2;//8x16x16x2:7种模式，以像素为单位
float **mb_directions;//4x4:7种模式，以4x4块为单位
//float **mb_directions2;//4x4:7种模式，以4x4块为单位
float **pic_directions2;//heightxwidth:7种模式，以像素为单位

FILE *pInfo;
FILE *pAvg;

#if _MB_RDDATA_ 
float **mb_rddata;
FILE *pmbrddata;
#define MB_DATA_NUM 16
#endif

static const float orient[FLTNUM]={0,PI/4,PI/2,PI*3/4};
static const float gauss8[8][8] = {
	{0.0003f,    0.0012f,    0.0029f,    0.0045f,    0.0045f,    0.0029f,    0.0012f,    0.0003f,},
	{0.0012f,    0.0045f,    0.0108f,    0.0169f,    0.0169f,    0.0108f,    0.0045f,    0.0012f,},
	{0.0029f,    0.0108f,    0.0264f,    0.0411f,    0.0411f,    0.0264f,    0.0108f,    0.0029f,},
	{0.0045f,    0.0169f,    0.0411f,    0.0641f,    0.0641f,    0.0411f,    0.0169f,    0.0045f,},
	{0.0045f,    0.0169f,    0.0411f,    0.0641f,    0.0641f,    0.0411f,    0.0169f,    0.0045f,},
	{0.0029f,    0.0108f,    0.0264f,    0.0411f,    0.0411f,    0.0264f,    0.0108f,    0.0029f,},
	{0.0012f,    0.0045f,    0.0108f,    0.0169f,    0.0169f,    0.0108f,    0.0045f,    0.0012f,},
	{0.0003f,    0.0012f,    0.0029f,    0.0045f,    0.0045f,    0.0029f,    0.0012f,    0.0003f,},
};
static const float gauss4[4][4] = {   
	 {0.0382f,    0.0595f,    0.0595f,    0.0382f,},
	 {0.0595f,    0.0928f,    0.0928f,    0.0595f,},
	 {0.0595f,    0.0928f,    0.0928f,    0.0595f,},
	 {0.0382f,    0.0595f,    0.0595f,    0.0382f,},
};
static float orients[ORIENTS]={0,PI/32,PI*2/32,PI*3/32,PI*4/32,PI*5/32,PI*6/32,PI*7/32,
PI*8/32,PI*9/32,PI*10/32,PI*11/32,PI*12/32,PI*13/32,PI*14/32,PI*15/32,
PI*16/32,PI*17/32,PI*18/32,PI*19/32,PI*20/32,PI*21/32,PI*22/32,PI*23/32,
PI*24/32,PI*25/32,PI*26/32,PI*27/32,PI*28/32,PI*29/32,PI*30/32,PI*31/32,
};
static float orients2[ORIENTS2]={0,PI/16,PI*2/16,PI*3/16,PI*4/16,PI*5/16,PI*6/16,PI*7/16,
PI*8/16,PI*9/16,PI*10/16,PI*11/16,PI*12/16,PI*13/16,PI*14/16,PI*15/16,};


float vFilter(int beta,int y,int x,float wa,float wb);
float hFilter(int beta,int y,int x,float wa,float wb);
float lFilter(int beta,int y,int x,float wa,float wb);
float rFilter(int beta,int y,int x,float wa,float wb);
float vFilter_4(int beta,int y,int x,float wa,float wb);
float hFilter_4(int beta,int y,int x,float wa,float wb);
float lFilter_4(int beta,int y,int x,float wa,float wb);
float rFilter_4(int beta,int y,int x,float wa,float wb);
void calOrit(float tmporient,short *orit);
void storeRefAndEncFrames();

float compute_SSIM (imgpel **imgRef, imgpel **imgEnc,int yRef,int yEnc, int xRef, int xEnc, int height, int width, int wint, int comp);
double distortionSSIM(Macroblock *currMB) ;
float compute_stVSSIM(imgpel**imgRef,imgpel**imgEnc,int yRef,int yEnc, int xRef, int xEnc, 
					  int height, int width,int wint,int gama,int comp,
					  float *ssimValue,float*ssim3dValue,float*stvssimValue);
double distortionstVSSIM(Macroblock *currMB);
float compute_SSIM3D(imgpel**imgRef,imgpel**imgEnc,int yRef,int yEnc, int xRef, int xEnc, int height, int width,int wint,int gama,int comp);
double distortionSSIM3D(Macroblock *currMB); 
float compute_stVSSIM2(imgpel**imgRef,imgpel**imgEnc,int yRef,int yEnc, int xRef, int xEnc, 
					  int height, int width,int wint,int gama,int comp,
					  float *ssimValue,float*ssim3dValue,float*stvssimValue);
double distortionstVSSIM2(Macroblock *currMB);
void InitialRDOVar();
void ReleaseRDOVar();
void getMV_macroblock(short****mb_mv,int mode,int block,int bestref);
short  getOrientation(short x,short y);
float chooseOrient(short*orit);
void getDirection_macroblock();
float compute_PSNR(imgpel**imgRef,imgpel**imgEnc,int yRef,int yEnc, int xRef, int xEnc, int height, int width);
void get_mb_rddata(Macroblock *currMB);
void writeout_slice_rddata();
void find_stvssim(DistMetric *mssim,DistMetric*mssim3d,DistMetric*mstvssim);
double adjust_lambda(double lambda,double eta);
#if _RDO_SSIM_ || _RDO_SSIM3D_ || _RDO_STVSSIM_ || _RDO_STVSSIM2_
double lambda_poly(int qp);
double lambda_expon(int qp); 
double lambda_gauss(int qp); 
double lambda_1(int qp); 
double lambda_2(int qp); 
#endif

#endif