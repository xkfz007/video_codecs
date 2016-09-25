#ifndef _ATT_STV_
#define _ATT_STV_

#define _RDO_STVSSIM2_ 0
#define _RDO_STVSSIM_ 1
#define _RDO_SSIM3D_ 0
#define _RDO_SSIM_ 0
#define _SSIM_WEIGHTED_ 0
#define _MB_RDDATA_ 0


#define _ATT_STATIC_ 1
#define _ATT_MOTION_ 1
#define _ATT_NOVELTY_ 1
#define _ATT_ALL_  ((_ATT_STATIC_)&&(_ATT_MOTION_)&&(_ATT_NOVELTY_))
#define _ATTENTION_  ((_ATT_STATIC_)||(_ATT_MOTION_)||(_ATT_NOVELTY_))

#define _M1_   0  //在计算ssim的时候进行加权，（像素级）
#define _M2_   1   //利用attention 调整宏块的lambda
#define _ATT_M1_ ((_M1_)&&(_ATTENTION_))
#define _ATT_M2_ ((_M2_)&&(_ATTENTION_))

#define _L1_ 0
#define _L2_ 1
#define _ADJUST_L1_ ((_L1_)&&(_ATT_M2_))
#define _ADJUST_L2_ ((_L2_)&&(_ATT_M2_))

#define _ADJUST_ ((_ADJUST_L1_)||(_ADJUST_L2_))


int  get_mem2Dfloat(float ***array2D, int dim0, int dim1);
void free_mem2Dfloat(float      **array2D);

//#define  REFNUM 30
#define  REFNUM 26
#define  SSIM3D_WGT 0.6f
#define ALPHA 0.9f
//#define ALPHA_0_9

#ifdef ALPHA_0_1
//alpha=0.1
#define _A1_ 0.001987
#define _B2_ -0.02347
#define _B1_ 0.1064
#endif

#ifdef ALPHA_0_3
//alpha=0.3
#define _A1_ 0.002637
#define _B2_ -0.02411 
#define _B1_ 0.09921
#endif

#ifdef ALPHA_0_5
//alpha=0.5
#define _A1_ 0.002748
#define _B2_ -0.02492
#define _B1_ 0.09777
#endif

#ifdef ALPHA_0_7
//alpha=0.7
#define _A1_ 0.002786 
#define _B2_  -0.02591 
#define _B1_ 0.09734 
#endif

#ifdef ALPHA_0_9
//alpha=0.9
#define _A1_ 0.002755 
#define _B2_ -0.02677
#define _B1_ 0.09781
#endif














#endif

