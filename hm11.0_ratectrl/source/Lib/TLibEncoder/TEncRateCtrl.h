/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.  
 *
 * Copyright (c) 2010-2013, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/** \file     TEncRateCtrl.h
    \brief    Rate control manager class
*/

#ifndef _HM_TENCRATECTRL_H_
#define _HM_TENCRATECTRL_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "../TLibCommon/CommonDef.h"
#include "../TLibCommon/TComDataCU.h"

//#include <vector>
//#include <algorithm>
#include <stdint.h>

//using namespace std;

//! \ingroup TLibEncoder
//! \{

#if RATE_CONTROL_LAMBDA_DOMAIN
#include "../TLibEncoder/TEncCfg.h"
#include <list>
#include <cassert>

const Int g_RCInvalidQPValue = -999;
const Int g_RCSmoothWindowSize = 40;
const Int g_RCMaxPicListSize = 32;
const Double g_RCWeightPicTargetBitInGOP    = 0.9;
const Double g_RCWeightPicRargetBitInBuffer = 1.0 - g_RCWeightPicTargetBitInGOP;
#if M0036_RC_IMPROVEMENT
const Int g_RCIterationNum = 20;
const Double g_RCWeightHistoryLambda = 0.5;
const Double g_RCWeightCurrentLambda = 1.0 - g_RCWeightHistoryLambda;
const Int g_RCLCUSmoothWindowSize = 4;
const Double g_RCAlphaMinValue = 0.05;
const Double g_RCAlphaMaxValue = 500.0;
const Double g_RCBetaMinValue  = -3.0;
const Double g_RCBetaMaxValue  = -0.1;
#endif

#if RATE_CONTROL_INTRA
#define ALPHA     6.7542;
#define BETA1     1.2517
#define BETA2     1.7860
#endif

struct TRCLCU
{
  Int m_actualBits;
  Int m_QP;     // QP of skip mode is set to g_RCInvalidQPValue
  Int m_targetBits;
  Double m_lambda;
#if M0036_RC_IMPROVEMENT
  Double m_bitWeight;
#else
  Double m_MAD;
#endif
  Int m_numberOfPixel;
#if RATE_CONTROL_INTRA
  Double m_costIntra;
  Int m_targetBitsLeft;
#endif
};

struct TRCParameter
{
  Double m_alpha;
  Double m_beta;
};

class TEncRCSeq
{
public:
  TEncRCSeq();
  ~TEncRCSeq();

public:
#if M0036_RC_IMPROVEMENT
  Void create( Int totalFrames, Int targetBitrate, Int frameRate, Int GOPSize, Int picWidth, Int picHeight, Int LCUWidth, Int LCUHeight, Int numberOfLevel, Bool useLCUSeparateModel, Int adaptiveBit );
#else
  Void create( Int totalFrames, Int targetBitrate, Int frameRate, Int GOPSize, Int picWidth, Int picHeight, Int LCUWidth, Int LCUHeight, Int numberOfLevel, Bool useLCUSeparateModel );
#endif
  Void destroy();
  Void initBitsRatio( Int bitsRatio[] );
  Void initGOPID2Level( Int GOPID2Level[] );
  Void initPicPara( TRCParameter* picPara  = NULL );    // NULL to initial with default value
  Void initLCUPara( TRCParameter** LCUPara = NULL );    // NULL to initial with default value
  Void updateAfterPic ( Int bits );
#if !RATE_CONTROL_INTRA
  Int  getRefineBitsForIntra( Int orgBits );
#endif
#if M0036_RC_IMPROVEMENT
  Void setAllBitRatio( Double basicLambda, Double* equaCoeffA, Double* equaCoeffB );
#endif

public:
  Int  getTotalFrames()                 { return m_totalFrames; }
  Int  getTargetRate()                  { return m_targetRate; }
  Int  getFrameRate()                   { return m_frameRate; }
  Int  getGOPSize()                     { return m_GOPSize; }
  Int  getPicWidth()                    { return m_picWidth; }
  Int  getPicHeight()                   { return m_picHeight; } 
  Int  getLCUWidth()                    { return m_LCUWidth; }
  Int  getLCUHeight()                   { return m_LCUHeight; }
  Int  getNumberOfLevel()               { return m_numberOfLevel; }
  Int  getAverageBits()                 { return m_averageBits; }
  Int  getLeftAverageBits()             { assert( m_framesLeft > 0 ); return (Int)(m_bitsLeft / m_framesLeft); }
  Bool getUseLCUSeparateModel()         { return m_useLCUSeparateModel; }

  Int  getNumPixel()                    { return m_numberOfPixel; }
  Int64  getTargetBits()                { return m_targetBits; }
  Int  getNumberOfLCU()                 { return m_numberOfLCU; }
  Int* getBitRatio()                    { return m_bitsRatio; }
  Int  getBitRatio( Int idx )           { assert( idx<m_GOPSize); return m_bitsRatio[idx]; }
  Int* getGOPID2Level()                 { return m_GOPID2Level; }
  Int  getGOPID2Level( Int ID )         { assert( ID < m_GOPSize ); return m_GOPID2Level[ID]; }
  TRCParameter*  getPicPara()                                   { return m_picPara; }
  TRCParameter   getPicPara( Int level )                        { assert( level < m_numberOfLevel ); return m_picPara[level]; }
  Void           setPicPara( Int level, TRCParameter para )     { assert( level < m_numberOfLevel ); m_picPara[level] = para; }
  TRCParameter** getLCUPara()                                   { return m_LCUPara; }
  TRCParameter*  getLCUPara( Int level )                        { assert( level < m_numberOfLevel ); return m_LCUPara[level]; }
  TRCParameter   getLCUPara( Int level, Int LCUIdx )            { assert( LCUIdx  < m_numberOfLCU ); return getLCUPara(level)[LCUIdx]; }
  Void           setLCUPara( Int level, Int LCUIdx, TRCParameter para ) { assert( level < m_numberOfLevel ); assert( LCUIdx  < m_numberOfLCU ); m_LCUPara[level][LCUIdx] = para; }

  Int  getFramesLeft()                  { return m_framesLeft; }
  Int64  getBitsLeft()                  { return m_bitsLeft; }

  Double getSeqBpp()                    { return m_seqTargetBpp; }
  Double getAlphaUpdate()               { return m_alphaUpdate; }
  Double getBetaUpdate()                { return m_betaUpdate; }

#if M0036_RC_IMPROVEMENT
  Int    getAdaptiveBits()              { return m_adaptiveBit;  }
  Double getLastLambda()                { return m_lastLambda;   }
  Void   setLastLambda( Double lamdba ) { m_lastLambda = lamdba; }
#endif

private:
  Int m_totalFrames;
  Int m_targetRate;
  Int m_frameRate; 
  Int m_GOPSize;
  Int m_picWidth;
  Int m_picHeight;
  Int m_LCUWidth;
  Int m_LCUHeight;
  Int m_numberOfLevel;
  Int m_averageBits;

  Int m_numberOfPixel;
  Int64 m_targetBits;
  Int m_numberOfLCU;
  Int* m_bitsRatio;
  Int* m_GOPID2Level;
  TRCParameter*  m_picPara;
  TRCParameter** m_LCUPara;

  Int m_framesLeft;
  Int64 m_bitsLeft;
  Double m_seqTargetBpp;
  Double m_alphaUpdate;
  Double m_betaUpdate;
  Bool m_useLCUSeparateModel;

#if M0036_RC_IMPROVEMENT
  Int m_adaptiveBit;
  Double m_lastLambda;
#endif
};

class TEncRCGOP
{
public:
  TEncRCGOP();
  ~TEncRCGOP();

public:
  Void create( TEncRCSeq* encRCSeq, Int numPic );
  Void destroy();
  Void updateAfterPicture( Int bitsCost );

private:
  Int  xEstGOPTargetBits( TEncRCSeq* encRCSeq, Int GOPSize );
#if M0036_RC_IMPROVEMENT
  Void   xCalEquaCoeff( TEncRCSeq* encRCSeq, Double* lambdaRatio, Double* equaCoeffA, Double* equaCoeffB, Int GOPSize );
  Double xSolveEqua( Double targetBpp, Double* equaCoeffA, Double* equaCoeffB, Int GOPSize );
#endif

public:
  TEncRCSeq* getEncRCSeq()        { return m_encRCSeq; }
  Int  getNumPic()                { return m_numPic;}
  Int  getTargetBits()            { return m_targetBits; }
  Int  getPicLeft()               { return m_picLeft; }
  Int  getBitsLeft()              { return m_bitsLeft; }
  Int  getTargetBitInGOP( Int i ) { return m_picTargetBitInGOP[i]; }

private:
  TEncRCSeq* m_encRCSeq;
  Int* m_picTargetBitInGOP;
  Int m_numPic;
  Int m_targetBits;
  Int m_picLeft;
  Int m_bitsLeft;
};

class TEncRCPic
{
public:
  TEncRCPic();
  ~TEncRCPic();

public:
  Void create( TEncRCSeq* encRCSeq, TEncRCGOP* encRCGOP, Int frameLevel, list<TEncRCPic*>& listPreviousPictures );
  Void destroy();

#if !RATE_CONTROL_INTRA
  Double estimatePicLambda( list<TEncRCPic*>& listPreviousPictures );
#endif
  Int    estimatePicQP    ( Double lambda, list<TEncRCPic*>& listPreviousPictures );
#if RATE_CONTROL_INTRA
  Int    getRefineBitsForIntra(Int orgBits);
  Double calculateLambdaIntra(double alpha, double beta, double MADPerPixel, double bitsPerPixel);
  Double estimatePicLambda( list<TEncRCPic*>& listPreviousPictures, SliceType eSliceType);

  Void   updateAlphaBetaIntra(double *alpha, double *beta);

  Double getLCUTargetBpp(SliceType eSliceType);
  Double getLCUEstLambdaAndQP(Double bpp, Int clipPicQP, Int *estQP);
#else
  Double getLCUTargetBpp();
#endif
  Double getLCUEstLambda( Double bpp );
  Int    getLCUEstQP( Double lambda, Int clipPicQP );

  Void updateAfterLCU( Int LCUIdx, Int bits, Int QP, Double lambda, Bool updateLCUParameter = true );
#if M0036_RC_IMPROVEMENT
#if RATE_CONTROL_INTRA
  Void updateAfterPicture( Int actualHeaderBits, Int actualTotalBits, Double averageQP, Double averageLambda, SliceType eSliceType);
#else
  Void updateAfterPicture( Int actualHeaderBits, Int actualTotalBits, Double averageQP, Double averageLambda );
#endif
#else
  Void updateAfterPicture( Int actualHeaderBits, Int actualTotalBits, Double averageQP, Double averageLambda, Double effectivePercentage );
#endif

  Void addToPictureLsit( list<TEncRCPic*>& listPreviousPictures );
#if !M0036_RC_IMPROVEMENT
  Double getEffectivePercentage();
#endif
  Double calAverageQP();
  Double calAverageLambda();

private:
  Int xEstPicTargetBits( TEncRCSeq* encRCSeq, TEncRCGOP* encRCGOP );
  Int xEstPicHeaderBits( list<TEncRCPic*>& listPreviousPictures, Int frameLevel );

public:
  TEncRCSeq*      getRCSequence()                         { return m_encRCSeq; }
  TEncRCGOP*      getRCGOP()                              { return m_encRCGOP; }

  Int  getFrameLevel()                                    { return m_frameLevel; }
  Int  getNumberOfPixel()                                 { return m_numberOfPixel; }
  Int  getNumberOfLCU()                                   { return m_numberOfLCU; }
  Int  getTargetBits()                                    { return m_targetBits; }
#if !RATE_CONTROL_INTRA 
  Void setTargetBits( Int bits )                          { m_targetBits = bits; }
#endif
  Int  getEstHeaderBits()                                 { return m_estHeaderBits; }
  Int  getLCULeft()                                       { return m_LCULeft; }
  Int  getBitsLeft()                                      { return m_bitsLeft; }
  Int  getPixelsLeft()                                    { return m_pixelsLeft; }
  Int  getBitsCoded()                                     { return m_targetBits - m_estHeaderBits - m_bitsLeft; }
  Int  getLCUCoded()                                      { return m_numberOfLCU - m_LCULeft; }
  TRCLCU* getLCU()                                        { return m_LCUs; }
  TRCLCU& getLCU( Int LCUIdx )                            { return m_LCUs[LCUIdx]; }
  Int  getPicActualHeaderBits()                           { return m_picActualHeaderBits; }
#if !M0036_RC_IMPROVEMENT
  Double getTotalMAD()                                    { return m_totalMAD; }
  Void   setTotalMAD( Double MAD )                        { m_totalMAD = MAD; }
#endif
#if RATE_CONTROL_INTRA
  Void setTargetBits( Int bits )                          { m_targetBits = bits; m_bitsLeft = bits;}
  Void setTotalIntraCost(Double cost)                     { m_totalCostIntra = cost; }
  Void getLCUInitTargetBits();
#endif

  Int  getPicActualBits()                                 { return m_picActualBits; }
  Int  getPicActualQP()                                   { return m_picQP; }
  Double getPicActualLambda()                             { return m_picLambda; }
  Int  getPicEstQP()                                      { return m_estPicQP; }
  Void setPicEstQP( Int QP )                              { m_estPicQP = QP; }
  Double getPicEstLambda()                                { return m_estPicLambda; }
  Void setPicEstLambda( Double lambda )                   { m_picLambda = lambda; }

private:
  TEncRCSeq* m_encRCSeq;
  TEncRCGOP* m_encRCGOP;

  Int m_frameLevel;
  Int m_numberOfPixel;
  Int m_numberOfLCU;
  Int m_targetBits;
  Int m_estHeaderBits;
  Int m_estPicQP;
  Double m_estPicLambda;

  Int m_LCULeft;
  Int m_bitsLeft;
  Int m_pixelsLeft;

  TRCLCU* m_LCUs;
  Int m_picActualHeaderBits;    // only SH and potential APS
#if !M0036_RC_IMPROVEMENT
  Double m_totalMAD;
#endif
#if RATE_CONTROL_INTRA
  Double m_totalCostIntra; 
  Double m_remainingCostIntra;
#endif
  Int m_picActualBits;          // the whole picture, including header
  Int m_picQP;                  // in integer form
  Double m_picLambda;
#if !M0036_RC_IMPROVEMENT
  TEncRCPic* m_lastPicture;
#endif
};

class TEncRateCtrl
{
public:
  TEncRateCtrl();
  ~TEncRateCtrl();

public:
#if M0036_RC_IMPROVEMENT
  Void init( Int totalFrames, Int targetBitrate, Int frameRate, Int GOPSize, Int picWidth, Int picHeight, Int LCUWidth, Int LCUHeight, Int keepHierBits, Bool useLCUSeparateModel, GOPEntry GOPList[MAX_GOP] );
#else
  Void init( Int totalFrames, Int targetBitrate, Int frameRate, Int GOPSize, Int picWidth, Int picHeight, Int LCUWidth, Int LCUHeight, Bool keepHierBits, Bool useLCUSeparateModel, GOPEntry GOPList[MAX_GOP] );
#endif
  Void destroy();
  Void initRCPic( Int frameLevel );
  Void initRCGOP( Int numberOfPictures );
  Void destroyRCGOP();

public:
  Void       setRCQP ( Int QP ) { m_RCQP = QP;   }
  Int        getRCQP ()         { return m_RCQP; }
  TEncRCSeq* getRCSeq()          { assert ( m_encRCSeq != NULL ); return m_encRCSeq; }
  TEncRCGOP* getRCGOP()          { assert ( m_encRCGOP != NULL ); return m_encRCGOP; }
  TEncRCPic* getRCPic()          { assert ( m_encRCPic != NULL ); return m_encRCPic; }
  list<TEncRCPic*>& getPicList() { return m_listRCPictures; }

private:
  TEncRCSeq* m_encRCSeq;
  TEncRCGOP* m_encRCGOP;
  TEncRCPic* m_encRCPic;
  list<TEncRCPic*> m_listRCPictures;
  Int        m_RCQP;
};
#elif defined(X264_RATECONTROL_2006)

enum slice_type_e
{
	SLICE_TYPE_P  = 0,
	SLICE_TYPE_B  = 1,
	SLICE_TYPE_I  = 2,
};

/* CPU flags
*/
#define X264_CPU_MMX        0x000001    /* mmx */
#define X264_CPU_MMXEXT     0x000002    /* mmx-ext*/
#define X264_CPU_SSE        0x000004    /* sse */
#define X264_CPU_SSE2       0x000008    /* sse 2 */
#define X264_CPU_3DNOW      0x000010    /* 3dnow! */
#define X264_CPU_3DNOWEXT   0x000020    /* 3dnow! ext */
#define X264_CPU_ALTIVEC    0x000040    /* altivec */

/* Analyse flags
*/
#define X264_ANALYSE_I4x4       0x0001  /* Analyse i4x4 */
#define X264_ANALYSE_I8x8       0x0002  /* Analyse i8x8 (requires 8x8 transform) */
#define X264_ANALYSE_PSUB16x16  0x0010  /* Analyse p16x8, p8x16 and p8x8 */
#define X264_ANALYSE_PSUB8x8    0x0020  /* Analyse p8x4, p4x8, p4x4 */
#define X264_ANALYSE_BSUB16x16  0x0100  /* Analyse b16x8, b8x16 and b8x8 */
#define X264_DIRECT_PRED_NONE        0
#define X264_DIRECT_PRED_SPATIAL     1
#define X264_DIRECT_PRED_TEMPORAL    2
#define X264_DIRECT_PRED_AUTO        3
#define X264_ME_DIA                  0
#define X264_ME_HEX                  1
#define X264_ME_UMH                  2
#define X264_ME_ESA                  3
#define X264_CQM_FLAT                0
#define X264_CQM_JVT                 1
#define X264_CQM_CUSTOM              2

static const char * const x264_direct_pred_names[] = { "none", "spatial", "temporal", "auto", 0 };
static const char * const x264_motion_est_names[] = { "dia", "hex", "umh", "esa", 0 };

/* Colorspace type
*/
#define X264_CSP_MASK           0x00ff  /* */
#define X264_CSP_NONE           0x0000  /* Invalid mode     */
#define X264_CSP_I420           0x0001  /* yuv 4:2:0 planar */
#define X264_CSP_I422           0x0002  /* yuv 4:2:2 planar */
#define X264_CSP_I444           0x0003  /* yuv 4:4:4 planar */
#define X264_CSP_YV12           0x0004  /* yuv 4:2:0 planar */
#define X264_CSP_YUYV           0x0005  /* yuv 4:2:2 packed */
#define X264_CSP_RGB            0x0006  /* rgb 24bits       */
#define X264_CSP_BGR            0x0007  /* bgr 24bits       */
#define X264_CSP_BGRA           0x0008  /* bgr 32bits       */
#define X264_CSP_VFLIP          0x1000  /* */

/* Slice type
*/
#define X264_TYPE_AUTO          0x0000  /* Let x264 choose the right type */
#define X264_TYPE_IDR           0x0001
#define X264_TYPE_I             0x0002
#define X264_TYPE_P             0x0003
#define X264_TYPE_BREF          0x0004  /* Non-disposable B-frame */
#define X264_TYPE_B             0x0005
#define IS_X264_TYPE_I(x) ((x)==X264_TYPE_I || (x)==X264_TYPE_IDR)
#define IS_X264_TYPE_B(x) ((x)==X264_TYPE_B || (x)==X264_TYPE_BREF)

/* Log level
*/
#define X264_LOG_NONE          (-1)
#define X264_LOG_ERROR          0
#define X264_LOG_WARNING        1
#define X264_LOG_INFO           2
#define X264_LOG_DEBUG          3

#define _NEW_SATD_EST_ 1
#define RC_P_WINDOW 8
#define _OVERFLOW_ADJUST_ 0
#define _USE_STD_PRED_ 1
#define _USE_QPOFFSET_ 1
#define _USE_SMALL_BIG_P_ 0
#define _USE_LCU_ 0


#define X264_RC_CQP                  0 //Constant quantizer, the QPs are simply based on whether the frame is P,I or B frame.
#define X264_RC_CRF                  1 //Constant rate factor, one pass mode that is optimal if the user doesn't desire a specific bitrate,but specify quality instead.
#define X264_RC_ABR                  2 //Average bitrate, One pass scheme

typedef struct
{
	int i_start, i_end;
	int b_force_qp;
	int i_qp;
	float f_bitrate_factor;
} x264_zone_t;

typedef struct
{
	/* Video Properties */
	int         i_width;
	int         i_height;

	int         i_fps_num;
	int         i_fps_den;

	/* Bitstream parameters */
	int         i_frame_reference;  /* Maximum number of reference frames */
	int         i_keyint_max;       /* Force an IDR keyframe at this interval */
	int         i_keyint_min;       /* Scenecuts closer together than this are coded as I, not IDR. */
	int         i_scenecut_threshold; /* how aggressively to insert extra I frames */
	int         i_bframe;   /* how many b-frame between 2 references pictures */
	int         b_bframe_adaptive;
	int         i_bframe_bias;
	int         b_bframe_pyramid;   /* Keep some B-frames as references */

	int         b_deblocking_filter;
	int         i_deblocking_filter_alphac0;    /* [-6, 6] -6 light filter, 6 strong */
	int         i_deblocking_filter_beta;       /* [-6, 6]  idem */

	int         b_cabac;
	int         i_cabac_init_idc;

#if 0
	int         i_cqm_preset;
	char        *psz_cqm_file;      /* JM format */
	uint8_t     cqm_4iy[16];        /* used only if i_cqm_preset == X264_CQM_CUSTOM */
	uint8_t     cqm_4ic[16];
	uint8_t     cqm_4py[16];
	uint8_t     cqm_4pc[16];
	uint8_t     cqm_8iy[64];
	uint8_t     cqm_8py[64];
#endif

	/* Log */
	void        (*pf_log)( void *, int i_level, const char *psz, va_list );
	void        *p_log_private;
	int         i_log_level;
	int         b_visualize;

	/* Encoder analyser parameters */
	struct
	{
		unsigned int intra;     /* intra partitions */
		unsigned int inter;     /* inter partitions */

		int          b_transform_8x8;
		int          b_weighted_bipred; /* implicit weighting for B-frames */
		int          i_direct_mv_pred; /* spatial vs temporal mv prediction */
		int          i_chroma_qp_offset;

		int          i_me_method; /* motion estimation algorithm to use (X264_ME_*) */
		int          i_me_range; /* integer pixel motion estimation search range (from predicted mv) */
		int          i_mv_range; /* maximum length of a mv (in pixels) */
		int          i_subpel_refine; /* subpixel motion estimation quality */
		int          b_bidir_me; /* jointly optimize both MVs in B-frames */
		int          b_chroma_me; /* chroma ME for subpel and mode decision in P-frames */
		int          b_bframe_rdo; /* RD based mode decision for B-frames */
		int          b_mixed_references; /* allow each mb partition in P-frames to have it's own reference number */
		int          i_trellis;  /* trellis RD quantization */
		int          b_fast_pskip; /* early SKIP detection on P-frames */
		int          i_noise_reduction; /* adaptive pseudo-deadzone */

		int          b_psnr;    /* Do we compute PSNR stats (save a few % of cpu) */
	} analyse;

	/* Rate control parameters */
	struct
	{
        int         i_rc_method;    /* X264_RC_* */
		int         i_qp_constant;  /* 0-51 */
		int         i_qp_min;       /* min allowed QP value */
		int         i_qp_max;       /* max allowed QP value */
		int         i_qp_step;      /* max QP step between frames */

		int         b_cbr;          /* use bitrate instead of CQP */
		int b_lcurc;
		int         i_bitrate;
		int         i_rf_constant;  /* 1pass VBR, nominal QP */
		float       f_rate_tolerance;
		int         i_vbv_max_bitrate;
		int         i_vbv_buffer_size;
		float       f_vbv_buffer_init;
		float       f_ip_factor;
		float       f_pb_factor;

		/* 2pass params (same as ffmpeg ones) */
		float       f_qcompress;    /* 0.0 => cbr, 1.0 => constant qp */
		float       f_qblur;        /* temporally blur quants */
		float       f_complexity_blur; /* temporally blur complexity */
	} rc;

	int b_aud;                  /* generate access unit delimiters */
	int b_repeat_headers;       /* put SPS/PPS before each keyframe */
	int m_numberOfLCU;
	int b_variable_qp;
	int picWidthInBU;
	int picHeightInBU;
} x264_param_t;

typedef struct
{
	int pict_type;
    int frame_type;
	int kept_as_ref;
	float qscale;
	int mv_bits;
	int i_tex_bits;
	int p_tex_bits;
	int tex_bits;
	int misc_bits;
    uint64_t expected_bits; /*total expected bits up to the current frame (current one excluded)*/
    double expected_vbv;
    double new_qscale;
	int new_qp;
	int i_count;
	int p_count;
	int s_count;
	float blurred_complexity;
	char direct_mode;
    int16_t weight[3][2];
    int16_t i_weight_denom[2];
    int refcount[16];
    int refs;
    int64_t i_duration;
    int64_t i_cpb_duration;
} ratecontrol_entry_t;

typedef struct
{
	double coeff_min;
	double coeff;
	double count;
	double decay;
	double offset;
} predictor_t;

struct x264_ratecontrol_t
{
	/* constants */
	int b_abr;
	int b_2pass;
	int b_vbv;
	int b_vbv_min_rate;
	double fps;
	double bitrate;
	double rate_tolerance;
	double qcompress;
	int nmb;                    /* number of macroblocks in a frame */
	int qp_constant[3];
	int gop_id;

	/* current frame */
	ratecontrol_entry_t *rce;
	int qp;                     /* qp for current frame */
	float qpm;                    /* qp for current macroblock */
	float qpa;                  /* average of macroblocks' qp */
	float qpa_prev;
	float qp_novbv;

	int slice_type;
	int qp_offset;
	double qp_factor;

	/* VBV stuff */
    double buffer_size;
    int64_t buffer_fill_final;
    double buffer_fill;         /* planned buffer, if all in-progress frames hit their bit budget */
    double buffer_rate;         /* # of bits added to buffer_fill after each frame */
    double vbv_max_rate;        /* # of bits added to buffer_fill per second */

	predictor_t pred[3];        /* predict frame size from satd */
	int last_satd_for[3];
	int bits_for[3];
	int ftype;

    int single_frame_vbv;
    double rate_factor_max_increment; /* Don't allow RF above (CRF + this value). */

	/* ABR stuff */
	int    last_satd;
	double last_rceq;
	double last_rceq2;
	double cplxr_sum;           /* sum of bits*qscale/rceq */
	double expected_bits_sum;   /* sum of qscale2bits after rceq, ratefactor, and overflow */
    int64_t filler_bits_sum;    /* sum in bits of finished frames' filler data */
	double wanted_bits_window;  /* target bitrate * window */
//	double wanted_bits_window2;
//	double cplxr_sum2;
	double cbr_decay;
	double short_term_cplxsum;
	double short_term_cplxcount;
	double rate_factor_constant;
	double ip_offset;
	double pb_offset;

	int num_entries;            /* number of ratecontrol_entry_ts */
	ratecontrol_entry_t *entry; /* FIXME: copy needed data and free this once init is done */
	double last_qscale;
//	double last_qscale2;
	double last_qscale_for[3];  /* last qscale for a specific pict type, used for max_diff & ipb factor stuff  */
	int last_non_b_pict_type;
	double accum_p_qp;          /* for determining I-frame quant */
	double accum_p_norm;
	double last_accum_p_norm;
	double lmin[3];             /* min qscale by frame type */
	double lmax[3];
	double lstep;               /* max change (multiply) in qscale per frame */

	/* MBRC stuff */
    float frame_size_estimated; /* Access to this variable must be atomic: double is
                                 * not atomic on all arches we care about */
    double frame_size_maximum;  /* Maximum frame size due to MinCR */
    double frame_size_planned;
    double slice_size_planned;
	int first_row, last_row;    /* region of the frame to be encoded by this thread */
	predictor_t *row_pred;//[2];
	predictor_t row_preds[3];//[2];
	predictor_t pred_b_from_p;  /* predict B-frame size from P-frame satd */
	int bframes;                /* # consecutive B-frames before this P-frame */
	int bframe_bits;            /* total cost of those frames */

//add
	int64_t bitcost;
	int i_frame;
	int i_mb_x;
	int i_mb_y;
	int i_type_last;
	int     i_slice_count[3];
	int     *i_row_satd;
	int i_row_satd_tmp;
	int     *i_row_bits;
	float *i_row_qp;
	int     *i_row_satd_last;
	int     *i_row_bits_last;
	float *i_row_qp_last;

	int *lcu_sad;
	int lcu_sad_avg;
	int lcu_sad_sum;
	int lcu_idx;
#if _NEW_SATD_EST_
	double framesad_Pavg;
	double sad_Ilast;
	double sad_Plast;
	double sad_Pfrm[RC_P_WINDOW];
	int i_statvalid_Pfrm;
	int i_index_Pfrm;
	int i_numvalid_Pfrm;
#endif
#if _USE_STD_PRED_
	double std_val;
#endif

};

#define MAX_DELTA_QP    2
#define MAX_CUDQP_DEPTH 0 

int x264_ratecontrol_new( x264_ratecontrol_t* rc, x264_param_t* pParam, int lcuwidth, int lcuheight);
void x264_ratecontrol_delete( x264_ratecontrol_t* rc,x264_param_t* pParam );
void x264_ratecontrol_start( x264_ratecontrol_t *rc, x264_param_t* pParam, int i_slice_type, int i_force_qp );
void x264_ratecontrol_end( x264_ratecontrol_t *rc, x264_param_t* pParam,int bits, int cost);
int x264_ratecontrol_qp( x264_ratecontrol_t *rc );
void x264_ratecontrol_mb( x264_ratecontrol_t *rc, x264_param_t* pParam, int bits, int cost );

#else

// ====================================================================================================================
// Class definition
// ====================================================================================================================
#define MAX_DELTA_QP    2
#define MAX_CUDQP_DEPTH 0 

typedef struct FrameData
{
  Bool       m_isReferenced;
  Int        m_qp;
  Int        m_bits;
  Double     m_costMAD;
}FrameData;

typedef struct LCUData
{
  Int     m_qp;                ///<  coded QP
  Int     m_bits;              ///<  actually generated bits
  Int     m_pixels;            ///<  number of pixels for a unit
  Int     m_widthInPixel;      ///<  number of pixels for width
  Int     m_heightInPixel;     ///<  number of pixels for height
  Double  m_costMAD;           ///<  texture complexity for a unit
}LCUData;

class MADLinearModel
{
private:
  Bool   m_activeOn;
  Double m_paramY1;
  Double m_paramY2;
  Double m_costMADs[3];

public:
  MADLinearModel ()   {};
  ~MADLinearModel()   {};
  
  Void    initMADLinearModel      ();
  Double  getMAD                  ();
  Void    updateMADLiearModel     ();
  Void    updateMADHistory        (Double costMAD);
  Bool    IsUpdateAvailable       ()              { return m_activeOn; }
};

class PixelBaseURQQuadraticModel
{
private:
  Double m_paramHighX1;
  Double m_paramHighX2;
  Double m_paramLowX1;
  Double m_paramLowX2;
public:
  PixelBaseURQQuadraticModel () {};
  ~PixelBaseURQQuadraticModel() {};

  Void    initPixelBaseQuadraticModel       ();
  Int     getQP                             (Int qp, Int targetBits, Int numberOfPixels, Double costPredMAD);
  Void    updatePixelBasedURQQuadraticModel (Int qp, Int bits, Int numberOfPixels, Double costMAD);
  Bool    checkUpdateAvailable              (Int qpReference );
  Double  xConvertQP2QStep                  (Int qp );
  Int     xConvertQStep2QP                  (Double qStep );
};

class TEncRateCtrl
{
private:
  Bool            m_isLowdelay;
  Int             m_prevBitrate;
  Int             m_currBitrate;
  Int             m_frameRate;
  Int             m_refFrameNum;
  Int             m_nonRefFrameNum;
  Int             m_numOfPixels;
  Int             m_sourceWidthInLCU;
  Int             m_sourceHeightInLCU;      
  Int             m_sizeGOP;
  Int             m_indexGOP;
  Int             m_indexFrame;
  Int             m_indexLCU;
  Int             m_indexUnit;
  Int             m_indexRefFrame;
  Int             m_indexNonRefFrame;
  Int             m_indexPOCInGOP;
  Int             m_indexPrevPOCInGOP;
  Int             m_occupancyVB;
  Int             m_initialOVB;
  Int             m_targetBufLevel;
  Int             m_initialTBL;
  Int             m_remainingBitsInGOP;
  Int             m_remainingBitsInFrame;
  Int             m_occupancyVBInFrame;
  Int             m_targetBits;
  Int             m_numUnitInFrame;
  Int             m_codedPixels;
  Bool            m_activeUnitLevelOn;
  Double          m_costNonRefAvgWeighting;
  Double          m_costRefAvgWeighting;
  Double          m_costAvgbpp;         
  
  FrameData*      m_pcFrameData;
  LCUData*        m_pcLCUData;

  MADLinearModel              m_cMADLinearModel;
  PixelBaseURQQuadraticModel  m_cPixelURQQuadraticModel;
  
public:
  TEncRateCtrl         () {};
  virtual ~TEncRateCtrl() {};

  Void          create                (Int sizeIntraPeriod, Int sizeGOP, Int frameRate, Int targetKbps, Int qp, Int numLCUInBasicUnit, Int sourceWidth, Int sourceHeight, Int maxCUWidth, Int maxCUHeight);
  Void          destroy               ();

  Void          initFrameData         (Int qp = 0);
  Void          initUnitData          (Int qp = 0);
  Int           getFrameQP            (Bool isReferenced, Int POC);
  Bool          calculateUnitQP       ();
  Int           getUnitQP             ()                                          { return m_pcLCUData[m_indexLCU].m_qp;  }
  Void          updateRCGOPStatus     ();
  Void          updataRCFrameStatus   (Int frameBits, SliceType eSliceType);
  Void          updataRCUnitStatus    ();
  Void          updateLCUData         (TComDataCU* pcCU, UInt64 actualLCUBits, Int qp);
  Void          updateFrameData       (UInt64 actualFrameBits);
  Double        xAdjustmentBits       (Int& reductionBits, Int& compensationBits);
  Int           getGOPId              ()                                          { return m_indexFrame; }
};
#endif

#endif


