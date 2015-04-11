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

/** \file     encmain.cpp
    \brief    Encoder application main
*/

#include <time.h>
#include <iostream>
#include "TAppEncTop.h"
#include "TAppCommon/program_options_lite.h"

using namespace std;
namespace po = df::program_options_lite;
#if _HFZ_CABAC_
ContextModel         g_contextModels[MAX_NUM_CTX_MOD];
Int                  g_numContextModels=0;
ContextModel3DBuffer  g_cCUSplitFlagSCModel       ( 1,             1, NUM_SPLIT_FLAG_CTX            , g_contextModels + g_numContextModels, g_numContextModels );
ContextModel3DBuffer  g_cCUSkipFlagSCModel        ( 1,             1, NUM_SKIP_FLAG_CTX             , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUMergeFlagExtSCModel    ( 1,             1, NUM_MERGE_FLAG_EXT_CTX        , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUMergeIdxExtSCModel     ( 1,             1, NUM_MERGE_IDX_EXT_CTX         , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUPartSizeSCModel        ( 1,             1, NUM_PART_SIZE_CTX             , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUPredModeSCModel        ( 1,             1, NUM_PRED_MODE_CTX             , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUIntraPredSCModel       ( 1,             1, NUM_ADI_CTX                   , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUChromaPredSCModel      ( 1,             1, NUM_CHROMA_PRED_CTX           , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUDeltaQpSCModel         ( 1,             1, NUM_DELTA_QP_CTX              , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUInterDirSCModel        ( 1,             1, NUM_INTER_DIR_CTX             , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCURefPicSCModel          ( 1,             1, NUM_REF_NO_CTX                , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUMvdSCModel             ( 1,             1, NUM_MV_RES_CTX                , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUQtCbfSCModel           ( 1,             2, NUM_QT_CBF_CTX                , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUTransSubdivFlagSCModel ( 1,             1, NUM_TRANS_SUBDIV_FLAG_CTX     , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUQtRootCbfSCModel       ( 1,             1, NUM_QT_ROOT_CBF_CTX           , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUSigCoeffGroupSCModel   ( 1,             2, NUM_SIG_CG_FLAG_CTX           , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUSigSCModel             ( 1,             1, NUM_SIG_FLAG_CTX              , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCuCtxLastX               ( 1,             2, NUM_CTX_LAST_FLAG_XY          , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCuCtxLastY               ( 1,             2, NUM_CTX_LAST_FLAG_XY          , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUOneSCModel             ( 1,             1, NUM_ONE_FLAG_CTX              , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUAbsSCModel             ( 1,             1, NUM_ABS_FLAG_CTX              , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cMVPIdxSCModel            ( 1,             1, NUM_MVP_IDX_CTX               , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cCUAMPSCModel             ( 1,             1, NUM_CU_AMP_CTX                , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cSaoMergeSCModel          ( 1,             1, NUM_SAO_MERGE_FLAG_CTX   , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cSaoTypeIdxSCModel        ( 1,             1, NUM_SAO_TYPE_IDX_CTX          , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_cTransformSkipSCModel     ( 1,             2, NUM_TRANSFORMSKIP_FLAG_CTX    , g_contextModels + g_numContextModels, g_numContextModels);
ContextModel3DBuffer  g_CUTransquantBypassFlagSCModel( 1,          1, NUM_CU_TRANSQUANT_BYPASS_FLAG_CTX, g_contextModels + g_numContextModels, g_numContextModels);
#endif
//! \ingroup TAppEncoder
//! \{

// ====================================================================================================================
// Main function
// ====================================================================================================================

int main(int argc, char* argv[])
{
  TAppEncTop  cTAppEncTop;

  // print information
  fprintf( stdout, "\n" );
  fprintf( stdout, "HM software: Encoder Version [%s]", NV_VERSION );
  fprintf( stdout, NVM_ONOS );
  fprintf( stdout, NVM_COMPILEDBY );
  fprintf( stdout, NVM_BITS );
  fprintf( stdout, "\n" );

  // create application encoder class
  cTAppEncTop.create();

  // parse configuration
  try
  {
    if(!cTAppEncTop.parseCfg( argc, argv ))
    {
      cTAppEncTop.destroy();
      return 1;
    }
  }
  catch (po::ParseFailure& e)
  {
    cerr << "Error parsing option \""<< e.arg <<"\" with argument \""<< e.val <<"\"." << endl;
    return 1;
  }

  // starting time
  double dResult;
  long lBefore = clock();

  // call encoding function
  cTAppEncTop.encode();

  // ending time
  dResult = (double)(clock()-lBefore) / CLOCKS_PER_SEC;
  printf("\n Total Time: %12.3f sec.\n", dResult);

  // destroy application encoder class
  cTAppEncTop.destroy();

  return 0;
}

//! \}
