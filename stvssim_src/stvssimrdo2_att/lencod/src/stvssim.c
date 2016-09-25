#include "stvssim.h"
#include "global.h"
#include "defines.h"
#include "memalloc.h"
#include "mbuffer.h"
#include "enc_statistics.h"

extern InputParameters *params;
extern ImageParameters *img;
extern StorablePicture *enc_picture;
extern StatParameters *stats;
void InitialRDOVar()
{
int i;
	int height=params->source.height;
	int width=params->source.width;
for(i=0;i<REFNUM;i++)
{
  get_mem2Dpel(&refPicsData[i][0],height,width);
  get_mem2Dpel(&refPicsData[i][1],height/2,width/2);
  get_mem2Dpel(&refPicsData[i][2],height/2,width/2);

  get_mem2Dpel(&srcPicsData[i][0],height,width);
  get_mem2Dpel(&srcPicsData[i][1],height/2,width/2);
  get_mem2Dpel(&srcPicsData[i][2],height/2,width/2);
}
	get_mem4Dshort(&mb_mv,MODENUM,4,4,2);
	get_mem2Dfloat(&mb_directions,4,4);
//	get_mem2Dfloat(&mb_directions2,16,16);
    get_mem2Dfloat(&pic_directions2,height,width);
#if _MB_RDDATA_
	if (strlen (params->MbRDdataFile) > 0 && (pmbrddata=fopen(params->MbRDdataFile,"w"))==NULL)
	{
		snprintf(errortext, ET_SIZE, "Error open file %s", params->MbRDdataFile);
		error (errortext, 500);
	}
	get_mem2Dfloat(&mb_rddata,height/16*width/16,MB_DATA_NUM);
#endif

   if(params->Distortion[SSIM]==1)
   {

  if (strlen (params->InfoFile) > 0 && (pInfo=fopen(params->InfoFile,"w"))==NULL)
  {
    snprintf(errortext, ET_SIZE, "Error open file %s", params->InfoFile);
    error (errortext, 500);
  }
 
  if (strlen (params->AvgFile) > 0 && (pAvg=fopen(params->AvgFile,"w"))==NULL)
  {
    snprintf(errortext, ET_SIZE, "Error open file %s", params->AvgFile);
    error (errortext, 500);
  }
   }

}
void ReleaseRDOVar()
{
	int i;
for(i=0;i<REFNUM;i++)
{
  free_mem2Dpel(refPicsData[i][0]);
  free_mem2Dpel(refPicsData[i][1]);
  free_mem2Dpel(refPicsData[i][2]);

  free_mem2Dpel(srcPicsData[i][0]);
  free_mem2Dpel(srcPicsData[i][1]);
  free_mem2Dpel(srcPicsData[i][2]);
}
	free_mem2Dfloat(pic_directions2);

	if(params->Distortion[SSIM]==1)
	{
		fclose(pInfo);
		fclose(pAvg);
	}
#if _MB_RDDATA_
	fclose(pmbrddata);
	free_mem2Dfloat(mb_rddata);
#endif

	free_mem4Dshort(mb_mv);
	free_mem2Dfloat(mb_directions);
//	free_mem2Dfloat(mb_directions2);
}
int get_mem2Dfloat(float ***array2D, int dim0, int dim1)
{
  int i;

  if((  *array2D  = (float**)calloc(dim0,       sizeof(float*))) == NULL)
    no_mem_exit("get_mem2Dfloat: array2D");
  if((*(*array2D) = (float* )calloc(dim0 * dim1,sizeof(float ))) == NULL)
    no_mem_exit("get_mem2Dfloat: array2D");

  for(i = 1; i < dim0; i++)
    (*array2D)[i] = (*array2D)[i-1] + dim1;

  return dim0 * dim1 * sizeof(float);
}
void free_mem2Dfloat(float **array2D)
{
  if (array2D)
  {
    if (*array2D)
      free (*array2D);
    else error ("free_mem2Dfloat: trying to free unused memory",100);

    free (array2D);
  } 
  else
  {
    error ("free_mem2Dfloat: trying to free unused memory",100);
  }
}

float vFilter_4(int beta,int y,int x,float wa,float wb)
{
	float weight;
	if(wa<0)
	{
		wa=1.0f;
	}
	if(wb<0)
	{
		wb=1.0f;
	}
	if(wa<wb)
	{
		float c=wb;
		wb=wa;
		wa=c;
	}
	if(y==beta/2-1)
	{
		weight=wa;
	}
	else
	{
		weight=wb;
	}
	return weight;
}

float hFilter_4(int beta,int y,int x,float wa,float wb)
{
	float weight;
	if(wa<0)
	{
		wa=1.0f;
	}
	if(wb<0)
	{
		wb=1.0f;
	}
	if(wa<wb)
	{
		float c=wb;
		wb=wa;
		wa=c;
	}
	if(x==beta/2-1)
	{
		weight=wa;
	}
	else
	{
		weight=wb;
	}
	return weight;
}
float lFilter_4(int beta,int y,int x,float wa,float wb)
{
	float weight;
	if(wa<0)
	{
		wa=1.0f;
	}
	if(wb<0)
	{
		wb=1.0f;
	}
	if(wa<wb)
	{
		float c=wb;
		wb=wa;
		wa=c;
	}
	if(x==y)
	{
		weight=wa;
	}
	else
	{
		weight=wb;
	}
	return weight;
}
float rFilter_4(int beta,int y,int x,float wa,float wb)
{
	float weight;
	if(wa<0)
	{
		wa=1.0f;
	}
	if(wb<0)
	{
		wb=1.0f;
	}
	if(wa<wb)
	{
		float c=wb;
		wb=wa;
		wa=c;
	}
	if(x+y==beta-1)
	{
		weight=wa;
	}
	else
	{
		weight=wb;
	}
	return weight;
}

float vFilter(int beta,int y,int x,float wa,float wb)
{
	float weight;
	if(wa<0)
	{
		wa=1.0f;
	}
	if(wb<0)
	{
		wb=1.0f;
	}
	if(wa<wb)
	{
		float c=wb;
		wb=wa;
		wa=c;
	}
	if(y>=beta/2-1&&y<=beta/2+1)
	{
		weight=wa;
	}
	else
	{
		weight=wb;
	}
	return weight;
}

float hFilter(int beta,int y,int x,float wa,float wb)
{
	float weight;
	if(wa<0)
	{
		wa=1.0f;
	}
	if(wb<0)
	{
		wb=1.0f;
	}
	if(wa<wb)
	{
		float c=wb;
		wb=wa;
		wa=c;
	}
	if(x>=beta/2-1&&x<=beta/2+1)
	{
		weight=wa;
	}
	else
	{
		weight=wb;
	}
	return weight;
}
float lFilter(int beta,int y,int x,float wa,float wb)
{
	float weight;
	if(wa<0)
	{
		wa=1.0f;
	}
	if(wb<0)
	{
		wb=1.0f;
	}
	if(wa<wb)
	{
		float c=wb;
		wb=wa;
		wa=c;
	}
	if(abs(x-y)<=1)
	{
		weight=wa;
	}
	else
	{
		weight=wb;
	}
	return weight;
}
float rFilter(int beta,int y,int x,float wa,float wb)
{
	float weight;
	if(wa<0)
	{
		wa=1.0f;
	}
	if(wb<0)
	{
		wb=1.0f;
	}
	if(wa<wb)
	{
		float c=wb;
		wb=wa;
		wa=c;
	}
	if(x+y-beta>=-2&&x+y-beta<=0)
	{
		weight=wa;
	}
	else
	{
		weight=wb;
	}
	return weight;
}

void calOrit(float tmporient,short *orit)
{
	short k;
	short inx;
	short j;
  float distn[FLTNUM]={0.0f};
  float distx=10000.0f;
	for(k=0;k<FLTNUM;++k)
	{
		distn[k]=(float)fabs(tmporient-orient[k]);
		if(distn[k]<distx)
		{
			distx=distn[k];
			inx=k;
		}
	}
	j=0;
	for(k=0;k<FLTNUM;k++)
	{
		if(fabs(distx-distn[k])<0.01f)
		{
			orit[k]++;
		}
	}
}

void storeRefAndEncFrames()
{
	int i,j;

	if(img->number<REFNUM)
	{
		int index=img->number;
		imgpel**tmpY=refPicsData[index][0];
		imgpel**tmpU=refPicsData[index][1];
		imgpel**tmpV=refPicsData[index][2];

		imgpel**tmpY1=srcPicsData[index][0];
		imgpel**tmpU1=srcPicsData[index][1];
		imgpel**tmpV1=srcPicsData[index][2];
		for(i=index;i>=1;i--)
		{
			refPicsData[i][0]=refPicsData[i-1][0];
			refPicsData[i][1]=refPicsData[i-1][1];
			refPicsData[i][2]=refPicsData[i-1][2];
		}
		refPicsData[0][0]=tmpY;
		refPicsData[0][1]=tmpU;
		refPicsData[0][2]=tmpV;
		for(j=0;j<img->height;j++)
		{
			for(i=0;i<img->width;i++)
			{
				refPicsData[0][0][j][i]=img_org_frm[0][j][i];
			}
		}

		for(j=0;j<img->height_cr;j++)
		{
			for(i=0;i<img->width_cr;i++)
			{
				refPicsData[0][1][j][i]=img_org_frm[1][j][i];
				refPicsData[0][2][j][i]=img_org_frm[2][j][i];
			}
		}

		for(i=index;i>=1;i--)
		{
			srcPicsData[i][0]=srcPicsData[i-1][0];
			srcPicsData[i][1]=srcPicsData[i-1][1];
			srcPicsData[i][2]=srcPicsData[i-1][2];
		}
		srcPicsData[0][0]=tmpY1;
		srcPicsData[0][1]=tmpU1;
		srcPicsData[0][2]=tmpV1;
		for(j=0;j<img->height;j++)
		{
			for(i=0;i<img->width;i++)
			{
				srcPicsData[0][0][j][i]=enc_picture->imgY[j][i];
			}
		}

		for(j=0;j<img->height_cr;j++)
		{
			for(i=0;i<img->width_cr;i++)
			{
				srcPicsData[0][1][j][i]=enc_picture->imgUV[0][j][i];
				srcPicsData[0][2][j][i]=enc_picture->imgUV[1][j][i];
			}
		}
	}
	else //if(img->number>=5)
	{
		imgpel**tmpY=refPicsData[REFNUM-1][0];
		imgpel**tmpU=refPicsData[REFNUM-1][1];
		imgpel**tmpV=refPicsData[REFNUM-1][2];

		imgpel**tmpY1=srcPicsData[REFNUM-1][0];
		imgpel**tmpU1=srcPicsData[REFNUM-1][1];
		imgpel**tmpV1=srcPicsData[REFNUM-1][2];
		for(i=REFNUM-1;i>=1;i--)
		{
			refPicsData[i][0]=refPicsData[i-1][0];
			refPicsData[i][1]=refPicsData[i-1][1];
			refPicsData[i][2]=refPicsData[i-1][2];
		}
		refPicsData[0][0]=tmpY;
		refPicsData[0][1]=tmpU;
		refPicsData[0][2]=tmpV;
		for(j=0;j<img->height;j++)
		{
			for(i=0;i<img->width;i++)
			{
				refPicsData[0][0][j][i]=img_org_frm[0][j][i];
			}
		}

		for(j=0;j<img->height_cr;j++)
		{
			for(i=0;i<img->width_cr;i++)
			{
				refPicsData[0][1][j][i]=img_org_frm[1][j][i];
				refPicsData[0][2][j][i]=img_org_frm[2][j][i];
			}
		}

		for(i=REFNUM-1;i>=1;i--)
		{
			srcPicsData[i][0]=srcPicsData[i-1][0];
			srcPicsData[i][1]=srcPicsData[i-1][1];
			srcPicsData[i][2]=srcPicsData[i-1][2];
		}
		srcPicsData[0][0]=tmpY1;
		srcPicsData[0][1]=tmpU1;
		srcPicsData[0][2]=tmpV1;
		for(j=0;j<img->height;j++)
		{
			for(i=0;i<img->width;i++)
			{
				srcPicsData[0][0][j][i]=enc_picture->imgY[j][i];
			}
		}

		for(j=0;j<img->height_cr;j++)
		{
			for(i=0;i<img->width_cr;i++)
			{
				srcPicsData[0][1][j][i]=enc_picture->imgUV[0][j][i];
				srcPicsData[0][2][j][i]=enc_picture->imgUV[1][j][i];
			}
		}
	}
}

float compute_SSIM (imgpel **imgRef, imgpel **imgEnc,int yRef,int yEnc, int xRef, int xEnc, int height, int width, int wint, int comp)
{
	static const float K1 = 0.01f, K2 = 0.03f;
	static float max_pix_value_sqd;
	float C1, C2;
	float win_pixels = (float) (wint * wint);
	float wgt;
	float mb_ssim,varOrg, varEnc, covOrgEnc;
	float imeanOrg, imeanEnc, ivarOrg, ivarEnc, icovOrgEnc;
	float cur_distortion = 0.0;
	int i, j, n, m, win_cnt = 0;
	int overlapSize = params->SSIMOverlapSize;
	imgpel**refImg,**encImg;
	imgpel pixRef,pixEnc;

	refImg=imgRef;
	encImg=imgEnc;

	max_pix_value_sqd = (float) (img->max_imgpel_value_comp[comp] * img->max_imgpel_value_comp[comp]);
	C1 = K1 * K1 * max_pix_value_sqd;
	C2 = K2 * K2 * max_pix_value_sqd;

	for (j = 0; j <= height - wint; j += overlapSize)
	{
		for (i = 0; i <= width - wint; i += overlapSize)
		{
			imeanOrg = 0;
			imeanEnc = 0; 
			ivarOrg  = 0;
			ivarEnc  = 0;
			icovOrgEnc = 0;

			for ( n = j;n < j + wint;n ++)
			{
				for (m = i;m < i + wint;m ++)
				{
#if _SSIM_WEIGHTED_
					if(wint==4)
					wgt=gauss4[n-j][m-i];
					else
					wgt=gauss8[n-j][m-i];
#else
	                wgt=1.0f/win_pixels;
#endif

						pixRef=refImg[yRef+n][xRef+m];
						pixEnc=encImg[yEnc+n][xEnc+m];
						imeanOrg   += wgt*pixRef;
						imeanEnc   += wgt*pixEnc;
						ivarOrg    += wgt*pixRef * pixRef;
						ivarEnc    += wgt*pixEnc * pixEnc;
						icovOrgEnc += wgt*pixRef * pixEnc;
				}
			}

			varOrg    = fabs(ivarOrg - imeanOrg*imeanOrg);/// win_pixels;
			varEnc    = fabs(ivarEnc - imeanEnc*imeanEnc) ;/// win_pixels;
			covOrgEnc = fabs(icovOrgEnc - imeanOrg*imeanEnc) ;/// win_pixels;


			mb_ssim  = (float) ((2.0 * imeanOrg *imeanEnc + C1) * (2.0 * covOrgEnc + C2));
			mb_ssim /= (float) (imeanOrg *imeanOrg +imeanEnc *imeanEnc + C1) * (varOrg + varEnc + C2);

			cur_distortion += mb_ssim;
			win_cnt++;
		}
	}

	cur_distortion /= (float) win_cnt;

	if (cur_distortion >= 1.0 && cur_distortion < 1.01) // avoid float accuracy problem at very low QP(e.g.2)
		cur_distortion = 1.0;


	return cur_distortion;
}
double distortionSSIM(Macroblock *currMB) 
{
  float distortionY = 0;
  float distortionCr[2] = {0};

  // LUMA
 distortionY = 1.0f-compute_SSIM(pCurImg,enc_picture->p_curr_img,img->opix_y,img->pix_y, img->opix_x, img->pix_x, MB_BLOCK_SIZE, MB_BLOCK_SIZE,8,0);

  //if (img->yuv_format != YUV400 && )
  if ((img->yuv_format != YUV400) && !IS_INDEPENDENT(params))
  {
    // CHROMA
    distortionCr[0] = 1.0f-compute_SSIM(pImgOrg[1], enc_picture->imgUV[0],img->opix_c_y,img->pix_c_y, img->opix_c_x, img->pix_c_x, img->mb_cr_size_y, img->mb_cr_size_x,8,1);
    distortionCr[1] = 1.0f-compute_SSIM(pImgOrg[2], enc_picture->imgUV[1],img->opix_c_y,img->pix_c_y, img->opix_c_x, img->pix_c_x, img->mb_cr_size_y, img->mb_cr_size_x,8,2);
  }

  return ( distortionY * params->WeightY + distortionCr[0] * params->WeightCb + distortionCr[1] * params->WeightCr );
}

//extern float **att_wgt;
float compute_stVSSIM(imgpel**imgRef,imgpel**imgEnc,int yRef,int yEnc, int xRef, int xEnc, 
					  int height, int width,int wint,int gama,int comp,
					  float *ssimValue,float*ssim3dValue,float*stvssimValue)
{
	static const float K1 = 0.01f, K2 = 0.03f;
	static float max_pix_value_sqd;
	float C1, C2;
	int frameused=min(gama,REFNUM);
	float win_pixels = (float) (wint * wint*frameused);

	float mb_ssim,varOrg, varEnc, covOrgEnc;
	float imeanOrg, imeanEnc, ivarOrg, ivarEnc, icovOrgEnc;
	float ssim3d[FLTNUM],ssim3dx,ssimx,stvssimx;
	float tmpssim3d,tmpssim,tmpstvssim;
	int i, j, n, m,o,k, win_cnt = 0;
	int overlapSize = params->SSIMOverlapSize;

	imgpel **refImg;
	imgpel **encImg;
//	int yRefn,xRefm,yEncn,xEncm;
	imgpel pixRef,pixEnc;
	int uv=comp>0?2:1;

	float (*filter[FLTNUM])(int,int,int,float,float);
	short orit[FLTNUM]={0};

	short tmporit,inx;

	float wa,wb,wgta[4],wgtb[4],wgt;

	wa=SSIM3D_WGT;
	wb=1.0f-wa;

	if(wint==4)
	{	
		filter[0]=hFilter_4;
		filter[1]=rFilter_4;
		filter[2]=vFilter_4;
		filter[3]=lFilter_4;
	wgta[0]=wgta[2]=wgta[1]=wgta[3]=wa/(wint*(frameused));
	wgtb[0]=wgtb[2]=wgtb[1]=wgtb[3]=wb/((wint*wint-wint)*(frameused));
	}
	else
	{
		filter[0]=hFilter;
		filter[1]=rFilter;
		filter[2]=vFilter;
		filter[3]=lFilter;

	wgta[0]=wgta[2]=wa/(3*wint*(frameused));
	wgta[1]=wgta[3]=wa/((3*wint-2)*(frameused));
	wgtb[0]=wgtb[2]=wb/((wint*wint-3*wint)*(frameused));
	wgtb[1]=wgtb[3]=wb/((wint*wint-3*wint+2)*(frameused));
	}


	max_pix_value_sqd = (float) (img->max_imgpel_value_comp[comp] * img->max_imgpel_value_comp[comp]);
	C1 = K1 * K1 * max_pix_value_sqd;
	C2 = K2 * K2 * max_pix_value_sqd;

	ssim3dx=0.0f;
	ssimx=0.0f;
	stvssimx=0.0f;

	for (j = 0; j <= height - wint; j += overlapSize)
	{
		for (i = 0; i <= width - wint; i += overlapSize)
		{
			//compute the ssim3d value
			for(k=0;k<FLTNUM;k++)
			{
//					FILE *pa=fopen("a.txt","w");
//					FILE *pb=fopen("b.txt","w");
				imeanOrg = 0;
				imeanEnc = 0; 
				ivarOrg  = 0;
				ivarEnc  = 0;
				icovOrgEnc = 0;
				for(o=0;o<frameused;o++)
				{	
					if(o!=frameused-1)
					{
						refImg=refPicsData[o][comp];
						encImg=srcPicsData[o][comp];
					}
					else
					{
						refImg=imgRef;
						encImg=imgEnc;
					}
					for ( n = j;n < j + wint;n ++)
					{
						for (m = i;m < i + wint;m ++)
						{
							wgt=filter[k](wint,n-j,m-i,wgta[k],wgtb[k]);

						pixRef=refImg[yRef+n][xRef+m];
						pixEnc=encImg[yEnc+n][xEnc+m];
//							fprintf(pa,"%f ",wgt);
//							fprintf(pb,"%d ",pixRef);
						imeanOrg   += wgt*pixRef;
						imeanEnc   += wgt*pixEnc;
						ivarOrg    += wgt*pixRef * pixRef;
						ivarEnc    += wgt*pixEnc * pixEnc;
						icovOrgEnc += wgt*pixRef * pixEnc;

						}
//						fprintf(pa,"\n");
//						fprintf(pb,"\n");
					}
//						fprintf(pa,"\n");
//						fprintf(pb,"\n");
				}
//					fclose(pa);
//					fclose(pb);

				varOrg    = fabs(ivarOrg - imeanOrg*imeanOrg);// / win_pixels_bias;
				assert(varOrg>=0);
				varEnc    = fabs(ivarEnc - imeanEnc*imeanEnc);// / win_pixels_bias;
				covOrgEnc = fabs(icovOrgEnc - imeanOrg*imeanEnc);// / win_pixels_bias;

				mb_ssim  = (float) ((2.0 *imeanOrg *imeanEnc + C1) * (2.0 * covOrgEnc + C2));
				mb_ssim /= (float) (imeanOrg *imeanOrg +imeanEnc *imeanEnc + C1) * (varOrg + varEnc + C2);

//				assert(mb_ssim>0);
				// stvssimx += mb_ssim[k];
				ssim3d[k] = mb_ssim;
				if (ssim3d[k]>= 1.0 && ssim3d[k]< 1.01) // avoid float accuracy problem at very low QP(e.g.2)
					ssim3d[k]= 1.0;
			}
			//统计方向，求出最佳的方向

			orit[0]=orit[1]=orit[2]=orit[3]=0;
			for ( n = j;n < j + wint;n ++)
			{
				for (m = i;m < i + wint;m ++)
				{
//					yRefn=(yRef+n)*uv;
//					xRefm=(xRef+m)*uv;
					calOrit(pic_directions2[(yRef+n)*uv][(xRef+m)*uv],orit);
				}
			}

			tmporit=0;
			inx=0;
			for(n=0;n<4;n++)
			{
				if(orit[n]>tmporit)
				{
					tmporit=orit[n];
					inx=n;
				}
			}
			for(n=0;n<4;++n)
			{
				if((tmporit-orit[n])<10&&inx!=n)
					break;
			}

			if(n==4)
			{
				ssim3dx+=ssim3d[inx];
				tmpssim3d=ssim3d[inx];
			}
			else
			{
				ssim3dx+=(ssim3d[inx]+ssim3d[n])/2;
				tmpssim3d=(ssim3d[inx]+ssim3d[n])/2;
			}
			//end of computing ssim3d
			//compute the ssim value
			{
				imeanOrg = 0;
				imeanEnc = 0; 
				ivarOrg  = 0;
				ivarEnc  = 0;
				icovOrgEnc = 0;

				refImg=imgRef;
				encImg=imgEnc;
				
				wgt=1.0f/(wint*wint);

				for ( n = j;n < j + wint;n ++)
				{
					for (m = i;m < i + wint;m ++)
					{
#if _SSIM_WEIGHTED_
						if(wint==4)
							wgt=gauss4[n-j][m-i];
						else
							wgt=gauss8[n-j][m-i];
#endif
//						yRefn=yRef+n;
//						yEncn=yEnc+n;
//						xRefm=xRef+m;
//						xEncm=xEnc+m;
						pixRef=refImg[yRef+n][xRef+m];
						pixEnc=encImg[yEnc+n][xEnc+m];
//#if _ATT_M1_
//						wgt=att_wgt[yRef+n][xRef+m];
//#endif
						imeanOrg   += wgt*pixRef;
						imeanEnc   += wgt*pixEnc;
						ivarOrg    += wgt*pixRef * pixRef;
						ivarEnc    += wgt*pixEnc * pixEnc;
						icovOrgEnc += wgt*pixRef * pixEnc;
					}
				}

				varOrg    = fabs(ivarOrg - imeanOrg*imeanOrg);/// win_pixels;
				varEnc    = fabs(ivarEnc - imeanEnc*imeanEnc);/// win_pixels;
				covOrgEnc = fabs(icovOrgEnc - imeanOrg*imeanEnc);/// win_pixels;

				mb_ssim  = (float) ((2.0 *imeanOrg *imeanEnc + C1) * (2.0 * covOrgEnc + C2));
				mb_ssim /= (float) (imeanOrg *imeanOrg +imeanEnc *imeanEnc + C1) * (varOrg + varEnc + C2);

				ssimx += mb_ssim;
				tmpssim=mb_ssim;

			}
			//end of computing ssim
			//compute the stvssim
			stvssimx+=tmpssim*tmpssim3d;
			win_cnt++;
		}
	}
	assert(ssim3dx>0);
	assert(ssimx>0);

	// stvssimx /= (float) win_cnt;
	ssim3dx /= (float) win_cnt;
	ssimx /= (float) win_cnt;
	stvssimx /= (float) win_cnt;

	tmpstvssim=ssimx*ssim3dx;//(ssim3d[0]+ssim3d[1]+ssim3d[2]+ssim3d[3])/4;
	if (stvssimx >= 1.0 && stvssimx < 1.01) // avoid float accuracy problem at very low QP(e.g.2)
		stvssimx = 1.0;
	*ssimValue=ssimx;
	*ssim3dValue=ssim3dx;
	*stvssimValue=stvssimx;
	return tmpstvssim;

}
double distortionstVSSIM(Macroblock *currMB) 
{
	float distortionY = 0;
	float distortionCr[2] = {0};
	float ssimValue,ssim3dValue,stvssimValue;

	// LUMA
	compute_stVSSIM(pCurImg, enc_picture->p_curr_img,img->opix_y,img->pix_y,img->opix_x, img->pix_x, 
		MB_BLOCK_SIZE, MB_BLOCK_SIZE,8,img->number+1,0,&ssimValue,&ssim3dValue,&stvssimValue);
	distortionY = 1.0f-stvssimValue;
	//if (img->yuv_format != YUV400 && )
	if ((img->yuv_format != YUV400) && !IS_INDEPENDENT(params))
	{
		// CHROMA
		compute_stVSSIM(pImgOrg[1], enc_picture->imgUV[0],img->opix_c_y,img->pix_c_y,img->opix_c_x, img->pix_c_x,
			img->mb_cr_size_y, img->mb_cr_size_x,8,img->number+1,1,&ssimValue,&ssim3dValue,&stvssimValue);
		distortionCr[0] = 1.0f-stvssimValue; 

		compute_stVSSIM(pImgOrg[2], enc_picture->imgUV[1],img->opix_c_y,img->pix_c_y,img->opix_c_x, img->pix_c_x,
			img->mb_cr_size_y, img->mb_cr_size_x,8,img->number+1,1,&ssimValue,&ssim3dValue,&stvssimValue);
		distortionCr[1] = 1.0f-stvssimValue; 
	}

	return ( distortionY * params->WeightY + distortionCr[0] * params->WeightCb + distortionCr[1] * params->WeightCr );
}



float compute_stVSSIM2(imgpel**imgRef,imgpel**imgEnc,int yRef,int yEnc, int xRef, int xEnc, 
					  int height, int width,int wint,int gama,int comp,
					  float *ssimValue,float*ssim3dValue,float*stvssimValue)
{
	static const float K1 = 0.01f, K2 = 0.03f;
	static float max_pix_value_sqd;
	float C1, C2;
	int frameused=min(gama,REFNUM);
	float win_pixels = (float) (wint * wint*frameused);

	float mb_ssim,varOrg, varEnc, covOrgEnc;
	float imeanOrg, imeanEnc, ivarOrg, ivarEnc, icovOrgEnc;
	float ssim3d[FLTNUM],ssim3dx,ssimx,stvssimx;
	float tmpssim3d,tmpssim,tmpstvssim;
	int i, j, n, m,o,k, win_cnt = 0;
	int overlapSize = params->SSIMOverlapSize;

	imgpel **refImg;
	imgpel **encImg;
	imgpel pixRef,pixEnc;
	int uv=comp>0?2:1;

	float (*filter[FLTNUM])(int,int,int,float,float)={hFilter,rFilter,vFilter,lFilter};
	short orit[FLTNUM]={0};

	short tmporit,inx;

	float wa,wb,wgta[4],wgtb[4],wgt;

	wa=SSIM3D_WGT;
	wb=1.0f-wa;

	wgta[0]=wgta[2]=wa/(3*wint*(frameused));
	wgta[1]=wgta[3]=wa/((3*wint-2)*(frameused));
	wgtb[0]=wgtb[2]=wb/((wint*wint-3*wint)*(frameused));
	wgtb[1]=wgtb[3]=wb/((wint*wint-3*wint+2)*(frameused));


	max_pix_value_sqd = (float) (img->max_imgpel_value_comp[comp] * img->max_imgpel_value_comp[comp]);
	C1 = K1 * K1 * max_pix_value_sqd;
	C2 = K2 * K2 * max_pix_value_sqd;

	ssim3dx=0.0f;
	ssimx=0.0f;
	stvssimx=0.0f;

	for (j = 0; j <= height - wint; j += overlapSize)
	{
		for (i = 0; i <= width - wint; i += overlapSize)
		{
			//compute the ssim3d value
			for(k=0;k<FLTNUM;k++)
			{
				imeanOrg = 0;
				imeanEnc = 0; 
				ivarOrg  = 0;
				ivarEnc  = 0;
				icovOrgEnc = 0;
				for(o=0;o<frameused;o++)
				{	
					if(o!=frameused-1)
					{
						refImg=refPicsData[o][comp];
						encImg=srcPicsData[o][comp];
					}
					else
					{
						refImg=imgRef;
						encImg=imgEnc;
					}
					for ( n = j;n < j + wint;n ++)
					{
						for (m = i;m < i + wint;m ++)
						{
							wgt=filter[k](wint,n-j,m-i,wgta[k],wgtb[k]);

						pixRef=refImg[yRef+n][xRef+m];
						pixEnc=encImg[yEnc+n][xEnc+m];
						imeanOrg   += wgt*pixRef;
						imeanEnc   += wgt*pixEnc;
						ivarOrg    += wgt*pixRef * pixRef;
						ivarEnc    += wgt*pixEnc * pixEnc;
						icovOrgEnc += wgt*pixRef * pixEnc;

						}
					}
				}

				varOrg    = fabs(ivarOrg - imeanOrg*imeanOrg);// / win_pixels_bias;
				varEnc    = fabs(ivarEnc - imeanEnc*imeanEnc);// / win_pixels_bias;
				covOrgEnc = fabs(icovOrgEnc - imeanOrg*imeanEnc);// / win_pixels_bias;

				mb_ssim  = (float) ((2.0 *imeanOrg *imeanEnc + C1) * (2.0 * covOrgEnc + C2));
				mb_ssim /= (float) (imeanOrg *imeanOrg +imeanEnc *imeanEnc + C1) * (varOrg + varEnc + C2);

				// stvssimx += mb_ssim[k];
				ssim3d[k] = mb_ssim;
				if (ssim3d[k]>= 1.0 && ssim3d[k]< 1.01) // avoid float accuracy problem at very low QP(e.g.2)
					ssim3d[k]= 1.0;
			}
			//统计方向，求出最佳的方向

			orit[0]=orit[1]=orit[2]=orit[3]=0;
			for ( n = j;n < j + wint;n ++)
			{
				for (m = i;m < i + wint;m ++)
				{
					calOrit(pic_directions2[(yRef+n)*uv][(xRef+m)*uv],orit);
				}
			}

			tmporit=0;
			inx=0;
			for(n=0;n<4;n++)
			{
				if(orit[n]>tmporit)
				{
					tmporit=orit[n];
					inx=n;
				}
			}
			for(n=0;n<4;++n)
			{
				if((tmporit-orit[n])<10&&inx!=n)
					break;
			}

			if(n==4)
			{
				ssim3dx+=ssim3d[inx];
				tmpssim3d=ssim3d[inx];
			}
			else
			{
				ssim3dx+=(ssim3d[inx]+ssim3d[n])/2;
				tmpssim3d=(ssim3d[inx]+ssim3d[n])/2;
			}
			//end of computing ssim3d
			//compute the ssim value
			{
				imeanOrg = 0;
				imeanEnc = 0; 
				ivarOrg  = 0;
				ivarEnc  = 0;
				icovOrgEnc = 0;

				refImg=imgRef;
				encImg=imgEnc;
				
				wgt=1.0f/(wint*wint);

				for ( n = j;n < j + wint;n ++)
				{
					for (m = i;m < i + wint;m ++)
					{
#if _SSIM_WEIGHTED_
						if(wint==4)
							wgt=gauss4[n-j][m-i];
						else
							wgt=gauss8[n-j][m-i];
#endif
						pixRef=refImg[yRef+n][xRef+m];
						pixEnc=encImg[yEnc+n][xEnc+m];
//#if _ATT_M1_
//						wgt=att_wgt[yRef+n][xRef+m];
//#endif
						imeanOrg   += wgt*pixRef;
						imeanEnc   += wgt*pixEnc;
						ivarOrg    += wgt*pixRef * pixRef;
						ivarEnc    += wgt*pixEnc * pixEnc;
						icovOrgEnc += wgt*pixRef * pixEnc;
					}
				}

				varOrg    = fabs(ivarOrg - imeanOrg*imeanOrg);/// win_pixels;
				varEnc    = fabs(ivarEnc - imeanEnc*imeanEnc);/// win_pixels;
				covOrgEnc = fabs(icovOrgEnc - imeanOrg*imeanEnc);/// win_pixels;

				mb_ssim  = (float) ((2.0 *imeanOrg *imeanEnc + C1) * (2.0 * covOrgEnc + C2));
				mb_ssim /= (float) (imeanOrg *imeanOrg +imeanEnc *imeanEnc + C1) * (varOrg + varEnc + C2);

				ssimx += mb_ssim;
				tmpssim=mb_ssim;

			}
			//end of computing ssim
			//compute the stvssim
			stvssimx+=(ALPHA*tmpssim3d+(1-ALPHA)*tmpssim);
			win_cnt++;
		}
	}
	assert(ssim3dx>0);
	assert(ssimx>0);

	// stvssimx /= (float) win_cnt;
	ssim3dx /= (float) win_cnt;
	ssimx /= (float) win_cnt;
	stvssimx /= (float) win_cnt;

	tmpstvssim=ALPHA*ssim3dx+(1-ALPHA)*ssimx;//(ssim3d[0]+ssim3d[1]+ssim3d[2]+ssim3d[3])/4;
	if (stvssimx >= 1.0 && stvssimx < 1.01) // avoid float accuracy problem at very low QP(e.g.2)
		stvssimx = 1.0;
	*ssimValue=ssimx;
	*ssim3dValue=ssim3dx;
	*stvssimValue=stvssimx;
	return tmpstvssim;

}
double distortionstVSSIM2(Macroblock *currMB) 
{
	float distortionY = 0;
	float distortionCr[2] = {0};
	float ssimValue,ssim3dValue,stvssimValue;

	// LUMA
	compute_stVSSIM2(pCurImg, enc_picture->p_curr_img,img->opix_y,img->pix_y,img->opix_x, img->pix_x, 
		MB_BLOCK_SIZE, MB_BLOCK_SIZE,8,img->number+1,0,&ssimValue,&ssim3dValue,&stvssimValue);
	distortionY = 1.0f-stvssimValue;
	//if (img->yuv_format != YUV400 && )
	if ((img->yuv_format != YUV400) && !IS_INDEPENDENT(params))
	{
		// CHROMA
		compute_stVSSIM2(pImgOrg[1], enc_picture->imgUV[0],img->opix_c_y,img->pix_c_y,img->opix_c_x, img->pix_c_x,
			img->mb_cr_size_y, img->mb_cr_size_x,8,img->number+1,1,&ssimValue,&ssim3dValue,&stvssimValue);
		distortionCr[0] = 1.0f-stvssimValue; 

		compute_stVSSIM2(pImgOrg[2], enc_picture->imgUV[1],img->opix_c_y,img->pix_c_y,img->opix_c_x, img->pix_c_x,
			img->mb_cr_size_y, img->mb_cr_size_x,8,img->number+1,1,&ssimValue,&ssim3dValue,&stvssimValue);
		distortionCr[1] = 1.0f-stvssimValue; 
	}

	return ( distortionY * params->WeightY + distortionCr[0] * params->WeightCb + distortionCr[1] * params->WeightCr );
}

float compute_SSIM3D(imgpel**imgRef,imgpel**imgEnc,int yRef,int yEnc, int xRef, int xEnc, 
					  int height, int width,int wint,int gama,int comp)
{
	static const float K1 = 0.01f, K2 = 0.03f;
	static float max_pix_value_sqd;
	float C1, C2;
	int frameused=min(gama,REFNUM);
	float win_pixels = (float) (wint * wint*frameused);

	float mb_ssim,varOrg, varEnc, covOrgEnc;
	float imeanOrg, imeanEnc, ivarOrg, ivarEnc, icovOrgEnc;
	float ssim3d[FLTNUM],ssim3dx;
	int i, j, n, m,o,k, win_cnt = 0;
	int overlapSize = params->SSIMOverlapSize;

	imgpel **refImg;
	imgpel **encImg;
	imgpel pixRef,pixEnc;
	int uv=comp>0?2:1;

	float (*filter[FLTNUM])(int,int,int,float,float)={hFilter,rFilter,vFilter,lFilter};
	short orit[FLTNUM]={0};

	short tmporit,inx;

	float wa,wb,wgta[4],wgtb[4],wgt;

	wa=SSIM3D_WGT;
	wb=1.0f-wa;

	wgta[0]=wgta[2]=wa/(3*wint*(frameused));
	wgta[1]=wgta[3]=wa/((3*wint-2)*(frameused));
	wgtb[0]=wgtb[2]=wb/((wint*wint-3*wint)*(frameused));
	wgtb[1]=wgtb[3]=wb/((wint*wint-3*wint+2)*(frameused));


	max_pix_value_sqd = (float) (img->max_imgpel_value_comp[comp] * img->max_imgpel_value_comp[comp]);
	C1 = K1 * K1 * max_pix_value_sqd;
	C2 = K2 * K2 * max_pix_value_sqd;

	ssim3dx=0.0f;

	for (j = 0; j <= height - wint; j += overlapSize)
	{
		for (i = 0; i <= width - wint; i += overlapSize)
		{
			//compute the ssim3d value
			for(k=0;k<FLTNUM;k++)
			{
				imeanOrg = 0;
				imeanEnc = 0; 
				ivarOrg  = 0;
				ivarEnc  = 0;
				icovOrgEnc = 0;
				for(o=0;o<frameused;o++)
				{	
					if(o!=frameused-1)
					{
						refImg=refPicsData[o][comp];
						encImg=srcPicsData[o][comp];
					}
					else
					{
						refImg=imgRef;
						encImg=imgEnc;
					}
					for ( n = j;n < j + wint;n ++)
					{
						for (m = i;m < i + wint;m ++)
						{
							wgt=filter[k](wint,n-j,m-i,wgta[k],wgtb[k]);

						pixRef=refImg[yRef+n][xRef+m];
						pixEnc=encImg[yEnc+n][xEnc+m];
						imeanOrg   += wgt*pixRef;
						imeanEnc   += wgt*pixEnc;
						ivarOrg    += wgt*pixRef * pixRef;
						ivarEnc    += wgt*pixEnc * pixEnc;
						icovOrgEnc += wgt*pixRef * pixEnc;

						}
					}
				}

				varOrg    = fabs(ivarOrg -imeanOrg*imeanOrg);// / win_pixels_bias;
				varEnc    = fabs(ivarEnc -imeanEnc*imeanEnc);// / win_pixels_bias;
				covOrgEnc = fabs(icovOrgEnc -imeanOrg*imeanEnc);// / win_pixels_bias;

				mb_ssim  = (float) ((2.0 *imeanOrg *imeanEnc + C1) * (2.0 * covOrgEnc + C2));
				mb_ssim /= (float) (imeanOrg * imeanOrg +imeanEnc *imeanEnc + C1) * (varOrg + varEnc + C2);

				// stvssimx += mb_ssim[k];
				ssim3d[k] = mb_ssim;
				if (ssim3d[k]>= 1.0 && ssim3d[k]< 1.01) // avoid float accuracy problem at very low QP(e.g.2)
					ssim3d[k]= 1.0;
			}
			//统计方向，求出最佳的方向

			orit[0]=orit[1]=orit[2]=orit[3]=0;
			for ( n = j;n < j + wint;n ++)
			{
				for (m = i;m < i + wint;m ++)
				{
					calOrit(pic_directions2[(yRef+n)*uv][(xRef+m)*uv],orit);
				}
			}

			tmporit=0;
			inx=0;
			for(n=0;n<4;n++)
			{
				if(orit[n]>tmporit)
				{
					tmporit=orit[n];
					inx=n;
				}
			}
			for(n=0;n<4;++n)
			{
				if((tmporit-orit[n])<10&&inx!=n)
					break;
			}

			if(n==4)
			{
				ssim3dx+=ssim3d[inx];
			}
			else
			{
				ssim3dx+=(ssim3d[inx]+ssim3d[n])/2;
			}
			//end of computing ssim3d
			win_cnt++;
		}
	}
	assert(ssim3dx>0);
	// stvssimx /= (float) win_cnt;
	ssim3dx /= (float) win_cnt;

	if (ssim3dx >= 1.0 && ssim3dx < 1.01) // avoid float accuracy problem at very low QP(e.g.2)
		ssim3dx = 1.0;
	return ssim3dx;

}
double distortionSSIM3D(Macroblock *currMB) 
{
	float distortionY = 0;
	float distortionCr[2] = {0};
	float ssim3dValue;

	// LUMA
	ssim3dValue=compute_SSIM3D(pCurImg, enc_picture->p_curr_img,img->opix_y,img->pix_y,img->opix_x, img->pix_x, 
		MB_BLOCK_SIZE, MB_BLOCK_SIZE,8,img->number+1,0);
	distortionY = 1.0f-ssim3dValue;
	//if (img->yuv_format != YUV400 && )
	if ((img->yuv_format != YUV400) && !IS_INDEPENDENT(params))
	{
		// CHROMA
		ssim3dValue=compute_SSIM3D(pImgOrg[1], enc_picture->imgUV[0],img->opix_c_y,img->pix_c_y,img->opix_c_x, img->pix_c_x,
			img->mb_cr_size_y, img->mb_cr_size_x,8,img->number+1,1);
		distortionCr[0] = 1.0f-ssim3dValue; 

		ssim3dValue=compute_SSIM3D(pImgOrg[2], enc_picture->imgUV[1],img->opix_c_y,img->pix_c_y,img->opix_c_x, img->pix_c_x,
			img->mb_cr_size_y, img->mb_cr_size_x,8,img->number+1,1);
		distortionCr[1] = 1.0f-ssim3dValue; 
	}

	return ( distortionY * params->WeightY + distortionCr[0] * params->WeightCb + distortionCr[1] * params->WeightCr );
}


void getMV_macroblock(short****mb_mv,int mode,int block,int bestref)
{
	int i,j;
	//	int block_x, block_y;
	short ***curr_mv = NULL;
	int64 curr_ref_idx = 0;

	int start_x = 0, start_y = 0, end_x = BLOCK_MULTIPLE, end_y = BLOCK_MULTIPLE; 
	switch (mode)
	{
	case 1:
		start_x = 0;
		start_y = 0;
		end_x   = BLOCK_MULTIPLE;
		end_y   = BLOCK_MULTIPLE;
		break;
	case 2:
		start_x = 0;
		start_y = block * 2;
		end_x   = BLOCK_MULTIPLE;
		end_y   = (block + 1) * 2;
		break;
	case 3:
		start_x = block * 2;
		start_y = 0;
		end_x   = (block + 1) * 2;
		end_y   = BLOCK_MULTIPLE;
		break;
	case 4:
	case 5:
	case 6:
	case 7:
		start_x=(block%2)<<1;
		end_x=start_x+2;
		start_y=(block&0x2);
		end_y=start_y+2;
		break;
	default:
		break;
	}
	curr_mv = img->all_mv[0][bestref][mode];	

	for (j = start_y; j < end_y; j++)
	{
		//	 block_y = img->block_y + j;
		for (i = start_x; i < end_x; i++)
		{
			mb_mv[mode][j][i][0]=curr_mv[j][i][0];
			mb_mv[mode][j][i][1]=curr_mv[j][i][1];
		}
	}

}
short  getOrientation(short x,short y)
{
	float ddd,dge;
	float distx;
	float dist[ORIENTS];
	short i,index;
	if(x==0)
		ddd=PI/2;
	else if(y==0)
		ddd=0;
	else
	{
		dge=(float)y/x;
		ddd=(float)atan(dge);
		if(ddd<0)
			ddd+=PI;
	}
	//	printf("ddd=%7.4f\n",ddd*180/PI);
	distx=10000.0f;
	for(i=0;i<ORIENTS;i++)
	{
		dist[i]=(float)fabs(ddd-orients[i]);
		if(dist[i]<distx)
		{
			distx=dist[i];
			index=i;
		}
	}
	return index;
}
float chooseOrient(short*orit)
{
	int i;
	short tmp;
	short index;
	short orit2[ORIENTS2]={0};
	for(i=1;i<MODENUM;++i)
	{
		orit2[orit[i]/2]++;
	}
	tmp=0;
	index=0;
	for(i=0;i<ORIENTS2;++i)
	{
		if(orit2[i]>tmp)
		{
			tmp=orit2[i];
			index=i;
		}
	}
	return orients2[index];
}
void getDirection_macroblock()
{
	int i,j,k;
	short orientation[MODENUM];
	short *mv;

	for(j=0;j<4;++j)
	{
		for(i=0;i<4;++i)
		{
			for(k=1;k<MODENUM;k++)
			{
				mv=mb_mv[k][j][i];
				orientation[k]=getOrientation(mv[0],mv[1]);
			}
			mb_directions[j][i]=chooseOrient(orientation);

		}
	}

	  for(j=0;j<16;++j)
	  {
		  for(i=0;i<16;++i)
		  {
			  pic_directions2[img->pix_y+j][img->pix_x+i]=mb_directions[j/4][i/4];
//			  mb_directions2[j][i]=mb_directions[j/4][i/4];
		  }
	  }

}

float compute_PSNR(imgpel**imgRef,imgpel**imgEnc,int yRef,int yEnc, int xRef, int xEnc, int height, int width)
{
	float sse_distortion;
	float dpsnr;
	int max_value=(1<<8)-1;
	int max_sample_sq=max_value*max_value;
	int samples=height*width;
	sse_distortion=(float)compute_SSE(&imgRef[yRef],&imgEnc[yEnc],xRef,xEnc,height,width);

	dpsnr=(float) (10.0 * log10(max_sample_sq * (double) ((double) samples / (sse_distortion == 0.0 ? 1.0 : sse_distortion))));

	return dpsnr;
}
float compute_VAR(imgpel**imgRef,int yRef,int xRef,int height, int width)
{
	int pixNum=height*width;
	int i,j;
	int yn,xm;
	imgpel  pixel;
	int imean,ivar;
	float mean,var;
	imean=0;
	ivar=0;

	for (j = 0; j < height; j++)
	{
		yn=j+yRef;
		for(i=0;i<width;i++)
		{
			xm=i+xRef;
			pixel=imgRef[yn][xm];
			imean+=pixel;
			ivar+=pixel*pixel;
		}
	}
	mean=((float)imean)/pixNum;
	var=((float)ivar-((float)imean)*mean)/pixNum;

	return var;
}

#if _MB_RDDATA_
void get_mb_rddata(Macroblock *currMB)
{
//	float psnrY,psnrU,psnrV;
	float varY,varU,varV;
	float ssimY,ssimU,ssimV;
	float ssim3dY,ssim3dU,ssim3dV;
	float stvssimY,stvssimU,stvssimV;
	float tmpstvssimY,tmpstvssimU,tmpstvssimV;
	int rate=currMB->bitcounter[BITS_TOTAL_MB];

//	psnrY=compute_PSNR(pCurImg,enc_picture->imgY,img->opix_y,0,img->opix_x,0,16,16);
//	psnrU=compute_PSNR(pImgOrg[1],enc_picture->imgUV[0],img->opix_c_y,0,img->opix_c_x,0,8,8);
//	psnrV=compute_PSNR(pImgOrg[2],enc_picture->imgUV[1],img->opix_c_y,0,img->opix_c_x,0,8,8);
	varY=compute_VAR(pCurImg,img->opix_y,img->opix_x,16,16);
	varU=compute_VAR(pImgOrg[1],img->opix_c_y,img->opix_c_x,8,8);
	varV=compute_VAR(pImgOrg[2],img->opix_c_y,img->opix_c_x,8,8);

	tmpstvssimY=compute_stVSSIM(pCurImg,enc_picture->imgY,img->opix_y,img->pix_y,img->opix_x,img->pix_x,
		16,16,8,img->number+1,0,&ssimY,&ssim3dY,&stvssimY);
	tmpstvssimU=compute_stVSSIM(pImgOrg[1],enc_picture->imgUV[0],img->opix_c_y,img->pix_c_y,img->opix_c_x,img->pix_c_x,
		8,8,8,img->number+1,1,&ssimU,&ssim3dU,&stvssimU);
	tmpstvssimV=compute_stVSSIM(pImgOrg[2],enc_picture->imgUV[1],img->opix_c_y,img->pix_c_y,img->opix_c_x,img->pix_c_x,
		8,8,8,img->number+1,2,&ssimV,&ssim3dV,&stvssimV);
//	mb_rddata[img->current_mb_nr][1]=psnrY;
//	mb_rddata[img->current_mb_nr][2]=psnrU;
//	mb_rddata[img->current_mb_nr][3]=psnrV;
	mb_rddata[img->current_mb_nr][1]=varY;
	mb_rddata[img->current_mb_nr][2]=varU;
	mb_rddata[img->current_mb_nr][3]=varV;
	mb_rddata[img->current_mb_nr][4]=ssimY;
	mb_rddata[img->current_mb_nr][5]=ssimU;
	mb_rddata[img->current_mb_nr][6]=ssimV;
	mb_rddata[img->current_mb_nr][7]=ssim3dY;
	mb_rddata[img->current_mb_nr][8]=ssim3dU;
	mb_rddata[img->current_mb_nr][9]=ssim3dV;
	mb_rddata[img->current_mb_nr][10]=stvssimY;
	mb_rddata[img->current_mb_nr][11]=stvssimU;
	mb_rddata[img->current_mb_nr][12]=stvssimV;
	mb_rddata[img->current_mb_nr][13]=tmpstvssimY;
	mb_rddata[img->current_mb_nr][14]=tmpstvssimU;
	mb_rddata[img->current_mb_nr][15]=tmpstvssimV;
	mb_rddata[img->current_mb_nr][0]=rate*1.0f;
//	mb_rddata[img->current_mb_nr][1]=currMB->mb_type;
//	mb_rddata[img->current_mb_nr][2]=currMB->skip_flag;
//	mb_rddata[img->current_mb_nr][3]=currMB->cbp;
}
void writeout_slice_rddata()
{
	int i,j;
	int height=img->height;
	int width=img->width;
	for(j=0;j<height/16*width/16;j++)
	{
		fprintf(pmbrddata,"%2d ",j);
		fprintf(pmbrddata,"%8.0f ",mb_rddata[j][0]);
		fprintf(pmbrddata,"%9.4f ",mb_rddata[j][1]);
//		fprintf(pmbrddata,"%2.0f ",mb_rddata[j][2]);
//		fprintf(pmbrddata,"%4.0f ",mb_rddata[j][3]);

		for(i=2;i<MB_DATA_NUM;i++)
			fprintf(pmbrddata,"%8.4f ",mb_rddata[j][i]);
		fprintf(pmbrddata,"\n");
	}

}
#endif
extern void accumulate_metric(float *ave_metric, float cur_metric, int frames);
void find_stvssim(DistMetric *mssim,DistMetric*mssim3d,DistMetric*mstvssim)
{
	float ssimY,ssimU,ssimV;
	float ssim3dY,ssim3dU,ssim3dV;
	float stvssimY,stvssimU,stvssimV;
	int frames=img->number;
	int slice_type=img->type;
	int frames2=stats->frame_ctr[img->type];

#if _RDO_STVSSIM2_
	compute_stVSSIM2(img_org_frm[0], enc_picture->imgY,0,0,0,0,img->height,img->width,8,img->number+1,0,&ssimY,&ssim3dY,&stvssimY);
#else
	compute_stVSSIM(img_org_frm[0], enc_picture->imgY,0,0,0,0,img->height,img->width,8,img->number+1,0,&ssimY,&ssim3dY,&stvssimY);
#endif
	mssim->value[0]=ssimY;
	mssim3d->value[0]=ssim3dY;
	mstvssim->value[0]=stvssimY;
  accumulate_metric(&mssim->average[0],  mssim->value[0],  frames);
  accumulate_metric(&mssim3d->average[0],  mssim3d->value[0],  frames);
  accumulate_metric(&mstvssim->average[0],  mstvssim->value[0],  frames);

	accumulate_metric(&mssim->avslice[slice_type][0], mssim->value[0],  frames2);
	accumulate_metric(&mssim3d->avslice[slice_type][0], mssim3d->value[0],  frames2);
	accumulate_metric(&mstvssim->avslice[slice_type][0], mstvssim->value[0],  frames2);
	if(params->yuv_format!=YUV400)
	{
#if _RDO_STVSSIM2_
	compute_stVSSIM2(img_org_frm[1],enc_picture->imgUV[0],0,0,0,0,img->height_cr,img->width_cr,8,img->number+1,0,&ssimU,&ssim3dU,&stvssimU);
	compute_stVSSIM2(img_org_frm[2],enc_picture->imgUV[1],0,0,0,0,img->height_cr,img->width_cr,8,img->number+1,0,&ssimV,&ssim3dV,&stvssimV);
#else
	compute_stVSSIM(img_org_frm[1],enc_picture->imgUV[0],0,0,0,0,img->height_cr,img->width_cr,8,img->number+1,0,&ssimU,&ssim3dU,&stvssimU);
	compute_stVSSIM(img_org_frm[2],enc_picture->imgUV[1],0,0,0,0,img->height_cr,img->width_cr,8,img->number+1,0,&ssimV,&ssim3dV,&stvssimV);
#endif
	mssim->value[1]=ssimU;
	mssim->value[2]=ssimV;
	mssim3d->value[1]=ssim3dU;
	mssim3d->value[2]=ssim3dV;
	mstvssim->value[1]=stvssimU;
	mstvssim->value[2]=stvssimV;
  accumulate_metric(&mssim->average[1],  mssim->value[1],  frames);
  accumulate_metric(&mssim->average[2],  mssim->value[2],  frames);
  accumulate_metric(&mssim3d->average[1],  mssim3d->value[1],  frames);
  accumulate_metric(&mssim3d->average[2],  mssim3d->value[2],  frames);
  accumulate_metric(&mstvssim->average[1],  mstvssim->value[1],  frames);
  accumulate_metric(&mstvssim->average[2],  mstvssim->value[2],  frames);

	accumulate_metric(&mssim->avslice[slice_type][1], mssim->value[1],  frames2);
	accumulate_metric(&mssim->avslice[slice_type][2], mssim->value[2],  frames2);
	accumulate_metric(&mssim3d->avslice[slice_type][1], mssim3d->value[1],  frames2);
	accumulate_metric(&mssim3d->avslice[slice_type][2], mssim3d->value[2],  frames2);
	accumulate_metric(&mstvssim->avslice[slice_type][1], mstvssim->value[1],  frames2);
	accumulate_metric(&mstvssim->avslice[slice_type][2], mstvssim->value[2],  frames2);
	}
	
}

double adjust_lambda(double lambda,double eta)
{
	double C=0.00000001;
	double lambda2;
	lambda2=lambda;
#if _ADJUST_L1_
//	lambda2=lambda+eta/100;//a1 和调整前比较：重合
//	lambda2=lambda+eta/1000;//a2 和调整前比较：重合 
	lambda2=lambda+eta/10;//a3  效果更差
#endif
#if _ADJUST_L2_
//	lambda2=lambda*exp(eta);//a1 基本重合
//	lambda2=lambda*exp(eta/5);//a2基本重合 
//	lambda2=lambda*exp(eta*2);//a3基本重合
//	lambda2=lambda*exp(eta*10);//a4
//	lambda2=lambda*exp(eta*20);//a5 一些有非常小非常小的提升，大部分还是重合
//	lambda2=lambda*exp(eta*30);//a6
//	lambda2=lambda*exp(eta*50);//a7
//	lambda2=lambda*exp(eta*100);//a8
//	lambda2=lambda*exp(eta*200);//a9
//	lambda2=lambda*exp(eta*300);//a10
//	lambda2=lambda*exp(eta*500);//a11
//	lambda2=lambda*exp(eta*800);//a12
//	lambda2=lambda*exp(eta*1000);//a13
//	lambda2=lambda*exp(-eta);//b1 ==
//	lambda2=lambda*exp(-eta*10);//b2
//	lambda2=lambda*exp(-eta*20);//b3 ==
//	lambda2=lambda*exp(-eta/10);//b4
//	lambda2=lambda*exp(-eta/20);//b5 ==
//	lambda2=lambda*exp(-eta*40);//b6  ==
//	lambda2=lambda*exp(-eta*80);//b7
//	lambda2=lambda*exp(-eta*100);//b8
//	lambda2=lambda*exp(-eta*200);//b9
//	lambda2=lambda*exp(-eta*400);//b10
//	lambda2=lambda*exp(-eta*45);//b11 
//	lambda2=lambda*exp(-eta*50);//b12
//	lambda2=lambda*exp(-eta*55);//b13
//	lambda2=lambda*exp(-eta*60);//b14
//	lambda2=lambda*exp(-eta*65);//b15
//	lambda2=lambda*exp(-eta*70);//b16
//	lambda2=lambda*exp(-eta*75);//b17
//	lambda2=lambda*exp(-eta*23);//c1 
//	lambda2=lambda*exp(-eta*25);//c2
//	lambda2=lambda*exp(-eta*28);//c3
//	lambda2=lambda*exp(-eta*30);//c4
//	lambda2=lambda*exp(-eta*33);//c5
//	lambda2=lambda*exp(-eta*35);//c6
//	lambda2=lambda*exp(-eta*38);//c7
//	lambda2=lambda*exp(-eta*45);//c8
//	lambda2=lambda*exp(-eta*55);//c9
//	lambda2=lambda*exp(-eta*55);//c10
//	lambda2=lambda*exp(-eta*20.5);//d1 
//	lambda2=lambda*exp(-eta*21);//d2 
//	lambda2=lambda*exp(-eta*21.5);//d3 
//	lambda2=lambda*exp(-eta*22);//d4 
//	lambda2=lambda*exp(-eta*22.5);//d5 
//	lambda2=lambda*exp(-eta*20.05);//e1 
//	lambda2=lambda*exp(-eta*20.1);//e2 
//	lambda2=lambda*exp(-eta*20.15);//e3 
//	lambda2=lambda*exp(-eta*20.2);//e4 
//	lambda2=lambda*exp(-eta*20.25);//e5 
//	lambda2=lambda*exp(-eta*20.3);//e6 
//	lambda2=lambda*exp(-eta*20.35);//e7 
//	lambda2=lambda*exp(-eta*20.4);//e8 
//	lambda2=lambda*exp(-eta*20.01);//f1 
//	lambda2=lambda*exp(-eta*20.02);//f2 
//	lambda2=lambda*exp(-eta*20.03);//f3 
//	lambda2=lambda*exp(-eta*20.001);//f4 
//	lambda2=lambda*exp(-eta*20.0001);//f5 
//	lambda2=lambda*exp(-eta*20.00001);//f6 
//	lambda2=lambda*exp(-eta*20.000001);//f7 
//	lambda2=lambda*exp(-eta*20.0000001);//f8 
//	lambda2=lambda*exp(-eta*20.00000001);//f9 
//	lambda2=lambda*exp(-eta*19);//f10
//	lambda2=lambda*exp(-eta*18);//f11
//	lambda2=lambda*exp(-eta*12);//f12
//	lambda2=lambda*exp(-eta*14);//f13
//	lambda2=lambda*exp(-eta*16);//f14
//	lambda2=lambda*exp(-eta*22);//f15
//	lambda2=lambda*exp(-eta*20.00000000001);//f16
//	lambda2=lambda*exp(-eta*11);//g1
//	lambda2=lambda*exp(-eta*11.5);//g2
//	lambda2=lambda*exp(-eta*9);//g3
//	lambda2=lambda*exp(-eta*8);//g4
//	lambda2=lambda*exp(-eta*7);//g5
//	lambda2=lambda*exp(-eta*6);//g6
//	lambda2=lambda*exp(-eta*5);//g7
//	lambda2=lambda*exp(-eta*4);//g8
//	lambda2=lambda*exp(-eta*3);//g9
//	lambda2=lambda*exp(-eta*2);//g10
//	lambda2=lambda*exp(-eta*1.5);//g11
//	lambda2=lambda*exp(-eta*1.3);//g12
//	lambda2=lambda*exp(-eta*1);//g13
//	lambda2=lambda*exp(-eta*0.5);//g14
//	lambda2=lambda*exp(-eta*0.1);//g15
//	lambda2=lambda*exp(-eta*0);//g16
//	lambda2=lambda*exp(-eta*(-1));//g17
//	lambda2=lambda*exp(-eta*(-10));//g18
//	lambda2=lambda*exp(-eta*(20));//g19
//	lambda2=lambda/(eta+C);//g20
//	lambda2=lambda/(eta*100+C);//h1
//	lambda2=lambda/(eta*90+C);//h2
//	lambda2=lambda/(eta*80+C);//h3
//	lambda2=lambda/(eta*100+0.01);//h4
//	lambda2=lambda/(eta*100+0.02);//h5
//	lambda2=lambda/(eta*100+0.03);//h6
//	lambda2=lambda/(eta*100+0.04);//h7
//	lambda2=lambda/(eta*100+0.06);//h8
//	lambda2=lambda/(eta*100+0.08);//h9
//	lambda2=lambda/(eta*100+0.10);//h10
//	lambda2=pow(lambda,eta)/1000;//i1
//	lambda2=pow(lambda,eta)/1500;//i2
//	lambda2=pow(lambda,eta)/2000;//i3
//	lambda2=lambda*exp(-eta*1);//i4
//	lambda2=lambda*exp(-eta*0.08);//i5
//	lambda2=lambda*exp(-eta*0.06);//i6
//	lambda2=lambda*exp(-eta*0.02);//i7
//	lambda2=lambda*exp(-eta*0.009);//i8

//	lambda2=lambda*pow(eta,1);//j1
//	lambda2=lambda*pow(eta,2);//j2 xxxx
//	lambda2=lambda*pow(eta,0.5);//j3
//	lambda2=lambda*pow(eta,0.25);//j4 xxxx
//	lambda2=lambda*pow(eta,0.1);//j5 xxx
//	lambda2=lambda*pow(eta,0.05);//j6 xx
//	lambda2=lambda*pow(eta,3);//j7 xxxx
//	lambda2=lambda*pow(eta,0.2);//j8 xx
//	lambda2=lambda*pow(eta,0.3);//j9 xx
//	lambda2=lambda*pow(eta,0.8);//j10 x
//	lambda2=lambda*pow(eta,1.2);//j11 xx
//	lambda2=lambda*pow(eta,1.5);//j12 xx
//	lambda2=lambda*pow(eta,0.4);//j13 
//	lambda2=lambda*pow(eta,0.5);//j14 
//	lambda2=lambda*pow(eta,0.6);//j15 
//	lambda2=lambda*pow(eta,0.7);//j16 
//	lambda2=lambda*pow(eta,0.89);//j17 x
//	lambda2=lambda*pow(eta,0.1);//j18 xx
//	lambda2=lambda*pow(eta,1.3);//k1 
//	lambda2=lambda*pow(eta,1.4);//k2 
//	lambda2=lambda*pow(eta,1.7);//k3 
//	lambda2=lambda*pow(eta,0.65);//k4
//	lambda2=lambda*pow(eta,0.75);//k5
	lambda2=lambda*pow(eta,0.85);//k6


#endif
	return lambda2;
}

#if _RDO_STVSSIM2_
double lambda_poly(int qp)
{
	double lambda=0;
	return lambda;
}
double lambda_expon(int qp)
{
	double lambda=0;
	return lambda;
}
double lambda_gauss(int qp)
{
	double lambda=0;
	return lambda;
}
double lambda_1(int qp)
{
	double lambda=0;
	return lambda;
}
double lambda_2(int qp)
{
	double a1,b2,b1,lambda;
	a1=  _A1_;
	b2= _B2_;
	b1= _B1_;
	
	lambda=-a1*b2*exp(b1*(qp));//lam2_1
	return lambda;
}
#elif _RDO_STVSSIM_
double lambda_poly(int qp)
{
	double lambda;
	double p1,p2;
	double x=qp*1.0;
	p1=1.0113*pow(10,-5)*pow(x,2)+2.764*pow(10,-4)*x-0.003128;
	p2=-0.06096*pow(x,2)+4.806*x+1245;
	lambda=p1/p2;
	return lambda;
}
double lambda_expon(int qp)
{
	double lambda;
	double x=qp*1.0;
//	lambda=7.5866*pow(10,-7)*exp(0.1789*x);
	lambda=7.5866*pow(10,-6)*exp(0.1789*x);
	return lambda;
}
double lambda_gauss(int qp)
{
	double lambda;
	double x=qp*1.0;
	double tmp;
	tmp=pow((x-60.16),2)/507.6009-pow(x+19.29,2)/560.7424;
	lambda=-3.4223*pow(10,-4)*(x-60.16)/(x+19.29)*exp(-tmp);
	return lambda;

}
double lambda_1(int qp)
{
	double a1=-1.119812399977852e-03;
	double b1=  7.880237189597844e-02;
	double lambda;
	lambda=-a1*exp(b1*qp*1.05)/8;
	return lambda;
}
double lambda_2(int qp)
{
	double a1=  5.883060266548170e-03;
	double b2= -2.229472265847692e-02;
	double b1= 9.279543980380707e-02;
	double lambda;
//	lambda=-a1*b2*exp(b1*qp)/1.8;//lambda_2
//	lambda=-a1*b2*exp(b1*qp*0.9)/1.8;//lam1 A
//	lambda=-a1*b2*exp(b1*qp*0.8)/1.8;//lam2 A
//	lambda=-a1*b2*exp(b1*qp*0.7)/1.8;//lam3
//	lambda=-a1*b2*exp(b1*qp*0.6)/1.8;//lam4
//	lambda=-a1*b2*exp(b1*qp*1.8)/18;//lam5
//	lambda=-a1*b2*exp(b1*qp*2)/18;//lam6
//	lambda=-a1*b2*exp(b1*qp*2.2)/18;//lam7
//	lambda=-a1*b2*exp(b1*qp*3)/250;//lam8
//	lambda=-a1*b2*exp(b1*qp*3)/300;//lam9
//	lambda=-a1*b2*exp(b1*qp*4)/3000;//lam10
//	lambda=-a1*b2*exp(b1*(qp+1)*0.9)/1.8;//
//	lambda=-a1*b2*exp(b1*(qp)*0.9)/1.8-0.0001;//
//	lambda=-a1*b2*exp(b1*(qp)*0.5)/2;//lam1_1
//	lambda=-a1*b2*exp(b1*(qp-12));//lam1_2 A
//	lambda=-a1*b2*exp(b1*(qp-13));//lam2_1 A-
//	lambda=-a1*b2*exp(b1*(qp-14));//lam2_2 A
	lambda=-a1*b2*exp(b1*(qp-15));//lam2_3 A
//	lambda=-a1*b2*exp(b1*(qp-16));//lam2_4 A-
//	lambda=-a1*b2*exp(b1*(qp-15))/1.3;//lam2_5
//	lambda=-a1*b2*exp(b1*(qp-12))*3;//lam3_1
//	lambda=-a1*b2*exp(b1*(qp-12))*2.7;//lam3_2
//	lambda=-a1*b2*exp(b1*(qp-12))*2.5;//lam3_3
//	lambda=-a1*b2*exp(b1*(qp-12))*2;//lam3_4
//	lambda=-a1*b2*exp(b1*(qp-12))*1.5;//lam3_5 +
//	lambda=-a1*b2*exp(b1*(qp-12))*1.8;//lam3_6
//	lambda=-a1*b2*exp(b1*(qp-10))*1.5;//lam3_7
//	lambda=-a1*b2*exp(b1*(qp-15))*3;//lam3_8
//	lambda=-a1*b2*exp(b1*(qp-15))*2;//lam3_9
//	lambda=-a1*b2*exp(b1*(qp-20))*3;//lam3_10
//	lambda=-a1*b2*exp(b1*(qp-20))*7;//lam3_11
//	lambda=-a1*b2*exp(b1*(qp-21))*7;//lam3_12
//	lambda=-a1*b2*exp(b1*(qp-25))*11;//lam3_13
//	lambda=-a1*b2*exp(b1*(qp-12))*1.1;//lam3_14
//	lambda=-a1*b2*exp(b1*(qp-12))*1.3;//lam3_15
//	lambda=-a1*b2*exp(b1*(qp-10))*2.3;//lam4_1
//	lambda=-a1*b2*exp(b1*(qp-10))*2.5;//lam4_2
//	lambda=-a1*b2*exp(b1*(qp-10))*2.7;//lam4_3
//	lambda=-a1*b2*exp(b1*(qp-11))*2.7;//lam4_4
//	lambda=-a1*b2*exp(b1*(qp-11))*2.9;//lam4_5
//	lambda=-a1*b2*exp(b1*(qp-11))*3.1;//lam4_6
//	lambda=-a1*b2*exp(b1*(qp-9))*2.1;//lam4_7
//	lambda=-a1*b2*exp(b1*(qp-9))*2.3;//lam4_8
	return lambda;
}
#elif _RDO_SSIM3D_
double lambda_poly(int qp)
{
	double lambda=0;
	return lambda;
}
double lambda_expon(int qp)
{
	double lambda=0;
	return lambda;
}
double lambda_gauss(int qp)
{
	double lambda=0;
	return lambda;
}
double lambda_1(int qp)
{
	double lambda=0;
	return lambda;
}

double lambda_2(int qp)
{
	double a1= 0.002212;
	double b2= -0.02707;
	double b1=  0.1042;
	double lambda;
//	lambda=-a1*b2*exp(b1*(qp-12));//lam2_1
//	lambda=-a1*b2*exp(b1*(qp-15));//lam2_2
//	lambda=-a1*b2*exp(b1*(qp-19));//lam2_3
//	lambda=-a1*b2*exp(b1*(qp-21));//lam2_4
	lambda=-a1*b2*exp(b1*(qp-23));//lam2_5

	return lambda;
}
#elif _RDO_SSIM_
/*
*/
double lambda_poly(int qp)
{
	double lambda;
	double p1,p2;
	double x=qp*1.0;
	p1=1.4583*pow(10,-5)*pow(x,2)-1.6906*pow(10,-4)*x+0.0004329;
	p2=-0.06096*pow(x,2)+4.806*x+1245;
	lambda=p1/p2;
	return lambda;
}
double lambda_expon(int qp)
{
	double lambda;
	double x=qp*1.0;
//	lambda=3.5470*pow(10,-6)*exp(0.18821*x);
	lambda=3.5470*pow(10,-6)*exp(0.18821*x);
	return lambda;
}
double lambda_gauss(int qp)
{
	double lambda;
	double x=qp*1.0;
	double tmp;
	tmp=pow((x-66.91),2)/587.0929-pow(x+19.29,2)/560.7424;
	lambda=-3.518*pow(10,-5)*(x-66.91)/(x+19.29)*exp(-tmp);
	return lambda;

}
double lambda_1(int qp)
{
	double a1=-5.531041695192410e-04;
	double b1= 8.811336654398498e-02;
	double lambda;
	lambda=-a1*exp(b1*qp)/10;
	return lambda;
}

double lambda_2(int qp)
{
	double a1=  2.001017510889751e-03;
	double b2= -2.324189844454854e-02;
	double b1=   1.061835914044467e-01;
	double lambda;
//	lambda=-a1*b2*exp(b1*qp)/2;
//	lambda=-a1*b2*exp(b1*qp*1.1)/2;//lam1_1
//	lambda=-a1*b2*exp(b1*(qp-12));//lam1_2 A
//	lambda=-a1*b2*exp(b1*(qp-12)*1.3)/2;//lam1_3
//	lambda=-a1*b2*exp(b1*(qp-13));//lam2_1 A
//	lambda=-a1*b2*exp(b1*(qp-14));//lam2_2 A
	lambda=-a1*b2*exp(b1*(qp-14))*8;//lam3_1
//	lambda=-a1*b2*exp(b1*(qp-20))*15;//lam3_2

	return lambda;
}
#endif