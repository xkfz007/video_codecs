#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "saptialattention.h"

int main(int argc,char*argv[])
{
#ifdef WIN32
    char* yuvpath="e:/sequences/qcif/suzie_qcif.yuv";
	char* savepath="d:/kk/";
	
	int xs=176;
	int ys=144;
#else
	if(argc<5)
	{
		printf("Usage:attention \"yuvpath\" \"savepath\" width height\n");
		exit(0);
	}
	char*yuvpath=argv[1];
	char*savepath=argv[2];
	int xs=atoi(argv[3]);
	int ys=atoi(argv[4]);

#endif

	seqSpatialAttention(yuvpath,xs,ys,savepath);


//	free(avipath);
//	free(name);

	//	testStaticNovelty(yuvpath,xs,ys);
	//SpatialAttention(yuvpath,xs,ys);
	return 0;
}