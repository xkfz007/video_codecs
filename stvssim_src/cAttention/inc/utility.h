#ifndef UTILITY_H
#define UTILITY_H
#include <stdio.h>
#include <stdlib.h>
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
//#define PI 3.1415926
typedef unsigned char byte;
typedef struct
{
	// Size info
	int x_size, y_framesize;  
	byte *yf, *uf, *vf;                    //!< frame representation
	//	char *yt, *ut, *vt;                    //!< top field
	//	char *yb, *ub, *vb;                    //!< bottom field
} Sourceframe;
void no_mem_exit(char*str);
int get_mem2D(byte ***array2D, int rows, int columns);
void free_mem2D(byte **array2D);
Sourceframe *AllocSourceframe (int xs, int ys);
void FreeSourceframe (Sourceframe *sf);
int ReadOneFrame (FILE*p_in,int FrameNoInFile, int start_frame,int HeaderSize, int xs, int ys, Sourceframe *sf);
void CopyFrameToOldImg (Sourceframe *sf,byte**imgY,byte**imgU,byte**imgV);

#endif