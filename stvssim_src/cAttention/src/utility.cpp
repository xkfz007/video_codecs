#include "utility.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

void no_mem_exit(char*str)
{
	printf("%s\n",str);
	exit(-1);
}
int get_mem2D(byte ***array2D, int rows, int columns)
{
	int i;

	if((*array2D      = (byte**)calloc(rows,        sizeof(byte*))) == NULL)
		no_mem_exit("get_mem2D: array2D");
	if(((*array2D)[0] = (byte* )calloc(columns*rows,sizeof(byte ))) == NULL)
		no_mem_exit("get_mem2D: array2D");

	for(i=1;i<rows;i++)
		(*array2D)[i] = (*array2D)[i-1] + columns ;

	return rows*columns;
}
void free_mem2D(byte **array2D)
{
	if (array2D)
	{
		if (array2D[0])
			free (array2D[0]);
		else 
		{
			printf("free_mem2D: trying to free unused memory");
			exit(-1);
		}

		free (array2D);
	} else
	{
		printf ("free_mem2D: trying to free unused memory");
		exit(-1);
	}
}
Sourceframe *AllocSourceframe (int xs, int ys)
{
	Sourceframe *sf = NULL;
	const unsigned int bytes_y = xs*ys;
	const unsigned int bytes_uv = (xs*ys)/4;

	if ((sf = (Sourceframe*)calloc (1, sizeof (Sourceframe))) == NULL) 
		no_mem_exit ("ReadOneFrame: sf");
	if (sf->yf == NULL) if ((sf->yf = (byte*)calloc (1, bytes_y)) == NULL) no_mem_exit ("ReadOneFrame: sf->yf");
	if (sf->uf == NULL) if ((sf->uf =  (byte*)calloc (1, bytes_uv)) == NULL) no_mem_exit ("ReadOneFrame: sf->uf");
	if (sf->vf == NULL) if ((sf->vf =  (byte*)calloc (1, bytes_uv)) == NULL) no_mem_exit ("ReadOneFrame: sf->vf");

	sf->x_size = xs;
	sf->y_framesize = ys;

	return sf;
}
void FreeSourceframe (Sourceframe *sf)
{
	if (sf!=NULL) 
	{
		if (sf->yf != NULL) free (sf->yf);
		if (sf->uf != NULL) free (sf->uf);
		if (sf->vf != NULL) free (sf->vf);
		free (sf);
	}
}
int ReadOneFrame (FILE*p_in,int FrameNoInFile, int start_frame,int HeaderSize, int xs, int ys, Sourceframe *sf)
{
	int i;

	const unsigned int bytes_y = xs*ys;
	const unsigned int bytes_uv = (xs*ys)/4;
	const int framesize_in_bytes = bytes_y + 2*bytes_uv;

	// 	 assert (xs % MB_BLOCK_SIZE == 0);
	// 	 assert (ys % MB_BLOCK_SIZE == 0);
	assert (p_in != NULL);
	assert (sf != NULL);
	assert (sf->yf != NULL);

	if (fseek (p_in, HeaderSize, SEEK_SET) != 0)
	{
		printf ("ReadOneFrame: cannot fseek to (Header size) in p_in");
		return -1;
	}

	// the reason for the following loop is to support source files bigger than
	// MAXINT.  In most operating systems, including Windows, it is possible to
	// fseek to file positions bigger than MAXINT by using this relative seeking
	// technique.  StW, 12/30/02
	// Skip starting frames
	for (i=0; i<start_frame; i++)
		if (fseek (p_in, framesize_in_bytes, SEEK_CUR) != 0) 
		{
			printf ("ReadOneFrame: cannot advance file pointer in p_in beyond frame %d, looping to picture zero\n", i);
			if (fseek (p_in, HeaderSize, SEEK_SET) != 0)
			{
				printf("fseek Error\n");
			}
			return (-1);
		} 

		for (i=0; i<FrameNoInFile; i++)
			if (fseek (p_in, framesize_in_bytes, SEEK_CUR) != 0) 
			{
				printf ("ReadOneFrame: cannot advance file pointer in p_in beyond frame %d, looping to picture zero\n", i);
				if (fseek (p_in, HeaderSize, SEEK_SET) != 0)
					printf ("ReadOneFrame: cannot fseek to (Header size) in p_in");
				return(-1);
			}

			// Here we are at the correct position for the source frame in the file.  Now
			// read it.
			if (fread (sf->yf, 1, bytes_y, p_in) != bytes_y)
			{
		//		printf ("ReadOneFrame: cannot read %d bytes from input file, unexpected EOF?, exiting", bytes_y);
				//		 report_stats_on_error();
				return (-1);
			}
			if (fread (sf->uf, 1, bytes_uv, p_in) != bytes_uv)
			{
		//		printf ("ReadOneFrame: cannot read %d bytes from input file, unexpected EOF?, exiting", bytes_uv);
				//	 report_stats_on_error();
				return (-1);
			}
			if (fread (sf->vf, 1, bytes_uv, p_in) != bytes_uv)
			{
			//	printf ("ReadOneFrame: cannot read %d bytes from input file, unexpected EOF?, exiting", bytes_uv);
				///	 report_stats_on_error();
				return (-1);
			}
			return 1;

}
void CopyFrameToOldImg (Sourceframe *sf,byte**imgY,byte**imgU,byte**imgV)
{
	int x, y;


	for (y=0; y<sf->y_framesize; y++)
	{
		for (x=0; x<sf->x_size; x++)
		{
			imgY[y][x] = sf->yf[y*sf->x_size+x];
		}
	}

	for (y=0; y<sf->y_framesize/2; y++)
		for (x=0; x<sf->x_size/2; x++)
		{
			imgU[y][x] = sf->uf[y*sf->x_size/2+x];
			imgV[y][x] = sf->vf[y*sf->x_size/2+x];
		}
}
