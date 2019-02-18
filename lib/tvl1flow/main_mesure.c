// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2011, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#ifndef DISABLE_OMP
#include <omp.h>
#endif//DISABLE_OMP

#include "iio.h"

#include "tvl1flow_lib.c"


#define PAR_DEFAULT_OUTFLOW "flow.flo"
#define PAR_DEFAULT_NPROC   8
#define PAR_DEFAULT_TAU     0.25
#define PAR_DEFAULT_LAMBDA  0.15
#define PAR_DEFAULT_THETA   0.3
#define PAR_DEFAULT_NSCALES 100
#define PAR_DEFAULT_ZFACTOR 0.5
#define PAR_DEFAULT_NWARPS  5
#define PAR_DEFAULT_EPSILON 0.01
#define PAR_DEFAULT_VERBOSE 0


/**
 *
 *  Function to read images using the iio library
 *  It always returns an allocated the image.
 *
 */
static float **read_video(const char *filename, int *w, int *h, int f, int l)
{
	float **output = xmalloc((l-f+1) * sizeof(*output));
	for(int i = f,j = 0; i <= l; ++i,++j)
	{
		char gen_filename[1024];
		sprintf(gen_filename, filename, i);
		float *f = iio_read_image_float(gen_filename, w, h);
		if (!f)
			fprintf(stderr, "ERROR: could not read image from file "
					"\"%s\"\n", gen_filename);
		output[j] = f;
	}
	return output;
}

static float **read_flow(const char *filename, int *w, int *h, int f, int l)
{
	int c;
	float **output = xmalloc((l-f) * sizeof(*output));
	for(int i = f,j = 0; i < l; ++i,++j)
	{
		char gen_filename[1024];
		sprintf(gen_filename, filename, i);
		float *f = iio_read_image_float_split(gen_filename, w, h, &c);
		if (!f)
			fprintf(stderr, "ERROR: could not read flow from file "
					"\"%s\"\n", gen_filename);
		output[j] = f;
	}
	return output;
}


/**
 *
 *  Main program:
 *   This program reads the following parameters from the console and
 *   then computes the optical flow:
 *   -nprocs      number of threads to use (OpenMP library)
 *   -I0          first image
 *   -I1          second image
 *   -tau         time step in the numerical scheme
 *   -lambda      data term weight parameter
 *   -theta       tightness parameter
 *   -nscales     number of scales in the pyramidal structure
 *   -zfactor     downsampling factor for creating the scales
 *   -nwarps      number of warps per scales
 *   -epsilon     stopping criterion threshold for the iterative process
 *   -out         name of the output flow field
 *   -verbose     switch on/off messages
 *
 */
int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: %s I0 I1 [out "
		//                       0 1  2   3
		"nproc tau lambda theta nscales zfactor nwarps epsilon "
		//  4  5   6      7     8       9       10     11
		"verbose]\n", *argv);
		// 12
		return EXIT_FAILURE;
	}

	//read the parameters
	int i = 1;
	char* vid_name  = argv[i]; i++;
	int   f         = atoi(argv[i]);   i++;
	int   l         = atoi(argv[i]);   i++;
	char* of_name   = argv[i]; i++;
	char* outpath   = argv[i]; i++;
	float lambda  = (argc>i)? atof(argv[i]): PAR_DEFAULT_LAMBDA;  i++;
	int verbose = PAR_DEFAULT_VERBOSE;
	int nproc = PAR_DEFAULT_NPROC;

	//check parameters
	if (lambda <= 0) {
		lambda = PAR_DEFAULT_LAMBDA;
		if (verbose) fprintf(stderr, "warning: "
				"lambda changed to %g\n", lambda);
	}

#ifndef DISABLE_OMP
	if (nproc > 0)
		omp_set_num_threads(nproc);
#endif//DISABLE_OMP

	// read the input images
	int    nx, ny;
	float **vid  = read_video(vid_name, &nx, &ny, f, l);
	float **flow = read_flow(of_name, &nx, &ny, f, l);
	
	char gen_outpath[1024];
	float error = 0.;
	//allocate memory for the diff
	float *diff = xmalloc(nx*ny * sizeof*diff);
	for(int i = 0, t=f; i < (l-f); ++i,++t)
	{	
		sprintf(gen_outpath, outpath, t);
		//compute the optical flow
		error += energy_optic_flow(
				vid[i], vid[i+1], flow[i], flow[i]+(nx*ny), diff, nx, ny, lambda
				);

		//save the difference
		iio_save_image_float(gen_outpath, diff, nx, ny);
	}
	error /= (l-f);
	printf("Energy error: %f\n", error);

	//delete allocated memory
	for(int t = 0; t < (l-f); ++t)
	{
		free(vid[t]);
		free(flow[t]);
	}
	free(vid[l-f]);
	free(vid);
	free(flow);
	free(diff);

	return EXIT_SUCCESS;
}
