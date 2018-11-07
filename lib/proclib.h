#ifndef PROCLIB_H
#define PROCLIB_H
#pragma once

#include <fitsio.h>
#include "fitting.h"

typedef enum {
	COORD_USE=0, COORD_REJECT
} COORDTYPE;

typedef struct {
	COORDTYPE type;
	double x,y; // centroid
	double range;
	double psf_fwhm, psf_max, psf_sky; // used if PSF measured
} COORD;

#define	MAXLEN	1000

#define	sqr(x) 	((x) * (x))

#define nint(x) ((int) floor(x+0.5))

typedef enum {
	FILT_MEANMIN, FILT_MEDIAN, FILT_SIGCLIP
} FILTTYPE;

typedef enum {
	APER_NOCENTROID, APER_CENTROID
} PHOTOM_TYPE;

typedef enum {
	FILT_GAUSSIAN
} FILTERTYPE;

typedef struct {
	PHOTOM_TYPE type;
	int centroid;
	long centre[2];
	double r1, r2;
	double counts;
	double sky;
} PHOTOM_APER;

typedef struct {
	int separable; // set to indicate that 2D filtering can be done as two passes of 1D filtering
	int n1, n2;
	double *coef;
} FILTER;

typedef struct {
	int bitpix;
	const char *obstime;
	const char *comment;
} FITSDESC;

typedef enum  {
	COMBINE_CLIP, COMBINE_MAX, COMBINE_MIN, COMBINE_MEDIAN, COMBINE_MEAN, COMBINE_SUM
} COMBINE_TYPE;

void create_new_img_from_src(fitsfile *tgt, fitsfile *src, long axes[2], int bitpix);

void dumpline(float *buf, long axes[], long line);
float subtract_sky(float *buf, long axes[], float *rms);
double aperture(float *buf, long axes[], double centroid[], int radius);
void find_centroid(float *buf, long axes[], long start[], int r);
void dumparea(float *buf, long axes[], int start[], int r);

// return the mean of an n-point array excluding outliers. Done in place
float mean_removeoutliers(int n, float buf[], int noremove, int verbose);
float combine_minpix(int n, float buf[]);
float combine_maxpix(int n, float buf[]);
float centroid1d(float *buf, float sky, int rows, int wid,  double *cent, double r, int debug);
float centroid2d(float *buf, float sky, long axes[2], double centre[2], double r);
void cliplimits(float *buf, long axes[], float clip, float *min, float *max, int cliptype);
float median(float vals[], int n);
float rmsneg(float inp[], int n, float offset);
int stats(float inp[], long axes[2], float level[], int steps);
float skylevel(float vals[], long axes[2], float nsigma);
float *get_pixel(float *image, long x, long y, long axes[2]);
float rms_sky(float *img, long axes[2], float background, int border);
void add_noise(float *img, long axes[2], float sd);

int add_const(float *a, float constval, long axes[2]);
int normalize(float *buf, long axes[2]);
int image_add(float *a, float *b, long axes[2], int norm);
int image_sub(float *a, float *b, long axes[2], int norm);
int image_div(float *a, float *b, long axes[2]);
int image_mult(float *a, float *b, long axes[2]);
int image_limit(float *a, float limit, long axes[2]);
int image_flip(float *img, long axes[2]);
int image_invert(float *buf, long axes[2]);
int image_ped(float *a, int pedestal, long axes[2]);
int handle_status(int *status, int warning, const char *msg);

float *img_filt(float *img, long axes[], int rad, FILTTYPE type);
float *img_medtile(float *img, long axes[], int n);

int checksize(long img1[2], long img2[2]);

int XYZ_to_sRGB(float *r, float *g, float *b, float *x, float *y, float *z, long axes[2]);

int aperture_photom(float *img, long axes[2], PHOTOM_APER *ap);
float *filter(float *buf, long axes[], FILTER *filt, int inplace);
int asinh_stretch(float *buf, long axes[2], double a, double Q, float *min, float *max);
int logsqrt_stretch(float *buf, long axes[2], float *min, float *max);

int fit_gradient(float *buf, long axes[], int order, int iterations);
FILTER *gen_filter(FILTERTYPE ft, double param1);

float *denoise(float *buf, long axes[], float skylev, float skyrms, int inplace, float sigma);

float *binning(float *buf, long axes[2], int bin);
int write_as_new_fits(float *buf, long axes[2], const char *name, FITSDESC *fd);

float *comet_stack(int n, float *images[], long axes[2], COORD starshift[], COORD cometshift[], int iterations, int verbose);
float *stack_images(int n, float *images[], long axes[2], int drizzle, COORD shiftlist[], COMBINE_TYPE type, int verbose);
int track_object(int n, float *images[], long axes[2], float skylevel[], COORD start, double r, COORD shiftlist[], GAUSSFIT *psf);
COORD centre_alignlist(int n, COORD list[]);
int mirror_image(float *buf, long axes[2]);

#endif // PROCLIB_H

