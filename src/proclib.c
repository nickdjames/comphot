#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"
#include "proclib.h"

// return non-zero if images are not the same size
int checksize(long img1[2], long img2[2])
{
	if (img1[0] != img2[0])
		return 1;
	if (img1[1] != img2[1])
		return 2;
	return 0;
}

// create a new image based on the source header but with potentially different size and bitpix
void create_new_img_from_src(fitsfile *tgt, fitsfile *src, long axes[2], int bitpix)
{
	int status = 0;
	int i, nkeys;
	char card[FLEN_CARD];

	fits_create_img(tgt, bitpix, 2, axes, &status);
	if (status)
		fits_report_error(stderr, status);

	fits_get_hdrspace(src, &nkeys, NULL, &status);

	for (i = 1; i <= nkeys; i++)  {
		fits_read_record(src, i, card, &status); /* read keyword */
		if (fits_get_keyclass(card) > TYP_CMPRS_KEY)
			fits_write_record(tgt, card, &status);
	}
}


// add a constant value to an image
int add_const(float *a, float constval, long axes[2])
{
	long i;
	for (i = 0; i < (axes[0]*axes[1]); i++)
		a[i] += constval;

	return 0;
}

// Normalize image to have a mean level of 1.0
int normalize(float *buf, long axes[2])
{
	double sum;
	long i;
	sum = 0;
	for (i = 0; i < (axes[0]*axes[1]); i++)
		sum += buf[i];
	sum /= axes[0]*axes[1];
	for (i = 0; i < (axes[0]*axes[1]); i++)
		buf[i] /= sum;
	return 0;
}

// subtract image b from a with optional normalisation of the median
int image_sub(float *a, float *b, long axes[2], int norm)
{
	long i;
	float med;

	if (norm)
		med = median(b, axes[0]*axes[1]);
	else
		med = 0;


	for (i = 0; i < (axes[0]*axes[1]); i++)
		a[i] -= (b[i]-med);
	return 0;
}

// divide image a by image b
int image_div(float *a, float *b, long axes[2])
{
	long i;
	for (i = 0; i < (axes[0]*axes[1]); i++) {
		if (b[i] > 0)
			a[i] /= b[i];
		else
			a[i] = 0;
	}
	return 0;
}

// Limit image pixles to 0..limit
int image_limit(float *a, float limit, long axes[2])
{
	long i;
	for (i = 0; i < (axes[0]*axes[1]); i++) {
		if (a[i] < 0)
			a[i] = 0;
		if (a[i] > limit)
			a[i] = limit;
	}
	return 0;
}



// Handle CFITSIO status reports
int handle_status(int *status, int warning, const char *msg)
{
	char s[FLEN_ERRMSG];
	int stat;
	stat = *status;
	if (stat) {
		while (fits_read_errmsg(s))
			fprintf(stderr, "%s (%d): %s\n", warning ? "WARNING" : "ERROR", stat, s);
		if (!warning) {
			fprintf(stderr, "FATAL: %s (status = %d), exiting...\n", msg, stat);
			exit(stat);
		}
	}
	*status = 0;
	return stat;
}


// return pointer to pixel at coordinate x,y (starting from 0,0) in top left of image
float *get_pixel(float *image, long x, long y, long axes[2])
{
	if ((x < 0) || (x >= axes[0]))
		return NULL;
	if ((y < 0) || (y >= axes[1]))
		return NULL;
	return &(image[x+y*axes[0]]);
}


int fsortfunc(const void *x1, const void *x2)
{
	return ( *( (float *) x1) > *( (float *) x2) ) ? 1 : -1;
}


// determine median of array of size n (non-destructive)
float median(float inp[], int n)
{
	float *vals, med;
	int i;
	vals = (float *) malloc(sizeof(float)*n);
	for (i = 0; i < n; i++)
		vals[i] = inp[i];
	qsort(vals, n, sizeof(float), fsortfunc); // sort to ascending order
	med = vals[n/2];
	free(vals);
	return med;
}

// determine mean of the lowest x %  of an array of size n (non-destructive)
float mean_min(float inp[], int n, float x)
{
	float *vals;
	float mean;
	int i, ct;
	vals = (float *) malloc(sizeof(float)*n);
	for (i = 0; i < n; i++)
		vals[i] = inp[i];
	qsort(vals, n, sizeof(float), fsortfunc); // sort to ascending order
	ct = (int) floor(n * x + 0.5);
	mean = 0;
	for (i = 0; i < ct; i++)
		mean += vals[i];
	free(vals);
	return mean/ct;
}



// determine image stats
int stats(float inp[], long axes[2], float level[], int steps)
{
	float *vals;
	int i;
	long size;
	int lev;
	size = axes[0] * axes[1];
	vals = (float *) malloc(sizeof(float)*size);
	for (i = 0; i < size; i++)
		vals[i] = inp[i];
	qsort(vals, size, sizeof(float), fsortfunc); // sort to ascending order

	for (i = 0; i < steps; i++) {
		lev = i * size / (steps-1);
		if (lev >= size)
			lev = size-1;
		level[i] = vals[lev];
	}
	free(vals);
	return 0;
}


// determine sky level of array (non-destructive)
float skylevel(float inp[], long axes[2], float nsigma)
{
	float *vals, med, sky;
	float rms;
	int i;
	int ct;
	int n = axes[0] * axes[1];

	vals = (float *) malloc(sizeof(float)*n);
	for (i = 0; i < n; i++)
		vals[i] = inp[i];

	// first determine the median
	qsort(vals, n, sizeof(float), fsortfunc);
	med = vals[n/2];

	// then estimate the RMS using pixels below the median only
	rms = 0;
	for (i = 0; i < n/2; i++)
		rms += powf((vals[i]-med),2);
	rms = sqrtf(rms/(n/2));

	// then calculate the  mean using pixels within nsigma * RMS of median
	// FIXME - This is not a good approach since 1-sigma below the median is asymmetric relative to above
	// For now just use the median
	ct = sky = 0;
	for (i = 0; i < n/2; i++) {
		if (fabsf(vals[i]-med) < (nsigma*rms)) {
			sky += vals[i];
			ct++;
		}
	}
	free(vals);
	// return sky/ct;
	return med;
}

// estimate the image noise level by calculting the RMS of negative pixels. This assumes the sky has been removed
float rms_sky(float *img, long axes[2])
{
	double sumsq = 0;
	unsigned i;
	unsigned ct = 0;
	for (i = 0; i < axes[0] * axes[1]; i++) {
		if (img[i] <= 0) {
			sumsq += img[i] * img[i];
			ct++;
		}
	}
	if (ct > 0)
		return (float) sqrt(sumsq/ct);
	else
		return -1;
}

// generate a Gaussian noise sample
float noise(float sd)
{
	float x;
	float noise, noiseph;
	do {
			x = (float) rand() / RAND_MAX;
	} while (x == 0);
	noise = sd * sqrtf(-logf(x)); /* noise amplitude */
	noiseph = 2 * M_PI * (float) rand() / RAND_MAX; /* noise phase */
	return noise * cosf(noiseph);
}

// add zero mean Gaussian noise of SD to an image
void add_noise(float *img, long axes[2], float sd)
{
	unsigned i;
	for (i = 0; i < axes[0] * axes[1]; i++)
		img[i] += noise(sd);
}

// print a row from an image
void dumpline(float *buf, long axes[], long line)
{
	unsigned i;
	unsigned start, end;

	start = axes[0] * line;
	end = start + axes[0];

	for (i = start; i < end; i++)
		printf("LINE %3ld: %8.2f\n", line, buf[i]);
}

// print an area from an image
void dumparea(float *buf, long axes[], int start[], int r)
{
	long x, y;

	for (y = start[1] - r; y <= start[1] + r ; y++) {
		for (x = start[0] - r; x <= start[0] + r ; x++) {
			printf("%5.0f ", buf[y*axes[0]+x]);
		}
		printf("\n");
	}
}

// determine the 1 dimensional centroid
float centroid1d(float *buf, float sky, int rows, int wid, double *cent, double r, int debug_cen)
{
	// float sum;
	int x, y, idx;
	float sum1, sum2;
	float dx;
	float val;
	int icent;

	icent = (int) floor(*cent+0.5); // nearest integer for centroid

	sum1 = sum2 = 0;
	for (y = 0; y < rows; y++) {
		for (x = -2*r; x <= 2*r; x++) {
			idx = icent + x;
			val = buf[y*wid+idx] - sky;
			if (val < 0)
				val = 0;
			if (debug_cen)
				printf("CEN1 %4d %4d %6.1f %6.1f %4d\n", y, idx, val, sky, icent);
			sum1 += val;
			sum2 += x * val;
		}
	}
	if (sum1 > 0)
		dx = sum2/sum1;
	else
		dx = 0;
	*cent = icent + dx;
	icent = (int) floor(*cent+0.5);

	if (debug_cen)
		printf("CEN2 %7.1f %7.1f %5.2f %6.2f\n", sum1, sum2, dx, *cent);

	sum1 = 0;
	for (y = 0; y < rows; y++) {
		for (x = -r; x <= r; x++) {
			idx = icent + x;
			val = buf[y*wid+idx] - sky;
			sum1 += val;
		}
	}

	return sum1;
}



// determine the 2-D centroid
void find_centroid(float *buf, long axes[2], int start[2], int r)
{
	float sum, sumy, sumx;
	float ycent, xcent;
	int x, y;
	float *pix;

	sum = sumx = sumy = 0;
	for (y = start[1] - r; y <= start[1] + r ; y++) {
		for (x = start[0] - r; x <= start[0] + r ; x++) {
			pix = get_pixel(buf, x, y, axes);
			if (pix) {
				sum += *pix;
				sumy += (y - start[1]) * (*pix);
				sumx += (x - start[0]) * (*pix);
			}
		}
	}
	ycent = sumy/sum;
	xcent = sumx/sum;

	start[0] += nint(xcent);
	start[1] += nint(ycent);

	// printf("find_centroid():  %.2f %.2f %4d %4d %.2f %.2f %.2f\n", xcent, ycent, start[0], start[1], sumx, sumy, sum);

#if 0
	for (y = -r; y <= r ; y++)
		for (x = -r; x <= r ; x++)
			printf("PIX %7.3f %7.3f %.3f\n", x-xcent,  y-ycent, buf[(y+start[1])*axes[0]+(x+start[0])]);
#endif

}

// determine level clip limits depending on proportion of pixels clipped
void cliplimits(float *buf, long axes[], float clip, float *min, float *max, int cliptype)
{
	float *tmp;
	int i;
	long size;

	size = axes[0] * axes[1];
	// allocate temporary image buffer and copy image
	if ((tmp = (float *) malloc(sizeof(float) * size)) == NULL) {
		fprintf(stderr, "Couldn't allocate temporary buffer\n")	;
		exit(5);
	}
	for (i = 0; i < size; i++)
		tmp[i] = buf[i];

	qsort(tmp, size, sizeof(float), fsortfunc); // sort to ascending order
	switch (cliptype) {
	case 0:
	default:
		*min = tmp[(unsigned) floor(clip*size)];
		*max = tmp[(unsigned) floor((1-clip)*size)-1];
		break;
	case 1:
		*max = tmp[(unsigned) floor(0.99*size)-1];
		*min = tmp[(unsigned) floor((0.5-clip)*size)];
		break;
	case 2:
		*min = tmp[(unsigned) floor(6*clip*size)];
		*max = tmp[(unsigned) floor((1-clip)*size)-1];
		break;
	}
	free(tmp);
}



float subtract_sky(float *buf, long axes[], float *rms)
{
	float *tmp;
	int i;
	long size;
	double mean, sum;
	int used;

	size = axes[0] * axes[1];

	// allocate temporary image buffer and copy image
	if ((tmp = (float *) malloc(sizeof(float) * size)) == NULL) {
		fprintf(stderr, "Couldn't allocate temporary buffer\n")	;
		exit(5);
	}
	for (i = 0; i < size; i++)
		tmp[i] = buf[i];

	// first determine the median
	qsort(tmp, size, sizeof(float), fsortfunc); 

	// Then the mean assuming that 10% of the image pixels around the median are sky
	sum = 0;
	used = 0;
	for (i = -size/20; i < size/20; i++) {
		sum += tmp[i + size/2];
		used++;
	}
	mean = sum/used;

	for (i = 0; i < size; i++)
		buf[i] -= mean; // subtract estimated sky level from the image

	// use negative pixels to estimate RMS sky noise
	if (rms) {
		sum = 0;
		used = 0;
		for (i = 0; i < size; i++) {
			if (buf[i] < 0) {
				sum += sqr(buf[i]);
				used++;
			}
		}
		*rms = sqrt(sum/used);
	}

	free(tmp);

	return mean;
}


double aperture(float *buf, long axes[], double centroid[], int radius)
{
	long x, y;
	double dx, dy, sum, sumx, sumy;
	int iteration;
	double count;

	for (iteration = 0; iteration < 6; iteration++) {
		// printf("Iteration %d, centroid %.2f %.2f\n", iteration, centroid[0], centroid[1]);
		sum = sumx = sumy = 0;
		count = 0;
		for (y = floor(centroid[1] - 1.5*radius); y <= floor(centroid[1]+1.5*radius); y++) {
			for (x = floor(centroid[0] - 1.5*radius); x <= floor(centroid[0]+1.5*radius); x++) {
				dx = x - centroid[0];
				dy = y - centroid[1];
				sum += buf[y*axes[0]+x];
				sumx += dx * buf[y*axes[0]+x];
				sumy += dy * buf[y*axes[0]+x];
				if (sqrt(sqr(dx)+sqr(dy)) <= radius) {
					//printf("%6.1f ", buf[y*axes[0]+x]);
					count += buf[y*axes[0]+x];
				} else {
					//printf("       ", buf[y*axes[0]+x]);
				}
				//printf("%6.1f ", buf[y*axes[0]+x]);
			}
			//printf("   PIX\n");
		}
		centroid[0] += sumx/sum;
		centroid[1] += sumy/sum;

	}
	for (y = floor(centroid[1] - 1.5*radius); y <= floor(centroid[1]+1.5*radius); y++) {
		for (x = floor(centroid[0] - 1.5*radius); x <= floor(centroid[0]+1.5*radius); x++) {
			dx = x - centroid[0];
			dy = y - centroid[1];
	//		printf("PROF %.3f %.3f\n", sqrt(sqr(dx)+sqr(dy)), buf[y*axes[0]+x]);
		}
	}

//	printf("Centroid: %.2f %.2f %.2f\n", centroid[0], centroid[1], count);
	return count;
}


static int compfunc(const void *e1, const void *e2)
{
	float *p1, *p2;
	p1 = (float *) e1;
	p2 = (float *) e2;
	return (*p2 < *p1) ? 1 : -1;
}

float combine_minpix(int n, float buf[])
{
	qsort(buf, n, sizeof(float), compfunc);
	return buf[0];
}

float combine_maxpix(int n, float buf[])
{
	qsort(buf, n, sizeof(float), compfunc);
	return buf[n-1];
}

// return the mean of an n-point array excluding outliers. Done in place
float mean_removeoutliers(int n, float buf[], int noremove, int verbose)
{
	int i;
	double m, s, std, mean, m0, m1, s1, s2, err;
	double t0, t1, c = 1.5;

	mean = 0;
	for (i = 0; i < n; i++)
		mean += buf[i];
	m0 = mean/n;

	if (noremove)
		return m0;

	if (verbose > 1) {
		printf("A %7.2f ", m0);
		for (i = 0; i < n; i++)
			printf("%5.0f ", buf[i]);
		printf("\n");
	}

	// sort set and determine median
	qsort(buf, n, sizeof(float), compfunc);
	m = buf[n/2];

	if (verbose > 1) {
		printf("B %7.2f ", m);
		for (i = 0; i < n; i++)
			printf("%5.0f ", buf[i]);
		printf("\n");
	}

	// calculate std based on median
	std = 0;
	for (i = 0; i < n; i++)
		std += sqr(buf[i]-m);
	s = sqrt(std/n);

	do {

		t0 = m - c*s;
		t1 = m + c*s; // limits
#if 1
		// winsorize the set
		for (i = 0; i < n; i++) {
			if (buf[i] < t0)
				buf[i] = t0;
			if (buf[i] > t1)
				buf[i] = t1;
		}
#endif

#if 0
		for (i = 0; i < n; i++) {
			if (buf[i] < t0)
				buf[i] = m;
			if (buf[i] > t1)
				buf[i] = m;
		}
#endif

		// calculate mean and std of new set (m1, s1)
		mean = 0;
		for (i = 0; i < n; i++)
			mean += buf[i];
		m1 = mean/n;
		std = 0;
		for (i = 0; i < n; i++)
			std += sqr(buf[i]-m1);
		s1 = sqrt(std/n);

		s2 = 1.134*s1; // assume normal distribution and c = 1.5

		err = fabs(s2 - s)/s;

		// printf("%.2f %.2f %.2f %.2f %.2f %.2f %.5f\n", m, s, t0, t1, m1, s1, err);

		m = m1;
		s = s2;

	} while (err > 0.002);

	if (verbose > 1) {
		printf("C %7.2f ", m1);
		for (i = 0; i < n; i++)
			printf("%5.0f ", buf[i]);
		printf("=====\n");
	}

	return m;
}

// Apply a filter to an image. Return a pointer to the result
float *img_filt(float *img, long axes[], int rad)
{
	long x,y;
	long i, j;
	int n;

	float *res;
	float *list;
	float *ptr;

	res = (float *) malloc(sizeof(float) * axes[0] * axes[1]);
	list = (float *) malloc(sizeof(float) * (2*rad+1) * (2*rad+1));

	for (x = 0; x < axes[0]; x++) {
		for (y = 0; y < axes[1]; y++) {
			n = 0;
			for (i = -rad; i <= rad; i++) {
				for (j = -rad; j <= rad; j++) {
					ptr = get_pixel(img, x+i, y+j, axes);
					if (ptr)
						list[n++] = *ptr;

				}
			}
			res[y*axes[0]+x] = mean_min(list, n, 0.1);
		}
	}

	free(list);
	return res;
}

// Generate a median tile image with tiles of size n
float *img_medtile(float *img, long axes[], int n)
{
	int x, y, x1, y1;
	int i;
	float *res, *tile, *ptr;
	float val;
	res = (float *) malloc(sizeof(float) * axes[0] * axes[1]);
	tile = (float *) malloc(sizeof(float) * n * n);

	for (y = 0; y < axes[1]; y += n) {
		for (x = 0; x < axes[0]; x += n) {
			i = 0;
			for (y1 = 0; y1 < n; y1++) {
				for (x1 = 0; x1 < n; x1++) {
					ptr = get_pixel(img, x+x1, y+y1, axes);
					if (ptr)
						tile[i++] = *ptr;
				}
			}
			val = median(tile, i);
			for (y1 = 0; y1 < n; y1++) {
				for (x1 = 0; x1 < n; x1++) {
					ptr = get_pixel(res, x+x1, y+y1, axes);
					if (ptr)
						*ptr = val;
				}
			}
		}
	}
	free (tile);
	return res;
}
