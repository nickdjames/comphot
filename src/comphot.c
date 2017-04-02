//
// Do comet photometry on FITS images
// Nick James
//
// Uses libgd (http://www.boutell.com/gd/manual2.0.33.html) for graphics file generation
// Uses cfitsio (http://heasarc.gsfc.nasa.gov/fitsio/) for FITS file manipulation
//
//
//
//

#include "comphot.h"
#include "proclib.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gd.h>
#include <gdfontmb.h>
#include <fitsio.h>
#include <time.h>
#include <assert.h>

#define	MAXBUF	1000
#define	APLIMIT 0.8 // sky background limit for aperture termination (in sigma)
#define SKYRING		40 // arcsec
#define SEARCHRING	20 // arcsec

void strrep(char *s, char c1, char c2)
{
	while (*s != '\0') {
		if (*s == c1)
			*s = c2;
		s++;
	}
}



// rotate the image by 0=>0, 1=>90, 2=>180, 3=>270 but only if the image is square
int rotate(long x, long y, long axes[2], int rot, long *xr, long *yr)
{
	long sz;
	if (axes[1] != axes[0]) // non square image
		rot = ((rot>>1)<<1); // only allow 0 or 180
	sz = axes[0];

	switch (rot) {
	case 0: // 0 deg
	default:
		*xr = x;
		*yr = y;
		break;
	case 1: // 90 deg
		*xr = sz - 1 - y;
		*yr = x;
		break;
	case 2: // 180 deg
		*xr = sz-1-x;
		*yr = sz-1-y;
		break;
	case 3: // 270 deg
		*xr = y;
		*yr = sz-1-x;
		break;
	}
	return 0;
}




int norm_subtract(float *img, float *grad, long axes[2])
{
	float meanpix;
	int i;
	meanpix = 0;
	for (i = 0; i < axes[0]*axes[1]; i++)
		meanpix += grad[i];
	meanpix /= axes[0]*axes[1];
	for (i = 0; i < axes[0]*axes[1]; i++)
		img[i] = img[i] - grad[i] + meanpix;
	return 0;
}



// Calculate the sky level using annular apertures of increasing inner radius until a termination condition is met
float sky_annulus(float *buf, long axes[2], int cent[2], float scale, float start, float width, float rms, float *maxr)
{
	int r, pixused, errors;
	float rmin, rmax, dist;
	float *ptr, *pixels;
	int x, y;
	float sky;

	pixels = (float *) malloc(sizeof(float)*axes[0]*axes[1]); // buffer for aperture pixels
	// First do increasing annuli to determine the coma size.
	if (*maxr > 0)
		rmin = *maxr;
	else
		rmin = start;
	do {
		rmax = rmin + width;
		pixused = 0;
		errors = 0;
		r = (int) ceil(rmax/scale); // Max aperture radius in pixels
		for (y = -r; y <= r; y++) {
			for (x = -r; x <= r; x++) {
				dist = scale * sqrtf(x*x + y*y); // distance from centroid in arcsec
				if ((dist >= rmin) && (dist < rmax)) { // inside annulus
					ptr = get_pixel(buf, cent[0]+x, cent[1]+y, axes);
					if (ptr)
						pixels[pixused++] = *ptr;
					else
						errors++;
				}
			}
		}
		if ((pixused > 100) && (pixused > errors)) {
			sky = median(pixels, pixused);
			// printf("# SKY: %.1f (%0.f -> %0.f) %d %d %.1f\n", rmin, ceil(rmin/scale), ceil(rmax/scale),  pixused, errors, sky);
			if (sky <= APLIMIT * rms)
				break;
			if (*maxr > 0)
				break;
		} else
			break;
		rmin += width/4;
	} while (1);
	*maxr = rmin;
	printf("# Annulus inner radius %.1f arcsec, ringlevel %.1f\n", rmin, sky);
	free(pixels);
	return sky;
}

void plot_circle(gdImagePtr img, int cent[2], double r)
{
	gdImageArc(img, cent[0], cent[1], (int) floor(2*r+0.5), (int) floor(2*r+0.5), 0, 360, gdTrueColor(255, 255, 128));
}

// Generate a check image which shows the comet's extent and the selected aperture sizes
int generate_falsecolour_image(const char *name, float *buf, long axes[2], float rms, float max, int cent[2], float scale, float coma_outer, float sky_inner, int rot)
{
	gdImagePtr image;
	FILE *out;
	float *ptr;
	long x, y, xr, yr;
	float pix;
	// float val;
	int green, blue;
	char txtbuf[MAXBUF];
	int cen[2];

	image = gdImageCreateTrueColor(axes[0], axes[1]);

	for (y = 0; y < axes[1]; y++) {
		for (x = 0; x < axes[0]; x++) {
			ptr = get_pixel(buf, x, y, axes);
			if (ptr)
				pix = *ptr;
			else
				pix = 0;

			green = (int) floorf(255 * logf(2*pix/rms) / fabsf(logf(2*max/rms)));
			if (green < 0) green = 0;
			if (green > 255) green = 255;

			blue = (int) floorf(200 * pix/rms);
			if (blue > 200) blue = 200;
			if (blue < 0) blue = 0;

			rotate(x, y, axes, rot, &xr, &yr);

			gdImageTrueColorPixel(image, xr, yr) = gdTrueColor(0, green, blue);
		}
	}

	rotate(cent[0], cent[1], axes, rot, &xr, &yr);
	cen[0] = xr;
	cen[1] = yr;

	plot_circle(image, cen, ceil(5.6/scale));
	plot_circle(image, cen, ceil(coma_outer/scale));
	plot_circle(image, cen, ceil(sky_inner/scale));
	plot_circle(image, cen, ceil((sky_inner+SKYRING)/scale));
	sprintf(txtbuf, "%.0f\"", coma_outer);
	gdImageString(image, gdFontMediumBold, cen[0], cen[1]+ (int) ceil(coma_outer/scale), (unsigned char *) txtbuf, 0x00FFFFFF);

	out = fopen(name, "wb");
	if (out) {
		gdImageJpeg(image, out, 90);
		fclose(out);
	} else {
		fprintf(stderr, "Could not open output file %s\n", name);
	}
	gdImageDestroy(image);
	return 0;
}


// Generate a check image which shows the comet's extent and the selected aperture sizes
int generate_mono_image(const char *name, float *buf, long axes[2], float rms, float max, int cent[2], float scale, float apradius, int rot)
{
	gdImagePtr image;
	FILE *out;
	float *ptr;
	long x, y, xr, yr;
	float pix;
	int green;
	char txtbuf[MAXBUF];
	int cen[2];

	image = gdImageCreateTrueColor(axes[0], axes[1]);

	for (y = 0; y < axes[1]; y++) {
		for (x = 0; x < axes[0]; x++) {
			ptr = get_pixel(buf, x, y, axes);
			if (ptr)
				pix = *ptr;
			else
				pix = 0;

			green = (int) floorf(255 * logf(2*pix/rms) / fabsf(logf(3*max/rms)));
			if (green < 0) green = 0;
			if (green > 255) green = 255;

			rotate(x, y, axes, rot, &xr, &yr);

			gdImageTrueColorPixel(image, xr, yr) = gdTrueColor(green, green, green);
		}
	}

	rotate(cent[0], cent[1], axes, rot, &xr, &yr);
	cen[0] = xr;
	cen[1] = yr;

	plot_circle(image, cen, ceil(apradius/scale));
	sprintf(txtbuf, "%.0f\"", apradius);
	gdImageString(image, gdFontMediumBold, cen[0], cen[1]+ceil(apradius/scale), (unsigned char *) txtbuf, 0x00FFFFFF);

	out = fopen(name, "wb");
	if (out) {
		gdImageJpeg(image, out, 90);
		fclose(out);
	} else {
		fprintf(stderr, "Could not open output file %s\n", name);
	}
	gdImageDestroy(image);
	return 0;
}

int generate_profile_plot(const char *name, int points, float mag[])
{
	gdImagePtr image;
	FILE *out;
	const int sizex=640, sizey=480;
	float scalex, scaley;
	int i;

	scaley = sizey / (mag[0] - mag[points-1]);
	scalex = sizex / points;

	image = gdImageCreateTrueColor(sizex, sizey);
	for (i = 0; i<points-1; i++) {
		gdImageLine(image,
			(int) floor(0.5+scalex*i),
			(int) floor(0.5*scaley*(mag[0]-mag[i])),
			(int) floor(0.5+scalex*(i+1)),
			(int) floor(0.5*scaley*(mag[0]-mag[i+1])),
			gdTrueColor(255, 255, 128));
	}

	out = fopen(name, "wb");
	if (out) {
		gdImageJpeg(image, out, 90);
		fclose(out);
	}
	gdImageDestroy(image);
	return 0;
}

// generate photometry check
void generate_photom_check(char *name, float *img, long axes[], float rms, int cent[2], float scale, int points, float rad[], float med[])
{
	gdImagePtr image;
	FILE *out;
	float *ptr;
	float pix, pixpos, pixneg;
	long x,y,i;
	int green, blue, red;
	float *checkimg;
	float max, min;
	float logratio;
	int r, r1, r2;

	image = gdImageCreateTrueColor(axes[0], axes[1]);

	checkimg = (float *) malloc(sizeof(float) * axes[0] * axes[1]);
	memcpy(checkimg, img, sizeof(float) * axes[0] * axes[1]);
	ptr = get_pixel(checkimg, cent[0], cent[1], axes);
	if (ptr)
		max = *ptr;
	else
		max = 1e5;
	min = APLIMIT * rms / max;
	logratio = log10f(max/min);

	for (i = points-1; i >= 0; i--) {
		r1 = floor(rad[i]); // outer
		r2 = i ? floor(rad[i-1]) : 0; // inner 
		for (y = -r1; y <= r1; y++) {
			for (x = -r1; x <= r1; x++) {
				r = sqrt(x*x + y*y);
				if ((r <= r1) && (r > r2)) {
					ptr = get_pixel(checkimg, cent[0]+x, cent[1]+y, axes);
					if (ptr)
						*ptr = *ptr - med[i];
				}
			}
		}
	}

	for (y = 0; y < axes[1]; y++) {
		for (x = 0; x < axes[0]; x++) {
			ptr = get_pixel(checkimg, x, y, axes);
			if (ptr)
				pix = *ptr;
			else
				pix = 0;
			if (pix > min)
				pixpos = log10f(pix/min) / logratio;
			else
				pixpos = 0;
			if (pixpos > 1)
				pixpos = 1;
			if (pix < -min)
				pixneg = log10f(-pix/min) / logratio;
			else
				pixneg = 0;
			if (pixneg > 1)
				pixneg = 1;

			blue = (int) floor(255 * pixpos);
			red = (int) floor(255 * pixneg);
			green = 0;

			gdImageTrueColorPixel(image, x, y) = gdTrueColor(red, green, blue);
		}
	}
		plot_circle(image, cent, rad[points-1]);

	out = fopen(name, "wb");
	if (out) {
		gdImageJpeg(image, out, 90);
		fclose(out);
	} else {
		fprintf(stderr, "Could not open output file %s\n", name);
	}
	gdImageDestroy(image);
	free(checkimg);
}


// dump a sky check image
void generate_sky_check(char *name, float *img, long axes[], float skyval, float rms, int cent[2], float scale, float sky_inner)
{
	gdImagePtr image;
	FILE *out;
	float *checkimg, *ptr;
	float pix;
	int green, blue, red;
	long x, y;
	char txtbuf[MAXBUF];

	checkimg = img_medtile(img, axes, 16);
	image = gdImageCreateTrueColor(axes[0], axes[1]);

	for (y = 0; y < axes[1]; y++) {
		for (x = 0; x < axes[0]; x++) {
			ptr = get_pixel(checkimg, x, y, axes);
			if (ptr)
				pix = *ptr-skyval;
			else
				pix = 0;

			if (pix > 0) {
				blue = (int) floor(256 * pix/rms);
				red = 0;
			}
			if (pix < 0) {
				red = (int) floor(256 * (-pix)/rms);
				blue = 0;
			}

			ptr = get_pixel(img, x, y, axes);
			if (ptr)
				pix = *ptr-skyval;
			else
				pix = 0;
			if (pix > 0)
				green = (int) floor(64 * (pix/rms - APLIMIT));
			else
				green  = 0;


			if (green < 0) green = 0;
			if (green > 255) green = 255;

			if (blue > 255) blue = 255;
			if (blue < 0) blue = 0;
			
			if (red > 255) red = 255;
			if (red < 0) red = 0;

			gdImageTrueColorPixel(image, x, y) = gdTrueColor(red, green, blue);
		}
	}
	plot_circle(image, cent, ceil(sky_inner/scale));
	sprintf(txtbuf, "%.0f\"", sky_inner);
	gdImageString(image, gdFontMediumBold, cent[0], cent[1]+ (int) ceil(sky_inner/scale), (unsigned char *) txtbuf, 0x00FFFFFF);

	out = fopen(name, "wb");
	if (out) {
		gdImageJpeg(image, out, 90);
		fclose(out);
	} else {
		fprintf(stderr, "Could not open output file %s\n", name);
	}
	gdImageDestroy(image);
	free(checkimg);
}


// Extract magnitude data from the offset stack
float extract_magnitudes(float *buf, long axes[2], int cent[2], float scale, float step, float max,
		float zp, float background, float rms, int rot, const char *object, float *coma, int border)
{
	int x, y;
	// int ofs;
	float *pixels, pix;
	float *sum_mean, *sum_med, *mean, *med;
	float *rad;
	float aprad, dist;
	int i, r;
	// int rmax;
	int *n1, *n2;
	int point=0, points;
	//float sky;
	float *ptr;
	int error;
	float *mag1, *mag2;
	float maxpix;
	float vem;
	float skyinner;
	char fname[MAXBUF];

	for (i = 0; i < axes[0] * axes[1]; i++) // subtract the initial sky estimate from the image
		buf[i] -= background;

	printf("# Residual sky background in offset stack: %.1f. RMS sky from fixed %.1f\n", median(buf, axes[0] * axes[1]), rms);

	find_centroid(buf, axes, cent, 8);
	printf("# Centroid at %d %d, Max pixel is %.1f\n", cent[0], cent[1], maxpix = *(get_pixel(buf, cent[0], cent[1],  axes)));

	sky_annulus(buf, axes, cent, scale, 10, SEARCHRING, rms, &max); // search outwards from 10 arcsec to get limit of coma

	skyinner = 1.3*max;

	background = sky_annulus(buf, axes, cent, scale, skyinner, SKYRING, rms, &skyinner);
	for (i = 0; i < axes[0] * axes[1]; i++)
		buf[i] -= background;

	// rmax = (int) ceil(max/scale); // max radius in pixels
	// pixels = malloc(4*rmax*rmax*sizeof(float)); // buffer for aperture pixels
	pixels = (float *) malloc(sizeof(float)*axes[0]*axes[1]); // buffer for aperture pixels

	// allocate buffers for the various results
	points = (int) ceil(max / step);
	n1 = (int *) malloc(sizeof(int) * points);
	n2 = (int *) malloc(sizeof(int) * points);
	sum_mean =(float *)  malloc(sizeof(float) * points);
	sum_med = (float *) malloc(sizeof(float) * points);
	mean = (float *) malloc(sizeof(float) * points);
	med = (float *) malloc(sizeof(float) * points);
	mag1 = (float *) malloc(sizeof(float) * points);
	mag2 = (float *) malloc(sizeof(float) * points);
	rad = (float *) malloc(sizeof(float) * points);

	// output the check images
	sprintf(fname, "%s_dump.jpg", object);
	strrep(fname, '/', '_');
	generate_falsecolour_image(fname, buf, axes, rms, maxpix, cent, scale, max, skyinner, rot);
	sprintf(fname, "%s_mono.jpg", object);
	strrep(fname, '/', '_');
	generate_mono_image(fname, buf, axes, rms, maxpix, cent, scale, max, rot);
	sprintf(fname, "%s_skycheck.jpg", object);
	strrep(fname, '/', '_');
	generate_sky_check(fname, buf, axes, -background, rms, cent, scale, max);

	aprad = step;
	do { // measure for each aperture
		n1[point] = n2[point] = 0;
		sum_mean[point] = mean[point] = 0;
		error = 0;
		r = (int) ceil(aprad/scale); // max radius for this aperture in pixels
		rad[point] = r;
		for (y = -r; y <= r; y++) {
			for (x = -r; x <= r; x++) {
				ptr = get_pixel(buf, cent[0]+x, cent[1]+y, axes);
				if (ptr) {
					pix = *ptr;
					dist = scale * sqrtf(x*x + y*y);
					if (dist < aprad) { // accumulate pixels inside an aperture of radius aprad
						sum_mean[point] += pix;
						n1[point]++;
					}
					if ((dist < aprad) && (dist >= (aprad-step))) { // Use pixels in annulus (aprad-step) to aprad
						mean[point] += pix;
						pixels[n2[point]] = pix;
						n2[point]++;
					}
				} else {
					error++;
				}
			}
		}
		mean[point] /= n2[point];
		med[point] = median(pixels, n2[point]);
		aprad += step;
		point++;
	} while (aprad < max);
	printf("# Total pixels requested but outside frame: %d\n", error);

	for (i = 0; i < point; i++) { // go through each annulus and dump the results
		double r_outer, r_inner, area;
		r_outer = (i+1) * step;
		r_inner = i * step;
		area = M_PI * (r_outer * r_outer - r_inner*r_inner) / (scale*scale); // area of annulus in square pixels
		if (i == 0)
			sum_med[0] = sum_mean[0]; // start with simple sum for inner circle
		else
			sum_med[i] = (float) (sum_med[i-1] + med[i] * area); // accumulate rings of median samples
		mag1[i] = zp - 2.5f * log10f(sum_mean[i]);
		mag2[i] = zp - 2.5f * log10f(sum_med[i]);
		printf("# %5.1f | %6d %7.0f | %5.1f %5.1f %5d %6.0f %6.0f %6.0f | %5.2f %5.2f\n", (i+1)*step, n1[i], sum_mean[i],  mean[i], med[i], n2[i],mean[i]*n2[i],  med[i]*n2[i],  sum_med[i], mag1[i], mag2[i]);
	}
	sprintf(fname, "%s_photcheck.jpg", object);
	strrep(fname, '/', '_');
	generate_photom_check(fname, buf, axes, rms, cent, scale, point, rad, med);

	generate_profile_plot("profile.jpg", point, mag2);

	// output the multibox estimates
	printf(" 10x10  20x20  30x30  40x40  50x50  60x60\n");
	for (i = 0; i < 6; i++)
		if (point > i)
			printf("%6.2f ", mag1[i]);
		else
			printf("       ");
	printf("# Counts\n");
	for (i = 0; i < 6; i++)
		if (point > i)
			printf("%6.2f ", mag2[i]);
		else
			printf("       ");
	printf("# Median annuli\n");
	printf("Total integrated magnitude: %.2f (radius %.1f arcsec)\n", mag2[point-1], (point+1)*step);
	vem = mag2[point-1];
	if (coma)
		*coma = 2.0f * (point+1)*step / 60.0f;

	// release all of the allocated storage
	free(n1);
	free(n2);
	free(sum_mean);
	free(sum_med);
	free(mean);
	free(med);
	free(mag2);
	free(mag1);
	free(pixels);

	return vem;
}

void process( const ComphotConfig* config )
{
	fitsfile *offset_img, *fixed_img, *grad_img;
	// FILE *flat;
	int status;
	long axes[2], axeschk[2];
	int cent[2];
	// int bitpix;
	float *offset_buf, *fixed_buf, *grad_buf;
	float nulval = 0;
	int anynul;
	float step;
	double scalex, scaley, exposure, rotation;
	float scale;
	double zp;
	char date_obs[FLEN_CARD];
	char telescope[FLEN_CARD];
	char instrument[FLEN_CARD];
	char observer[FLEN_CARD];
	char object[FLEN_CARD];
	char creator[FLEN_CARD];
	char fcombine[FLEN_CARD];
	int isgrad = 0;
	int rot;
	float skyfix, skyofs, skymag, rms;
	float vem, coma;
	int obs_yr, obs_mn, obs_da, obs_hr, obs_min, obs_sec;

	status = 0;
	fits_clear_errmsg();

	assert(config);
	if (config->flatimage) { // optional flat image
		if (fits_open_file(&grad_img, config->flatimage, READONLY, &status))
			handle_status(&status, 0, "Flat image");
		isgrad = 1;
	}

	assert(config->offsetimage);
	if (fits_open_file(&offset_img, config->offsetimage, READONLY, &status)) {
		printf("File is %s\n", config->offsetimage);
		handle_status(&status, 0, "Offset image");
	}

	assert(config->fixedimage);
	if (fits_open_file(&fixed_img, config->fixedimage, READONLY, &status))
		handle_status(&status, 0, "Non offset image");

	cent[0] = config->cenx;
	cent[1] = config->ceny; // start point for centroid

	status = 0;
	fits_clear_errmsg();

	// get dimensions - both images need to be the same size
	fits_get_img_size(offset_img, 2, axes, &status);
	fits_get_img_size(offset_img, 2, axeschk, &status);
	handle_status(&status, 0, "Image size");
	if ((axes[0] != axeschk[0]) || (axes[1] != axeschk[1])) {
		fprintf(stderr, "Fixed and offset images are of different size! %ld %ld %ld %ld\n", axes[0], axes[1], axeschk[0], axeschk[1]);
		exit(1);
	}

	// Read essential keys
	fits_read_key(fixed_img, TDOUBLE, "CDELT1", &scalex, NULL, &status);
	handle_status(&status, 0, "CDELT1 tag missing - non-offset frame needs astrometric cal");
	fits_read_key(fixed_img, TDOUBLE, "CDELT2", &scaley, NULL, &status);
	handle_status(&status, 0, "CDELT2 tag missing - non-offset frame needs astrometric cal");
	fits_read_key(fixed_img, TDOUBLE, "CROTA2", &rotation, NULL, &status);
	handle_status(&status, 0, "CROTA2 tag missing - non-offset frame needs astrometric cal");
	fits_read_key(fixed_img, TDOUBLE, "MZERO", &zp, NULL, &status);
	handle_status(&status, 0, "MZERO tag missing - non-offset frame needs photometric cal");
	fits_read_key(fixed_img, TSTRING, "DATE-OBS", date_obs, NULL, &status);
	handle_status(&status, 0, "DATE-OBS tag missing");

	// read optional keys
	fits_read_key(fixed_img, TDOUBLE, "EXPTIME", &exposure, NULL, &status);
	if (handle_status(&status, 1, ""))
		exposure = 0;
	fits_read_key(fixed_img, TSTRING, "TELESCOP", telescope, NULL, &status);
	if (handle_status(&status, 1, ""))
		strcpy(telescope, "Unknown");
	fits_read_key(fixed_img, TSTRING, "INSTRUME", instrument, NULL, &status);
	if (handle_status(&status, 1, ""))
		strcpy(telescope, "Unknown");
	fits_read_key(fixed_img, TSTRING, "OBSERVER", observer, NULL, &status);
	if (handle_status(&status, 1, ""))
		strcpy(observer, "Unknown");
	fits_read_key(fixed_img, TSTRING, "OBJECT", object, NULL, &status);
	if (handle_status(&status, 1, ""))
		strcpy(object, "Unknown");
	fits_read_key(fixed_img, TSTRING, "CREATOR", creator, NULL, &status);
	if (handle_status(&status, 1, ""))
		strcpy(object, "Unknown");

	printf("# Comphot version %s\n", VERSION);

	// Look for fcombine comment and use time string if available
	fits_read_str(fixed_img, "Mid time:", fcombine, &status);
	if (handle_status(&status, 1, "") == 0)
		printf("%s\n", &fcombine[8]);

	if ((fabs(scalex/scaley) - 1) > 0.0001) { // check for square pixels
		fprintf(stderr, "Image does not have square pixels (%f x %f). Exiting...\n", scalex, scaley);
		exit(2);
	}

	if (scaley > 0) // negative scale implies image inverted
		rotation += 180;


#if 1
	rot = (int) floor(0.5 + (rotation/90));
	while (rot < 0)
		rot += 4;
	while (rot > 3)
		rot -= 4;
	rotation = fmod(360 + rotation - 90 * rot, 360);
#endif

	scale = (float) fabs(scalex * 3600.0); // convert scale to arcsec/pix

	if (scale < 3) // select aperture step depending on how big the pixels are
		step = 5.64f;
	else if (scale < 6)
		step = 2*5.64f;
	else
		step = 3*5.64f;

	printf("Date: %s, Exposure: %.1f s\n", date_obs, exposure);
	sscanf(date_obs, "%d-%d-%dT%d:%d:%d", &obs_yr, &obs_mn, &obs_da, &obs_hr, &obs_min, &obs_sec);

	printf("Telescope: %s, Camera: %s\n", telescope, instrument);
	printf("Observer: %s\n", observer);
	printf("Object: %s\n", object);
	printf("Scale: %.2f \"/pix, FoV %.1fx%.1f arcmin, PA: %.1f deg, ZP: %.2f mag\n", scale, 60*fabs(scalex)*axes[0], 60*fabs(scaley)*axes[1], rotation,  zp);

	// allocate image buffers
	offset_buf = (float *) malloc(sizeof(float) * axes[0] * axes[1]);
	fixed_buf = (float *) malloc(sizeof(float) * axes[0] * axes[1]);
	if (isgrad)
		grad_buf = (float *) malloc(sizeof(float) * axes[0] * axes[1]);

	// Read image data
	fits_read_img(offset_img, TFLOAT, 1, axes[0]*axes[1], &nulval, offset_buf, &anynul, &status);
	fits_read_img(fixed_img, TFLOAT, 1, axes[0]*axes[1], &nulval, fixed_buf, &anynul, &status);
	if (isgrad)
		fits_read_img(grad_img, TFLOAT, 1, axes[0]*axes[1], &nulval, grad_buf, &anynul, &status); // FIXME dimensions!!
	handle_status(&status, 0, "Read buffer");

	// If the image types are integer add a small amount of float Gaussian noise - this helps the median estimators converge
	add_noise(offset_buf, axes, 1.0);
	add_noise(fixed_buf, axes, 1.0);
	// sky = median(fixed_buf, axes[0]*axes[1]); // determine sky background for fixed image
	skyfix = skylevel(fixed_buf, axes, 0.2f);
	skymag = (float) zp - 2.5f * log10f(skyfix/powf(scale, 2));
	printf("Sky background before normalization (non-offset): %.1f  (%.2f mag/sqarcsec)\n", skyfix, skymag);
	skyofs = skylevel(offset_buf, axes, 0.2f);
	skymag = (float) zp - 2.5f * log10f(skyofs/powf(scale, 2));
	printf("# Sky background before normalization (offset): %.1f  (%.2f mag/sqarcsec)\n", skyofs, skymag);

	if (isgrad) {
		norm_subtract(offset_buf, grad_buf, axes);
		norm_subtract(fixed_buf, grad_buf, axes);
	}
	// skyfix = median(fixed_buf, axes[0]*axes[1]); // determine sky background for fixed image
	// skyofs = median(offset_buf, axes[0]*axes[1]); // determine sky background for fixed image
	skyfix = skylevel(fixed_buf, axes, 0.2f);
	skyofs = skylevel(offset_buf, axes, 0.2f);
	printf("# Sky background: %.1f (non-offset) %.1f (offset)\n", skyfix, skyofs);

	// do the main processing job
	rot = 0; // FIXME
	rms = rms_sky(fixed_buf, axes, skyfix, config->border); // calculate RMS noise from fixed stack
	vem = extract_magnitudes(offset_buf, axes, cent, scale, step, config->apradius, (float) zp, skyofs, rms, rot, object, &coma, config->border);

	printf("ICQ:  %4d %2d %5.2f    %4.1f   %4.1f\n",
		obs_yr, obs_mn, obs_da + (3600 * obs_hr + 60 * obs_min + obs_sec)/86400.0,
		vem, coma
	);

	printf("COMPHOT: %s %4d %02d %06.3f %6.2f %6.2f %6.2f %6.2f %6.2f %7.1f %6.2f %s %s %s %s\n",
		VERSION,
		obs_yr, obs_mn, obs_da + (3600 * obs_hr + 60 * obs_min + obs_sec)/86400.0,
		vem, coma, skymag, zp, rms, skyofs, scale, creator, observer, object, config->offsetimage);

	// then release storage and close image files
	free(offset_buf);
	free(fixed_buf);
	fits_close_file(offset_img, &status);
	fits_close_file(fixed_img, &status);
	if (isgrad) {
		free(grad_buf);
		fits_close_file(grad_img, &status);
	}
}
