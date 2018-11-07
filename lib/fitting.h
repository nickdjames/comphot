#ifndef FITTING_H

#define	FITTING_H

typedef	struct {
	int n;
	double A, SD;
	double sky;
	double reject;
	double centre[2];
} GAUSSFIT;


int gauss_psf(float *img, long axes[2], double centre[2], double r1, double r2, double r3, GAUSSFIT *fit, int verbose);

#endif
