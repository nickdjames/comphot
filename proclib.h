#define	MAXLEN	1000

typedef	struct {
	double x,y;
} POINT;

#define	sqr(x) 	((x) * (x))

#define nint(x) ((int) floor(x+0.5))

void create_new_img_from_src(fitsfile *tgt, fitsfile *src, long axes[2], int bitpix);

void dumpline(float *buf, long axes[], long line);
float subtract_sky(float *buf, long axes[]);
double aperture(float *buf, long axes[], double centroid[], int radius);
void find_centroid(float *buf, long axes[], int start[], int r);
void dumparea(float *buf, long axes[], int start[], int r);

// return the mean of an n-point array excluding outliers. Done in place
float mean_removeoutliers(int n, float buf[], int noremove, int verbose);
float combine_minpix(int n, float buf[], int noremove);
float combine_maxpix(int n, float buf[], int noremove);
float centroid1d(float *buf, float sky, int rows, int wid,  double *cent, double r, int debug);
void cliplimits(float *buf, long axes[], float clip, float *min, float *max, int cliptype);
float median(float vals[], int n);
int stats(float inp[], long axes[2], float level[], int steps);
float skylevel(float vals[], long axes[2]);
float *get_pixel(float *image, long x, long y, long axes[2]);
float rms_sky(float *img, long axes[2]);
void add_noise(float *img, long axes[2], float sd);

int add_const(float *a, float constval, long axes[2]);
int normalize(float *buf, long axes[2]);
int image_sub(float *a, float *b, long axes[2], int norm);
int image_div(float *a, float *b, long axes[2]);
int image_limit(float *a, float limit, long axes[2]);
int handle_status(int *status, int warning, const char *msg);

float *img_filt(float *img, long axes[], int rad);

int checksize(long img1[2], long img2[2]);
