#ifndef SHAPE_TOOLS
#define SHAPE_TOOLS

typedef
struct{
	double ambient;
	double max_diffuse;
	double specpow;
	double *location;
} LIGHT;

typedef
struct{
	double ustart;
	double uend;
	double vstart;
	double vend;
	int (*paint)(double u, double v, double rgb[3]);
	int (*point)(double u, double v,  double xy[3]);
	double  m[4][4];
	double  v[4][4];
} SHAPE;
 

int print_shape( SHAPE  *shape);
int init_light( LIGHT *light, double ambient, double diffuse, double specpow, double location[3]);

int init_shape( SHAPE *shape, double ustart, double uend, double vstart, double vend,
		int (*paint)(double u, double v, double rgb[3]), int (*point)(double u, double v,  double xy[3]),
		double m[4][4], 
		double v[4][4]
 		);

double * plot(double eye[3],  LIGHT *light,  SHAPE *shape, int Half_window_size);

int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3], LIGHT *light);



//shapes 
int sphere(double u, double v, double xy[3]);
int hyperboloid(double u, double v, double xy[3]);


//patterns 

int checkers(double u, double v, double rgb[3]);
int stripe(double u, double v, double rgb[3]);
int smooth(double u, double v, double rgb[3]);
int horizontalStripe(double u, double v, double rgb[3]);


//tools
void makemat(double m[4][4], double minv[4][4], double sx, double sy, double sz, 
        double rx,double ry,double rz, 
        double tx, double ty, double tz);


// double offest(double * array, double x, double y, double z);
void print_point(double p[3]);

#endif