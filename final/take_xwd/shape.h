

typedef void (*Operation)(double vals[4], double x, double y, double z, double rad, int pos, int total);


typedef struct shape
{
	double x;
	double y;
	double z;
	Operation calcChild;
	double radius;
	double color[3];
	double ref;
	int numChildren;

}SHAPE;

	
SHAPE new_shape(double x, double y, double z, double radius, Operation calcChild, double color[3], double ref, int numChildren);

SHAPE get_child(SHAPE shape, int childNumber );

void make_matrix(SHAPE shape, double m[4][4], double minv[4][4]);