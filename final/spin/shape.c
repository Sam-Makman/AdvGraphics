#include <math.h>
#include "shape.h"
#include "D3d_matrix.h"



SHAPE new_shape(double x, double y, double z, double radius, Operation calcChild, double color[3], double ref, int numChildren){
	SHAPE shape;
	shape.x = x;
	shape.y = y;
	shape.z = z;
	shape.radius = radius;
	shape.calcChild = calcChild;
	shape.color[0] = color[0];
	shape.color[1] = color[1];
	shape.color[2] = color[2];
	shape.ref = ref;
	shape.numChildren = numChildren;
	return shape;	
}

SHAPE get_child(SHAPE shape, int childNumber ){
	SHAPE child;
	childNumber = childNumber%shape.numChildren;
	double vals[4];

	//gets x,y,z,radius of child shape 
	shape.calcChild(vals, shape.x, shape.y, shape.z, shape.radius, childNumber, shape.numChildren);

	child.x = vals[0] ;
	child.y = vals[1] ;
	child.z = vals[2];
	child.color[0] = shape.color[1];
	child.color[1] = shape.color[2];
	child.color[2] = shape.color[0];
	child.ref = shape.ref;
	child.radius = vals[3];
	child.calcChild = shape.calcChild;
	child.numChildren = shape.numChildren;
	return child;

}

void make_matrix(SHAPE shape, double m[4][4], double minv[4][4]){
	int num = 0 ; 
	int tlist[6];
	double plist[6];

	tlist[num] = SX ; plist[num] =  shape.radius ; num++ ;
	tlist[num] = SY ; plist[num] =  shape.radius ; num++ ;
  	tlist[num] = SZ ; plist[num] =  shape.radius ; num++ ;

	tlist[num] = TX ; plist[num] = shape.x ; num++ ;
	tlist[num] = TY ; plist[num] = shape.y ; num++ ;
	tlist[num] = TZ ; plist[num] = shape.z ; num++ ;

	D3d_make_movement_sequence_matrix(m,minv, num, tlist, plist) ;

}

