#include <shape.h>
#include <math.h>

Shape new_shape(double x, double y, double radius, int numChildren){
	Shape shape;
	shape.x = x;
	shape.y = y;
	shape.radius = radius;
	shape.numChildren = numChildren;
}

Shape get_child(Shape shape, int childNumber ){
	Shape child;
	childNumber = childNumber%shape.numChildren;
	child.x = 
	child.y = 
	child.radius = shape.radius/2;
	child.numChildren = shape.numChildren;

}