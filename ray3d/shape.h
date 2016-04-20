typedef struct Shape
{
	double x;
	double y;
	double radius;
	int numChildren;

};

Shape new_shape(double x, double y, double radius, int numChildren);

Shape get_child(Shape shape, int childNumber );