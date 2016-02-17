#include <FPT.h>
#include <D3d_matrix.h>

int sgn(double val){
	if(val > 0) return 1;
	if(val < 0) return -1;
	return 0;
}

void makemat(double m[4][4], double minv[4][4], double sx, double sy, double rz, double tx, double ty){
	int num = 0 ; 
	int tlist[5];
	double plist[5];
	tlist[num] = SX ; plist[num] =  sx ; num++ ;
	tlist[num] = SY ; plist[num] =  sy ; num++ ;
	tlist[num] = RZ ; plist[num] = rz ; num++;
	tlist[num] = TX ; plist[num] = tx ; num++ ;
	tlist[num] = TY ; plist[num] = ty ; num++ ;

	D3d_make_movement_sequence_matrix(m,minv, num, tlist, plist) ;
}
double dist(double p[3] , double p2[3]){
	double dx = p[0] - p2[0];
	double dy = p[1] - p2[1];
	return sqrt((dx*dx) + (dy*dy));

}

double dot_product(double a[3], double b[3]){
	return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
}

int make_unit_vector(double v[3]){
	double length = sqrt((v[0]*v[0])+(v[1]*v[1])+(v[2]*v[2]));
			v[0] = v[0]/length;
			v[1] = v[1]/length;
			v[2] = v[2]/length;
}

void plot(double start, double end, double m[4][4], int (*func)(double rad, double xy[2])){
	double i;
	double ident[4][4];
	D3d_make_identity(ident);
	for(i=start; i<end;i = i + .001){
		double point[3];
		func(i,point);
		D3d_mat_mult_pt(point, m, point);
		G_point(point[0],point[1]);

		if(fmod(sgn(i)*i,(end-start)/1000) < .0001){
			double center[3];
			center[0] = 0; center[1] = 0; center[2] = 0;

			D3d_mat_mult_pt(center, m, center);
			double p[3];
			func(i+.001,p);
			D3d_mat_mult_pt(p, m, p);
			double slope[3];
			double z[3];
			z[0] = 0;
			z[1] = 0;
			z[2] = 1;
			p[0] = point[0] - p[0];
			p[1] = point[1] - p[1];
			p[2] = point[2] - p[2];

			make_unit_vector(p);

			D3d_x_product(slope, p, z);

			make_unit_vector(slope);

			G_line(point[0] , point[1] ,point[0] - (slope[0]*10), (point[1] ) - (slope[1]*10));
		}
	}	
} 

int circle(double rad, double xy[3]){
	xy[0] = cos(rad);
	xy[1] = sin(rad);
	xy[2] = 0;
}

int sum4(double rad, double xy[3]){
	xy[0] = sgn(cos(rad)) * sqrt(fabs(cos(rad)));
	xy[1] = sgn(sin(rad)) * sqrt(fabs(sin(rad)));
	xy[2] = 0;
}

int square(double rad, double xy[3]){
	xy[0] = sgn(cos(rad)) * pow(cos(rad),2);
	xy[1] = sgn(sin(rad)) * pow(sin(rad),2);
	xy[2] = 0;
}

int astroid(double rad, double xy[3]){
	xy[0] = sgn(cos(rad)) * pow(cos(rad),4);
	xy[1] = sgn(sin(rad)) * pow(sin(rad),4);
	xy[2] = 0;
}

int hyperbola(double rad, double xy[3]){
	xy[0] = cosh(rad);
	xy[1] = sinh(rad);
	xy[2] = 0;
}

int parabola(double rad, double xy[3]){
	xy[0] = rad;
	xy[1] = rad*rad;
	xy[2] = 0;
}

int lemon(double rad, double xy[3]){
	xy[0] = pow(cos(rad),3);
	xy[1] = sin(rad);
	xy[2] = 0;
}
void clear(){
	G_wait_key();
	G_rgb(0,0,0) ;
    G_clear();
    G_rgb(0,1,0);
    
}
int main(){
	G_init_graphics(700,700);
	G_rgb(0,0,0) ;
    G_clear();
    G_rgb(0,1,0);

    double m[4][4], minv[4][4];

    makemat(m, minv, 50, 100, 0, 300,500);
    plot(0 * M_PI, M_PI * 2, m, circle);
    // clear();

	makemat(m, minv, 30, 60, 0, 300,300);
    plot(.5 * M_PI, M_PI * 1.75, m, sum4);
    // clear();

    makemat(m, minv, 150, 100, 0, 500,500);
    plot(0, M_PI * 2, m, square);
    // clear();

    makemat(m, minv, 80, 40, 45, 500,300);
    plot(0, M_PI * 2, m, astroid);
    // clear();

    makemat(m, minv, 100, 100, 0, 250,250);
    plot(-1, 2.123  , m, hyperbola);
    // clear();

    makemat(m, minv, 150, 50, 60, 250,250);
    plot(-1, 2.123, m, parabola);
    // clear();

    makemat(m, minv, 150, 150, 60, 600,150);
    plot(0, M_PI * 2, m, lemon);
    // clear();


    G_wait_key();
}