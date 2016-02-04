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

void plot(double start, double end, double m[4][4], int (*func)(double rad, double xy[2])){
	double i;
	double ident[4][4];
	D3d_make_identity(ident);
	for(i=start; i<end;i = i + .001){
		double points[3];
		func(i,points);
		D3d_mat_mult_pt(points, m, points);
		G_point(points[0],points[1]);

		if(fmod(sgn(i)*i,(end-start)/1000) < .0001){
			double center[3];
			center[0] = 0; center[1] = 0; center[2] = 0;

			D3d_mat_mult_pt(center, m, center);
			double p[3];
			func(i+.001,p);
			D3d_mat_mult_pt(p, m, p);
			double slope = -1/((p[1] - points[1])/(p[0]-points[0]));
			double tp[3];
			D3d_mat_mult_pt(tp, ident, points);
			tp[0] +=  10;
			tp[1] = tp[1] + (slope*10);
			int sign =1;
			double dir = (points[0] - center[0]) / (points[1] - center[1]);
			if(dir * slope > 0){
				G_line(points[0] , points[1] ,points[0] + 10, (points[1] ) + (slope*10));
				// printf("x = %lf y = %lf\n", tp[0], tp[1]);
				// printf("further\n");
			}else{
				G_line(points[0] , points[1] ,points[0] - 10, points[1]  - (slope*10));
				// printf("x = %lf y = %lf\n", tp[0], tp[1]);
				// printf("closer\n");
			}
		}
	}	
} 

int circle(double rad, double xy[3]){
	xy[0] = cos(rad);
	xy[1] = sin(rad);
}

int sum4(double rad, double xy[3]){
	xy[0] = sgn(cos(rad)) * sqrt(fabs(cos(rad)));
	xy[1] = sgn(sin(rad)) * sqrt(fabs(sin(rad)));
}

int square(double rad, double xy[3]){
	xy[0] = sgn(cos(rad)) * pow(cos(rad),2);
	xy[1] = sgn(sin(rad)) * pow(sin(rad),2);
}

int astroid(double rad, double xy[3]){
	xy[0] = sgn(cos(rad)) * pow(cos(rad),4);
	xy[1] = sgn(sin(rad)) * pow(sin(rad),4);
}

int hyperbola(double rad, double xy[3]){
	xy[0] = cosh(rad);
	xy[1] = sinh(rad);
}

int parabola(double rad, double xy[3]){
	xy[0] = rad;
	xy[1] = rad*rad;
}

int lemon(double rad, double xy[3]){
	xy[0] = pow(cos(rad),3);
	xy[1] = sin(rad);
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
    plot(.24 * M_PI, M_PI * 1.5, m, circle);
    clear();

	makemat(m, minv, 30, 60, 0, 300,300);
    plot(.5 * M_PI, M_PI * 1.75, m, sum4);
    clear();

    makemat(m, minv, 150, 100, 0, 500,500);
    plot(0, M_PI * 2, m, square);
    clear();

    makemat(m, minv, 80, 40, 45, 500,300);
    plot(0, M_PI * 2, m, astroid);
    clear();

    makemat(m, minv, 100, 100, 0, 250,250);
    plot(-1, 2  , m, hyperbola);
    clear();

    makemat(m, minv, 150, 50, 60, 250,250);
    plot(-1, 2 , m, parabola);
    clear();

    makemat(m, minv, 150, 150, 60, 600,150);
    plot(0, M_PI * 2, m, lemon);
    clear();


    G_wait_key();
}