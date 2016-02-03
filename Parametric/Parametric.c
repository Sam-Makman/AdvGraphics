#include <FPT.h>
#include <D3d_matrix.h>

double absd(double d){
	if(d < 0){
		return d * -1;
	}
	return d;
}

double fixPi(double pi){
	if(pi >=  0){
		pi = fmod(pi, 2*  M_PI);
	}else{
		pi = fmod(-pi, M_PI);
		pi =  (2* M_PI) - pi;
	}
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

void findPoint(double point[3], double p[3], double slope, double u){
	double ident[4][4];
	D3d_make_identity(ident);
	D3d_mat_mult_pt(point,  ident,  p) ;
	double invslope = -1/slope;
	double c = point[1] - (invslope * point[0]);
	if(fixPi(u)==M_PI/2){
		point[1] +=100;
	}
	else if(fixPi(u) >= .5 * M_PI && fixPi(u) <= 1.5 * M_PI){
		point[0] = point[0] - 10;
		point[1] = (point[0] * invslope ) + c;
	}else{
		point[0] = point[0] + 10;
		point[1] = (point[0] * invslope ) + c;
	}
}

void circle(double points[][3], 
	int n,
	double start, double end,
	double sx, double sy, 
	double rz,
	double tx, double ty){
	int i;
	double inc  = (end - start )/ n;

	double m[4][4], minv[4][4], ident[4][4];
	D3d_make_identity(ident);
	makemat(m, minv, sx, sy, rz, tx, ty); 

	for(i = 0; i< n ; i++){
		double u = start + (i*inc);
		points[i][0] = cos(u);
		points[i][1] = sin(u);
		points[i][2] = 0;

		D3d_mat_mult_pt(points[i],  m,  points[i]) ;

		if(i%(n/50) == 0){
			double slope = (sy * cos(u))/(-sx * sin(u));
			double point[3];
			findPoint(point, points[i], slope, u);
			G_line(points[i][0], points[i][1], point[0], point[1]);
		}
	}
}

void sum4(double points[][3], 
	int n,
	double start, double end,
	double sx, double sy, 
	double rz, 
	double tx, double ty){

	int i;
	double inc  = (end - start )/ n;

	double m[4][4], minv[4][4];
	makemat(m, minv, sx, sy, rz, tx, ty); 
	
	for(i = 0; i< n ; i++){
		points[i][2] = 0;
		double u = start + (i*inc);
		double compare = fixPi(u);
		if(compare < M_PI/2 && compare > 0 ){
		points[i][0] = sqrt(absd(cos(u)));
		points[i][1] = sqrt(absd(sin(u)));
		}else if(compare < M_PI){
		points[i][0] = -sqrt(absd(cos(u)));
		points[i][1] = sqrt(absd(sin(u)));
		}else if(compare < (M_PI *3)/2 ){
		points[i][0] = -sqrt(absd(cos(u)));
		points[i][1] = -sqrt(absd(sin(u)));
		}else{
		points[i][0] = sqrt(absd(cos(u)));
		points[i][1] = -sqrt(absd(sin(u)));
		}
		D3d_mat_mult_pt(points[i],  m,  points[i]) ;

		if(i%(n/50) == 0){
			double dx = -sin(u) / (2 * sqrt(fabs(cos(u))));
			double dy = cos(u)/(2*sqrt(fabs(sin(u))));
			double slope = (sy*dy)/(sx*dx);
			double point[3];
			findPoint(point, points[i], slope, u);
			G_line(points[i][0], points[i][1], point[0], point[1]);
		}

		
	}
}

void square(double points[][3], 
	int n,
	double start, double end,
	double sx, double sy, 
	double rz,
	double tx, double ty){

	int i;
	double inc  = (end - start )/ n;

	double m[4][4], minv[4][4];
	makemat(m, minv, sx, sy, rz, tx, ty); 

	
	for(i = 0; i< n ; i++){
		double u = start + (i*inc);
		double compare = fixPi(u);
		double slope;
		if(compare < .5 * M_PI){
		points[i][0] = pow(cos(u),2);
		points[i][1] = pow(sin(u),2);
		slope = -sy/sx;
		}else if(compare < M_PI ){
		points[i][0] = -pow(cos(u),2);
		points[i][1] = pow(sin(u),2);
		slope = sy/sx;
		}else if(compare < (M_PI *3)/2 ){
		points[i][0] = -pow(cos(u),2);
		points[i][1] = -pow(sin(u),2);
		slope = -sy/sx;
		}else{
		points[i][0] = pow(cos(u),2);
		points[i][1] = -pow(sin(u),2);
		slope = sy/sx;
		}
		points[i][2] = 0;
		D3d_mat_mult_pt(points[i],  m,  points[i]) ;

		if(i%(n/50) == 0){
			double point[3];
			findPoint(point, points[i], slope, u);
			G_line(points[i][0], points[i][1], point[0], point[1]);
		}
	}
}

void astroid(double points[][3], 
	int n,
	double start, double end,
	double sx, double sy, 
	double rz,
	double tx, double ty){

	int i;
	double inc  = (end - start )/ n;

	double m[4][4], minv[4][4];
	makemat(m, minv, sx, sy, rz, tx, ty); 
	
	for(i = 0; i< n ; i++){
		double u = start + (i*inc);
		if(fixPi(u) < M_PI/2 ){
		points[i][0] = pow(cos(u),4);
		points[i][1] = pow(sin(u),4);
		}else if(fixPi(u) < M_PI ){
		points[i][0] = -pow(cos(u),4);
		points[i][1] = pow(sin(u),4);
		}else if(fixPi(u) < (M_PI *3)/2 ){
		points[i][0] = -pow(cos(u),4);
		points[i][1] = -pow(sin(u),4);
		}else{
		points[i][0] = pow(cos(u),4);
		points[i][1] = -pow(sin(u),4);
		}
		points[i][2] = 0;
		D3d_mat_mult_pt(points[i],  m,  points[i]) ;
	
	if(i%(n/50) == 0){
			double dx = -4*(pow(cos(u),3))*sin(u);
			double dy = 4*(pow(sin(u),3))*cos(u);
			double slope = (sy*dy)/(sx*dx);
			double point[3];
			findPoint(point, points[i], slope, u);
			G_line(points[i][0], points[i][1], point[0], point[1]);
		}
	}
}

void hyperbola(double points[][3], 
	int n,
	double start, double end,
	double sx, double sy, 
	double rz,
	double tx, double ty){

	int i;
	double inc  = (end - start )/ n;

	double m[4][4], minv[4][4];
	makemat(m, minv, sx, sy, rz, tx, ty); 
	
	for(i = 0; i< n ; i++){
		double u = start + (i*inc);
		points[i][0] = cosh(u);
		points[i][1] = sinh(u);
		points[i][2] = 0;
		D3d_mat_mult_pt(points[i],  m,  points[i]) ;
		if(i%(n/50) == 0){
			double slope = (sy*cosh(u))/(sx*sinh(u));
			double point[3];double ident[4][4];
			D3d_make_identity(ident);
			D3d_mat_mult_pt(point,  ident,  points[i]) ;
			double invslope = -1/slope;
			double c = point[1] - (invslope * point[0]);
			point[0] = point[0] + 10;
			point[1] = (point[0] * invslope ) + c;		
			G_line(points[i][0], points[i][1], point[0], point[1]);
			
		}
	}
}

void parabola(double points[][3], 
	int n,
	double start, double end,
	double sx, double sy, 
	double rz,
	double tx, double ty){

	int i;
	double inc  = (end - start )/ n;

	double m[4][4], minv[4][4];
	makemat(m, minv, sx, sy, rz, tx, ty); 
	
	for(i = 0; i< n ; i++){
		double u = start + (i*inc);
		points[i][0] = u;
		points[i][1] = pow(u, 2);
		points[i][2] = 0;
		D3d_mat_mult_pt(points[i],  m,  points[i]) ;

		if(i%(n/50) == 0){
			double slope = (sy*2*u)/sx;
			double point[3];
			findPoint(point, points[i], slope, u);
			G_line(points[i][0], points[i][1], point[0], point[1]);
		}
	}
}

void lemon(double points[][3], 
	int n,
	double start, double end,
	double sx, double sy, 
	double rz,
	double tx, double ty){

	int i;
	double inc  = (end - start )/ n;

	double m[4][4], minv[4][4];
	makemat(m, minv, sx, sy, rz, tx, ty); 
	
	for(i = 0; i< n ; i++){
		double u = start + (i*inc);
		points[i][0] = pow(cos(u),3);
		points[i][1] = sin(u);
		points[i][2] = 0;
		D3d_mat_mult_pt(points[i],  m,  points[i]) ;

		if(i%(n/50) == 0){
			double slope = ((cos(u))/(3*pow(cos(u),2)*-sin(u)));
			double point[3];
			findPoint(point, points[i], slope, u);
			G_line(points[i][0], points[i][1], point[0], point[1]);
		}
	}
}

void drawpoints(double points[][3], int n){
	int i;
	for(i = 0; i< n; i++){
		G_point(points[i][0], points[i][1]);
	}
}

// void plot(double start, double end, double (*func)(double)(*double)){
	
// }

int main(){

	G_init_graphics(700,700);
	G_rgb(0,0,0) ;
    G_clear();
    G_rgb(0,1,0);

    int n = 100000;
	double points[n][3];
	

	circle(points,n, .25 * M_PI , 1.5 * M_PI, 50,100, 0 , 300,500);
	drawpoints(points,n);
	G_wait_key();

	G_rgb(0,0,0) ;
    G_clear();
    G_rgb(0,1,0);

    sum4(points,n, .5 * M_PI, 1.75 * M_PI, 30,60, 0, 300,300);
	drawpoints(points,n);
	G_wait_key();

	G_rgb(0,0,0) ;
    G_clear();
    G_rgb(0,1,0);

    square(points,n, 0 , 2 * M_PI, 150,100, 0, 500,500);
	drawpoints(points,n);
	G_wait_key();

	G_rgb(0,0,0) ;
    G_clear();
    G_rgb(0,1,0);

    astroid(points,n, 0 , 2 * M_PI, 180,140,45, 500,300);
	drawpoints(points,n);
	G_wait_key();

	G_rgb(0,0,0) ;
    G_clear();
    G_rgb(0,1,0);

    hyperbola(points,n, -1 , 2 , 100,100, 0, 250,250);
	drawpoints(points,n);
	G_wait_key();

	G_rgb(0,0,0) ;
    G_clear();
    G_rgb(0,1,0);

    parabola(points,n, -1 , 2 , 150,50, 60, 250,250);
	drawpoints(points,n);
	G_wait_key();

	G_rgb(0,0,0) ;
    G_clear();
    G_rgb(0,1,0);

    lemon(points,n, -1 , 2 * M_PI, 150,150, 60, 600,150);
	drawpoints(points,n);
	G_wait_key();
}