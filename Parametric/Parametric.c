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

void findPoint(double m[4][4], double point[3], double p[3], double slope, double u){
	double ident[4][4];
	D3d_make_identity(ident);
	D3d_mat_mult_pt(point,  ident,  p) ;
	double invslope = -1/slope;
	double c = point[1] - (invslope * point[0]);
	if(fixPi(u)==M_PI/2){
		point[0] = point[0];
		point[1] +=100;
	}
	else if(fixPi(u) >= .5 * M_PI && fixPi(u) <= 1.5 * M_PI){
		point[0] = point[0] - 10;
		point[1] = (point[0] * invslope ) + c;
	}else{
		point[0] = point[0] + 10;
		point[1] = (point[0] * invslope ) + c;
	}
			
	// G_line(points[i][0], points[i][1], point[0], point[1]);
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

		if(i%(n/100) == 0){
			double slope = (sy * cos(u))/(-sx * sin(u));
			double point[3];
			findPoint(m,point, points[i], slope, u);
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
		double compare = fixPi(start + (i*inc));
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
		if(i%(n/100) == 0){
			double slope = ((sy/2)*(pow(sin(u),-.5))*cos(u))/(((sx/2)*(pow(cos(u),-.5))*-cos(u)));
			double point[3];
			findPoint(m,point, points[i], slope, u);
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
		double compare = fixPi(start + (i*inc));
		if(compare < .5 * M_PI){
		points[i][0] = pow(cos(start + (i*inc)),2);
		points[i][1] = pow(sin(start + (i*inc)),2);
		}else if(compare < M_PI ){
		points[i][0] = -pow(cos(start + (i*inc)),2);
		points[i][1] = pow(sin(start + (i*inc)),2);
		}else if(compare < (M_PI *3)/2 ){
		points[i][0] = -pow(cos(start + (i*inc)),2);
		points[i][1] = -pow(sin(start + (i*inc)),2);
		}else{
		points[i][0] = pow(cos(start + (i*inc)),2);
		points[i][1] = -pow(sin(start + (i*inc)),2);
		}
		points[i][2] = 0;
		D3d_mat_mult_pt(points[i],  m,  points[i]) ;
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
		if(fixPi(start + (i*inc)) < M_PI/2 ){
		points[i][0] = pow(cos(start + (i*inc)),4);
		points[i][1] = pow(sin(start + (i*inc)),4);
		}else if(fixPi(start + (i*inc)) < M_PI ){
		points[i][0] = -pow(cos(start + (i*inc)),4);
		points[i][1] = pow(sin(start + (i*inc)),4);
		}else if(fixPi(start + (i*inc)) < (M_PI *3)/2 ){
		points[i][0] = -pow(cos(start + (i*inc)),4);
		points[i][1] = -pow(sin(start + (i*inc)),4);
		}else{
		points[i][0] = pow(cos(start + (i*inc)),4);
		points[i][1] = -pow(sin(start + (i*inc)),4);
		}
		points[i][2] = 0;
		D3d_mat_mult_pt(points[i],  m,  points[i]) ;
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
	
	for(i = 0; i< n ; i+=2){
		points[i][0] = cosh(start + (i*inc));
		points[i][1] = sinh(start + (i*inc));
		points[i][2] = 0;
		D3d_mat_mult_pt(points[i],  m,  points[i]) ;
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
		points[i][0] = start + (i*inc);
		points[i][1] = pow(start + (i*inc), 2);
		points[i][2] = 0;
		D3d_mat_mult_pt(points[i],  m,  points[i]) ;
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
		points[i][0] = pow(cos(start + (i*inc)),3);
		points[i][1] = sin(start + (i*inc));
		points[i][2] = 0;
		D3d_mat_mult_pt(points[i],  m,  points[i]) ;
	}
}

void drawpoints(double points[][3], int n){
	int i;
	for(i = 0; i< n; i++){
		G_point(points[i][0], points[i][1]);
	}
}

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

    sum4(points,n, .5 * M_PI , 1.75 * M_PI, 30,60, 0, 300,300);
	drawpoints(points,n);
	G_wait_key();

	// G_rgb(0,0,0) ;
 //    G_clear();
 //    G_rgb(0,1,0);

 //    square(points,n, 0 , 2 * M_PI, 150,100, 0, 500,500);
	// drawpoints(points,n);
	// G_wait_key();

	// G_rgb(0,0,0) ;
 //    G_clear();
 //    G_rgb(0,1,0);

 //    astroid(points,n, 0 , 2 * M_PI, 80,40,45, 500,300);
	// drawpoints(points,n);
	// G_wait_key();

	// G_rgb(0,0,0) ;
 //    G_clear();
 //    G_rgb(0,1,0);

 //    hyperbola(points,n, -1 , 2 , 100,100, 0, 250,250);
	// drawpoints(points,n);
	// G_wait_key();

	// G_rgb(0,0,0) ;
 //    G_clear();
 //    G_rgb(0,1,0);

 //    parabola(points,n, -1 , 2 , 150,50, 60, 250,250);
	// drawpoints(points,n);
	// G_wait_key();

	// G_rgb(0,0,0) ;
 //    G_clear();
 //    G_rgb(0,1,0);

 //    lemon(points,n, -1 , 2 * M_PI, 150,150, 60, 600,150);
	// drawpoints(points,n);
	// G_wait_key();
}