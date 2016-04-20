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

		// if(fmod(sgn(i)*i,(end-start)/1000) < .0001){
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

		// 	G_line(point[0] , point[1] ,point[0] - (slope[0]*10), (point[1] ) - (slope[1]*10));
		// }
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

void getPoint(double point[3]){
	double coords[2];

		G_wait_click(coords);
		point[0] = coords[0];
		point[1] = coords[1];
		point[2] = 1;

		G_rectangle(coords[0] -2, coords[1]-2, 4,4);
}

void quadratic(double res[2], double a, double b, double c){
	res[0] = (-b + sqrt((b*b)-(4*a*c)))/(2*a);
	res[1] = (-b - sqrt((b*b)-(4*a*c)))/(2*a);
}

void copy(double *a , double * b , int n){
	int i;
	for(i=0;i<n;i++){
		a[i] = b[i];
	}
}

double reflection(double nor[3], double inc[3], double ref[3]){
  double m;
  m = 2*(nor[0]*inc[0] + nor[1]*inc[1] + nor[2]*inc[2]);
  
  ref[0] = m*nor[0] - inc[0];
  ref[1] = m*nor[1] - inc[1];
  ref[2] = m*nor[2] - inc[2];
}

double tracer(double m[5][4][4],double minv[5][4][4], double points[2][3], int n, int numRef){
    int i;
    int nfinal;
    double resFinal = 100000000;
    double pfinal[3];

    double tpoint[3], tpoint2[3];

	printf("reflections left = %d \n", numRef);
	if(numRef <=0){
		printf("no more relfections \n");
		return;
	}

   for(i=0; i < n; i++){

   		//convert points to object space
	    D3d_mat_mult_pt(tpoint, minv[i], points[0]);
	    D3d_mat_mult_pt(tpoint2, minv[i], points[1]);

	    //find slope of two points 
	    double p = tpoint2[0] - tpoint[0];
	    double q = tpoint2[1] - tpoint[1];

	    //find a,b,c for quadtratic equation
	    double a = pow(p,2) + pow(q,2);
	    double b = (2*tpoint[0]*p) + (2*tpoint[1] * q);
	    double c = pow(tpoint[0],2) + pow(tpoint[1],2) - 1;

	    double res[2];
	    quadratic(res, a, b ,c);

	    G_rgb(1,0,0);

	    printf("res 0 = %lf , res 1 = %lf\n",res[0], res[1] );
	    //find closest point 
	    if(res[0] < res[1] || res[0] <= 0){
	    	if(resFinal > res[0] && res[0] > 0){
			    pfinal[0] = tpoint[0] + (res[0]*p);
			    pfinal[1] = tpoint[1] + (res[0] * q);
			    pfinal[2] = 0;

			    nfinal = i;
			    resFinal = res[0];
			}

		}else{
		    if(resFinal > res[1] && res[1] > 0){
			    pfinal[0] = tpoint[0] + (res[1] * p);
			    pfinal[1] = tpoint[1] + (res[1] * q);
			    pfinal[2] = 0;

			    nfinal = i;
			    resFinal = res[1];
			}
		}
	}

	//return if not intersections
	if(resFinal > 100000){
		printf("no intersections\n");
		return;
	}
    double fnorm[3];

    //normal vector in world space based at orign
    fnorm[0] = 2*pfinal[0]*minv[nfinal][0][0] + 2*pfinal[1]*minv[nfinal][1][0];
    fnorm[1] = 2*minv[nfinal][0][1]*pfinal[0] + 2*minv[nfinal][1][1]*pfinal[1];
    fnorm[2] = 0 ;

    //draw rectangle where ray intersects circle
    D3d_mat_mult_pt(pfinal,m[nfinal],pfinal);

    //translate normal to intersection point and extend
    double tx = pfinal[0] + (1000*fnorm[0]);
    double ty = pfinal[1] + (1000*fnorm[1]);

    printf("tx = %lf , ty= %lf\n",tx,ty );

    //draw interseciton
    G_rectangle(pfinal[0],pfinal[1], 4, 4);
    //draw light vector
    G_line(points[0][0], points[0][1], pfinal[0], pfinal[1]);
    double inc[3], ref[3], nor[3];

    //caluclate light vector based at origin
    inc[0] = points[1][0] - pfinal[0];
    inc[1] = points[1][1] - pfinal[1];
    inc[2] = 0 ;
    make_unit_vector(inc);
    
    //calculate normal vector based at origin 
    nor[0] = tx - pfinal[0];
    nor[1] = ty - pfinal[1];
    nor[2] = 0 ;
    make_unit_vector(nor);

    //reflection 
    reflection(nor, inc,ref);

    //draw reclection vector
    G_line(100*ref[0]+pfinal[0], 100*ref[1]+pfinal[1], pfinal[0], pfinal[1]);
    
    points[0][0] = pfinal[0] + (0.001*ref[0]);
    points[0][1] = pfinal[1] + (0.001*ref[1]);
    points[0][2] = 0 ;
    

    points[1][0] = (pfinal[0] + 100*ref[0]);
    points[1][1] = (pfinal[1] + 100*ref[1]);
    points[1][2] = 0 ;

	tracer(m,minv,points, n, numRef-1);
}

int main(){
	G_init_graphics(700,700);
	G_rgb(0,0,0) ;
    G_clear();
    G_rgb(0,1,0);

    double m[5][4][4], minv[5][4][4];

    makemat(m[0], minv[0], 20, 350, 0, 20,350);
    plot(0 * M_PI, M_PI * 2, m[0], circle);

    makemat(m[1], minv[1], 350, 20, 0, 350,20);
    plot(0 * M_PI, M_PI * 2, m[1], circle);

    makemat(m[2], minv[2], 20, 350, 0, 680, 350);
    plot(0 * M_PI, M_PI * 2, m[2], circle);

    makemat(m[3], minv[3], 350, 20, 0, 350,680);
    plot(0 * M_PI, M_PI * 2, m[3], circle);

    makemat(m[4], minv[4], 70, 40, 0, 350,350);
    plot(0 * M_PI, M_PI * 2, m[4], circle);
   	
   	
   	double points[2][3];
   	getPoint(points[0]);
   	getPoint(points[1]);

    tracer(m,minv,points, 5,10);

    G_wait_key();
}