#include <FPT.h>
#include <D3d_matrix.h>


double Half_window_size = 350;
double Half_angle_degrees ;
double Tan_half_angle ;


double light_in_eye_space[3] ;
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 50 ;

double eye[3];
double zbuff[700][700];

/**===================================================================================================*/

int hyperboloid(double u, double v, double xyz[3]){
  xyz[0] = sqrt(1+(v*v))*cos(u);
  xyz[1] = v;
  xyz[2] = sqrt(1+(v*v))*sin(u);
}

int sgn(double val){
	if(val > 0) return 1;
	if(val < 0) return -1;
	return 0;
}


void makemat(double m[4][4], double minv[4][4], double sx, double sy, double sz, 
				double rz, double tx, double ty, double tz){
	int num = 0 ; 
	int tlist[7];
	double plist[7];
	tlist[num] = SX ; plist[num] =  sx ; num++ ;
	tlist[num] = SY ; plist[num] =  sy ; num++ ;
        tlist[num] = SZ ; plist[num] =  sz ; num++ ;

	tlist[num] = RZ ; plist[num] = rz ; num++;

	tlist[num] = TX ; plist[num] = tx ; num++ ;
	tlist[num] = TY ; plist[num] = ty ; num++ ;
	tlist[num] = TZ ; plist[num] = tz ; num++ ;

	D3d_make_movement_sequence_matrix(m,minv, num, tlist, plist) ;
}


/**===================================================================================================*/


int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3])
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// s = location of start of ray (probably the eye)
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]

// return 1 if successful, 0 if error
{

  double len ;
  double N[3] ; 
  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
  if (len == 0) return 0 ;
  N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;

  double E[3] ;
  E[0] = s[0] - p[0] ; 
  E[1] = s[1] - p[1] ; 
  E[2] = s[2] - p[2] ; 
  len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
  if (len == 0) return 0 ;
  E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
  double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;

  double L[3] ;
  L[0] = light_in_eye_space[0] - p[0] ; 
  L[1] = light_in_eye_space[1] - p[1] ; 
  L[2] = light_in_eye_space[2] - p[2] ; 
  len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
  if (len == 0) return 0 ;
  L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;





  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
     // this needs to occur BEFORE you possibly jump to LLL below




  double intensity ;
  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = AMBIENT ; 
    goto LLL ;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // eye and light on same side but normal pointing "wrong" way
    N[0] *= (-1.0) ;    N[1] *= (-1.0) ;    N[2] *= (-1.0) ; 
    NdotL *= (-1.0) ;
    NdotE *= (-1.0) ;   // don't use NdotE below, probably should eliminate this
  }


  // ignore Blinn's variant
  double R[3] ; // Reflection vector of incoming light
  R[0] = 2*NdotL*N[0] - L[0] ;
  R[1] = 2*NdotL*N[1] - L[1] ;
  R[2] = 2*NdotL*N[2] - L[2] ;

  double EdotR = E[0]*R[0] + E[1]*R[1] + E[2]*R[2] ;

  double diffuse ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = MAX_DIFFUSE*NdotL ; }

  double specular ;
  if (EdotR <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPECPOW) ;}

  // printf("%lf %lf\n",diffuse,specular) ;
  intensity = AMBIENT + diffuse + specular ;



 LLL : ;

  double f,g ;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse ;
    argb[0] = f * irgb[0] ;
    argb[1] = f * irgb[1] ;
    argb[2] = f * irgb[2] ;
  } else {
    f = (intensity - max_ambient_and_diffuse) / 
                           (1.0 - max_ambient_and_diffuse) ;
    g = 1.0 - f ;
    argb[0] = g * irgb[0] + f ;
    argb[1] = g * irgb[1] + f ;
    argb[2] = g * irgb[2] + f ;
  }

  return 1 ;
}


/**===================================================================================================*/






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

/**===================================================================================================*/


void plot(double ustart, double uend, double vstart, double vend, 
	double rgb[3], double m[4][4], double V[4][4] , int (*func)(double u, double v,  double xy[3])){
	
  double i,j;
  double point[3],  color[3];  
  double p[3], p2[3], normal[3];
  
	for(i=ustart; i<uend;i = i + .01){
		for(j=vstart; j<vend;j = j + .01){

			func(i,j,point);
			D3d_mat_mult_pt(point, m, point);
      double tempPoint[3];
      D3d_mat_mult_pt(tempPoint, V, point);

      if(tempPoint[2] < 0)  continue;
      else if(fabs(tempPoint[1]/tempPoint[2]) > Tan_half_angle)continue;
      else if(fabs(tempPoint[0]/tempPoint[2]) > Tan_half_angle)continue;

			func(i+.001,j,p);
			D3d_mat_mult_pt(p, m, p);

			func(i,j+.001,p2);
			D3d_mat_mult_pt(p2, m, p2);


      p[0] = p[0] - point[0];
      p[1] = p[1] - point[1];
      p[2] = p[2] - point[2];

      p2[0] = p2[0] - point[0];
      p2[1] = p2[1] - point[1];
      p2[2] = p2[2] - point[2];

			D3d_x_product(normal, p, p2);

			Light_Model (rgb, eye, point, normal, color);
  
			G_rgb(color[0], color[1], color[2]);

      D3d_mat_mult_pt(point, V, point);

      point[0] = ((Half_window_size * point[0])/(point[2]*Tan_half_angle)) + Half_window_size;
      point[1] = ((Half_window_size * point[1])/(point[2]*Tan_half_angle)) + Half_window_size;
      
      if(zbuff[(int)point[0]][(int)point[1]] > point[2] ){
        G_point(point[0],point[1]);
        zbuff[(int)point[0]][(int)point[1]] = point[2];
      }
			
		}
	}	
} 


/**===================================================================================================*/



int sphere(double u, double v, double xy[3]){
	xy[0] = sqrt(1-(v*v))*cos(u);
	xy[1] = v;
	xy[2] = sqrt(1-(v*v))*sin(u);
}

int cone(double u, double v, double xy[3]){
  xy[0] = (1-v)*cos(u);
  xy[1] = v;
  xy[2] = (1-v)*sin(u);
}

// 
void clear(){
	G_wait_key();
	G_rgb(0,0,0) ;
  G_clear();
  G_rgb(0,1,0);
    
}

int main(){

  char prefix[100];
  printf("Enter file extension\n");
  scanf("%s", prefix);
  int x,y;

	double rgb[3];
 
  AMBIENT = 0.2 ;
  MAX_DIFFUSE = 0.5 ;
  SPECPOW = 30 ;

	eye[0] = 0 ;
	eye[1] = 0 ;
	eye[2] = 0 ;

  light_in_eye_space[0] = 100;
  light_in_eye_space[1] = 100;
  light_in_eye_space[2] = -100;

	Half_window_size  = 300;
	Half_angle_degrees = 30;
	Tan_half_angle = tan(Half_angle_degrees*M_PI/180) ;


  G_init_graphics(Half_window_size*2,Half_window_size*2);
	G_rgb(0,0,0) ;
  G_clear();
  G_rgb(0,1,0);
  
  double t;
  double coi[3], up[3];
  double fnum = 0;
  double numframes = 1000;
//===========================================================================//
  int i = 0;
  // while (i < numframes) {
    for (x = 0; x < Half_window_size * 2; ++x)
      { 
        for (y = 0; y < Half_window_size *2; ++y)
          {
            zbuff[x][y] = 100000;
          }
      }
    i++;
    double V[4][4], Vi[4][4];

    t = 2.0 * M_PI *fnum/numframes ; // goes from 0 to 1

    eye[0] = 300*cos(2*M_PI*t) ; 
    eye[1] =  150*t ; 
    eye[2] =  150*sin(2*M_PI*t) ; 

    coi[0] =  0 ;
    coi[1] =  0 ; 
    coi[2] =  1 ;

    up[0]  = eye[0] ; 
    up[1]  = eye[1] + 1 ;
    up[2]  = eye[2] ; 

    D3d_make_identity (V) ;
    D3d_make_identity (Vi) ;

    D3d_view (V, Vi,  eye,coi,up) ;
//--------------------------------------------------------------------//

  double m[4][4], minv[4][4];

  rgb[0] = 0 ;
  rgb[1] = 1 ;
  rgb[2] = 1 ;

  // makemat(m, minv,50, 10, 10, 0, 0,50, 150);
  // plot(0 , 2*M_PI, -1,1, rgb,  m, V, sphere);

  D3d_make_identity(m);
  D3d_make_identity(minv);

  makemat(m, minv,1, 1, 1, 0, 0,0, 1);
  D3d_rotate_x(m,minv,M_PI/2);

  rgb[0] = 1;
  rgb[1] = 0;
  rgb[2] = 0;


    plot(0 , 2*M_PI, -1,1, rgb,  m,V, hyperboloid);
    G_wait_key();
//=============================================================================
  //   char filename[100];
  //   sprintf(filename, "%s%04d.xwd", prefix, i);
  //   G_save_image_to_file(filename);
  // fnum++;
  printf("image number %d\n", i);
  G_rgb(0,0,0);
  G_clear();
// }

//---------------------------------------------------------------------------

}
