#include <shape_tools.h>
#include <math.h>
#include <D3d_matrix.h>


//====================================================================================================
//patterns 

int checkers(double u, double v, double rgb[3]){
  double a = floor(7*u);
  double b = floor(7*v);

  rgb[0] = fmod(fabs((a+b)),2);
  rgb[1] = 0;
  rgb[2] = 1;
}

int stripe(double u, double v, double rgb[3]){
  rgb[0] = 1 - fmod(u,1);
  rgb[1] = 1 - fmod(v,1);
  rgb[2] = 1;

}

int smooth(double u, double v, double rgb[3]){
  rgb[0] = 1 - (u/2*M_PI);
  rgb[1] = 1 - (v+1/2);
  rgb[2] = 1;  
}

int horizontalStripe(double u, double v, double rgb[3]){
  double a = floor(7*v);

  rgb[0] =fmod(fabs(a),2) ;
  rgb[1] = 0;
  rgb[2] = 0;  
}

//========================================================================================================

//shapes 

int hyperboloid(double u, double v, double xyz[3]){
  xyz[0] = sqrt(fabs(1.0+(v*v)))*cos(u);
  xyz[1] = v;
  xyz[2] = sqrt(fabs(1.0+(v*v)))*sin(u);
}


int sphere(double u, double v, double xy[3]){
  xy[0] = sqrt(fabs(1.0-(v*v)))*cos(u);
  xy[1] = v;
  xy[2] = sqrt(fabs(1.0-(v*v)))*sin(u);
  // printf("u = %lf , v = %lf\n",u, v );
  // printf("%lf , %lf, %lf \n",xy[0], xy[1], xy[2] );
}

//========================================================================================================
//tools
void makemat(double m[4][4], double minv[4][4], double sx, double sy, double sz, 
        double rx,double ry,double rz, 
        double tx, double ty, double tz){

  int num = 0 ; 
  int tlist[9];
  double plist[9];
  tlist[num] = SX ; plist[num] =  sx ; num++ ;
  tlist[num] = SY ; plist[num] =  sy ; num++ ;
  tlist[num] = SZ ; plist[num] =  sz ; num++ ;

  tlist[num] = RX ; plist[num] = rx ; num++;
  tlist[num] = RY ; plist[num] = ry ; num++;
  tlist[num] = RZ ; plist[num] = rz ; num++;

  tlist[num] = TX ; plist[num] = tx ; num++ ;
  tlist[num] = TY ; plist[num] = ty ; num++ ;
  tlist[num] = TZ ; plist[num] = tz ; num++ ;

  D3d_make_movement_sequence_matrix(m,minv, num, tlist, plist) ;
}

int print_shape(SHAPE *shape){
  printf("Ustart %lf. Uend %lf\n", shape->ustart, shape->uend);
  printf("Vstart %lf. Vend %lf\n", shape->vstart, shape->vend);
  printf("m\n");
  D3d_print_mat(shape->m);
  printf("v\n");
  D3d_print_mat(shape->v);
}

void print_point(double p[3]){
  printf("x = %lf , y = %lf, z = %lf\n", p[0], p[1], p[2] );
}

//================================================================================================
int init_light( LIGHT *light, double ambient, double diffuse, double specpow, double location[3]){

	light->ambient = ambient;
	light->max_diffuse = diffuse;
	light->specpow = specpow;
	light->location = location;
}

int init_shape( SHAPE *shape, double ustart, double uend, double vstart, double vend,
		int (*paint)(double u, double v, double rgb[3]), int (*point)(double u, double v,  double xy[3]),
		double m[4][4], 
		double v[4][4]
 		){

	shape->ustart = ustart;
	shape->uend = uend;
	shape->vstart = vstart;
	shape->vend = vend;
	shape->paint = paint;
	shape->point = point;
  int i,j;
  for(i=0; i<4;i++){
    for(j=0;j<4;j++){
      shape->m[i][j] = m[i][j];
      shape->v[i][j] = v[i][j];
    }
  }

}



//returns an array  [half_window_size * 2][half_window_size * 2][4]
//
//
double * plot(double eye[3],  LIGHT *light,  SHAPE *shape, int Half_window_size){
	double i,j;
  	double point[3],  color[3];  
  	double p[3], p2[3], normal[3];
    int length = Half_window_size*2 * Half_window_size*2 * 4;
  	double frame[length];
  	double Tan_half_angle = tan((M_PI*25)/180);

    int f = 0;
    int g = 0;
    for(f=0; f< Half_window_size*2; f++){
      for(g=0; g< Half_window_size*2; g++){
          frame[(int)((f*(Half_window_size*2)*4) + (g*4) + 3 )] = 100000;
      }
    }

  	for(i=shape->ustart; i < shape->uend;i = i + .005){
    	for(j=shape->vstart; j<shape->vend;j = j + .005){

      shape->point(i,j,point);
      D3d_mat_mult_pt(point, shape->m, point);

      double tempPoint[3];
      D3d_mat_mult_pt(tempPoint, shape->v, point);

      if(tempPoint[2] < 0)  continue;
      else if(fabs(tempPoint[1]/tempPoint[2]) > Tan_half_angle)continue;
      else if(fabs(tempPoint[0]/tempPoint[2]) > Tan_half_angle)continue;
      shape->point(i+.001,j,p);
      D3d_mat_mult_pt(p, shape->m, p);

      shape->point(i,j+.001,p2);
      D3d_mat_mult_pt(p2, shape->m, p2);


      p[0] = p[0] - point[0];
      p[1] = p[1] - point[1];
      p[2] = p[2] - point[2];

      p2[0] = p2[0] - point[0];
      p2[1] = p2[1] - point[1];
      p2[2] = p2[2] - point[2];

      D3d_x_product(normal, p, p2);

      double rgb[3];
      shape->paint(i,j,rgb);
      Light_Model (rgb, eye, point, normal, color, light);

      D3d_mat_mult_pt(point, shape->v, point);

      point[0] = ((Half_window_size * point[0])/(point[2]*Tan_half_angle)) + Half_window_size;
      point[1] = ((Half_window_size * point[1])/(point[2]*Tan_half_angle)) + Half_window_size;
      if(frame[(int)((point[0]* Half_window_size*2*4) + (point[1] * 4 ) + 3)] > point[2]){
      	frame[(int)((point[0] * Half_window_size*2*4) + (point[1] * 4 ) + 0)] = rgb[0];
      	frame[(int)((point[0] * Half_window_size*2*4) + (point[1] * 4 ) + 1)] = rgb[1];
      	frame[(int)((point[0] * Half_window_size*2*4) + (point[1] * 4 ) + 2)] = rgb[2];
        frame[(int)((point[0] * Half_window_size*2*4) + (point[1] * 4 ) + 3)] = point[2];
      }
      
    }
  } 
  return frame;
}

int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3],  LIGHT *light)
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// s = location of start of ray (probably the eye)
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light[3]

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
  L[0] = light->location[0] - p[0] ; 
  L[1] = light->location[1] - p[1] ; 
  L[2] = light->location[2] - p[2] ; 
  len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
  if (len == 0) return 0 ;
  L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;





  double max_ambient_and_diffuse = light->ambient + light->max_diffuse ;
     // this needs to occur BEFORE you possibly jump to LLL below




  double intensity ;
  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = light->ambient ; 
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
  else { diffuse = light->max_diffuse*NdotL ; }

  double specular ;
  if (EdotR <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,light->specpow) ;}

  // printf("%lf %lf\n",diffuse,specular) ;
  intensity = light->ambient + diffuse + specular ;



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


