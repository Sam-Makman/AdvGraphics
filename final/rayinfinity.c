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

double colors[10][3];
double reflectivity[10];

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

int quadratic(double res[2], double a, double b, double c){
  if((b*b)-(4*a*c) < 0){
    return -1;
  }
  res[0] = (-b + sqrt((b*b)-(4*a*c)))/(2*a);
  res[1] = (-b - sqrt((b*b)-(4*a*c)))/(2*a);
  return 1;
}



double reflection(double nor[3], double inc[3], double ref[3]){
  double m;
  m = 2*(nor[0]*inc[0] + nor[1]*inc[1] + nor[2]*inc[2]);
  
  ref[0] = m*nor[0] - inc[0];
  ref[1] = m*nor[1] - inc[1];
  ref[2] = m*nor[2] - inc[2];
}
  
  void print_pt(double point[3]){
    printf("%lf , %lf , %lf\n",point[0], point[1], point[2] );
  }

/**======================================================================================================*/
int tracer(double m[10][4][4],double minv[10][4][4], double points[2][3], int n, int numRef,
          double intersect[3], double normal[3], double color[3]){
    int i;
    int nfinal = -1;
    double pfinal[3];

    double tpoint[3], tpoint2[3];
    double tsave[100];
    int osave[100];
    double mint;
    int mino;
    int numsave = 0;

  if(numRef <=0){ return -1; }

   for(i=0; i < n; i++){

      //convert points to object space
      D3d_mat_mult_pt(tpoint, minv[i], points[0]);
      D3d_mat_mult_pt(tpoint2, minv[i], points[1]);

      //find slope of two points 
      double p = tpoint2[0] - tpoint[0];
      double q = tpoint2[1] - tpoint[1];
      double r = tpoint2[2] - tpoint[2];

      //find a,b,c for quadtratic equation
      double a = pow(p,2) + pow(q,2) + pow(r,2);
      double b = 2*tpoint[0]*p + 2*tpoint[1]*q + 2*tpoint[2]*r;
      double c = pow(tpoint[0],2) + pow(tpoint[1],2) + pow(tpoint[2],2) - 1;

      double res[2];
      int result = quadratic(res, a, b ,c);

      if(result < 0){ continue; }

      //find closest point 
      
      if(res[0] > 0){

        tsave[numsave] = res[0]; osave[numsave] = i; numsave++;
        }
      if(res[1] > 0){
        tsave[numsave] = res[1]; osave[numsave] = i; numsave++;
      }  
  }//end loop

  if(numsave == 0){
      color[0] = 0;
      color[1] = 0;
      color[2] = 0;
    return -1;
  }


  mint=tsave[0]; mino = osave[0];
  int j;
  for(j=1; j<numsave; j++){
    if(tsave[j] < mint){
      mint = tsave[j];
      mino=osave[j];
    }
  }


    double fnorm[3];

    D3d_mat_mult_pt(tpoint, minv[mino], points[0]);
    D3d_mat_mult_pt(tpoint2, minv[mino], points[1]);

    pfinal[0] = tpoint[0] + (mint * (tpoint2[0] - tpoint[0]));
    pfinal[1] = tpoint[1] + (mint * (tpoint2[1] - tpoint[1]));
    pfinal[2] = tpoint[2] + (mint * (tpoint2[2] - tpoint[2]));

    nfinal = mino;
    //normal vector in world space based at orign
    fnorm[0] = 2*minv[nfinal][0][0]*pfinal[0] + 2*minv[nfinal][1][0]*pfinal[1] + 2*pfinal[2]*minv[nfinal][2][0];
    fnorm[1] = 2*minv[nfinal][0][1]*pfinal[0] + 2*minv[nfinal][1][1]*pfinal[1] + 2*pfinal[2]*minv[nfinal][2][1];
    fnorm[2] = 2*minv[nfinal][0][2]*pfinal[0] + 2*minv[nfinal][1][2]*pfinal[1] + 2*pfinal[2]*minv[nfinal][2][2] ;

    //draw rectangle where ray intersects circle
   
     D3d_mat_mult_pt(pfinal,m[mino],pfinal);

    //translate normal to intersection point and extend
    double tx = pfinal[0] + (0.001*fnorm[0]);
    double ty = pfinal[1] + (0.001*fnorm[1]);
    double tz = pfinal[2] + (0.001*fnorm[2]);

    double inc[3], ref[3], nor[3];

    //caluclate light vector based at origin
    inc[0] = points[1][0] - pfinal[0];
    inc[1] = points[1][1] - pfinal[1];
    inc[2] = points[1][2] - pfinal[2] ;
    make_unit_vector(inc);
    
    //calculate normal vector based at origin 
    nor[0] = tx - pfinal[0];
    nor[1] = ty - pfinal[1];
    nor[2] = tz - pfinal[2] ;
    make_unit_vector(nor);

    //reflection 
    reflection(nor, inc,ref);
    
    points[0][0] = pfinal[0] + (0.001*ref[0]);
    points[0][1] = pfinal[1] + (0.001*ref[1]);
    points[0][2] = pfinal[2] + (0.001*ref[2]) ;
    
    points[1][0] = (pfinal[0] + ref[0]);
    points[1][1] = (pfinal[1] + ref[1]);
    points[1][2] = (pfinal[2] + ref[2]) ;

    intersect[0] = pfinal[0];
    intersect[1] = pfinal[1];
    intersect[2] = pfinal[2];

    normal[0] = fnorm[0];
    normal[1] = fnorm[1];
    normal[2] = fnorm[2];

  if(reflectivity[mino] == 0){
    Light_Model (colors[nfinal], inc, intersect, normal, color);
    return mino;
  }

  int flag = tracer(m,minv,points, n, numRef-1, pfinal, nor, color);

  if(flag == -1){
    Light_Model (colors[nfinal], inc, intersect, normal, color);
  }
  return nfinal;
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

//**======================================================================**/

void draw(double m[10][4][4], double minv[10][4][4],  int n){
  int i,j;
  for(i = -Half_window_size; i< Half_window_size; i++){
    for(j = -Half_window_size; j< Half_window_size; j++){
      double points[2][3];
      
      points[0][0] = 0;
      points[0][1] = 0;
      points[0][2] = 0;

      points[1][0] = i;
      points[1][1] = j;
      points[1][2] = 300/Tan_half_angle;

      double intersect[3], normal[3];

      double color[3];
      int sphere =  tracer(m, minv, points, n, 2, intersect, normal, color);

      if(sphere == -1){
        G_rgb(0,0,0);
        G_point(i + Half_window_size, j + Half_window_size);
      }else{


      double eye[3];
      eye[0] = 0; eye[1] = 0; eye[2] = 0;


      // Light_Model (colors[sphere], eye, intersect, normal, color);

      G_rgb(color[0], color[1], color[2]);
      G_point(i + Half_window_size,j + Half_window_size);
      }
    }   
  }
  printf("Done\n");
}



//**======================================================================**/

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
  int x,y,i;


	double rgb[3];

  colors[0][0] = 1;
  colors[0][1] = 0;
  colors[0][2] = 0;

  colors[1][0] = 0;
  colors[1][1] = 0;
  colors[1][2] = 1;

  colors[2][0] = 0;
  colors[2][1] = 1;
  colors[2][2] = 0;
  
  reflectivity[0] = 1;
  reflectivity[1] = 0;
  reflectivity[2] = 0;

  AMBIENT = 0.2 ;
  MAX_DIFFUSE = 0.5 ;
  SPECPOW = 30 ;

  light_in_eye_space[0] = 100;
  light_in_eye_space[1] = 100;
  light_in_eye_space[2] = 0;


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
  double numframes = 20;

   for(i=0;i<1;i++){
//===========================================================================//

    double V[4][4], Vi[4][4];

    t = 2.0 * M_PI *fnum/numframes ; // goes from 0 to 1
    printf("%lf\n",t );

    // eye[0] = 300*cos(2*M_PI*t) ; 
    // eye[1] =  150*t ; 
    // eye[2] =  150*sin(2*M_PI*t) ; 

    eye[0] = 0 ; 
    eye[1] = 0 ; 
    eye[2] = 0 ; 

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

  double m[10][4][4], minv[10][4][4];


  D3d_make_identity(m[0]);
  D3d_make_identity(minv[0]);

  makemat(m[0], minv[0],30, 30, 30, 0, 0,0, 150);
  makemat(m[1], minv[1],30, 30, 30, 0, 50,0, 150);
  makemat(m[2], minv[2],30, 30, 30, 0, -50,0, 150);

  
  draw(m, minv, 3);
  G_wait_key();
//=============================================================================

  //   char filename[100];
  //   sprintf(filename, "%s%04d.xwd", prefix, i);
  //   G_save_image_to_file(filename);
  fnum++;

  G_rgb(0,0,0);
  G_clear();

//---------------------------------------------------------------------------
}
}
