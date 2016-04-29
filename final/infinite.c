#include <FPT.h>
#include <D3d_matrix.h>
#include <xwd_tools.h>
#include "shape.h"


#define AMBIENT       0.2 
#define MAX_DIFFUSE   0.5 
#define SPECPOW       50 
#define SHOWN         100000
#define DEPTH         1
#define CHILDREN      16
#define REFLECTIONS   2

double Half_window_size = 350;
double Half_angle_degrees ;
double Tan_half_angle ;
double light_in_eye_space[3] ;

double eye[3];

double ms[SHOWN][4][4];
double minvs[SHOWN][4][4];
int shapeNum = 0;
double offset = 0;
double colors[SHOWN][3];
double reflectivity[SHOWN];

/**===================================================================================================*/

void halfed(double vals[4], double x, double y, double z, double rad, int pos, int total){
  double pi = (2*M_PI*pos)/total;
  vals[0] = x + cos(pi)*rad*1.5;
  vals[1] = y + sin(pi)*rad*1.5;
  vals[2] = z ;
  vals[3] = rad/3;
}

void equa(double vals[4], double x, double y, double z, double rad, int pos, int total){
 
  double pi = (2*M_PI*pos)/total;
  vals[0] = x + cos(pi)*rad*1.5;
  vals[1] = y + sin(pi)*rad*1.5;
  vals[2] = z + (((double)pos/total)*6 -3)*rad;
  vals[3] = rad/3;

}

void spin(double vals[4], double x, double y, double z, double rad, int pos, int total){
  double  temp = offset/10.0;
  // printf("%lf \n",temp );
  double pi = (2*M_PI*pos)/total;
  vals[0] = x + cos(pi)*rad*1.5;
  vals[1] = y + sin(pi)*rad*1.5;
  vals[2] = z + (((double)pos/total)*(2*temp) -temp)*rad;
  vals[3] = rad/3;
}

void spiral(double vals[4], double x, double y, double z, double rad, int pos, int total){
  double u = (2.0*M_PI*pos)/(double)total;
  double v = ((2.0*pos)/(double)total) - 1.0;
  vals[3] = rad/3;
  vals[0] = x + rad*1.5*sqrt(1 -(v*v))*cos(u);
  vals[1] = y + rad*1.5*v;
  vals[2] = z + rad*1.5*sqrt(1-(v*v))*sin(u);
  printf("%lf , %lf , %lf \n",vals[0],vals[1],vals[2] );

}

void spaced(double vals[4], double x, double y, double z, double rad, int pos, int total){

  int square = sqrt(total);

  double u = (2*M_PI/square)*(pos%square);
  double v = (2.0/square)*(pos/square) - 1;

  // u = sqrt(1 - (u*u));
  printf("\n");
  printf("position = %d \n", pos);
  printf("v = %lf , u = %lf \n", v, u);

  vals[3] = rad/3;
  vals[0] = x + rad*1.5*sqrt(1 -(v*v))*cos(u);
  vals[1] = y + rad*1.5*v;
  vals[2] = z + rad*1.5*sqrt(1-(v*v))*sin(u);
  printf("%lf , %lf , %lf \n",vals[0],vals[1],vals[2] );

}

void test(double vals[4], double x, double y, double z, double rad, int pos, int total){

  double area = (4*M_PI*(rad*rad))/(double)total;
  double radius = sqrt(area/(2*M_PI));
  double increment = radius/(M_PI*rad);
  
  int square = sqrt(total);

  double u = increment*pos;
  double v = 0;//(2.0)*((u/2*M_PI)/(2*M_PI)) - 1;

  printf("increment = %lf , u = %lf , v = %lf \n",increment, u, v );
  
  vals[3] = rad/3;
  vals[0] = x + rad*1.5*sqrt(1 -(v*v))*cos(u);
  vals[1] = y + rad*1.5*v;  
  vals[2] = z + rad*1.5*sqrt(1-(v*v))*sin(u);
}

void npoints(double vals[4], double x, double y, double z, double rad, int pos, int total){
  double a = (4.0*M_PI*rad*rad)/(double)total;
  double d = sqrt(a);
  double Mo = ceil(M_PI/d);
  double dt = M_PI/Mo;
  double dy = a/dt;
  double m = (double)(pos%total);
  double V = (M_PI*(m+.5))/Mo;
  double My = ceil((2*M_PI*sin(V))/dy);
  double n = ceil(pos/total);
  double U = (2*M_PI*n)/My;

  printf("u = %lf , v = %lf\n", U,V );
  vals[3] = rad/3;
  vals[0] = x + rad*1.5*sqrt(1 -(V*V))*cos(U);
  vals[1] = y + rad*1.5*V;  
  vals[2] = z + rad*1.5*sqrt(1-(V*V))*sin(U);
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
int tracer(double points[2][3], int numRef,
          double intersect[3], double normal[3], double color[3]){
    int i;
    double pfinal[3];
    double tpoint[3], tpoint2[3];
    double tsave[100];
    int osave[100];
    double mint;
    int mino;
    int numsave = 0;

  if(numRef <=0){ return -1; }
   for(i=0; i < shapeNum; i++){
      //convert points to object space
      D3d_mat_mult_pt(tpoint, minvs[i], points[0]);
      D3d_mat_mult_pt(tpoint2, minvs[i], points[1]);

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
      if(res[0] > 0 && res[0] == res[0]){
        tsave[numsave] = res[0]; osave[numsave] = i; numsave++;
        }
      if(res[1] > 0 && res[1] == res[1]){
        tsave[numsave] = res[1]; osave[numsave] = i; numsave++;
      }  
  }//end loop

  if(numsave == 0){
      color[0] = 0;
      color[1] = 0;
      color[2] = 0;
    return -1;
  }

  mint= tsave[0]; mino = osave[0];
  int j;
  for(j=1; j<numsave; j++){
    if(tsave[j] < mint){
      mint = tsave[j];
      mino = osave[j];
    }
  }


    double fnorm[3];

    D3d_mat_mult_pt(tpoint, minvs[mino], points[0]);
    D3d_mat_mult_pt(tpoint2, minvs[mino], points[1]);

    pfinal[0] = tpoint[0] + (mint * (tpoint2[0] - tpoint[0]));
    pfinal[1] = tpoint[1] + (mint * (tpoint2[1] - tpoint[1]));
    pfinal[2] = tpoint[2] + (mint * (tpoint2[2] - tpoint[2]));

    //normal vector in world space based at orign
    fnorm[0] = 2*minvs[mino][0][0]*pfinal[0] + 2*minvs[mino][1][0]*pfinal[1] + 2*pfinal[2]*minvs[mino][2][0];
    fnorm[1] = 2*minvs[mino][0][1]*pfinal[0] + 2*minvs[mino][1][1]*pfinal[1] + 2*pfinal[2]*minvs[mino][2][1];
    fnorm[2] = 2*minvs[mino][0][2]*pfinal[0] + 2*minvs[mino][1][2]*pfinal[1] + 2*pfinal[2]*minvs[mino][2][2] ;

    //draw rectangle where ray intersects circle
   
     D3d_mat_mult_pt(pfinal,ms[mino],pfinal);

    //translate normal to intersection point and extend
    double tx = pfinal[0] + (0.001*fnorm[0]);
    double ty = pfinal[1] + (0.001*fnorm[1]);
    double tz = pfinal[2] + (0.001*fnorm[2]);

    double inc[3], ref[3], nor[3];

    //caluclate light vector based at origin
    inc[0] = points[0][0] - points[1][0];
    inc[1] = points[0][1] - points[1][1];
    inc[2] = points[0][2] - points[1][2];
    make_unit_vector(inc);
    
    //calculate normal vector based at origin 
    nor[0] = tx - pfinal[0];
    nor[1] = ty - pfinal[1];
    nor[2] = tz - pfinal[2] ;
    make_unit_vector(nor);

    //reflection 
    reflection(nor, inc,ref);
    
    points[0][0] = pfinal[0] + (0.00001*ref[0]);
    points[0][1] = pfinal[1] + (0.00001*ref[1]);
    points[0][2] = pfinal[2] + (0.00001*ref[2]) ;
    
    points[1][0] = (pfinal[0] + (.01*ref[0]));
    points[1][1] = (pfinal[1] + (.01*ref[1]));
    points[1][2] = (pfinal[2] + (.01*ref[2]));

  if(reflectivity[mino] == 0){
    Light_Model (colors[mino], inc, pfinal, normal, color);
    return mino;
  }

  int flag = tracer(points, numRef-1, pfinal, nor, color);
  double tcolor[3];
  Light_Model (colors[mino], inc, pfinal, fnorm, tcolor);
  
  color[0] = (1-reflectivity[mino])*tcolor[0] + reflectivity[mino] * color[0];
  color[1] = (1-reflectivity[mino])*tcolor[1] + reflectivity[mino] * color[1];
  color[2] = (1-reflectivity[mino])*tcolor[2] + reflectivity[mino] * color[2];

  if(flag == -1){
    Light_Model (colors[mino], inc, pfinal, fnorm, color);
  }
  return mino;
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



int make_unit_vector(double v[3]){
	double length = sqrt((v[0]*v[0])+(v[1]*v[1])+(v[2]*v[2]));
			v[0] = v[0]/length;
			v[1] = v[1]/length;
			v[2] = v[2]/length;
}

//**======================================================================**/

void draw(char filename[100]){
  int i,j;

  int map = create_new_xwd_map (Half_window_size*2,Half_window_size*2);
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
      int sphere =  tracer(points, REFLECTIONS,  intersect, normal, color);

      if(sphere == -1){

      set_xwd_map_color(map,i + Half_window_size,j + Half_window_size, 0, 0, 0);
      }else{


      double eye[3];
      eye[0] = 0; eye[1] = 0; eye[2] = 0;


      set_xwd_map_color(map,i + Half_window_size,j + Half_window_size,color[0], color[1], color[2]);


      // G_rgb(color[0], color[1], color[2]);
      // G_point(i + Half_window_size,j + Half_window_size);
      }
    }   
  }

  xwd_map_to_named_xwd_file(map, filename) ;
}

void fractal(SHAPE seed, int levels){
  if(levels <=0){return;}
  int i;
  double tcol[seed.numChildren+1][3];
  double tref[seed.numChildren+1];

  colors[shapeNum][0] = seed.color[0];
  colors[shapeNum][1] = seed.color[1];
  colors[shapeNum][2] = seed.color[2];

  reflectivity[shapeNum] = seed.ref;

  make_matrix(seed, ms[shapeNum], minvs[shapeNum]);
  shapeNum++;

  for(i=0; i<seed.numChildren; i++){
    if(i >= SHOWN){ return; }

      SHAPE child = get_child(seed,i);
      make_matrix(child, ms[shapeNum], minvs[shapeNum] );
    



      reflectivity[0] = seed.ref;
      colors[shapeNum][0] = child.color[0];
      colors[shapeNum][1] = child.color[1];
      colors[shapeNum][2] = child.color[2];

      reflectivity[shapeNum] = child.ref;
      shapeNum++;
      // draw(m, minv, seed.numChildren + 1);
      fractal(child, levels - 1);
  }


  // draw(m, minv, seed.numChildren + 1);
  // G_wait_key();
}


//**======================================================================**/

int main(){

  char prefix[100];

prefix[0] = 't';
prefix[1] = 'e';
prefix[2] = 's';
prefix[3] = 't';

  int i;


  light_in_eye_space[0] = 100;
  light_in_eye_space[1] = 100;
  light_in_eye_space[2] = 0;

	Half_window_size  = 300;
	Half_angle_degrees = 30;
	Tan_half_angle = tan(Half_angle_degrees*M_PI/180) ;

  eye[0] = 0;
  eye[1] = 0;
  eye[2] = 0;

  double color[3];
  color[0] = 0;
  color[1] = 0;
  color[2] = 1;

  // calc fractal here if only once is needed
    SHAPE seed = new_shape(0,0,300 ,30,npoints , color, .5,CHILDREN);
   fractal(seed,DEPTH);

  int numframes = 100;
   for(i=0;i<numframes;i++){
    offset = (double)i;



    //create fractal here if shape num is used in calculation 
        // shapeNum = 0;
   //  SHAPE seed = new_shape(0,0,300 ,30,spiral , color, .5,CHILDREN);
   // fractal(seed,DEPTH);


    int j;
    for(j=0;j<shapeNum; j++){
      D3d_translate(ms[j], minvs[j],0,0,-300);
      D3d_rotate_y(ms[j], minvs[j],.2);
      D3d_translate(ms[j], minvs[j],0,0,300);
    }
 
  //=============================================================================

    char filename[100];
    sprintf(filename, "%s%04d.xwd", prefix, i);
    draw(filename);
    printf("Image %d \n", i);

//---------------------------------------------------------------------------
}
}
