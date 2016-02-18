#include <FPT.h>
#include <D3d_matrix.h>


double light_in_eye_space[3];
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 50 ;
double zbuff[700][700];
double Half_window_size = 350;
double Tan_half_angle;


//========================================================================================================

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



//========================================================================================================

//shapes 

int hyperboloid(double u, double v, double xyz[3]){
  xyz[0] = sqrt(1+(v*v))*cos(u);
  xyz[1] = v;
  xyz[2] = sqrt(1+(v*v))*sin(u);
}


int sphere(double u, double v, double xy[3]){
  xy[0] = sqrt(1-(v*v))*cos(u);
  xy[1] = v;
  xy[2] = sqrt(1-(v*v))*sin(u);
}


//==========================================================================================================

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



/**===================================================================================================*/


void plot(double ustart, double uend, double vstart, double vend, 
  double rgb[3], double m[4][4], double V[4][4] , double eye[3], int (*func)(double u, double v,  double xy[3])){
  
  double i,j;
  double point[3],  color[3];  
  double p[3], p2[3], normal[3];
  
  for(i=ustart; i<uend;i = i + .001){
    for(j=vstart; j<vend;j = j + .001){

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



// tetrahedron model

void path (int frame_number, double path_xyz[3])
{
  double u,v,r ;
  double x,y,z ;

  u = 5*frame_number*M_PI/180 ;
  v = 0.3*u ;
  r = 2.0 + 1.4*sin(u) ;

  x = r*cos(u)*cos(v) ;
  y = r*sin(u) ;
  z = r*cos(u)*sin(v) ;
 
  path_xyz[0] = x ;
  path_xyz[1] = y ;
  path_xyz[2] = z ;

}




int init_scene (int frame_number)
{
  // model variables
  double xcen[4],ycen[4],zcen[4],brad ; // four nodes of tetrahedron
  double ccx,ccy,ccz,ccr ; // location of center of center sphere and radius

  double degrees_of_half_angle ;
  double eye[3],coi[3],up[3] ;
  double light_position[3], amb, diff, spow ;
  
  int i, j;
  for(i =0; i<Half_window_size*2;i++){
    for(j = 0; j< Half_window_size*2; j++){
        zbuff[i][j] = 1000000;
    }
  }

  //////////////////////////////////////////////
  //////////////////////////////////////////////
  // build a ball and stick model of a tetrahedron
  //////////////////////////////////////////////
  //////////////////////////////////////////////
  // 3 equally spaced pts around unit circle in the xz-plane 
  // form the base
  
  int k;
  for (k = 0 ; k < 3 ; k++) {
    double theta = 2*M_PI*k/3 ;
    xcen[k] = cos(theta) ;
    ycen[k] = 0 ;
    zcen[k] = sin(theta) ;
  }

  // you figure where the 4th node of the regular tetrahedron
  xcen[3] = 0 ; ycen[3] = sqrt(2) ; zcen[3] = 0 ;

  // also, figure out location of the 5th node of the model
  // which is at the center of mass of the tetrahedron
  ccx = (xcen[0] + xcen[1] + xcen[2] + xcen[3] )/4 ;
  ccy = (ycen[0] + ycen[1] + ycen[2] + ycen[3] )/4 ;
  ccz = (zcen[0] + zcen[1] + zcen[2] + zcen[3] )/4 ;
  //average of other points 

  brad = 0.08 ; // radius of the 4 verts of the tetrahedron
  ccr  = 0.20 ; // the radius of the center node of the model


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  degrees_of_half_angle = 25 ;

  path (frame_number, eye) ;

  coi[0] = ccx ;
  coi[1] = ccy ;
  coi[2] = ccz ;

  path (frame_number + 1, up) ;
  // printf("eye = %lf %lf %lf\n",eye[0],eye[1],eye[2]) ;
  // printf("coi = %lf %lf %lf\n",coi[0],coi[1],coi[2]) ;
  // printf("up  = %lf %lf %lf\n",up[0],up[1],up[2]) ;


  //////////////////////////////////////////////
  //////////////////////////////////////////////


  path (frame_number + 10, light_in_eye_space) ;
  amb  = 0.2 ;
  diff = 0.5 ;
  spow = 80 ;

  double nodeRgb[3];
  nodeRgb[0] = 1;
  nodeRgb[1] = 0;
  nodeRgb[2] = 0;
  
  double m[4][4], minv[4][4];
  D3d_make_identity(m);
  D3d_make_identity(minv);

  double V[4][4], Vi[4][4];
  D3d_make_identity(V);
  D3d_make_identity(Vi);

  D3d_view(V,Vi, eye, coi, up);
  
  // make points

  for(k=0; k<4;k++){
    makemat(m, minv, brad, brad, brad, 0,0,0,xcen[k],ycen[k], zcen[k]);
    plot(0, 2 * M_PI, -1, 1, nodeRgb, m, V, eye, sphere );
    D3d_make_identity(m);
    D3d_make_identity(minv);
  }
  //make center of mass
  nodeRgb[2]=1;
  makemat(m, minv, ccr, ccr, ccr, 0,0,0,ccx,ccy, ccz);
  plot(0, 2 * M_PI, -1, 1, nodeRgb, m, V, eye, sphere );
nodeRgb[0]=0;
int f;
for(k=0;k<4;k++){
    nodeRgb[0] = 0;
  nodeRgb[1] = 1;
  nodeRgb[2] = 0;
  for(f=0;f<4;f++){
    if(f==k){
      continue;
    }
    D3d_make_identity(m);
    D3d_make_identity(minv);
    
    double teye[3], tcoi[3], tup[3];

    teye[0] = xcen[k];
    teye[1] = ycen[k];
    teye[2] = zcen[k];

    tcoi[0] = xcen[f];
    tcoi[1] = ycen[f];
    tcoi[2] = zcen[f];

    tup[0] = xcen[k];
    tup[1] = ycen[k]+1;
    tup[2] = zcen[k];

    D3d_scale(m,minv, .03,.8,.03);
    D3d_rotate_x(m,minv,M_PI/2);
    D3d_translate(m,minv,0,0,1);
    // makemat(m, minv, brad, brad, brad, 0,0,0, 0,0,1 );

    double n[4][4], ninv[4][4];

    D3d_view(n,ninv, teye, tcoi, tup);

    D3d_mat_mult(m,ninv,m);

    plot(0, 2 * M_PI, -1, 1, nodeRgb, m, V, eye, hyperboloid );
}
  nodeRgb[0] = 1;
  nodeRgb[1] = 1;
  nodeRgb[2] = 1;
    D3d_make_identity(m);
    D3d_make_identity(minv);
    
    double teye[3], tcoi[3], tup[3];

    teye[0] = ccx;
    teye[1] = ccy;
    teye[2] = ccz;

    tcoi[0] = xcen[k];
    tcoi[1] = ycen[k];
    tcoi[2] = zcen[k];

    tup[0] = ccx;
    tup[1] = ccy+1;
    tup[2] = ccz;

    D3d_scale(m,minv, .03,.5,.03);
    D3d_rotate_x(m,minv,M_PI/2);
    D3d_translate(m,minv,0,0,.5);
    // makemat(m, minv, brad, brad, brad, 0,0,0, 0,0,1 );

    double n[4][4], ninv[4][4];

    D3d_view(n,ninv, teye, tcoi, tup);

    D3d_mat_mult(m,ninv,m);

    plot(0, 2 * M_PI, -1, 1, nodeRgb, m, V, eye, hyperboloid );
}
printf("frame number %d \n",frame_number );
}


int main(){

  char prefix[100];
  printf("Enter file extension\n");
  scanf("%s", prefix);


  Tan_half_angle = tan((M_PI*25)/180);
  G_init_graphics(Half_window_size*2, Half_window_size*2);
  G_rgb(0,0,0);
  G_clear();


  int i;
  for(i= 0; i < 73; i++){
    init_scene(i);
    char filename[100];
    sprintf(filename, "%s%04d.xwd", prefix, i);
    G_save_image_to_file(filename);
  }


}