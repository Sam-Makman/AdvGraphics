#include <shape_tools.h>
#include <math.h>
#include <stdio.h>


int main(){
	SHAPE shape;
	LIGHT light;
	int half_window = 250;
	double location[3];
	location[0] = 100; location[1] = 200; location[2] = -100;
	init_light(&light,.3,.5,70,location);

 	double V[4][4], m[4][4], minv[4][4],eye[3];
    D3d_make_identity(V);
    D3d_make_identity(m);
    D3d_make_identity(minv);
    eye[0]=0;
    eye[1]=0;
    eye[2]=0;
    makemat(m, minv, 1,1,1, 0,0,0,0,3,20);


	init_shape(&shape, 0, 2, 0, 2, checkers, hyperboloid,m,V);
	print_shape(&shape);
	double * result;

	result = plot(eye,&light, &shape, half_window);
	int map = create_new_xwd_map (half_window*2,half_window*2);
	int i,j;
	for(i=0; i< half_window*2; i++){
		for(j=0;j<half_window*2;j++){
			double r =result[(i*(half_window*2)*(4)) + (j*4) + 0];
			double g = result[(i*(half_window*2)*(4)) + (j*4) + 1];
			double b = result[(i*(half_window*2)*(4)) + (j*4) + 2];
			// printf("r = %lf , g = %lf, b = %lf\n", r, g , b );
			set_xwd_map_color(map,i,j,r,g,b);
			}
		}
  	xwd_map_to_named_xwd_file(map, "file0000.xwd") ;
	
}
