#include <stdio.h>
#include <stdlib.h>
//  gcc -cpp -fPIC -shared calc_adv.c -lm -o calc_adv.so -O3

void one(double *charge_g,double *charge_h,double *charge_i,double *charge_s, int time_idx, int z_idx, int lat_idx, int lon_idx, int time, int time_step,double *u,double *v,double *w,double dx, double dy,double *z) {
  int h, i, j, k,dt;
  double dz;    

  #define charge_g(h,i,j,k)   charge_g[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define charge_h(h,i,j,k)   charge_h[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define charge_i(h,i,j,k)   charge_i[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define charge_s(h,i,j,k)   charge_s[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define u(h,i,j,k)                 u[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define v(h,i,j,k)                 v[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define w(h,i,j,k)                 w[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define z(h,i,j,k)                 z[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]

  // #時間ごとのループはpythonで回す
  // #両側差分を採用しているので、縁の格子点はゼロのままになる。
  dt=time_step;
  // printf("%d\n",dt);
    for(i=1; i<z_idx-1; i++){
      for(j=1; j<lat_idx-1; j++){
        for(k=1; k<lon_idx-1; k++){
          dz=(z(time,i+1,j,k)-z(time,i-1,j,k))/2;
          charge_g(time,i,j,k)=charge_g(time,i,j,k)
          -u(time,i,j,k)*(charge_g(time,i,j,k+1)-charge_g(time,i,j,k-1))/(2*dx)*dt
          -v(time,i,j,k)*(charge_g(time,i,j+1,k)-charge_g(time,i,j-1,k))/(2*dy)*dt
          -(w(time,i,j,k)-2)*(charge_g(time,i+1,j,k)-charge_g(time,i-1,j,k))/(2*dz)*dt;
          charge_h(time,i,j,k)=charge_h(time,i,j,k)
          -u(time,i,j,k)*(charge_h(time,i,j,k+1)-charge_h(time,i,j,k-1))/(2*dx)*dt
          -v(time,i,j,k)*(charge_h(time,i,j+1,k)-charge_h(time,i,j-1,k))/(2*dy)*dt
          -(w(time,i,j,k)-3)*(charge_h(time,i+1,j,k)-charge_h(time,i-1,j,k))/(2*dz)*dt;
          charge_i(time,i,j,k)=charge_i(time,i,j,k)
          -u(time,i,j,k)*(charge_i(time,i,j,k+1)-charge_i(time,i,j,k-1))/(2*dx)*dt
          -v(time,i,j,k)*(charge_i(time,i,j+1,k)-charge_i(time,i,j-1,k))/(2*dy)*dt
          -w(time,i,j,k)*(charge_i(time,i+1,j,k)-charge_i(time,i-1,j,k))/(2*dz)*dt;
          charge_s(time,i,j,k)=charge_s(time,i,j,k)
          -u(time,i,j,k)*(charge_s(time,i,j,k+1)-charge_s(time,i,j,k-1))/(2*dx)*dt
          -v(time,i,j,k)*(charge_s(time,i,j+1,k)-charge_s(time,i,j-1,k))/(2*dy)*dt
          -w(time,i,j,k)*(charge_s(time,i+1,j,k)-charge_s(time,i-1,j,k))/(2*dz)*dt;
        }
      }
    }
    // 地表面の帯電を最下層として計算する。
    for(j=1; j<lat_idx-1; j++){
      for(k=1; k<lon_idx-1; k++){
        dz=(z(time,1,j,k)-z(time,0,j,k));
        if ((w(time,0,j,k)-2)<0){
          charge_g(time,0,j,k)=charge_g(time,0,j,k)
          -(w(time,0,j,k)-2)*(charge_g(time,1,j,k)-charge_g(time,0,j,k))/dz*dt;
          }
        if ((w(time,0,j,k)-3)<0){
          charge_h(time,0,j,k)=charge_h(time,0,j,k)
          -(w(time,0,j,k)-3)*(charge_h(time,1,j,k)-charge_h(time,0,j,k))/dz*dt;
          }
        if (w(time,0,j,k)<0){
          charge_i(time,0,j,k)=charge_i(time,0,j,k)
          -w(time,0,j,k)*(charge_i(time,1,j,k)-charge_i(time,0,j,k))/dz*dt;
          charge_s(time,0,j,k)=charge_s(time,0,j,k)
          -w(time,0,j,k)*(charge_s(time,1,j,k)-charge_s(time,0,j,k))/dz*dt;
          }
      }
    }
    // for(i=1; i<z_idx-1; i++){
    //   for(j=1; j<lat_idx-1; j++){
    //     for(k=1; k<lon_idx-1; k++){
    //       dz=(z(time,i+1,j,k)-z(time,i-1,j,k))/2;
    //       charge_g(time,i,j,k)=charge_g(time,i,j,k)
    //       -u(time,i,j,k)*(charge_g(time,i,j,k+1)-charge_g(time,i,j,k-1))/(2*dx)*dt
    //       -v(time,i,j,k)*(charge_g(time,i,j+1,k)-charge_g(time,i,j-1,k))/(2*dy)*dt
    //       -(w(time,i,j,k)-2)*(charge_g(time,i+1,j,k)-charge_g(time,i-1,j,k))/(2*dz)*dt;
    //       charge_h(time,i,j,k)=charge_h(time,i,j,k)
    //       -u(time,i,j,k)*(charge_h(time,i,j,k+1)-charge_h(time,i,j,k-1))/(2*dx)*dt
    //       -v(time,i,j,k)*(charge_h(time,i,j+1,k)-charge_h(time,i,j-1,k))/(2*dy)*dt
    //       -(w(time,i,j,k)-3)*(charge_h(time,i+1,j,k)-charge_h(time,i-1,j,k))/(2*dz)*dt;
    //       charge_i(time,i,j,k)=charge_i(time,i,j,k)
    //       -u(time,i,j,k)*(charge_i(time,i,j,k+1)-charge_i(time,i,j,k-1))/(2*dx)*dt
    //       -v(time,i,j,k)*(charge_i(time,i,j+1,k)-charge_i(time,i,j-1,k))/(2*dy)*dt
    //       -w(time,i,j,k)*(charge_i(time,i+1,j,k)-charge_i(time,i-1,j,k))/(2*dz)*dt;
    //       charge_s(time,i,j,k)=charge_s(time,i,j,k)
    //       -u(time,i,j,k)*(charge_s(time,i,j,k+1)-charge_s(time,i,j,k-1))/(2*dx)*dt
    //       -v(time,i,j,k)*(charge_s(time,i,j+1,k)-charge_s(time,i,j-1,k))/(2*dy)*dt
    //       -w(time,i,j,k)*(charge_s(time,i+1,j,k)-charge_s(time,i-1,j,k))/(2*dz)*dt;
    //     }
    //   }
    // }
  }