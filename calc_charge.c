#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//  gcc -cpp -fPIC -shared calc_charge.c -lm -o calc_charge.so -O3

void one(double *charge_g,double *charge_h,double *charge_i,double *charge_s,double *fc,int time_idx,int z_idx,int lat_idx,int lon_idx,int time,int time_step,double *qt_phaci,double *qt_phacs,double *qt_pgaci,double *qt_pgacs,double *qc,double *temp,double *mass_g,double *mass_h,double *mass_i) {
  int h, i, j, k;
  double loop_a;
  double loop_b;
  double loop_c;
  double loop_d;

  #define charge_g(h,i,j,k)   charge_g[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define charge_h(h,i,j,k)   charge_h[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define charge_i(h,i,j,k)   charge_i[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define charge_s(h,i,j,k)   charge_s[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define fc(h,i,j,k)               fc[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define qt_phaci(h,i,j,k)   qt_phaci[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define qt_phacs(h,i,j,k)   qt_phacs[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define qt_pgaci(h,i,j,k)   qt_pgaci[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define qt_pgacs(h,i,j,k)   qt_pgacs[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define qc(h,i,j,k)               qc[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define temp(h,i,j,k)           temp[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define mass_g(h,i,j,k)       mass_g[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define mass_h(h,i,j,k)       mass_h[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define mass_i(h,i,j,k)       mass_i[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]

  // #時間ごとのループはpythonで回す
    for(i=0; i<z_idx; i++){
      for(j=0; j<lat_idx; j++){
        for(k=0; k<lon_idx; k++){
          // loop_a=i*lat_idx*lon_idx+j*lon_idx+k;
          loop_a=fabs(qt_phaci(time,i,j,k))*time_step/1e-9*1.0e-15*fc(time,i,j,k);
          loop_b=fabs(qt_pgaci(time,i,j,k))*time_step/1e-9*1.0e-15*fc(time,i,j,k);
          loop_c=fabs(qt_phacs(time,i,j,k))*time_step/1e-9*1.0e-15*fc(time,i,j,k);
          loop_d=fabs(qt_pgacs(time,i,j,k))*time_step/1e-9*1.0e-15*fc(time,i,j,k);
          charge_h(time,i,j,k)=charge_h(time,i,j,k)+loop_a+loop_c;
          charge_g(time,i,j,k)=charge_g(time,i,j,k)+loop_b+loop_d;
          charge_i(time,i,j,k)=charge_i(time,i,j,k)-loop_a-loop_b;
          charge_s(time,i,j,k)=charge_s(time,i,j,k)-loop_c-loop_d;
        }
      }
    }
  }