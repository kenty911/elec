#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//  gcc -cpp -fPIC -shared flash.c -lm -o flash.so -O3

void one(double *rho_a,double *rho_g,double *rho_h,double *rho_i,double *rho_s,double *phi,double *fod, int z_idx, int lat_idx, int lon_idx, int time,double dx,double dy,double *z) {
  int h, i, j, k, p, q;
  double dz;    
  double dphi;    
  double sgn_rho;    
  double u_rho;    
  double d_rho;    
  int flag1;
  int count1;
  double dum_rho;

  const double lim  = 3.00e6 ; /* 閾値*/
  const double radi = 6000;  /*中和を行う半径*/
  

  #define rho_a(h,i,j,k)   rho_a[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define rho_g(h,i,j,k)   rho_g[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define rho_h(h,i,j,k)   rho_h[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define rho_i(h,i,j,k)   rho_i[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define rho_s(h,i,j,k)   rho_s[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define phi(h,i,j,k)       phi[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define fod(h,i,j,k)       fod[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
  #define z(h,i,j,k)           z[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]

  // 時間ごとのループはpythonで回す
  // 今回だけはｚ方向のループを最後に回す
  for(j=1; j<lat_idx-1; j++){
    for(k=1; k<lon_idx-1; k++){
      flag1=0;
      dum_rho=0;
      for(i=1; i<z_idx-1; i++){
          dz=0.0;
          dphi=0.0;
          dz=(z(time,i+1,j,k)-z(time,i-1,j,k))/2;
          dphi=(phi(time,i+1,j,k)-phi(time,i-1,j,k))/2;  
          if (fabs(dphi/dz)>lim){
              sgn_rho=0.0;
              d_rho=rho_a(time,i+1,j,k);
              u_rho=rho_a(time,i-1,j,k);
              sgn_rho=d_rho*u_rho;
              // (上下で電荷の符号が違う・1nC以上)が条件
              if (sgn_rho<-1e-18){
                // printf("%e\n",sgn_rho);
                flag1=1;
                }
              }
          }
      // if (flag1==1 && fod(time,0,j,k)<=fod(time,0,j-1,k) && fod(time,0,j,k)<=fod(time,0,j,k-1) && fod(time,0,j,k)<=fod(time,0,j+1,k) &&  fod(time,0,j,k)<=fod(time,0,j,k+1)){
      if (flag1==1){
        count1=0;
        //　雷の回数のカウント
        fod(time,0,j,k) = fod(time,0,j,k)+1;
        // 放電後に電荷の中和を行う。中和は4クラスで等分に分けているが、質量重みしたほうがいい？！
        // 中和を行う際は、横方向にも広がりを持たせる。
        // for(i=1; i<z_idx-1; i++){
        for(i=0; i<z_idx-1; i++){
          for(p=-floor(radi/dy); p<radi/dy; p++){
            for(q=-floor(pow(pow(radi,2)-pow(p*dy,2),0.5)/dx); q<pow(pow(radi,2)-pow(p*dy,2),0.5)/dx; q++){
              dum_rho=dum_rho+rho_a(time,i,j+p,k+q);
              count1++;    
            }
          }
        } 
        // for(i=1; i<z_idx-1; i++){
        for(i=0; i<z_idx-1; i++){
          for(p=-floor(radi/dy); p<radi/dy; p++){
            for(q=-floor(pow(pow(radi,2)-pow(p*dy,2),0.5)/dx); q<pow(pow(radi,2)-pow(p*dy,2),0.5)/dx; q++){       
              // rho_g(time,i,j+p,k+q) = 0;
              // rho_h(time,i,j+p,k+q) = 0;
              // rho_i(time,i,j+p,k+q) = 0;
              // rho_s(time,i,j+p,k+q) = 0;
              // rho_g(time,i,j+p,k+q) = rho_g(time,i,j+p,k+q)*0.5 + dum_rho/4.0/count1*0.5;
              // rho_h(time,i,j+p,k+q) = rho_h(time,i,j+p,k+q)*0.5 + dum_rho/4.0/count1*0.5;
              // rho_i(time,i,j+p,k+q) = rho_i(time,i,j+p,k+q)*0.5 + dum_rho/4.0/count1*0.5;
              // rho_s(time,i,j+p,k+q) = rho_s(time,i,j+p,k+q)*0.5 + dum_rho/4.0/count1*0.5;
              rho_g(time,i,j+p,k+q) = dum_rho/4.0/count1*0.5;
              rho_h(time,i,j+p,k+q) = dum_rho/4.0/count1*0.5;
              rho_i(time,i,j+p,k+q) = dum_rho/4.0/count1*0.5;
              rho_s(time,i,j+p,k+q) = dum_rho/4.0/count1*0.5;
              // rho_g(time,i,j,k) = 0;
              // rho_h(time,i,j,k) = 0;
              // rho_i(time,i,j,k) = 0;
              // rho_s(time,i,j,k) = 0;
              phi(time,i,j+p,k+q) = 0;
              }
            }
          }
        }
      }
    }
  }
// }