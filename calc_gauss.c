#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//  gcc -cpp -fPIC -shared calc_gauss.c -lm -o calc_gauss.so -O3

void calc_gauss(double *rho, double *phi, int time_idx, int z_idx, int lat_idx, int lon_idx, int time,double dx, double dy,double *z) {
    const double Conv  = 1.0e-6    ;     /* 収束と判定する差*/
    const double e0    = 8.85e-12  ;

    int h, i, j, k;
    int loop;
    double Prev_phi;    
    // double dx_test;    
    // double dy_test;    
    double dz;    
    double Ex,Ey,Ez; 
    double MaxPhi;     /* 最大電位*/
    double MaxErr;     /* 最大のエラー*/
    double CurErr;     /* 現在のエラー*/

    #define rho(h,i,j,k) rho[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
    #define phi(h,i,j,k) phi[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
    #define z(h,i,j,k)     z[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]

    loop=0;
    MaxPhi=1.0e-10; /* ゼロ割り防止*/

    // printf("%e\n",dx);
    // printf("%e\n",dy);
    // dx_test=1000;
    // dy_test=1000;
    do{
        // if(!(loop%1000))printf("%5d %e\n", loop, MaxPhi); /* 1000ループ毎に経過表示 */
        MaxErr = CurErr = 0.0;
        for(i=1; i<z_idx-1; i++){
            for(j=1; j<lat_idx-1; j++){
                for(k=1; k<lon_idx-1; k++){
                        Prev_phi = phi(time,i,j,k);
                        // dzを求める過程を追加する。
                        dz=(z(time,i+1,j,k)-z(time,i-1,j,k))/2;
                        // printf("%e\n",dz);
                        phi(time,i,j,k) = (0.166666666666)*(rho(time,i,j,k)*dx*dy*dz/e0+phi(time,i+1,j,k)+phi(time,i-1,j,k)+phi(time,i+1,j,k)+phi(time,i,j-1,k)+phi(time,i,j,k+1)+phi(time,i,j,k-1));
                        if (MaxPhi<fabs(phi(time,i,j,k)))MaxPhi=phi(time,i,j,k);
                        CurErr = (fabs(phi(time,i,j,k) - Prev_phi))/MaxPhi;
                        if (MaxErr < CurErr) MaxErr = CurErr;       
                }
            }
        }
        loop++;
    }while (MaxErr>Conv); /* 領域全ての点の誤差がConvを下回ったらおしまい*/
    // printf("%5d %e\n", loop, MaxPhi);
}

void calc_sor(double *rho, double *phi, int time_idx, int z_idx, int lat_idx, int lon_idx, int time,double dx, double dy,double *z) {
    const double Conv  = 1.0e-6    ;     /* 収束と判定する差*/
    const double e0    = 8.85e-12  ;
    const double omega_SOR = 1.50  ;   
    // const double omega_SOR = 1.01  ;   

    int h, i, j, k;
    int loop;
    double Prev_phi;    
    // double dx_test;    
    // double dy_test;    
    double dz;    
    double Ex,Ey,Ez; 
    double MaxPhi;     /* 最大電位*/
    double MaxErr;     /* 最大のエラー*/
    double CurErr;     /* 現在のエラー*/
    double def_phi;     /* 現在のエラー*/

    #define rho(h,i,j,k) rho[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
    #define phi(h,i,j,k) phi[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
    #define z(h,i,j,k)     z[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]

    loop=0;
    MaxPhi=1.0e-10; /* ゼロ割り防止*/

    // printf("%e\n",dx);
    // printf("%e\n",dy);
    // dx_test=1000;
    // dy_test=1000;

    do{
        // if(!(loop%1000))printf("%5d %e\n", loop, MaxPhi); /* 1000ループ毎に経過表示 */
        MaxErr = CurErr = 0.0;
        for(i=1; i<z_idx-1; i++){
            for(j=1; j<lat_idx-1; j++){
                for(k=1; k<lon_idx-1; k++){
                        Prev_phi = phi(time,i,j,k);
                        // dzを求める過程を追加する。
                        dz=(z(time,i+1,j,k)-z(time,i-1,j,k))/2;
                        // printf("%e\n",dz);
                        // phi(time,i,j,k) = (0.16666666666)*(rho(time,i,j,k)*dx*dy*dz/e0+phi(time,i+1,j,k)+phi(time,i-1,j,k)+phi(time,i+1,j,k)+phi(time,i,j-1,k)+phi(time,i,j,k+1)+phi(time,i,j,k-1));
                        def_phi=(0.166666666666)*(rho(time,i,j,k)*dx*dy*dz/e0+phi(time,i+1,j,k)+phi(time,i-1,j,k)+phi(time,i+1,j,k)+phi(time,i,j-1,k)+phi(time,i,j,k+1)+phi(time,i,j,k-1))-Prev_phi;
                        if (loop==0){
                            phi(time,i,j,k) = Prev_phi + def_phi;
                            }
                        else{
                            phi(time,i,j,k) = Prev_phi + omega_SOR*def_phi;
                            }
                        if (MaxPhi<fabs(phi(time,i,j,k)))MaxPhi=phi(time,i,j,k);
                        CurErr = (fabs(phi(time,i,j,k) - Prev_phi))/MaxPhi;
                        if (MaxErr < CurErr) MaxErr = CurErr;       
                }
            }
        }
        loop++;
    }while (MaxErr>Conv); /* 領域全ての点の誤差がConvを下回ったらおしまい*/
    // printf("%5d %e\n", loop, MaxPhi);
}