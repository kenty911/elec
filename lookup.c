#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//  gcc -cpp -fPIC -shared lookup.c -lm -o lookup.so -O3

void takahashi(double *cw,double *temp, double *fc, int z_idx, int lat_idx, int lon_idx, int time) {
    // double *table;    
    int h, i, j, k, p, q;
    int cw_idx, temp_idx;

    #define temp(h,i,j,k) temp[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
    #define cw(h,i,j,k)     cw[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
    #define fc(h,i,j,k)     fc[(h)*z_idx*lat_idx*lon_idx+(i)*lat_idx*lon_idx+(j)*lon_idx+(k)]
    
    // 30(cw)*31(temp)
    static double table[30][31]={
        {31.35,31.02,30.69,30.36,30.03,29.7,29.37,29.04,28.71,28.38,28.05,27.72,27.39,27.06,26.73,26.73,26.4,26.07,25.41,24.75,24.42,24.42,24.09,24.09,23.76,23.76,23.76,23.76,23.43,23.43,23.43},
        {31.35,31.35,31.02,30.69,30.03,29.7,29.04,28.71,28.38,28.05,27.72,27.39,27.06,26.73,26.4,26.4,25.74,25.41,25.08,24.75,24.42,24.09,23.76,23.43,23.1,23.1,22.77,22.44,22.11,21.78,21.45},
        {31.35,31.35,30.36,30.03,30.03,29.04,28.71,28.05,27.39,27.06,26.73,26.4,25.41,24.75,24.42,24.09,23.43,22.77,22.11,21.45,20.46,20.13,19.8,19.8,19.47,19.47,19.14,18.81,18.15,18.15,18.15},
        {31.35,31.35,31.02,31.02,30.36,30.03,29.37,29.04,28.38,27.39,26.73,25.74,25.41,24.42,23.76,23.43,22.44,21.45,20.79,20.46,19.8,19.47,19.47,18.81,18.48,18.15,17.82,17.49,17.49,17.16,17.16},
        {31.68,31.35,31.02,30.36,30.03,29.04,28.71,28.38,27.39,27.06,26.73,25.41,24.42,23.76,23.1,22.44,21.12,20.13,19.8,19.47,18.81,18.48,18.15,17.82,17.16,16.83,16.5,16.5,16.5,14.85,14.85},
        {32.01,32.01,31.68,30.69,30.36,29.37,28.71,28.05,27.39,27.06,26.07,25.08,23.76,23.1,21.78,20.79,19.47,19.14,18.81,18.48,17.82,17.16,16.83,16.5,16.5,14.85,14.85,13.2,11.55,9.9,9.9},
        {32.01,32.01,31.68,31.35,30.69,29.37,28.71,27.72,27.06,26.4,25.08,24.09,23.1,21.45,19.8,19.47,19.14,17.82,17.49,17.16,16.5,14.85,13.2,11.55,8.25,6.6,6.6,6.6,6.6,6.6,6.6},
        {32.34,32.34,32.01,31.68,30.36,29.37,28.38,27.39,26.73,25.41,23.76,22.77,20.13,19.14,18.48,17.49,16.83,16.5,13.2,6.6,4.95,4.29,3.3,3.3,2.97,2.97,2.64,2.64,2.31,2.31,2.31},
        {36.3,33,32.01,31.35,30.69,29.7,29.71,27.39,25.74,23.43,21.45,19.8,18.15,16.83,14.85,4.95,3.3,2.64,1.65,1.32,0.66,0,0,0,0,0,0,0,0,0,0},
        {39.6,36.3,33,31.35,30.69,30.03,28.05,26.4,23.1,20.46,18.15,16.5,3.3,0,-0.33,-1.32,-2.97,-3.3,-3.3,-3.3,-3.3,-3.3,-6.6,-9.9,-6.6,-6.6,-6.6,-6.6,-4.95,-4.95,-3.3},
        {49.5,49.5,42.9,42.9,33,31.35,28.05,24.75,16.5,0.99,-6.6,-16.5,-19.8,-26.4,-33,-33,-29.7,-28.1,-26.4,-24.8,-23.1,-21.5,-19.8,-18.2,-17.5,-16.5,-16.5,-16.5,-15.8,-9.99,-9.99},
        {52.8,52.8,56.1,59.4,52.8,49.5,42.9,33,19.8,9.9,0,-4.95,-11.55,-13.2,-14.85,-14.85,-14.85,-13.2,-13.2,-13.2,-13.2,-13.2,-13.2,-11.55,-11.55,-9.9,-9.9,-8.25,-4.95,-4.95,-3.3},
        {49.5,52.8,56.1,59.4,62.7,66,52.8,46.2,39.6,26.4,3.3,-0.99,-9.9,-9.9,-9.9,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.99,-9.9,-9.9,-8.25,-4.95,-3.63,-2.97},
        {46.2,52.8,56.1,59.4,62.7,56.1,42.9,39.6,33,19.8,9.9,0.99,-1.98,-4.29,-6.6,-8.25,-9.99,-9.99,-9.99,-9.99,-9.9,-9.9,-9.9,-9.9,-9.9,-9.9,-8.25,-6.6,-3.3,-2.97,-2.97},
        {42.9,49.5,52.8,56.1,59.4,56.1,46.2,42.9,36.3,23.1,16.5,3.3,0,-1.98,-3.3,-3.3,-4.95,-6.6,-6.6,-6.6,-6.6,-6.6,-6.6,-6.6,-6.6,-4.95,-4.29,-3.3,-2.97,-2.64,-2.64},
        {33,36.3,39.6,46.2,52.8,56.1,52.8,46.2,39.6,29.7,21.45,16.5,2.97,0,-1.32,-2.31,-2.97,-3.3,-3.3,-3.3,-3.3,-3.3,-3.3,-3.3,-3.3,-3.3,-3.3,-3.3,-2.97,-2.97,-2.64},
        {31.02,33,36.3,39.6,47.2,49.5,49.5,42.9,39.6,33,26.4,21.45,16.5,3.3,0.33,-0.33,-0.99,-1.65,-1.98,-2.31,-2.64,-2.64,-2.64,-2.64,-2.64,-2.64,-2.64,-2.64,-2.64,-2.31,-2.31},
        {26.4,29.7,32.34,36.3,39.6,39.6,42.9,42.9,39.6,36.3,33,28.71,26.4,22.44,19.47,9.9,0,-0.33,-0.99,-1.32,-1.98,-2.31,-2.31,-2.31,-1.98,-1.98,-1.98,-1.98,-1.98,-1.98,-1.65},
        {19.8,23.1,26.4,28.71,31.35,33,36.3,36.3,33,31.02,28.71,27.39,26.4,23.43,20.79,17.82,9.9,3.3,1.65,0,-0.66,-0.66,-0.66,-0.99,-0.99,-1.32,-1.32,-1.32,-1.32,-1.32,-1.32},
        {13.2,15.18,17.16,19.8,22.44,24.09,25.08,26.4,27.39,27.39,27.06,26.73,25.74,24.42,23.1,21.45,19.47,17.16,9.9,6.6,3.3,2.64,0.99,0,0,-0.33,-0.66,-0.66,-0.66,-0.66,-0.66},
        {4.95,6.93,9.24,10.56,11.88,13.2,14.85,16.5,16.5,18.15,19.14,18.81,18.81,18.48,17.82,16.83,16.17,15.51,14.85,13.86,13.2,11.88,10.56,9.57,8.58,7.59,6.6,3.63,2.97,2.31,1.65},
        {3.3,6.6,7.59,9.9,11.22,12.21,13.2,14.85,15.84,16.5,16.5,16.83,16.83,16.83,16.5,16.5,15.84,15.18,14.85,14.19,13.2,12.21,10.89,9.9,8.91,7.92,6.93,4.95,3.3,2.31,1.98},
        {1.65,4.95,6.6,8.25,9.9,11.55,13.2,14.85,15.18,15.51,15.84,15.84,16.17,15.84,15.84,15.51,15.18,14.85,14.52,13.86,12.87,11.88,10.89,10.23,9.24,7.59,7.26,6.6,4.62,2.31,1.65},
        {0,2.64,4.95,7.26,9.24,10.23,11.22,11.22,11.88,12.54,13.2,13.86,14.85,15.18,14.85,14.85,14.52,14.19,13.86,13.2,12.54,11.88,10.89,10.23,9.24,8.25,7.59,6.6,4.29,2.64,0.99},
        {0,0.33,3.3,5.61,7.59,8.91,9.9,10.23,12.21,12.87,13.53,13.53,13.53,13.86,13.86,13.53,13.53,13.2,13.2,12.54,11.88,11.55,10.89,9.9,9.24,8.25,7.26,6.6,4.29,3.3,0.66},
        {0,0,0.33,2.97,5.61,7.26,8.91,9.9,10.56,11.22,11.88,11.88,12.21,12.54,12.54,12.54,12.21,12.21,11.88,11.55,11.22,10.56,9.9,9.24,8.91,7.59,6.93,5.61,3.3,0.99,0},
        {0,0,0,0.33,2.64,4.29,6.6,7.59,9.24,9.57,10.23,10.89,11.22,11.22,11.22,11.22,10.89,10.89,10.89,10.56,9.9,9.57,8.91,8.25,7.59,6.6,5.28,3.63,0.99,0,0},
        {0,0,0,0,0.33,0.33,2.31,4.29,6.27,6.93,7.59,8.25,8.91,9.24,9.57,9.57,9.57,9.57,8.91,8.58,8.25,7.59,6.93,6.27,5.61,4.29,2.64,0.33,0,0,0},
        {0,0,0,0,0,0,0,0,0,0.66,2.64,3.3,4.62,5.28,5.94,5.94,6.27,5.94,5.61,5.28,4.95,4.29,3.96,2.97,0.66,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
        };

    static double cw_pt[30]={
        30.00,20.00,10.00,9.00 ,8.00 ,7.00 ,6.00 ,5.00 ,4.00 ,3.00 ,2.00 ,1.00 ,0.90 ,0.80 ,0.70 ,0.60 ,0.50 ,0.40 ,0.30 ,0.20 ,0.10 ,0.09 ,0.08 ,0.07 ,0.06 ,0.05 ,0.04 ,0.03 ,0.02 ,0.01 
        };
    static double temp_pt[31]={
        0.00,-1.00,2.00,3.00,4.00,5.00,6.00,7.00,8.00,9.00,10.00,-11.00,-12.00,-13.00,-14.00,-15.00,-16.00,-17.00,-18.00,-19.00,-20.00,-21.00,-22.00,-23.00,-24.00,-25.00,-26.00,-27.00,-28.00,-29.00,-30.00
        };
    
    // #時間ごとのループはpythonで回す
    for(i=0; i<z_idx; i++){
      for(j=0; j<lat_idx; j++){
        for(k=0; k<lon_idx; k++){
            cw_idx=0;
            while(cw_pt[cw_idx]>=cw(time,i,j,k) && cw_idx<30) cw_idx++;
            // while(cw_pt[cw_idx]<=cw(time,i,j,k) && cw_idx<30) cw_idx++;
            temp_idx=0;
            while(temp_pt[temp_idx]>=temp(time,i,j,k) && temp_idx<31) temp_idx++;
            // while(temp_pt[temp_idx]<=temp(time,i,j,k) && temp_idx<31) temp_idx++;
          fc(time,i,j,k)=table[cw_idx][temp_idx];
        }
      }
    }
    }


