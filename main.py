from datetime import datetime
print(datetime.now().replace(microsecond = 0),"START",__file__)

import numpy as np
from netCDF4 import Dataset
import os
from wrf import getvar, ALL_TIMES
# import time
import math

#soファイルの読み込み
from ctypes import *
import ctypes.util  
import sys
import gc

# バックグランドで実行するためのオーバーライド
_original_print = print
def print(*args, **kwargs):
    kwargs.setdefault('flush', True)
    _original_print(*args, **kwargs)

#soファイルを指定
lookup      = np.ctypeslib.load_library("lookup.so",".")
calc_charge = np.ctypeslib.load_library("calc_charge.so",".")
calc_adv    = np.ctypeslib.load_library("calc_adv.so",".")
calc_poi    = np.ctypeslib.load_library("calc_gauss.so",".")
flash       = np.ctypeslib.load_library("flash.so",".")

#計算に用いる定数
den_g=500
den_h=912
# den_i=100 #数密度が1e6が基本だから、混合比をこれで割る
# den_s=100


#解析領域を設定する
# ana_s_time_idx,ana_e_time_idx=0,1000
ana_s_time_idx,ana_e_time_idx=72,144+72
ana_s_z_idx,   ana_e_z_idx=0,40
ana_s_lat_idx, ana_e_lat_idx=100,250
ana_s_lon_idx, ana_e_lon_idx=150,300


path=os.path.join("/home1/nakamura_kento/WRF/work/ndown_TR_Nkyushu_20170705/wrfout_d03_2017-07-04_12:00:00_origin")
nc = Dataset(path, "r")
nc_var=nc.variables.keys()

u=nc.variables["U"]
in_time_idx,in_z_idx,in_lat_idx,in_lon_idx,=u.shape

ana_e_time_idx=min(ana_e_time_idx,in_time_idx)
ana_e_z_idx   =min(ana_e_z_idx   ,in_z_idx)
ana_e_lat_idx =min(ana_e_lat_idx ,in_lat_idx)
ana_e_lon_idx_f =ana_e_lon_idx
ana_e_lon_idx =min(ana_e_lon_idx ,in_lon_idx)

time_idx=ana_e_time_idx-ana_s_time_idx
z_idx   =ana_e_z_idx-ana_s_z_idx
lat_idx =ana_e_lat_idx-ana_s_lat_idx

# #何故か引き渡しのときにエラーが起きるので場合分け
if ana_e_lon_idx_f>in_lon_idx:
    lon_idx =ana_e_lon_idx-ana_s_lon_idx-1
else:
    lon_idx =ana_e_lon_idx-ana_s_lon_idx


#begin_時間の取得-------------------------------------------------------------------------------------------------------------
wrf_s_time = "" 
times = nc.variables["Times"][ana_s_time_idx:ana_e_time_idx]
for i in times[0]:
    wrf_s_time += i.decode()
wrf_s_time=datetime.strptime(wrf_s_time, '%Y-%m-%d_%H:%M:%S')#datetime型に変える
print("wrf_s_time",wrf_s_time)
#2つ目の時刻も求めて、1つ目の時間との差からtime stepを求める
wrf_2_time = "" 
for i in times[1]:
    wrf_2_time += i.decode()
wrf_2_time=datetime.strptime(wrf_2_time, '%Y-%m-%d_%H:%M:%S')#datetime型に変える


wrf_datetime=np.array([])
for j in range(time_idx):
    hoge_time=""
    for i in times[j]:
        hoge_time+=i.decode()
    hoge_time=str(hoge_time)
    hoge_time=str(hoge_time[:4]+hoge_time[5:7]+hoge_time[8:10]+hoge_time[11:13]+hoge_time[14:16]+hoge_time[17:29])
    hoge_time=int(hoge_time)
    hoge_time=np.array([hoge_time])
    wrf_datetime=np.append(wrf_datetime,hoge_time)


time_step=(wrf_datetime[1]-wrf_datetime[0])//10000*3600+((wrf_datetime[1]-wrf_datetime[0])%10000)//100*60+((wrf_datetime[1]-wrf_datetime[0]))%100
time_step=int(time_step)
print("time_step_output is",time_step,"seconds")
#end_時間の取得----------------------------------------------------------------------------------------------------------



print("get ua")
ua = getvar(nc, "ua", timeidx=ALL_TIMES)[ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
print(ua.shape, ua.__sizeof__())
gc.collect()

print("get va")
va = getvar(nc, "va", timeidx=ALL_TIMES)[ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
print(va.shape, va.__sizeof__())
gc.collect()

print("get wa")
wa = getvar(nc, "wa", timeidx=ALL_TIMES)[ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
print(wa.shape, wa.__sizeof__())
gc.collect()

print("get z")
z = getvar(nc, "z", timeidx=ALL_TIMES)[ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
print(z.shape, z.__sizeof__())
gc.collect()

# """
#begin水平方向の格子間隔を求める#################################################################################
print("get lat lon")
wrf_lat = nc.variables["XLAT"][0,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]#緯度の取得
wrf_lon = nc.variables["XLONG"][0,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]#経度の取得
#ランベルト法で座標を取っているので、少し複雑なのでgeopyを使う。
#現段階では水平方向の格子間隔は全てのグリットにおいて等しいと近似している。
from geopy.distance import geodesic
dx=geodesic((wrf_lat[lat_idx//2,lon_idx//2],wrf_lon[lat_idx//2,lon_idx//2]),(wrf_lat[lat_idx//2,lon_idx//2],wrf_lon[lat_idx//2,lon_idx//2+1])).m
dy=geodesic((wrf_lat[lat_idx//2,lon_idx//2],wrf_lon[lat_idx//2,lon_idx//2]),(wrf_lat[lat_idx//2+1,lon_idx//2],wrf_lon[lat_idx//2,lon_idx//2])).m
print("dx is "+str('{:.5g}'.format(dx))+" m")
print("dy is "+str('{:.5g}'.format(dy))+" m")
#end水平方向の格子間隔を求める#################################################################################
#"""
# dx=1000.1
# dy=996.92


#Registryでalt(密度の逆数を出力している場合はこれを使う)#######################################################
if 'ALT' in nc_var:
    in_rho_air=nc.variables['ALT'][ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
    rho_air= in_rho_air**(-1)
else:
    rho_air=0.584
print("rho_air",np.mean(rho_air))
#############################################################################################################

# """
#wrfoutから値の読み取り
# qv=nc.variables['QVAPOR'][ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
# qr= nc.variables['QRAIN'][ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
qc=nc.variables['QCLOUD'][ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
qi=  nc.variables['QICE'][ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
qs= nc.variables['QSNOW'][ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
qg=nc.variables['QGRAUP'][ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
if 'QHAIL' in nc_var:
    qh=nc.variables['QHAIL'][ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
else:
    qh=np.zeros([time_idx,z_idx,lat_idx,lon_idx])
temp=nc.variables['T'][ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
# """

# """
#ラムダ(スロープパラメーター)を計算する
#'QGRAUP'########################################################################################################
print("LAMBDA_G")
LAMBDAgmax=6.e4
avtg,bvtg=330.,0.8
n0g=4.e6
qg=np.where(qg<1.e-9,1.e-9,qg)
LAMBDA_G=((math.pi*den_g*n0g)/(qg*rho_air))**0.25
LAMBDA_G=np.where(LAMBDA_G>LAMBDAgmax,LAMBDAgmax,LAMBDA_G)
LAMBDA_G=np.where(LAMBDA_G>0.1,LAMBDA_G,0.1)
# velocity=(math.gamma(4+bvtg)/6) * (avtg/LAMBDA_G**bvtg)
#1粒子あたりの質量を求める
mass_g=math.pi*qg/LAMBDA_G**3
print(str(np.mean(mass_g))+"(kg/one particle)")

#'QHAIL'########################################################################################################
print("LAMBDA_H")
if 'QHAIL' in nc_var:
    LAMBDAhmax=2.e4
    avth,bvth=285.,0.8
    n0h=4.e4
    qh=np.where(qh<1.e-9,1.e-9,qh)
    LAMBDA_H=((math.pi*den_h*n0h)/(qh*rho_air))**0.25
    LAMBDA_H=np.where(LAMBDA_H>LAMBDAhmax,LAMBDAhmax,LAMBDA_H)
    LAMBDA_H=np.where(LAMBDA_H>0.1,LAMBDA_H,0.1)
    # velocity=(math.gamma(4+bvth)/6) * (avth/LAMBDA_H**bvth)
    #1粒子あたりの質量を求める
    mass_h=math.pi*qh/LAMBDA_H**3
else:
    mass_h=np.ones([time_idx,z_idx,lat_idx,lon_idx])
print(str(np.mean(mass_h))+"(kg/one particle)")

#'QICE'########################################################################################################
print("ICE")
mass_i=np.where(qi<1e-6,1e6,qi)*rho_air/1e6
print(str(np.mean(mass_i))+"(kg/one particle)")
# """

#変換項の読み込み##############################################################################################################

# if 'QT_PCACT' in nc_var:
#     qt_pcact=nc.variables['QT_PCACT'][ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
# else:
#     qt_pcact=np.zeros([time_idx,z_idx,lat_idx,lon_idx])
# # qt_pcact=np.random.rand(time_idx,z_idx,lat_idx,lon_idx)  #乱数

if 'QT_PHACI' in nc_var:
    qt_phaci=nc.variables['QT_PHACI'][ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
else:
    qt_phaci=np.zeros([time_idx,z_idx,lat_idx,lon_idx])

if 'QT_PHACS' in nc_var:
    qt_phacs=nc.variables['QT_PHACS'][ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
else:
    qt_phacs=np.zeros([time_idx,z_idx,lat_idx,lon_idx])

if 'QT_PGACI' in nc_var:
    qt_pgaci=nc.variables['QT_PGACI'][ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
else:
    qt_pgaci=np.zeros([time_idx,z_idx,lat_idx,lon_idx])

if 'QT_PGACS' in nc_var:
    qt_pgacs=nc.variables['QT_PGACS'][ana_s_time_idx:ana_e_time_idx,ana_s_z_idx:ana_e_z_idx,ana_s_lat_idx:ana_e_lat_idx,ana_s_lon_idx:ana_e_lon_idx]
else:
    qt_pgacs=np.zeros([time_idx,z_idx,lat_idx,lon_idx])

nc.close()
#################################################################################################################################

#長くなるので4次元のポインターの書式を定義しておく
dim4_cp=ctypes.POINTER(ctypes.c_double*time_idx*z_idx*lat_idx*lon_idx)

#電荷分布の初期化
charge_a=np.zeros([time_idx,z_idx,lat_idx,lon_idx]).astype(np.float64)
charge_a_c = charge_a.ctypes.data_as(dim4_cp).contents  

charge_g=np.zeros([time_idx,z_idx,lat_idx,lon_idx]).astype(np.float64)
charge_g_c = charge_g.ctypes.data_as(dim4_cp).contents  

charge_h=np.zeros([time_idx,z_idx,lat_idx,lon_idx]).astype(np.float64)
charge_h_c = charge_h.ctypes.data_as(dim4_cp).contents  

charge_i=np.zeros([time_idx,z_idx,lat_idx,lon_idx]).astype(np.float64)
charge_i_c = charge_i.ctypes.data_as(dim4_cp).contents  

charge_s=np.zeros([time_idx,z_idx,lat_idx,lon_idx]).astype(np.float64)
charge_s_c = charge_s.ctypes.data_as(dim4_cp).contents  


#静電ポテンシャルの初期化
phi=np.zeros([time_idx,z_idx,lat_idx,lon_idx]).astype(np.float64)
phi_c  =phi.ctypes.data_as(dim4_cp).contents  

#電荷分離のファクターの初期化
fc=np.zeros([time_idx,z_idx,lat_idx,lon_idx]).astype(np.float64)
fc_c  =fc.ctypes.data_as(dim4_cp).contents  

#発雷回数の初期化
fod=np.zeros([time_idx,z_idx,lat_idx,lon_idx]).astype(np.float64)
fod_c  =fod.ctypes.data_as(dim4_cp).contents  


# 体積当たりの雲水(g/m^3)
CWpV=qc*rho_air*1000
CWpV = CWpV.astype(np.float64)
CWpV_c  =CWpV.ctypes.data_as(dim4_cp).contents  

#int型もctypeのc_int型へ変換して渡す
time_idx_c = ctypes.c_int(time_idx)
z_idx_c    = ctypes.c_int(z_idx   )
lat_idx_c  = ctypes.c_int(lat_idx )
lon_idx_c  = ctypes.c_int(lon_idx )

time_step_c= ctypes.c_int(time_step)
dx_c= ctypes.c_double(dx)
dy_c= ctypes.c_double(dy)

#numpy配列は一度1次元の配列に戻して、ポインター渡しを行う
qc=qc.astype(np.float64)
qc_c= qc.ctypes.data_as(dim4_cp).contents 

qg=qg.astype(np.float64)
qg_c= qg.ctypes.data_as(dim4_cp).contents  

qi=qi.astype(np.float64)
qi_c= qi.ctypes.data_as(dim4_cp).contents  

qh=qh.astype(np.float64)
qh_c= qh.ctypes.data_as(dim4_cp).contents  

mass_g=mass_g.astype(np.float64)
mass_g_c= mass_g.ctypes.data_as(dim4_cp).contents  

mass_h=mass_h.astype(np.float64)
mass_h_c= mass_h.ctypes.data_as(dim4_cp).contents  

mass_i=mass_i.astype(np.float64)
mass_i_c= mass_i.ctypes.data_as(dim4_cp).contents  

temp=temp.astype(np.float64)
temp_c= temp.ctypes.data_as(dim4_cp).contents  

qt_phaci=qt_phaci.astype(np.float64)
qt_phaci_c= qt_phaci.ctypes.data_as(dim4_cp).contents  

qt_phacs=qt_phacs.astype(np.float64)
qt_phacs_c= qt_phacs.ctypes.data_as(dim4_cp).contents  

qt_pgaci=qt_pgaci.astype(np.float64)
qt_pgaci_c= qt_pgaci.ctypes.data_as(dim4_cp).contents  

qt_pgacs=qt_pgacs.astype(np.float64)
qt_pgacs_c= qt_pgacs.ctypes.data_as(dim4_cp).contents  

#wrf-pythonで得たデータは一度numpy配列に戻すことで同様に扱うことが出来る
ua=np.array(ua).astype(np.float64)
ua_c =ua.ctypes.data_as(dim4_cp).contents  

va=np.array(va).astype(np.float64)
va_c =va.ctypes.data_as(dim4_cp).contents   

wa=np.array(wa).astype(np.float64)
wa_c =wa.ctypes.data_as(dim4_cp).contents   

z=np.array(z).astype(np.float64) #m単位のデータ
z_c  =z.ctypes.data_as(dim4_cp).contents  


print(time_idx_c, z_idx_c, lat_idx_c, lon_idx_c,time_step_c)


#時間ごとのループはpythonで回す
for time_loop in range(time_idx):
# for time_loop in range(30,50):
# for time_loop in range(3):
    print("\ntime_step is "+str(time_loop+1)+"/"+str(time_idx))
    #時間の表示
    tmp = ""
    for i in times[time_loop]:
        tmp += i.decode()
    print(tmp)
    time_c=ctypes.c_int(time_loop)

    if time_loop>0:
        charge_g[time_loop,:,:,:]=charge_g[time_loop-1,:,:,:]
        charge_h[time_loop,:,:,:]=charge_h[time_loop-1,:,:,:]
        charge_i[time_loop,:,:,:]=charge_i[time_loop-1,:,:,:]
        charge_s[time_loop,:,:,:]=charge_s[time_loop-1,:,:,:]

    # """
    # begin電荷分離のファクターを求める####################################################################################################
    #降水粒子(氷晶と霰)および地面のみが電荷を持つことが出来るので、これで種類を分ける必要がある。
    #本来、電荷の総和は常に0にならなければならない
    #関数の引数の型を指定(ctypes)
    lookup.takahashi.argtypes = [
        dim4_cp,dim4_cp,dim4_cp,
        c_int32,c_int32,c_int32,c_int32
        ]
    #関数が返す値の型を指定(ctypes)
    lookup.takahashi.restype =  None
    #cのモジュールに投げる
    lookup.takahashi(CWpV_c,temp_c,fc_c, z_idx_c, lat_idx_c, lon_idx_c,time_c)
    #end電荷分離のファクターを求める######################################################################################################
    # """

    # mini_time_step=6 #計算に与えるタイムステップ
    mini_time_step=int(dx/1000*6) #計算に与えるタイムステップ
    #wrfoutのタイムステップを用いると計算不安定を引き起こすのでここでループ
    print("mini time step is "+str(mini_time_step))
    mini_time_step_c=ctypes.c_int(mini_time_step)
    for time_loop_mini in range(time_step//mini_time_step):


        # """
        # begin電荷分離####################################################################################################
        #降水粒子(氷晶と霰)および地面のみが電荷を持つことが出来るので、これで種類を分ける必要がある。
        #本来、電荷の総和は常に0にならなければならない
        #関数の引数の型を指定(ctypes)
        calc_charge.one.argtypes = [
            dim4_cp,dim4_cp,dim4_cp,dim4_cp,dim4_cp,
            c_int32,c_int32,c_int32,c_int32,c_int32,c_int32,
            dim4_cp,dim4_cp,dim4_cp,dim4_cp,dim4_cp,dim4_cp,dim4_cp,dim4_cp,dim4_cp
            ]
        #関数が返す値の型を指定(ctypes)
        calc_charge.one.restype =  None
        #cのモジュールに投げる
        calc_charge.one(charge_g_c,charge_h_c,charge_i_c,charge_s_c,fc_c, time_idx_c, z_idx_c, lat_idx_c, lon_idx_c,time_c,mini_time_step_c
        ,qt_phaci_c,qt_phacs_c,qt_pgaci_c,qt_pgacs_c
        ,qc_c,temp_c,mass_g_c,mass_h_c,mass_i_c)
        #end電荷分離######################################################################################################
        # """


        # """
        #begin移流#############################################################################################
        #関数の引数の型を指定(ctypes)
        calc_adv.one.argtypes = [
            dim4_cp,dim4_cp,dim4_cp,dim4_cp,
            c_int32,c_int32,c_int32,c_int32,c_int32,c_int32,
            c_double,c_double,
            dim4_cp,dim4_cp,dim4_cp,dim4_cp
            ]
        #関数が返す値の型を指定(ctypes)
        calc_adv.one.restype =  None
        #cのモジュールに投げる
        calc_adv.one(
            charge_g_c,charge_h_c,charge_i_c,charge_s_c,
            time_idx_c,z_idx_c,lat_idx_c,lon_idx_c,time_c,mini_time_step_c,
            dx_c,dy_c,
            ua_c,va_c,wa_c,z_c)
        #end移流##############################################################################################
        # """

        # """
        #begin静電ポテンシャルを求める
        ######################################################################################
        charge_a[time_loop,:,:,:]=charge_g[time_loop,:,:,:]+charge_h[time_loop,:,:,:]+charge_i[time_loop,:,:,:]+charge_s[time_loop,:,:,:]
        # charge_a=np.random.rand(time_idx,z_idx,lat_idx,lon_idx)  #乱数
        # charge_a_c = charge_a.ctypes.data_as(dim4_cp).contents  
        #関数の引数の型を指定(ctypes)　
        calc_poi.calc_gauss.argtypes = [
            dim4_cp,dim4_cp,
            c_int32,c_int32,c_int32,c_int32,c_int32,
            c_double,c_double,
            dim4_cp
            ]
        #関数が返す値の型を指定(ctypes)
        calc_poi.calc_gauss.restype =  None
        #cのモジュールに投げる
        calc_poi.calc_gauss(
            charge_a_c,phi_c,
            time_idx_c,z_idx_c,lat_idx_c,lon_idx_c,time_c,
            dx_c,dy_c,
            z_c)
        #end静電ポテンシャルを求める####################################################################################################
        # """

        # """
        #begin放電で中和を行う################################################################################################
        #関数の引数の型を指定(ctypes)
        charge_a[time_loop,:,:,:]=charge_g[time_loop,:,:,:]+charge_h[time_loop,:,:,:]+charge_i[time_loop,:,:,:]+charge_s[time_loop,:,:,:]
        flash.one.argtypes = [
            dim4_cp,dim4_cp,dim4_cp,dim4_cp,dim4_cp,dim4_cp,dim4_cp,
            c_int32,c_int32,c_int32,c_int32,
            c_double,c_double,
            dim4_cp
            ]
        #関数が返す値の型を指定(ctypes)
        flash.one.restype =  None
        #cのモジュールに投げる
        flash.one(
            charge_a_c,charge_g_c,charge_h_c,charge_i_c,charge_s_c,phi_c,fod_c,
            z_idx_c,lat_idx_c,lon_idx_c,time_c,
            dx_c,dy_c,
            z_c)
        #end放電で中和を行う####################################################################################################
        # """
    
    print("charge_g(C/m^3) "+f'{np.nanmin(charge_g[time_loop,:,:,:]):.2e}'+" "+f'{np.nanmax(charge_g[time_loop,:,:,:]):.2e}')
    print("charge_h(C/m^3) "+f'{np.nanmin(charge_h[time_loop,:,:,:]):.2e}'+" "+f'{np.nanmax(charge_h[time_loop,:,:,:]):.2e}')
    print("charge_i(C/m^3) "+f'{np.nanmin(charge_i[time_loop,:,:,:]):.2e}'+" "+f'{np.nanmax(charge_i[time_loop,:,:,:]):.2e}')
    print("charge_s(C/m^3) "+f'{np.nanmin(charge_s[time_loop,:,:,:]):.2e}'+" "+f'{np.nanmax(charge_s[time_loop,:,:,:]):.2e}')
    print("charge_a(C/m^3) "+f'{np.nanmin(charge_a[time_loop,:,:,:]):.2e}'+" "+f'{np.nanmax(charge_a[time_loop,:,:,:]):.2e}')
    print("charge_mean(C/m^3) "+f'{np.mean(charge_a[time_loop,:,:,:]):.2e}')
    print("charge_sum(C) " +f'{np.nansum(charge_a[time_loop,:,:,:]):.2e}')
    print("factor  (fC/) " +f'{np.nanmin(fc[time_loop,:,:,:]):.2e}'+" "+f'{np.nanmax(fc[time_loop,:,:,:]):.2e}')
    print("phi(V) "+f'{np.nanmin(phi[time_loop,:,:,:]):.2e}'+" "+f'{np.nanmax(phi[time_loop,:,:,:]):.2e}')
    print("fod(/grid/time_step) "+f'{np.max(fod[time_loop,:,:,:]):.2e}'+" "+f'{np.mean(fod[time_loop,:,:,:]):.2e}')
    print(datetime.now().replace(microsecond = 0))
    # time.sleep(0.3)




    # if (time_loop%30==0)|(time_loop==time_idx-1): #適宜上書きする。やりすぎると時間がかかる。
        # charge=np.array(charge_a).reshape(time_idx,z_idx,lat_idx,lon_idx) #多分不要

# """
#netCDFで書き出す####################################################################
#V3で書かれているのでV4のものにする
print("write netCDF4")
nc = Dataset('out.nc', 'w', format="NETCDF4")

nc.createDimension('Time', time_idx) 
nc.createDimension('south_north', lat_idx)           
nc.createDimension('west_east', lon_idx)           
nc.createDimension('bottom_top', z_idx)            

# Times = nc.createVariable('Times', 'S1', ('Time',19))
Times = nc.createVariable('Times', 'f8', ('Time'))
Times.long_name = 'time of test variable'
Times.units = 'days since 1968-05-23 00:00:00'

XLONG = nc.createVariable('XLONG',"f8", ('south_north','west_east'))
XLONG.long_name = 'east longitude'
XLONG.units = 'degree of east longitude'

XLAT = nc.createVariable('XLAT',"f8", ('south_north','west_east'))
XLAT.long_name = 'north latitude'
XLAT.units = 'degree of north latitude'

charge_nc = nc.createVariable('CHARGE',"f8", ('Time','bottom_top','south_north','west_east'))
charge_nc.long_name = 'charge'
charge_nc.units = 'C/m^3'

charge_g_nc = nc.createVariable('CHARGE_G',"f8", ('Time','bottom_top','south_north','west_east'))
charge_g_nc.long_name = 'charge_G'
charge_g_nc.units = 'C/m^3'

charge_h_nc = nc.createVariable('CHARGE_h',"f8", ('Time','bottom_top','south_north','west_east'))
charge_h_nc.long_name = 'charge_H'
charge_h_nc.units = 'C/m^3'

charge_i_nc = nc.createVariable('CHARGE_I',"f8", ('Time','bottom_top','south_north','west_east'))
charge_i_nc.long_name = 'charge_I'
charge_i_nc.units = 'C/m^3'

charge_s_nc = nc.createVariable('CHARGE_S',"f8", ('Time','bottom_top','south_north','west_east'))
charge_s_nc.long_name = 'charge_S'
charge_s_nc.units = 'C/m^3'


fod_nc = nc.createVariable('FOD',"f8", ('Time','bottom_top','south_north','west_east'))
fod_nc.long_name = 'FOD'
fod_nc.units = ''

fc_nc = nc.createVariable('FC',"f8", ('Time','bottom_top','south_north','west_east'))
fc_nc.long_name = 'FACTOR'
fc_nc.units = 'fC/collision'

# 予め np.ndarrayで作成しておいた値を代入する
Times[:] = wrf_datetime
XLONG[:,:] = wrf_lon
XLAT[:,:] = wrf_lat
charge_nc[:,:,:,:] =charge_a
charge_g_nc[:,:,:,:] =charge_g
charge_h_nc[:,:,:,:] =charge_h
charge_i_nc[:,:,:,:] =charge_i
charge_s_nc[:,:,:,:] =charge_s
fod_nc[:,:,:,:] = fod
fc_nc[:,:,:]=fc

nc.close()
#netCDFで書き出す####################################################################
#"""

from datetime import datetime
print(datetime.now().replace(microsecond = 0),"END",__file__)
