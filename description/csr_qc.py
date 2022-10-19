from datetime import datetime 
print(datetime.now().replace(microsecond = 0))

import cartopy.crs as ccrs
import cartopy.util as cutil
from cartopy.mpl.ticker import * #LongitudeFormatter, LatitudeFormatter
import cartopy.feature as feature #陸地や海洋を色づけ

from matplotlib.colors import *
from matplotlib.dates import *
from matplotlib.ticker import *
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.cm as cm
from netCDF4 import *
from wrf import *
import numpy as np
import sys#nohup実行時にnohup.outを随時更新するため
import os


path=os.path.join("out1.nc") #パスの指定

output_dir = os.path.join("./csr_qc")
os.makedirs(output_dir, exist_ok=True)

nc = Dataset(path)
nc_var=nc.variables.keys()
xlat=nc.variables["XLAT"][:]
xlong=nc.variables["XLONG"][:]

# rho=nc.variables["FC"][:]
rho=nc.variables["FO"][:]
idx_time_full,idx_z_full,idx_lat_full,idx_lon_full=rho.shape

path_a=os.path.join("/home1/nakamura_kento/WRF/work/TR_Nkyushu_20170705/WDM7/wrfout_d03_2017-07-04_12:00:00")
nc_a = Dataset(path_a, "r")
rho=nc_a.variables["QCLOUD"][:,:,100:250,150:300]
alt=nc_a.variables["ALT"][:,:,100:250,150:300]
rho=rho*1000/alt

d_lat_idx=70 #指定するlatのidx

lev=np.arange(0,16001,400) #描写の高度を指定する(最小値は地表面を描くために0にすべし)

lev_cont_rho=np.arange(0, 10, 0.5) #範囲の設定
# cmap = cm.bwr
# cmap_data = cmap(np.arange(cmap.N))
# cmap_data[127-5:127+5,0] = 0.999 # 0 のときのα値を0(透明)にする
# cmap_data[127-5:127+5,1] = 0.999 # 0 のときのα値を0(透明)にする
# customized_cmap = ListedColormap(cmap_data)

img_num=0
for time in range(idx_time_full):
    img_num+=1
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1,1,1)

    #時間の取得
    tmp = ""
    times = nc_a.variables["Times"]
    for i in times[time]:
        tmp += i.decode()
    print(tmp)


    #変数を得る
    ht =  getvar(nc_a, "z", timeidx=time)[:,100:250,150:300]
    ht=np.array(ht)
    # ter = getvar(nc_a, "ter", timeidx=time)
    # dbz = getvar(nc, "dbz", timeidx=time)

    rho_show= interplevel(rho[time,:,:,:],ht,lev)
    # ter_show= interplevel(ter[:,d_lat_idx,:],ht[:,d_lat_idx,:],lev)

    shade_plot = ax.contourf(np.arange(150),lev,rho_show[:,d_lat_idx,:],levels=lev_cont_rho, cmap=cm.PuRd, extend='max',alpha = 0.7)
    # ax.contourf(ter,lev,levels=lev_cont_rho, colors =["purple"], extend='max',alpha = 0.3)


    ax.set_ylim(0, 16000)
    ax.set_ylabel("Height (m)", fontsize=12)

    #カラーバーを描く
    cax = ax.inset_axes([1.04, 0.22, 0.03, 0.8], transform=ax.transAxes)#cbarの座標
    cbar = plt.colorbar(shade_plot, cax = cax, orientation='vertical',extendfrac = 'auto', ticks=lev_cont_rho )
    cbar.ax.tick_params(labelsize=8)

    fig.suptitle(tmp, fontsize=14)
    fig.savefig(output_dir+"/"+str(img_num).zfill(3)+".png", dpi=100, bbox_inches='tight', pad_inches=0.1)
    print('save_'+str(img_num).zfill(3)+".png"),sys.stdout.flush()
    plt.gca().clear()
    plt.close()

os.system("convert -delay 50 -loop 0 "+os.path.join(output_dir,"*.png")+" csr_qc.gif")