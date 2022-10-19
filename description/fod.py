#wrfoutからdbzを描くスクリプト
from datetime import datetime
print(datetime.now().replace(microsecond = 0),"START",__file__)

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


path=os.path.join("out.nc") #パスの指定

# output_dir = os.path.join("./charge")
output_dir = os.path.join("./fod")
os.makedirs(output_dir, exist_ok=True)

nc = Dataset(path)
nc_var=nc.variables.keys()
xlat=nc.variables["XLAT"][:]
xlong=nc.variables["XLONG"][:]

fod=nc.variables["FOD"][:,0,:,:]
# fod=nc.variables["CHARGE"][:,1,:,:]
# fod=fod*1e12
idx_time_full,idx_lat_full,idx_lon_full=fod.shape

path_a=os.path.join("/home1/nakamura_kento/WRF/work/TR_Nkyushu_20170705/WDM7/wrfout_d03_2017-07-04_12:00:00")
nc_a = Dataset(path_a, "r")

img_num=0
for time in range(idx_time_full):
    img_num+=1
    fig = plt.figure(figsize=(10,10))

    #時間の取得
    tmp = ""
    times = nc_a.variables["Times"][:]
    for i in times[time]:
        tmp += i.decode()
    print(tmp)


    # 描画部分---------------------------------------------------------------------------------------------------
    #各種指定
    fontsize = 8
    nl, sl, el, wl = np.min(xlat), np.max(xlat), np.min(xlong), np.max(xlong)
    proj = ccrs.PlateCarree()
    fig.subplots_adjust(right=0.85)#右に隙間を空ける(多分)
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    ax.coastlines(linewidths=0.3, zorder=0, resolution='10m',) 

    angle=0.25    #決めた角度ごとに線を引く
    xticks = np.arange(0, 360.1, angle)
    yticks = np.arange(-90, 90.1, angle)
    ax.set_xticks(xticks, crs=proj)
    ax.set_yticks(yticks, crs=proj)

    ax.xaxis.set_major_formatter(LongitudeFormatter(zero_direction_label=True))
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.xaxis.set_minor_locator(MultipleLocator(2.5))
    ax.yaxis.set_minor_locator(MultipleLocator(2.5))

    ax.set_extent([wl, el, sl, nl], crs=ccrs.PlateCarree())    #緯度経度の範囲設定

    ax.tick_params(direction='out', length=3, width=0.7, colors='k',grid_color='k', grid_alpha=0.5, labelsize=fontsize)
    ax.tick_params(which='minor', axis='x', length=2, color='k')
    ax.tick_params(which='minor', axis='y', length=2, color='k')

    gl = ax.gridlines(crs=proj, draw_labels=False, linewidth=0.5, linestyle='--', color='gray', alpha=0.5)
    gl.xlocator = FixedLocator(xticks) # 経度線を描く値
    gl.ylocator = FixedLocator(yticks) # 緯度線を描く値

    #シェイドを描く
    lev=np.arange(1, 10, 1) #範囲の設定
    shade_plot = ax.contourf(xlong,xlat,fod[time,:,:],transform=ccrs.PlateCarree(),levels=lev, cmap=cm.jet, extend='max')

    # cmap を弄って最小値の時のα値を0(透明)にする
    """
    lev=np.arange(-40, 41, 2) #範囲の設定
    cmap = cm.bwr
    cmap_data = cmap(np.arange(cmap.N))
    cmap_data[127-5:127+5,0] = 0.999 # 0 のときのα値を0(透明)にする
    cmap_data[127-5:127+5,1] = 0.999 # 0 のときのα値を0(透明)にする
    customized_cmap = ListedColormap(cmap_data)
    shade_plot = ax.contourf(xlong,xlat,fod[time,:,:],transform=ccrs.PlateCarree(),levels=lev, cmap=customized_cmap, extend='both',alpha=0.8)
    """


    #カラーバーを描く
    # cax = ax.inset_axes([1.04, 0.3, 0.03, 0.4], transform=ax.transAxes)#cbarの座標
    cax = ax.inset_axes([1.04, 0.22, 0.03, 0.8], transform=ax.transAxes)#cbarの座標
    cbar = plt.colorbar(shade_plot, cax = cax, orientation='vertical',extendfrac = 'auto', ticks=lev )
    cbar.ax.tick_params(labelsize=fontsize)

    fig.suptitle(tmp, fontsize=20)
    fig.savefig(output_dir+"/"+str(img_num).zfill(3)+".png", dpi=100, bbox_inches='tight', pad_inches=0.1)
    print('save_'+str(img_num).zfill(3)+".png"),sys.stdout.flush()
    plt.gca().clear()
    plt.close()

os.system("convert -delay 50 -loop 0 "+os.path.join(output_dir,"*.png")+" fod.gif")

print(datetime.now().replace(microsecond = 0),"END",__file__)