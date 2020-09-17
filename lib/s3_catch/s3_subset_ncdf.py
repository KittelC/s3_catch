# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 15:50:59 2019

@author: ceki
"""

import os
import numpy as np

import netCDF4
import geopandas as gpd


def subset_scihub_netcdf(download_dir, dest_dir, extent, file_id=r'enhanced_measurement.nc'):
    
    [ulx, uly, lrx, lry] = extent
    if ulx > lrx or lry > uly:
        raise ValueError('Study area extent incorrect - please check coordinates.')

    for folder in os.listdir(download_dir):
        if folder.endswith('.SEN3'):
            file = os.path.join(os.path.join(download_dir, folder), file_id)
            if os.path.isfile(file):
                nc_cropped = os.path.join(dest_dir, 'sub_' + os.path.split(folder)[-1].split('.SEN3')[0]+'.nc')
        
                if not os.path.isfile(nc_cropped):
                    nc = netCDF4.Dataset(file)
                    
                    # SciHub contains 1Hz and 20Hz variables
                    lat20 = nc['lat_20_ku'][:]
                    lon20 = nc['lon_20_ku'][:]
                    
                    lat01 = nc['lat_01'][:]
                    lon01 = nc['lon_01'][:]
            
                    
                    VS = {}
                    min_index20 = len(lat20)+1
                    max_index20 = 0
                    
                    min_index01 = len(lat01) + 1
                    max_index01 = 0 
                    
                    
                    selected = np.where((lon20 <= lrx) & (lon20 >= ulx) & (lat20 >= lry) & (lat20 <= uly))[0]
    
                    if len(selected) > 0:
                        if selected[0] < min_index20:                
                            min_index20 = selected[0]
                        if selected[-1] > max_index20:
                            max_index20 = selected[-1]
                    selected01 = np.where((lon01 <= lrx) & (lon01 >= ulx) & (lat01 >= lry) & (lat01 <= uly))[0]
              
                    if len(selected01) > 0:
                        if selected01[0] < min_index01:                
                            min_index01 = selected01[0]
                        if selected01[-1] > max_index01:
                            max_index01 = selected01[-1]
                        
                    if max_index20 > 0:
                        print(file)
                        with netCDF4.Dataset(file) as src, netCDF4.Dataset(nc_cropped, "w") as dst:
                            # copy global attributes all at once via dictionary
                            dst.setncatts(src.__dict__)
                            # copy dimensions
                            for name, dimension in src.dimensions.items():
                                if '20_' in name:
                                    dst.createDimension(
                                        name, (len(np.arange(min_index20,max_index20+1,1)) if not dimension.isunlimited() else None))
                                elif '01' in name:
                                    dst.createDimension(
                                        name, (len(np.arange(min_index01,max_index01+1,1)) if not dimension.isunlimited() else None))
                                else:
                                    dst.createDimension(
                                        name, (len(dimension) if not dimension.isunlimited() else None))
    
                            # copy all file data except for the excluded
                            for name, variable in src.variables.items():
                                if '20_ku' in name:
                                    x = dst.createVariable(name, variable.datatype, variable.dimensions)
                                    x = dst[name].setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
                                    dst[name][:] = src[name][min_index20:max_index20+1]
                                elif '01' in name:
                                    x = dst.createVariable(name, variable.datatype, variable.dimensions)
                                    x = dst[name].setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
                                    dst[name][:] = src[name][min_index01:max_index01+1]
        # else:
        #     print('No Sentinel-3 files in directory.')
        
if __name__ == '__main__':
    
    download_dir = r'..\..\..\test\SciHub\Full'
    dest_dir = os.path.join(download_dir, r'..\Subsets2')
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)
        
    extent = [19, -9, 36, -20] #upper left x, upper left y, lower right x, lower right y
    
    subset_scihub_netcdf(download_dir, dest_dir, extent, file_id=r'enhanced_measurement.nc')
