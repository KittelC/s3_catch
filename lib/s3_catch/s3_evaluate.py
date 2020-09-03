# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 13:30:14 2020

@author: Cecile M. M. Kittel
Copyright: Cecile M. M. Kittel and Technical University of Denmark
"""

import os
import numpy as np
import math
from math import sin, cos, sqrt, atan2, radians
from datetime import datetime, timedelta

import netCDF4
from osgeo import gdal
import geopandas as gpd

from s3_utils import dist_meters, create_vs, count_peaks, create_outliers, wm_outliers, outlier_filter, prep_vsd, year_fraction, l1b_filter


def sort_l2(vs, outliers, vsd, oltc=False, oltc_date=datetime(2019,3,1).date(),
            vs3=None, outliers3=None, vsd3=None):

    nodata = []
    nowater = []
    somedata = []
    mostdata = []
    postoltc = []
    windext = []
    
    for p in sorted(outliers.keys()):
        
        (x,y) = p
        if not p in vs.keys():
            nodata.append(p) # No data at VS
        elif len(vs[p]['lat']) == 0:
            nowater.append(p) # No water data at VS
        elif len(vsd[p]['lat'][~np.isnan(vsd[p]['height'])]) <= 0.8 * len(outliers[p]['Outliers_WM_DEM']):
            somedata.append(p) # Less than 80% of data at VS
        elif len(vsd[p]['lat'][~np.isnan(vsd[p]['height'])]) > 0.8 * len(outliers[p]['Outliers_WM_DEM']):
            mostdata.append(p) # At least 80% data at VS
    
        # For the remaining data, check window extension:
        if oltc:
            if isinstance(vs3, dict) and p in vsd3.keys():
                valid_obs3 = vsd3[p]['height'][~np.isnan(vsd3[p]['height'])]
                valid_dates3 = vsd3[p]['TAI'][~np.isnan(vsd3[p]['height'])]
                all_dates3 = np.array(list(outliers3[p]['Outliers_WM_DEM'].keys()))
        

            if p in vsd.keys() and p not in mostdata:
                valid_obs = vsd[p]['height'][~np.isnan(vsd[p]['height'])]
                valid_dates = vsd[p]['TAI'][~np.isnan(vsd[p]['height'])]
                
                all_dates = np.array(list(outliers[p]['Outliers_WM_DEM'].keys()))
            
                
                if (len(valid_obs[np.where(valid_dates >=  oltc_date)]) >=                   #Criteria 1: after OLTC, 80% data
                    0.8 * len(np.where(all_dates >= oltc_date)[0]) and
                    len(valid_obs[np.where(valid_dates < oltc_date)]) <                      #Criteria 2: before OLTC < 80% data
                                0.8 * len(np.where(all_dates < oltc_date)[0])):
                    
                    if isinstance(vs3, dict) and p in vsd3.keys():
                        if (len(valid_obs3[np.where(valid_dates3 >=  oltc_date)]) >=                 #Criteria 3: 3 bin window shows same improvement
                            0.8 * len(np.where(all_dates3 >= oltc_date)[0]) and
                            len(valid_obs3[np.where(valid_dates3 < oltc_date)]) <
                                        0.8 * len(np.where(all_dates3 < oltc_date)[0])):
                            postoltc.append(p)
                        
                    elif not isinstance(vs3, dict):
                        postoltc.append(p)
                
                elif isinstance(vs3, dict) and (len(valid_obs3) > len(valid_obs) and
                            len(valid_obs3) > 0.8 * len(all_dates3) and
                            len(valid_obs) <= 0.8 * len(all_dates)):
                            windext.append(p)

        
        if isinstance(vs3, dict):
            if (p in vsd3.keys()) and p not in windext and p not in mostdata:
                if len(vsd3[p]['lat'][~np.isnan(vsd3[p]['height'])]) > 0.8 * len(outliers3[p]['Outliers_WM_DEM']):
                    windext.append(p)
                
    tot = len(mostdata) + len(windext) + len(postoltc)
    print('Total stations: ', tot)

    rej = (1 - (len(mostdata) + len(windext) + len(postoltc))/len(outliers.keys()))*100
    print('Rejection rate, S3A: ', rej)
    
    return mostdata, windext, postoltc


def sort_l1b(vsd, mostdata, windext, postoltc, ext=None, oltc=None):
    
    spa, npa, ppa, mpa = l1b_filter(vsd, mostdata)
    s3_valid = npa
    if ext:
        spa3, npa3, ppa3, mpa3 = l1b_filter(ext, windext)
        s3_valid_3 = npa3
    else:
        s3_valid_3 = []
    if oltc:
        spao, npao, ppao, mpao = l1b_filter(vsd, postoltc)
        s3_valid_o = npao
    else:
        s3_valid_o = []

    return s3_valid, s3_valid_3, s3_valid_o




