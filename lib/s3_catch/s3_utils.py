# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 13:06:16 2020

@author: Cecile M. M. Kittel
Copyright: Cecile M. M. Kittel and Technical University of Denmark
"""

import glob
import os
import numpy as np
import math
from math import sin, cos, sqrt, atan2, radians
from datetime import datetime, timedelta

import netCDF4
from osgeo import gdal
import geopandas as gpd



def year_fraction(dt):
    # Calculate date in decimal year
    start = datetime(dt.year, 1, 1).toordinal()
    year_length = datetime(dt.year+1, 1, 1).toordinal() - start
    return dt.year + float(dt.toordinal() - start) / year_length

def dist(p1, p2):
    # Calculate degree distance between two points in lat-lon coordinates
    (x1, y1), (x2, y2) = p1, p2
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)


def dist_meters(p1, p2):
    # Calculate the approximate distance in m between two points using lat-lon coordinates
    (x1, y1) = p1
    (x2, y2) = p2
    # approximate radius of earth in km
    R = 6373.0
    
    lat1 = radians(y1)
    lon1 = radians(x1)
    lat2 = radians(y2)
    lon2 = radians(x2)
    
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    
    distance = R * c * 1000.
    
    return distance


def create_vs(vs, p):
    vs[p] = {}
    vs[p]['lat'] = np.array([])
    vs[p]['lon'] = np.array([])
    vs[p]['height'] = np.array([])
    vs[p]['epoch'] = np.array([])
    vs[p]['TAI'] = np.array([])
    vs[p]['wm'] = np.array([])
    vs[p]['DEM_S3'] = np.array([])
    vs[p]['orbit'] = np.array([])
    vs[p]['pass'] = np.array([])
    vs[p]['sat_path'] = np.array([])
    vs[p]['sigma0'] = np.array([])
    vs[p]['misfit'] = np.array([])
    vs[p]['RIP_vector'] = []
    vs[p]['WF_vector'] = []
 
    return vs

def create_outliers(outliers, p):
    outliers[p] = {}
    outliers[p]['Outliers_DEM'] = {}
    outliers[p]['Outliers_WM'] = {}
    outliers[p]['Outliers_WM_DEM'] = {}

    return outliers
 

def prep_vsd(vs, p):
    vsd = {}
    vsd['TAI'] = np.array(sorted(list(set(vs[p]['TAI']))))
    vsd['lat'] = np.ones(len(vsd['TAI']))*np.nan
    vsd['lon'] = np.ones(len(vsd['TAI']))*np.nan
    vsd['height'] = np.ones(len(vsd['TAI']))*np.nan
    vsd['DEM'] = np.ones(len(vsd['TAI']))*np.nan
    vsd['RIP'] = np.ones(len(vsd['TAI']))
    vsd['height_RIP'] = np.ones(len(vsd['TAI']))*np.nan
    vsd['orbit'] = np.ones(len(vsd['TAI']))*np.nan
    vsd['height_std'] = np.ones(len(vsd['TAI']))*np.nan

    if len(np.where(vs[p]['sat_path'] == 'descending')[0]) != len(vs[p]['sat_path']):
        vsd['TAI_asc'] = np.array(sorted(list(set(vs[p]['TAI'][np.where(vs[p]['sat_path'] == 'ascending')[0]]))))
        vsd['lat_asc'] = np.ones(len(vsd['TAI_asc']))*np.nan
        vsd['lon_asc'] = np.ones(len(vsd['TAI_asc']))*np.nan
        vsd['height_asc'] = np.ones(len(vsd['TAI_asc']))*np.nan
        vsd['DEM_asc'] = np.ones(len(vsd['TAI_asc']))*np.nan
        vsd['RIP_asc'] = np.ones(len(vsd['TAI_asc']))
        vsd['height_RIP_asc'] = np.ones(len(vsd['TAI_asc']))*np.nan
        vsd['orbit_asc'] = np.ones(len(vsd['TAI_asc']))*np.nan

        vsd['TAI_desc'] = np.array(sorted(list(set(vs[p]['TAI'][np.where(vs[p]['sat_path'] == 'descending')[0]]))))
        vsd['lat_desc'] = np.ones(len(vsd['TAI_desc']))*np.nan
        vsd['lon_desc'] = np.ones(len(vsd['TAI_desc']))*np.nan
        vsd['height_desc'] = np.ones(len(vsd['TAI_desc']))*np.nan
        vsd['DEM_desc'] = np.ones(len(vsd['TAI_desc']))*np.nan
        vsd['RIP_desc'] = np.ones(len(vsd['TAI_desc']))
        vsd['height_RIP_desc'] = np.ones(len(vsd['TAI_desc']))*np.nan
        vsd['orbit_desc'] = np.ones(len(vsd['TAI_desc']))*np.nan
            
    return vsd

def outlier_filter(p, nc, wm_folder, selected, dem_filter, lat, lon, sigma0, sigma_thresh, 
                   source, rip, rip_thresh):
    
    # Get water mask coordinates
    wm_coords = [(float(fn.split('_')[1]), float(fn.split('_')[2].split('.tif')[0])) for fn in os.listdir(wm_folder)]
    (x, y) = sorted([(dist(p1, p), p1) for p1 in wm_coords])[0][1]
    src_filename = glob.glob(os.path.join(wm_folder, '*' + '_' + str(np.round(x, 1)) + '_' + str(np.round(y, 1)) + '.tif'))[0]
    src_ds=gdal.Open(src_filename) 
    rb=src_ds.GetRasterBand(1)
    gt=src_ds.GetGeoTransform() #geotransform
    
    mask = np.ones(len(selected[dem_filter]))*np.nan
    for i in range(0,len(selected[dem_filter])):
        mx = nc.variables[lon][selected][dem_filter][i]
        my = nc.variables[lat][selected][dem_filter][i]
        
        px = int((mx - gt[0]) / gt[1]) #x pixel
        py = int((my - gt[3]) / gt[5]) #y pixel
    
        if px>0 and py>0 and px<src_ds.RasterXSize and py<src_ds.RasterYSize:
            intval = rb.ReadAsArray(px,py,1,1)[0]
            if intval >= 10: ## Not within river mask AND too wide
                mask[i] = intval            
            
            elif mask[i-1] >= 10:
                mask[i] = 1
    
            elif source == 'GPOD':
                if (((np.nanmax(nc.variables[rip][selected][dem_filter][i-1].filled()) >= rip_thresh) and 
                  (np.nanmax(nc.variables[rip][selected][dem_filter][i].filled()) >= rip_thresh)) or 
                  ((nc.variables[sigma0][selected][dem_filter][i-1] >= sigma_thresh) and 
                   (nc.variables[sigma0][selected][dem_filter][i] >= sigma_thresh))):
                    mask[i] = 1
            elif source == 'SciHub':
                if (((nc.variables[sigma0][selected][dem_filter][i-1] >= sigma_thresh) and 
                   (nc.variables[sigma0][selected][dem_filter][i] >= sigma_thresh))):
                    mask[i] = 1
                    
    return mask



def wm_outliers(p, nc, wm_folder, selected, dem_filter, lat, lon, sigma0, sigma_thresh, source, rip, rip_thresh):
    # Get water mask coordinates
    wm_coords = [(float(fn.split('_')[1]), float(fn.split('_')[2].split('.tif')[0])) for fn in os.listdir(wm_folder)]
    (x, y) = sorted([(dist(p1, p), p1) for p1 in wm_coords])[0][1]
        
    src_filename = glob.glob(os.path.join(wm_folder, '*' + '_' + str(np.round(x, 2)) + '_' + str(np.round(y, 2)) + '.tif'))[0]
    src_ds=gdal.Open(src_filename) 
    rb=src_ds.GetRasterBand(1)
    gt=src_ds.GetGeoTransform() #geotransform
    
    mask_wm = np.ones(len(selected))*np.nan
    for i in range(0,len(selected)):
        mx = nc.variables[lon][selected][i]
        my = nc.variables[lat][selected][i]
        
        px = int((mx - gt[0]) / gt[1]) #x pixel
        py = int((my - gt[3]) / gt[5]) #y pixel
    
        if px>0 and py>0 and px<src_ds.RasterXSize and py<src_ds.RasterYSize:
            intval = rb.ReadAsArray(px,py,1,1)[0]
            if intval >= 10: ## Not within river mask AND too wide
                mask_wm[i] = intval            
            
            elif mask_wm[i-1] >= 10:
                mask_wm[i] = 1

            elif source == 'GPOD':
                if (((np.nanmax(nc.variables[rip][selected][i-1].filled()) >= rip_thresh) and 
                  (np.nanmax(nc.variables[rip][selected][i].filled()) >= rip_thresh)) or 
                  ((nc.variables[sigma0][selected][i-1] >= sigma_thresh) and 
                   (nc.variables[sigma0][selected][i] >= sigma_thresh))):
                    mask_wm[i] = 1
                    
            elif source == 'SciHub':
                if (((nc.variables[sigma0][selected][i-1] >= sigma_thresh) and 
                   (nc.variables[sigma0][selected][i] >= sigma_thresh))):
                    mask_wm[i] = 1

    return mask_wm

def count_peaks(wf, threshold_per=0.25):
    # Count peaks exceeding 25% of waveform peak
    count = np.zeros(len(wf))
    for i in range(0,len(wf)):
        threshold = threshold_per * np.nanmax(wf[i])
        if wf[i,0] >= threshold:
            count[i] += 1
        for j in range(1,len(wf[i])):
            if wf[i,j] >= threshold and wf[i,j-1] < threshold:
                count[i] += 1
    return count

def wf_metrics(vs):
    # Calculate waveform metrics for all observations at vs
    for p in vs.keys():
        if 'RIP' in vs[p].keys():
            vs[p]['SP'] = np.array(np.nanmax(vs[p]['RIP'], axis = 1) /
                                  np.nansum(vs[p]['RIP'], axis = 1))
        vs[p]['MP'] = np.array(np.nanmax(vs[p]['WF'], axis = 1))
        vs[p]['PP'] = np.array(np.nanmax(vs[p]['WF'], axis = 1) /
                                np.nansum(vs[p]['WF'], axis = 1))
        vs[p]['NP'] = count_peaks(vs[p]['WF'])

    return vs

def wf_metrics_d(vsd):
    # Calculate waveform metrics for valid observations used to create WSE time series
    for p in vsd.keys():
        if 'RIP' in vsd[p].keys():
            vsd[p]['SP'] = np.array(np.nanmax(vsd[p]['RIP_used'], axis = 1) /
                                  np.nansum(vsd[p]['RIP_used'], axis = 1))
        vsd[p]['MP'] = np.array(np.nanmax(vsd[p]['WF_used'], axis = 1))
        vsd[p]['PP'] = np.array(np.nanmax(vsd[p]['WF_used'], axis = 1) /
                                np.nansum(vsd[p]['WF_used'], axis = 1))
        nps = count_peaks(vsd[p]['WF_used'])
        vsd[p]['NP'] = np.ones(len(vsd[p]['TAI']))*np.nan

        # Get median number of peaks at vsd
        for i, d in enumerate(vsd[p]['TAI']):
            if np.isfinite(vsd[p]['height'][np.where(vsd[p]['TAI'] == d)]):
                vsd[p]['NP'][i] = np.percentile(nps[np.where(vsd[p]['TAI_used'] == d)], 50)

    return vsd


def l1b_filter(vs, selection):
    # Sort stations depending on l1b statistics using NP criteria
    countSP = []
    countNP  = []
    countPP = []
    countMP = []

    for p in selection:
        if len(vs[p]['NP'][vs[p]['NP'] == 1]) >= 0.9*len(vs[p]['NP'][np.isfinite(vs[p]['NP'])]):
            countNP.append(p)
        if 'SP' in vs[p].keys():
            if np.percentile(vs[p]['SP'], 50) > 0.2:
                countSP.append(p)
        if 'RIP' in vs[p].keys():
            if np.percentile(vs[p]['MP'], 50) > 1e-15:
                countMP.append(p)
        if np.percentile(vs[p]['PP'], 50) > 0.1:
            countPP.append(p)
    return countSP, countNP, countPP, countMP


