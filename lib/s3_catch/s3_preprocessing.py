# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 13:12:58 2020

@author: Cecile M. M. Kittel
Copyright: Cecile M. M. Kittel and Technical University of Denmark
"""

import os
import numpy as np
from datetime import datetime, timedelta

import netCDF4
from osgeo import gdal
import geopandas as gpd

from s3_utils import create_vs, create_outliers, wm_outliers, outlier_filter, prep_vsd, year_fraction


def read_s3_nc(s3_folder, vs_coords, wm_folder, source, 
               dem_file=r'..\Zambezi\Model\RR_Model\MERIT\merit_egm2008.tif',
               sigma_thresh=30, dem_thresh=30, vs_buffer=0.015, rip_thresh=1e-13,
               stack=False, stack_folder=r'E:\L2_WDIR\L1b_RESDIR\STACK'):

    """
    
    Parameters
    ----------
    s3_folder : STRING
        Foldername and path for Level 2 netcdf files
    vs_coords : GEOPANDAS OBJECT
        Geopandas with virtual station coordinates
    wm_folder : STRING
        Foldername and path for water mask files at each vs
    source : STRING
        "GPOD" or "SciHub" - source of Sentinel-3 netcdf file
    dem_file : STRING,
        Path to DEM file. The default is r'..\Zambezi\Model\srtm_3arc_v4_clipped.tif'.

    sigma_thresh : Float, optional
        Backscatter coefficient threshold. The default is 30.
    dem_thresh : Float, optional
        DEM elevation threshold [m].
    vs_buffer : Float, optional
        Degree threshold in distance to vs. The default is 0.015.
    rip_thresh : Float, optional
        Maximum RIP threshold. The default is 1e-13.
    stack : Boolean, optional - only implemented for GPOD
        

    Returns
    -------
    vs : Dictionary
        Keys = vs coordinates, values = latitude, longitude, elevation, time...
        of all data points for each vs
    outliers: Dictionary
        Keys = vs Coordinates, values = DEM and water mask values for each outlier

    """
    
    if source == 'SciHub':
        lat = 'lat_20_ku'
        lon = 'lon_20_ku'
        elev = 'elevation_ocog_20_ku'
        wf = 'waveform_20_ku'
        sig0 = 'sig0_ocog_20_ku'
        tai = 'time_20_ku'
        rip=None
    
    elif source == 'GPOD':
        lat = 'latitude_20Hz'
        lon = 'longitude_20Hz'
        altitude = 'altitude_20Hz'
        range_unc = 'Range_Unc_20Hz'
        geo_corr = 'GEO_Corr_Land_20Hz'
        geoid = 'EGM_2008_20Hz'
        demid = 'Land_Alt_20Hz'
        rip = 'Substack_RIP_Data'
        wf = 'SAR_Echo_Data'
        sig0 = 'Sigma0_20Hz'
        tai = 'TAI_Time_20Hz'
        misfit = 'Misfit_20Hz'
    
    else: 
        raise ValueError('Source should be "GPOD" or "SciHub"') 

    vs = {}
    outliers = {}
    filetracker = {}
    count = 0
    start = datetime.now()
    for f in os.listdir(s3_folder):
        s3_file = os.path.join(s3_folder, f)
        if s3_file.endswith('.nc'):
            nc = netCDF4.Dataset(s3_file)
            for p in zip(vs_coords['xcoord'], vs_coords['ycoord']):
                (x, y) = p
                # Get all points at a distance of max 2.3 (0.15 degrees) from vs
                selected = np.where((nc.variables[lon][:] <= x + vs_buffer) &
                                    (nc.variables[lon][:] >= x - vs_buffer) &
                                    (nc.variables[lat][:] >= y - vs_buffer) &
                                    (nc.variables[lat][:] <= y + vs_buffer))[0]
                if len(selected) > 0:
                    if source == 'GPOD':
                        height = (nc.variables[altitude][selected].filled() -
                                  nc.variables[range_unc][selected].filled() -
                                  nc.variables[geo_corr][selected].filled() -
                                  nc.variables[geoid][selected].filled())
                        dem = (nc.variables[demid][selected].filled() -
                           nc.variables[geoid][selected].filled())

                    elif source == 'SciHub':
                        # Get geoid elevation at 20Hz resolution from 1Hz dataset 
                        lat_01 = nc['lat_01'][:].filled()
                        geoid_01 = nc['geoid_01'][:].filled()
                        geoid = np.interp(nc.variables[lat][:].filled(), lat_01, geoid_01)
                    
                        # Get retracked WSE 
                        height = (nc[elev][:].filled() - geoid)[selected]
                        src_ds=gdal.Open(dem_file) 
                        rb = src_ds.GetRasterBand(1)
                        gt=src_ds.GetGeoTransform() #geotransform
                        
                        n = len(selected)
                        dem = np.ones(n)*np.nan
                        for i in range(0,n):
                            mx = nc.variables[lon][selected][i]
                            my = nc.variables[lat][selected][i]
                            
                            px = int((mx - gt[0]) / gt[1]) #x pixel
                            py = int((my - gt[3]) / gt[5]) #y pixel
                        
                            if px>0 and py>0 and px<src_ds.RasterXSize and py<src_ds.RasterYSize:
                                intval = rb.ReadAsArray(px,py,1,1)[0]
                                dem[i] = intval

                    dem_filter = np.where((height <= dem + dem_thresh) &
                                          (height >= dem - dem_thresh))[0]
                    
                    if len(dem_filter) > 0:
                        mask = outlier_filter(p, nc, wm_folder, selected, dem_filter, lat, lon, sig0, sigma_thresh, source, rip, rip_thresh)
                                
                    mask_wm = wm_outliers(p, nc, wm_folder, selected, dem_filter, lat, lon, sig0, sigma_thresh, source, rip, rip_thresh)
                    
                    dem_filter_WM = np.where((height[~np.isnan(mask_wm)] <= dem[~np.isnan(mask_wm)] + dem_thresh) &
                                          (height[~np.isnan(mask_wm)] >= dem[~np.isnan(mask_wm)] - dem_thresh))[0]

                    if not p in outliers.keys():
                        outliers = create_outliers(outliers, p)
                        
                    sensing_date = [(datetime(2000,1,1) + timedelta(seconds = c2_time)).date()
                                    for c2_time in nc.variables[tai][selected].filled()][0]
                    outliers[p]['Outliers_WM_DEM'][sensing_date] = len(dem[~np.isnan(mask_wm)]) - len(dem[~np.isnan(mask_wm)][dem_filter_WM])

                    if len(dem_filter) == 0:
                        outliers[p]['Outliers_DEM'][sensing_date] = len(dem) 
                        

                    elif len(dem_filter) > 0:
                        if not p in vs.keys():
                            vs = create_vs(vs, p)
                        if stack and os.path.exists(stack_folder):
                            root = f.split('RES')[1].split('.nc')[0]
                            if not p in filetracker.keys():
                                filetracker[p] = {}
                            filetracker[p][root] = np.array([])
    
                        vs[p]['wm'] = np.concatenate([vs[p]['wm'], mask[~np.isnan(mask)]])
                        outliers[p]['Outliers_DEM'][sensing_date] = len(dem) - len(dem[dem_filter])
                        outliers[p]['Outliers_WM'][sensing_date] = len(dem[dem_filter]) - len(dem[dem_filter][~np.isnan(mask)])

                        vs[p]['lat'] = np.concatenate([vs[p]['lat'], nc.variables[lat][selected][dem_filter][~np.isnan(mask)].filled()])
                        vs[p]['lon'] = np.concatenate([vs[p]['lon'], nc.variables[lon][selected][dem_filter][~np.isnan(mask)].filled()])
#                            print(p, len(vs[p]['lat']))
                        vs[p]['height'] = np.concatenate([vs[p]['height'], height[dem_filter][~np.isnan(mask)]])
                        
                        
                        vs[p]['TAI'] = np.concatenate([vs[p]['TAI'], [(datetime(2000,1,1) +
                                                         timedelta(seconds = c2_time)).date() for c2_time in nc.variables[tai][selected][dem_filter][~np.isnan(mask)].filled()]])
                        vs[p]['DEM_S3'] = np.concatenate([vs[p]['DEM_S3'], dem[dem_filter][~np.isnan(mask)]])
                        vs[p]['sigma0'] = np.concatenate([vs[p]['sigma0'], nc.variables[sig0][selected][dem_filter][~np.isnan(mask)].filled()])
                        

                        vs[p]['WF_vector'].append(nc.variables[wf][selected][dem_filter][~np.isnan(mask)].filled())

                        if source == 'GPOD':
                            vs[p]['misfit'] = np.concatenate([vs[p]['misfit'], nc.variables[misfit][selected][dem_filter][~np.isnan(mask)].filled()])
                            
                            vs[p]['sat_path'] = np.concatenate([vs[p]['sat_path'], np.repeat(nc.getncattr('PASS_DIRECTION_START'),
                                 len(dem[dem_filter][~np.isnan(mask)]))])

                            vs[p]['RIP_vector'].append(nc.variables[rip][selected][dem_filter][~np.isnan(mask)].filled())
                            vs[p]['orbit'] = np.concatenate([vs[p]['orbit'], np.repeat(nc.getncattr('RELATIVE_ORBIT_NUMBER_START'),
                                 len(dem[dem_filter][~np.isnan(mask)]))])
                            vs[p]['pass'] = np.concatenate([vs[p]['pass'], np.repeat(nc.getncattr('RELATIVE_PASS_NUMBER_START'),
                                 len(dem[dem_filter][~np.isnan(mask)]))])
                            
                            vs[p]['epoch'] = np.concatenate([vs[p]['epoch'],
                                     (nc['Epoch_20Hz'][selected][dem_filter][~np.isnan(mask)].filled()*
                                              (320e6*nc['Echo_Oversampling_Factor'][:]) +
                                              nc['Epoch_Reference_Gate'][:]-1)])    
                            if stack:
                                filetracker[p][root] = np.concatenate([filetracker[p][root],
                                                                       nc.variables['Meas_Index_20Hz'][selected][dem_filter][~np.isnan(mask)].filled()])
                        if source == 'SciHub':
                            vs[p]['sat_path'] = np.concatenate([vs[p]['sat_path'], np.repeat('descending' if nc.getncattr('first_meas_lat')-nc.getncattr('last_meas_lat')<0 else 'ascending',
                                 len(dem[dem_filter][~np.isnan(mask)]))])
                            vs[p]['pass'] = np.concatenate([vs[p]['pass'], np.repeat(nc.getncattr('pass_number'),
                                 len(dem[dem_filter][~np.isnan(mask)]))])
                            vs[p]['epoch'] = np.concatenate([vs[p]['epoch'],
                                 (-nc['tracker_range_20_ku'][:][selected][dem_filter][~np.isnan(mask)].filled()
                                 +nc['range_ocog_20_ku'][:][selected][dem_filter][~np.isnan(mask)].filled()-0.55)*128/60+43])
        

            count += 1
            print(count)
            nc.close()
            
    for p in vs.keys():
        vs[p]['WF'] =  np.zeros([sum([len(vs[p]['WF_vector'][es]) for es in range(0,len(vs[p]['WF_vector']))]), 
                       max([np.shape(vs[p]['WF_vector'][es])[1] for es in range(0,len(vs[p]['WF_vector']))])])
        
        if source == 'GPOD':
            vs[p]['RIP'] =  np.zeros([sum([len(vs[p]['RIP_vector'][es]) for es in range(0,len(vs[p]['RIP_vector']))]), 
                   max([np.shape(vs[p]['RIP_vector'][es])[1] for es in range(0,len(vs[p]['RIP_vector']))])])
        count = 0    
        for t in range(0,len(vs[p]['WF_vector'])):
        #    for pt in range(0, len(wm[test])):
        #        if wm[test][pt] > 0:
        #            plt.plot(np.array(RIP[test][pt]).T)    
            vs[p]['WF'][count:count + np.shape(vs[p]['WF_vector'][t])[0],:np.shape(vs[p]['WF_vector'][t])[1]] = vs[p]['WF_vector'][t]
            if source == 'GPOD':
                vs[p]['RIP'][count:count + np.shape(vs[p]['RIP_vector'][t])[0],:np.shape(vs[p]['RIP_vector'][t])[1]] = vs[p]['RIP_vector'][t]

            count += len(vs[p]['WF_vector'][t])
        del vs[p]['WF_vector']
        if source == 'GPOD':
            del vs[p]['RIP_vector']
    if stack and source == 'GPOD':        
        for p in vs.keys():
            print(filetracker[p])
            filetracker[p]['SARStacks'] = []
            for root in filetracker[p].keys():
                for f in os.listdir(stack_folder):
                    stack_file = os.path.join(stack_folder, f)
                    if root in f:
                        starts = int(f.split('.nc')[0].split('_')[-3])
                        ends = int(f.split('.nc')[0].split('_')[-1])
                        if any(starts <= o <= ends for o in filetracker[p][root]):
                            nc = netCDF4.Dataset(stack_file)
                            ind = [int(o-starts) for o in filetracker[p][root][np.where((starts <= filetracker[p][root]) & (filetracker[p][root] <= ends))[0]]]
                            filetracker[p]['SARStacks'].append(nc.variables['SAR_Stack_Data'][ind].filled())
                            nc.close()
        for p in vs.keys():
            print([np.shape(filetracker[p]['SARStacks'][es])[1] for es in range(0,len(filetracker[p]['SARStacks']))])

            vs[p]['SARStacks'] =  np.zeros([sum([len(filetracker[p]['SARStacks'][es]) for es in range(0,len(filetracker[p]['SARStacks']))]), 
                            max([np.shape(filetracker[p]['SARStacks'][es])[1] for es in range(0,len(filetracker[p]['SARStacks']))]),
                            max([np.shape(filetracker[p]['SARStacks'][es])[2] for es in range(0,len(filetracker[p]['SARStacks']))])])
            count = 0
            for t in range(0,len(filetracker[p]['SARStacks'])):
                vs[p]['SARStacks'][count:count + np.shape(filetracker[p]['SARStacks'][t])[0],
                        :np.shape(filetracker[p]['SARStacks'][t])[1],
                        :np.shape(filetracker[p]['SARStacks'][t])[2]] = filetracker[p]['SARStacks'][t]
        
                count += len(filetracker[p]['SARStacks'][t])

    end = datetime.now()
    timer = end-start
    print(timer)
    
    return vs, outliers


def create_vs_ts(vs, subset_vs=False, sigma_thresh=30, source='GPOD'):
    """
    

    Parameters
    ----------
    vs : Dictionary
        Virtual stations to be processed for WSE time series.
    subset_vs : List, optional
        Coordinates of selected virtual stations. The default is False.
    sigma_thresh : string, optional
        DESCRIPTION. The default is 30.

    Returns
    -------
    vs_d

    """
    vs_d = {}
    if not subset_vs:
        subset_vs = list(vs.keys())
    ## Track average
    for p in vs.keys():
        if len(vs[p]['height'] > 0):
            if np.any(np.isclose(subset_vs, p, rtol=1e-10)):
                vs_d[p] = prep_vsd(vs, p)

                # Can be edited to filter out according to time series statistics
                filth = np.arange(0,len(vs[p]['height']))

                selected = np.array([])
                for ind, uniq_d in enumerate(sorted(list(set(vs[p]['TAI'])))):
                    day_ind = np.where(vs[p]['TAI'][filth] == uniq_d)[0]
                    if len(day_ind) > 0:
                        filt = np.where((vs[p]['sigma0'][filth][day_ind] >= sigma_thresh))[0]
    
                        if len(filt) >= 1:
                            vs_d[p]['lat'][ind] = np.mean(vs[p]['lat'][filth][day_ind[filt]])
                            vs_d[p]['lon'][ind] = np.mean(vs[p]['lon'][filth][day_ind[filt]])
                            vs_d[p]['height'][ind] = np.median(vs[p]['height'][filth][day_ind[filt]])
                            vs_d[p]['height_std'][ind] = np.std(vs[p]['height'][filth][day_ind[filt]])

                            vs_d[p]['DEM'][ind] = np.mean(vs[p]['DEM_S3'][filth][day_ind[filt]])
                            if source == 'GPOD':
                                vs_d[p]['RIP'][ind] = np.mean(np.array([np.nanmax(vs[p]['RIP'][t])
                                                                        for t in filth[day_ind[filt]]]))
                                vs_d[p]['height_RIP'][ind] = (vs[p]['height'][filth][day_ind[filt]]
                                                             [np.argmax(np.array([np.nanmax(vs[p]['RIP'][t])
                                                    for t in filth[day_ind[filt]]]))])
                                vs_d[p]['orbit'][ind] = int(vs[p]['orbit'][np.where(vs[p]['TAI'] == uniq_d)][0])
                            selected = np.append(selected, filth[day_ind[filt]])
    
        
                # RETRIEVE DAYS KNOWN TO A CERTAIN orbit
                if len(np.where(vs[p]['sat_path'] == 'descending')[0]) != len(vs[p]['sat_path']):
                    select = vs[p]['TAI'][np.where(vs[p]['sat_path'] == 'ascending')[0]]
                    for asc_ind, uniq_d in enumerate(sorted(list(set(select)))):
                        ind = np.where(np.array(vs_d[p]['TAI']) == uniq_d)
                        vs_d[p]['lat_asc'][asc_ind] = vs_d[p]['lat'][ind]
                        vs_d[p]['lon_asc'][asc_ind] =  vs_d[p]['lon'][ind]
                        vs_d[p]['height_asc'][asc_ind] =  vs_d[p]['height'][ind]
                        vs_d[p]['DEM_asc'][asc_ind] = vs_d[p]['DEM'][ind]
                        if source == 'GPOD':
                            vs_d[p]['RIP_asc'][asc_ind] = vs_d[p]['RIP'][ind]
                            vs_d[p]['height_RIP_asc'][asc_ind] = vs_d[p]['height_RIP'][ind]
                            vs_d[p]['orbit_asc'][asc_ind] = vs_d[p]['orbit'][ind]
            
                    select = vs[p]['TAI'][np.where(vs[p]['sat_path'] == 'descending')[0]]
                    for desc_ind, uniq_d in enumerate(sorted(list(set(select)))):
                        ind = np.where(np.array(vs_d[p]['TAI']) == uniq_d)
                        vs_d[p]['lat_desc'][desc_ind] = vs_d[p]['lat'][ind]
                        vs_d[p]['lon_desc'][desc_ind] =  vs_d[p]['lon'][ind]
                        vs_d[p]['height_desc'][desc_ind] =  vs_d[p]['height'][ind]
                        vs_d[p]['DEM_desc'][desc_ind] = vs_d[p]['DEM'][ind]
                        if source == 'GPOD':
                            vs_d[p]['RIP_desc'][desc_ind] = vs_d[p]['RIP'][ind]
                            vs_d[p]['height_RIP_desc'][desc_ind] = vs_d[p]['height_RIP'][ind]
                            vs_d[p]['orbit_desc'][desc_ind] = vs_d[p]['orbit'][ind]
            
                for k in vs[p].keys():
                    if len(vs[p][k]) > 0: 
                        print(k)
                        try:
                            vs_d[p][k + '_used'] = vs[p][k][selected.astype(int)]
                        except:
                            vs_d[p][k + '_used',:] = vs[p][k][selected.astype(int)]

    return vs_d


def read_s3_nc_floodplain(s3_folder, wm_file, source='GPOD',
                          tracks=['070','427','184','541','298'], dem_thresh=30,
                          sigma_thresh=30, rip_thresh=1e-13,
                          extent=[22.4,23.5,-16.4,-13.8]):
    """
    Extract observations over floodplain using track IDs to group WSE.

    Parameters
    ----------
    s3_folder : STRING - folder path
        Foldername and path for Level 2 netcdf files
    wm_file : STRING - path
        Water mask over floodplain
    source : STRING
        "GPOD" or "SciHub" - source of Sentinel-3 netcdf file, default is GPOD
    tracks : array of strings
        Unique id for tracks crossing the floodplain of interest
        The default is for Sentinel-3A over Kafue: ['070','427','184','541','298'].
    dem_thresh : Float, optional
        DEM elevation threshold [m].
    sigma_thresh : Float, optional
        Backscatter coefficient threshold. The default is 30.
    rip_thresh : Float, optional
        Maximum RIP threshold. The default is 1e-13.
    extent : array - [xmin, xmax, ymin, ymax]
        DESCRIPTION. The default is Kafue: [26.0, 28.2,-15.9,-15.4].

    Returns
    -------
    vs : Dictionary
        Keys = vs coordinates, values = latitude, longitude, elevation, time...
        of all data points for each track crossing the floodplain
    
    """

    if source == 'SciHub':
        lat = 'lat_20_ku'
        lon = 'lon_20_ku'
        elev = 'elevation_ocog_20_ku'
        wf = 'waveform_20_ku'
        sig0 = 'sig0_ocog_20_ku'
        tai = 'time_20_ku'
        passid = 'pass_number'
    
    elif source == 'GPOD':
        lat = 'latitude_20Hz'
        lon = 'longitude_20Hz'
        altitude = 'altitude_20Hz'
        range_unc = 'Range_Unc_20Hz'
        geo_corr = 'GEO_Corr_Land_20Hz'
        geoid = 'EGM_2008_20Hz'
        demid = 'Land_Alt_20Hz'
        rip = 'Substack_RIP_Data'
        wf = 'SAR_Echo_Data'
        sig0 = 'Sigma0_20Hz'
        tai = 'TAI_Time_20Hz'
        misfit = 'Misfit_20Hz'
        passid = 'RELATIVE_PASS_NUMBER_START'

    vs = {}
    count = 0
    start = datetime.now()
    for f in os.listdir(s3_folder):
        s3_file = os.path.join(s3_folder, f)  

        if s3_file.endswith('.nc'):
            nc = netCDF4.Dataset(s3_file)
            if nc.getncattr(passid) in tracks:
                selected = np.where((nc.variables[lon][:] <= extent[1]) &
                                    (nc.variables[lon][:] >= extent[0]) &
                                    (nc.variables[lat][:] >= extent[2]) &
                                    (nc.variables[lat][:] <= extent[3]))[0]
                print((nc.variables[lon][0]))
                # Get all points at a distance of max 2.3 (0.15 degrees) from VS
                if len(selected) > 0:
                    if source == 'GPOD':
                        height = (nc.variables[altitude][selected].filled() -
                                  nc.variables[range_unc][selected].filled() -
                                  nc.variables[geo_corr][selected].filled() -
                                  nc.variables[geoid][selected].filled())
                    elif source == 'SciHub':
                        # Get geoid elevation at 20Hz resolution from 1Hz dataset 
                        lat_01 = nc['lat_01'][:].filled()
                        geoid_01 = nc['geoid_01'][:].filled()
                        geoid = np.interp(nc.variables[lat][:].filled(), lat_01, geoid_01)
                    
                        # Get retracked WSE 
                        height = (nc[elev][:].filled() - geoid)[selected]

                    dem = (nc.variables[demid][selected].filled() -
                           nc.variables[geoid][selected].filled())
                                            
                    dem_filter = np.where((height <= dem + dem_thresh) &
                                          (height >= dem - dem_thresh))[0]
                
                    if len(dem_filter) > 0:
                        src_filename = wm_file
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
                        
                                if (intval >= 1): ## Not within river mask AND too wide
                                    mask[i] = intval            
                                
                                elif mask[i-1] > 1:
                                    mask[i] = 1
        
                                elif source == 'GPOD':
                                    if (((np.nanmax(nc.variables[rip][selected][dem_filter][i-1].filled()) > rip_thresh) and 
                                      (np.nanmax(nc.variables[rip][selected][dem_filter][i].filled()) > rip_thresh)) or 
                                      ((nc.variables[sig0][selected][dem_filter][i-1] > sigma_thresh) and 
                                       (nc.variables[sig0][selected][dem_filter][i] > sigma_thresh))):
                                        mask[i] = 1
                                elif source == 'SciHub':
                                    if ((nc.variables[sig0][selected][dem_filter][i-1] > sigma_thresh) and 
                                       (nc.variables[sig0][selected][dem_filter][i] > sigma_thresh)):
                                        mask[i] = 1
       
                                if nc.variables[sig0][selected][dem_filter][i] <= sigma_thresh:
                                    mask[i] = np.nan
                        p = nc.getncattr(passid)
                       
                        if not p in vs.keys():
                            vs = create_vs(vs, p)
                                    
                        vs[p]['wm'] = np.concatenate([vs[p]['wm'], mask[~np.isnan(mask)]])
                        vs[p]['lat'] = np.concatenate([vs[p]['lat'], nc.variables[lat][selected][dem_filter][~np.isnan(mask)].filled()])
                        vs[p]['lon'] = np.concatenate([vs[p]['lon'], nc.variables[lon][selected][dem_filter][~np.isnan(mask)].filled()])
    #                            print(p, len(vs[p]['lat']))
                        vs[p]['height'] = np.concatenate([vs[p]['height'], height[dem_filter][~np.isnan(mask)]])
                        
                        
                        vs[p]['TAI'] = np.concatenate([vs[p]['TAI'], [(datetime(2000,1,1) +
                                                         timedelta(seconds = c2_time)).date() for c2_time in nc.variables[tai][selected][dem_filter][~np.isnan(mask)].filled()]])
                        vs[p]['DEM_S3'] = np.concatenate([vs[p]['DEM_S3'], dem[dem_filter][~np.isnan(mask)]])
                        vs[p]['sigma0'] = np.concatenate([vs[p]['sigma0'], nc.variables[sig0][selected][dem_filter][~np.isnan(mask)].filled()])
                        
    
                        vs[p]['WF_vector'].append(nc.variables[wf][selected][dem_filter][~np.isnan(mask)].filled())
    
                        if source == 'GPOD':
                            vs[p]['misfit'] = np.concatenate([vs[p]['misfit'], nc.variables[misfit][selected][dem_filter][~np.isnan(mask)].filled()])
                            
                            vs[p]['sat_path'] = np.concatenate([vs[p]['sat_path'], np.repeat(nc.getncattr('PASS_DIRECTION_START'),
                                 len(dem[dem_filter][~np.isnan(mask)]))])
    
                            vs[p]['RIP_vector'].append(nc.variables[rip][selected][dem_filter][~np.isnan(mask)].filled())
                            vs[p]['orbit'] = np.concatenate([vs[p]['orbit'], np.repeat(nc.getncattr('RELATIVE_ORBIT_NUMBER_START'),
                                 len(dem[dem_filter][~np.isnan(mask)]))])
                            vs[p]['pass'] = np.concatenate([vs[p]['pass'], np.repeat(nc.getncattr('RELATIVE_PASS_NUMBER_START'),
                                 len(dem[dem_filter][~np.isnan(mask)]))])
                            
                            vs[p]['epoch'] = np.concatenate([vs[p]['epoch'],
                                     (nc['Epoch_20Hz'][selected][dem_filter][~np.isnan(mask)].filled()*
                                              (320e6*nc['Echo_Oversampling_Factor'][:]) +
                                              nc['Epoch_Reference_Gate'][:]-1)])    
    
                        if source == 'SciHub':
                            vs[p]['sat_path'] = np.concatenate([vs[p]['sat_path'], np.repeat('descending' if nc.getncattr('first_meas_lat')-nc.getncattr('last_meas_lat')<0 else 'ascending',
                                 len(dem[dem_filter][~np.isnan(mask)]))])
                            vs[p]['pass'] = np.concatenate([vs[p]['pass'], np.repeat(nc.getncattr('pass_number'),
                                 len(dem[dem_filter][~np.isnan(mask)]))])
                            vs[p]['epoch'] = np.concatenate([vs[p]['epoch'],
                                 (-nc['tracker_range_20_ku'][:][selected][dem_filter][~np.isnan(mask)].filled()
                                 +nc['range_ocog_20_ku'][:][selected][dem_filter][~np.isnan(mask)].filled()-0.55)*128/60+43])
                                
    
                count += 1
                print(count)
                nc.close()
            
            
    for p in vs.keys():
        vs[p]['WF'] =  np.zeros([sum([len(vs[p]['WF_vector'][es]) for es in range(0,len(vs[p]['WF_vector']))]), 
                       max([np.shape(vs[p]['WF_vector'][es])[1] for es in range(0,len(vs[p]['WF_vector']))])])
        
        if source == 'GPOD':
            vs[p]['RIP'] =  np.zeros([sum([len(vs[p]['RIP_vector'][es]) for es in range(0,len(vs[p]['RIP_vector']))]), 
                   max([np.shape(vs[p]['RIP_vector'][es])[1] for es in range(0,len(vs[p]['RIP_vector']))])])
        count = 0    
        for t in range(0,len(vs[p]['WF_vector'])):
        #    for pt in range(0, len(wm[test])):
        #        if wm[test][pt] > 0:
        #            plt.plot(np.array(RIP[test][pt]).T)    
            vs[p]['WF'][count:count + np.shape(vs[p]['WF_vector'][t])[0],:np.shape(vs[p]['WF_vector'][t])[1]] = vs[p]['WF_vector'][t]
            if source == 'GPOD':
                vs[p]['RIP'][count:count + np.shape(vs[p]['RIP_vector'][t])[0],:np.shape(vs[p]['RIP_vector'][t])[1]] = vs[p]['RIP_vector'][t]

            count += len(vs[p]['WF_vector'][t])
        del vs[p]['WF_vector']
        if source == 'GPOD':
            del vs[p]['RIP_vector']

    end = datetime.now()
    timer = end-start
    print(timer)
    
    return vs



def write_wse_files(vs, vs_all, selection=None, vs_to_write=None, folder='Time_Series', key='Zambezi_S3A_VS_'):

    """
    Write time series files
    """
    if selection == None:
        selection = vs.keys()
    if not os.path.exists(folder):
        os.mkdir(folder)
        print('Creating WSE time series folder: '+folder)
    for p in selection:
        (x, y) = p
        if len(vs[p]['height']) > 1:
            if vs_to_write:
                fname = key + vs_to_write[p]['ID']
            else:
                fname = (key+ str(round(vs_all[p]['lon'][0],2))+'_'+
                                               str(round(vs_all[p]['lat'][0], 2))+'.txt')
            with open(os.path.join(folder, fname), 'w') as outfile:
                if 'pass' in vs_all[p].keys():
                    outfile.writelines("# Lat = " + str(y) + ', ' + 'Lon = ' + str(x) + "\n"+
                        "# Data file format" + "\n"+
                                        "# Source: SciHub" + "\n" +
                    "#(1): decimal year"+ "\n"+
                    "#(2): date = Y/m/d"+ "\n"+
                    "#(3): height above geoid (EGM2008)"+ "\n"+
                    "#(4): standard deviation from height"+ "\n"+
                    "#(5): number of meas. in height"+ "\n"+
                    "#(6): Pass number" + "\n")
            
                    for ind, uniq_d in enumerate(vs[p]['TAI']):
                        day_ind = np.where(vs_all[p]['TAI'] == uniq_d)[0]
                        if not np.isnan(vs[p]['height'][ind]):
            #                    filth = np.where((VS[p]['height'][day_ind] > (math.floor(mean_h)-2*std_h)) &
            #                                     (VS[p]['height'][day_ind] < (math.ceil(mean_h)+2*std_h)))[0]
    
                            outfile.writelines((',').join([str(round(year_fraction(uniq_d),8)),
                                                    datetime.strftime(uniq_d,'%Y/%m/%d'),
                                                           str(vs[p]['height'][ind]),
                                                           str(np.std(vs_all[p]['height'][day_ind])),
                                                           str(len(vs_all[p]['height'][day_ind])),
                                                           str(int(vs_all[p]['pass'][0])).zfill(3)]))
                            outfile.writelines('\n')
    
                elif 'orbit' in vs[p].keys():
                    outfile.writelines("# Lat = " + str(y) + ', ' + 'Lon = ' + str(x) + "\n"+
                        "# Data file format" + "\n"+
                                       '# Source: GPOD' + "\n" +
                    "#(1): decimal year"+ "\n"+
                    "#(2): date = Y/m/d"+ "\n"+
                    "#(3): height above geoid (EGM2008)"+ "\n"+
                    "#(4): standard deviation from height"+ "\n"+
                    "#(5): number of meas. in height"+ "\n"+
                    "#(6): Orbit number" + "\n")
            
                    for ind, uniq_d in enumerate(vs[p]['TAI']):
                        day_ind = np.where(vs_all[p]['TAI'] == uniq_d)[0]
                        if not np.isnan(vs[p]['height'][ind]):
            #                    filth = np.where((VS[p]['height'][day_ind] > (math.floor(mean_h)-2*std_h)) &
            #                                     (VS[p]['height'][day_ind] < (math.ceil(mean_h)+2*std_h)))[0]
    
                                outfile.writelines((',').join([str(round(year_fraction(uniq_d),8)),
                                    datetime.strftime(uniq_d,'%Y/%m/%d'),
                                           str(vs[p]['height'][ind]),
                                           str(np.std(vs_all[p]['height'][day_ind])),
                                           str(len(vs_all[p]['height'][day_ind])),
                                           str(int(vs[p]['orbit'][ind])).zfill(3)]))
        
                                outfile.writelines('\n')

if __name__ == '__main__':
    
    # s3a_folder_g2 = r'..\..\test\GPOD_subset\s3a_2bin'
    # s3a_folder_g3 = r'..\..\test\GPOD_subset\s3a_3bin'
    # s3b_folder_g2 = r'..\..\test\GPOD_subset\s3b_2bin'
    # s3b_folder_g3 = r'..\..\test\GPOD_subset\s3b_3bin'
    

    s3a_folder_s = r'C:\test\SciHub_subset\s3a'
    
    s3a_folder_stack = r'..\..\..\test\GPOD_subset\s3a_stacks'
    s3a_folder_stackfolder = r'..\..\..\test\GPOD_subset\s3a_stacks\SARStacks'



    
    # s3a_folder_stack = r'..\..\test\GPOD_subset\s3a_stacks'
    s3a_folder_stackfolder = r'..\..\test\GPOD_subset\s3a_stacks\SARStacks'

    # vs_s3a = gpd.read_file(r'..\..\test\Zambezi_VS/S3A_VS_all_w.shp')
    # vs_s3b = gpd.read_file(r'..\..\test\Zambezi_VS/S3B_VS_all_w.shp')


    vs_s3a_Amur = gpd.read_file(r'h:\RACZIW\Amur\S3A_VS_Amur_coord.shp')
    vs_s3b_Amur = gpd.read_file(r'h:\RACZIW\Amur\S3B_VS_Amur_coord.shp')
    
    # wm_file = r'..\..\test\occurrence_zambezi_kafue.tif'
    
    # wm_folder_S3A = r'..\..\test\New_VS_Nature'
    # wm_folder_S3B = r'..\..\test\New_VS_Nature'
    wm_folder_Amur = r'h:\RACZIW\Amur\GSWE\submasks'
    
    # s3_folder2 = r'E:\Sentinel_3\L2_WDIR\L2__RESDIR'
    # s3_folder1 = r'E:\Sentinel_3\L2_WDIR\L1b_RESDIR\STACK'

    # s3_folder1 = r'E:\L2_WDIR\L1b_RESDIR\STACK'
    # s3_folder2 = r'E:\L2_WDIR'


    # vs_s3a_s, outliers_s3a_s = read_s3_nc(s3a_folder_s, vs_coords=vs_s3a,
    #                                       wm_folder=wm_folder_S3A, source='SciHub', 
    #                                         dem_file=r'C:\test\merit_egm2008.tif',
    #                                         sigma_thresh=30, dem_thresh=30, vs_buffer=0.015, rip_thresh=1e-13,
    #                                         stack=False, stack_folder=s3a_folder_stackfolder)
    
    
    vs_s3a_s_d = create_vs_ts(vs_s3a_s, subset_vs=False, sigma_thresh=30, source='SciHub')
    write_wse_files(VSA_d, VSA, s3a_valid_samosa+s3a_valid_samosa_o, VS_to_write_A, folder=r'..\..\test\Time_Series', key='Zambezi_S3A_GPOD_2bin_VS_')
    write_wse_files(VSAE_d, VSAE, s3a_valid_ocog_o+s3a_valid_ocog, VS_to_write_A, folder=r'..\..\test\Time_Series', key='Zambezi_S3A_SciHub_VS_')

    write_wse_files(VSB_d, VSB, s3b_valid_samosa, VS_to_write_B, folder=r'..\..\test\Time_Series', key='Zambezi_S3B_GPOD_2bin_VS_')
    write_wse_files(vs_s3a_s_d, vs_s3a_s, folder=r'C:\test\Time_Series', key='Zambezi_S3B_SciHub_VS_')

