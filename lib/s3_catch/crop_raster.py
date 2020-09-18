# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 09:34:37 2019

@author: Cecile M. M. Kittel
Copyright: Cecile M. M. Kittel and Technical University of Denmark
"""

from datetime import datetime, timedelta
import os
from osgeo import gdal
gdal.UseExceptions()  # Enable errors

import geopandas as gpd
import numpy as np


def vs_wm(full_raster, dest_folder, lrc, ulc, subset_size=0.5, name='water_mask_'):
    """
    

    Parameters
    ----------
    full_raster : string - file path
        Water mask or DEM covering full catchment.
    dest_folder : string - file path
        Folder to store subsets of full raster.
    llc : tuple (x, y)
        Lat lon coordinates of lower left corner.
    urc : tuple (x, y)
        Lat lon coordinates of upper right corner.

    Returns
    -------
    The destination folder contains the raster subsets covering the entire area of interest.

    """
    (x1, y0) = lrc
    (x0, y1) = ulc
    if x1 < x0 or y1 < y0:
        raise ValueError('Study area extent incorrect - please check coordinates.')
    
    xext = np.arange(x0, x1, 2*subset_size)
    yext = np.arange(y0, y1, 2*subset_size)
    
    
    ds = gdal.Open(full_raster)
    ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
    lrx = ulx + (ds.RasterXSize * xres)
    lry = uly + (ds.RasterYSize * yres)
    
    if x1 < ulx or x0 > lrx or y1 < lry or y0 > uly:
        raise ValueError('Study area extent does not overlap with full raster ' +
                         '- please check coordinates.')

    for x in xext:
        for y in yext:
            extent = " ".join([str(x-subset_size), str(y+subset_size), 
                               str(x+subset_size), str(y-subset_size)]) #ulx uly lrx lry
        
            translateoptions = gdal.TranslateOptions(gdal.ParseCommandLine("-of Gtiff -co COMPRESS=LZW -projwin " + extent))
        
            dst_filename = dest_folder + 'subset_' + str(np.round(x, 1)) + '_' + str(np.round(y, 1)) + '.tif'
            ds_out = gdal.Translate(dst_filename, ds, options = translateoptions)
            ds_out = None
    
    ds = None


