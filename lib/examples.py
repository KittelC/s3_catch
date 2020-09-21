# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 16:16:32 2020

@author: ceki
"""

from s3_catch import s3_utils, crop_raster, s3_subset_ncdf, s3_preprocessing
import geopandas as gpd


if __name__ == '__main__':
    
    """
    This script is intended as a step-by-step example of using the 
    Sentinel-3 code. Prior to running the script, you need:
        - Shapefiles containing the coordinates of the virtual stations in columns (xcoord and ycoord)
        - DEM of the region of interest
        - Water mask of the region of interest
    """
    
    # First preliminary step: Exctract water mask subsets
    # Define your raster for the entire AOI
    full_raster = r'h:\RACZIW\Amur\GSWE\occurrence_Amur_clip.tif'
    # where to save the subsets
    dest_folder = r'h:\RACZIW\Amur\GSWE\submasks\\'
    if not os.path.exists(dest_folder):
        os.mkdir(dest_folder)
    # Define extent and name for masks:  
    ulc = (107, 56)
    lrc = (141, 41)
    name='water_mask_'
    
    # Subset raster:
    crop_raster.vs_wm(full_raster, dest_folder, lrc, ulc, subset_size=0.5, name=name)   
    
    # Second preliminary step (ONLY SCIHUB): 
    # Directory with full netcdf files:
    download_dir = r'..\..\..\test\SciHub\Full'
    # Where to store subsets:
    dest_dir = os.path.join(download_dir, r'..\Subsets2')
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)
      
    # AOI:
    extent = [19, -9, 36, -20] #upper left x, upper left y, lower right x, lower right y
  
    # Subset netcdf files:
    s3_subset_ncdf.subset_scihub_netcdf(download_dir, dest_dir, extent, file_id=r'enhanced_measurement.nc')

        
    # First processing step: 
    # Define folders with relevant netcdf files to be processed:
    # GPOD folders:
    # s3a_folder_g2 = r'..\..\test\GPOD_subset\s3a_2bin'n
    # s3a_folder_g3 = r'..\..\test\GPOD_subset\s3a_3bin'
    # s3b_folder_g2 = r'..\..\test\GPOD_subset\s3b_2bin'
    # s3b_folder_g3 = r'..\..\test\GPOD_subset\s3b_3bin'
    
    # SciHub folders:
    s3a_folder_s = r'C:\test\SciHub_subset\s3a'
    
    # GPOD folders with Level-1b data.
    s3a_folder_stack = r'..\..\..\test\GPOD_subset\s3a_stacks'
    s3a_folder_stackfolder = r'..\..\..\test\GPOD_subset\s3a_stacks\SARStacks'

    # Define VS coordinates:
    # vs_s3a = gpd.read_file(r'..\..\test\Zambezi_VS/S3A_VS_all_w.shp')
    # vs_s3b = gpd.read_file(r'..\..\test\Zambezi_VS/S3B_VS_all_w.shp')

    vs_s3a_Amur = gpd.read_file(r'h:\RACZIW\Amur\S3A_VS_Amur_coord.shp')
    vs_s3b_Amur = gpd.read_file(r'h:\RACZIW\Amur\S3B_VS_Amur_coord.shp')
    
    # For floodplain processing, a single water mask is used - for VS processing, 
    # subsets in a folder (see above) are used:
    # wm_file = r'..\..\test\occurrence_zambezi_kafue.tif'
    
    # wm_folder_S3A = r'..\..\test\New_VS_Nature'
    # wm_folder_S3B = r'..\..\test\New_VS_Nature'
    wm_folder_Amur = r'h:\RACZIW\Amur\GSWE\submasks'
    

    # Second step: Process netcdf files - extracts all data for all VS
    # Returns two dictionaries: 
        # VS - contains the virtual stations in dictionary form
        # outliers - contains information about the removed points for each VS
    vs_s3a_s, outliers_s3a_s = s3_preprocessing.read_s3_nc(s3a_folder_s, vs_coords=vs_s3a,
                                          wm_folder=wm_folder_S3A, source='SciHub', 
                                            dem_file=r'C:\test\merit_egm2008.tif',
                                            sigma_thresh=30, dem_thresh=30, vs_buffer=0.015, rip_thresh=1e-13,
                                            stack=False, stack_folder=s3a_folder_stackfolder)
    
    # Third step: Create time series - calculates along-track means to produce
    # daily wSE observations
    vs_s3a_s_d = create_vs_ts(vs_s3a_s, subset_vs=False, sigma_thresh=30, source='SciHub')
    
    
    # Fourth step: Write text files with the observations at each VS
    # first argument: vs from step 2
    # second argument: vsd from step 3
    # optional arguments: 
        # selection: if only writing subset, list of coordinates (as keys of vs and vsd)
        # vs_to_write: if using different VS names, use to provide "translation" dictionary
        # folder and key: where to store and common name suffix
    write_wse_files(VSA_d, VSA, s3a_valid_samosa+s3a_valid_samosa_o, VS_to_write_A, folder=r'..\..\test\Time_Series', key='Zambezi_S3A_GPOD_2bin_VS_')
    write_wse_files(VSAE_d, VSAE, s3a_valid_ocog_o+s3a_valid_ocog, VS_to_write_A, folder=r'..\..\test\Time_Series', key='Zambezi_S3A_SciHub_VS_')

    write_wse_files(VSB_d, VSB, s3b_valid_samosa, VS_to_write_B, folder=r'..\..\test\Time_Series', key='Zambezi_S3B_GPOD_2bin_VS_')
    write_wse_files(vs_s3a_s_d, vs_s3a_s, folder=r'C:\test\Time_Series', key='Zambezi_S3B_SciHub_VS_')


    # Evaluate results (not against observations):
        # mostdata: at least 80% of expected observations
        # windext: improvement between 2x and 3x extension on GPOD
        # postoltc: improvement after S3A OLTC update
    mostdata, windext, postoltc = sort_l2(vs, outliers, vsd, oltc=False, oltc_date=datetime(2019,3,1).date(),
            vs3=None, outliers3=None, vsd3=None)
    # Level-1b evaluation.
        # also divided through total (s3-valid), 3x extension (s3_valid_3)
        # OLTC update (s3_valid_o)
    s3_valid, s3_valid_3, s3_valid_o = sort_l1b(vsd, mostdata, windext, postoltc, ext=None, oltc=None)