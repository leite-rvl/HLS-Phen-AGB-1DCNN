"""
Extracting values from a list images at GEDI footprint location. 
Inputs: 
gpkg with points 
list of images
grid gpkg for parallel processing  


"""

### EXporting outputs to
# /panfs/ccds02/nobackup/projects/hls/HLS/mwooten3/hls_agb/


##TODO
# Check warnings

import os
import re
import sys
import glob
import time
import multiprocessing as mp
from datetime import datetime, timedelta

from collections import OrderedDict

import fiona
import numpy as np
import pandas as pd
import geopandas as gpd

from pyproj import CRS
from shapely.geometry import box, shape

import rasterio
from rasterio.windows import from_bounds
from rasterstats import zonal_stats, point_query

# CHANGE
TEST = False

# Project directory
dir_prj = '/panfs/ccds02/nobackup/people/rvieiral/projects/GEDI_EdgeEffects/'
# Directory for shapefiles
dir_shp = f'{dir_prj}/data/shp/'
# Directory for tif
# dir_tif = f'{dir_prj}/data/tif/'
# tif directories set below nwo


##########################################
# pick a rasterstats type
extract_type = 'zonal'#'zonal' # zonal or point
product = 'HLS.30m'#'S1-SAR'# 'HLS.30m'## HLS.30m or S1-SAR 

ref_source = 'l2al4a_metrics'#'chm.all'#'chm.lp'#'chm.mcd'#'chm.all'#'chm.saojose'#'l2al4a_metrics'# for GEDI # chm.* for eg hls_agb/data/vector/chm/chm_1m_mcd_test6__randomPoints.gpkg

# for HLS - can leave as is for SAR
exclude_bands = ['Fmask','Edge1','Edge2','Edge3','NIR_Broad','ValidCount']

# for zonal extraction
zonal_stats_list = ['median','std']
if extract_type != 'zonal': 
    zonal_stats_list = None

##########################################
# #Input/output .gpkg files
# input_gpkg_path = '/explore/nobackup/people/rvieiral/projects/HLS/data/shp/gedi/class_qfilt_lsmetrics_l2al4a_mg.gpkg'
input_gpkg_path = '/explore/nobackup/people/rvieiral/projects/HLS/data/shp/gedi/class_qfilt_lsmetrics_l2al4a_mg__clean.gpkg'

grid_path = f'{dir_shp}/base/ls_unit20km.shp'

output_dir = '/explore/nobackup/projects/hls/HLS/mwooten3/hls_agb/data/vector/zonal_point_stats'
if TEST: output_dir = os.path.join(output_dir, '_tests')

#* 7/25
# if input ref is chm, get chm random points and other files specific to it
if ref_source.startswith('chm.'):
    site = ref_source.split('.')[1]
    input_gpkg_path = os.path.join('/explore/nobackup/projects/hls/HLS/mwooten3/hls_agb/data/vector/chm/', f'chm_1m_{site}_test6__randomPoints.gpkg')

    if ref_source == 'chm.all':
        input_gpkg_path = os.path.join('/explore/nobackup/projects/hls/HLS/mwooten3/hls_agb/data/vector/chm/', 'chm_1m_allSites_test6__randomPoints-zonalStats.gpkg')

    output_dir = os.path.join(output_dir, 'chm')
    
    # grid_path = '/explore/nobackup/projects/hls/HLS/mwooten3/hls_agb/data/grid/ls_unit20km_chmSub.gpkg' # old grid subset for the two filesls_unit20km
    grid_path = '/explore/nobackup/projects/hls/HLS/mwooten3/hls_agb/data/grid/ls_unit20km_chm.all-subset.gpkg' # new subset for all 3 (3 tiles)

os.makedirs(output_dir, exist_ok=True)

## INPUT TIFS TO EXTRACT
# Default output .gpkg proj - this can be overridden by passing different epsg to calculate_zonal_stats
default_output_epsg = 4326

if product == 'HLS.30m':

    prodtype = "seasonal"#"annual.2020.2022"#"2020.2022" # "2020.2022"# or "annual.2020.2022"# #"seasonal"
    
    # 2020-2022 composites
    if prodtype == "2020.2022":
        tifs_list = glob.glob('/explore/nobackup/projects/hls/HLS/mwooten3/Composites/HLS.30m/_vrt/??/*.vrt')
        # test one chunk, point: Processing completed in 0.39 minutes. ( 0.41 minutes for reproj test)
        # test one chunk, zonal: Processing completed in 9.92 minutes. ( ?? minutes for reproj test)
        
    # Monthly composites for entire time series
    elif prodtype.startswith("annual"): 
        tifs_list = glob.glob('/explore/nobackup/projects/hls/HLS/mwooten3/Composites/HLS.30m/_vrt/???????.???????/*.vrt')
        # test one chunk, point: Processing completed in 0.90 minutes (0.01 hours)
        # test one chunk, zonal: Processing completed in 33.00 minutes.
        # point all chunks (9-tile .vrt): Processing completed in 152.89 minutes (2.55 hours).

    # Seasonal composites, e.g. /explore/nobackup/projects/hls/HLS/mwooten3/Composites/HLS.30m/_seasonal_composites/2020.01-2020.04/HLS.M30.2020.01-2020.04.median.v2.0.Red.tif
    # 2020.01-2020.04, 2020.05-2020.10, 2020.11-2021.04, 2021.05-2021.10, 2021.11-2022.04, 2022.05-2022.10
    # these median composites are .tif not vrt
    elif prodtype == "seasonal":
        comp_stat = 'median'
        tifs_list = glob.glob(f'/explore/nobackup/projects/hls/HLS/mwooten3/Composites/HLS.30m/_seasonal_composites/????.??-????.??/*{comp_stat}*.tif')        

elif product == 'S1-SAR':

    # prodtype = "2020.2021" # "2020.2021" or "annual.2020.2022"#"2020.2021"#
    prodtype = "2020.2021"#"annual.2020.2022"#"2020.2021"#

    # vrt of vrts:
    # /S1-SAR.2020.2021.07.VV.MonthlyMedian.vrt
    if prodtype == "2020.2021":
        tifs_list = glob.glob('/explore/nobackup/projects/hls/HLS/mwooten3/hls_agb/data/SAR/_vrt/_2020.2021/??/*vrt')
        # test one chunk, point: Processing completed in 0.26 minutes.
        # test one chunk, zonal: Processing completed in 2.92 minutes.
        # zonal all chunks (9-tile .vrt): Processing completed in 578.24 minutes (9.64 hours)

    # individual vrts:
    elif prodtype.startswith("annual"): 
        tifs_list = glob.glob('/explore/nobackup/projects/hls/HLS/mwooten3/hls_agb/data/SAR/_vrt/????.??/*vrt')
        # test one chunk, point: Processing completed in 0.84 minutes.
        # test one chunk, zonal: Processing completed in 8.15 minutes.
        
output_gpkg = os.path.join(output_dir, f'{ref_source}_{product}-{prodtype}_{extract_type}.gpkg')
if TEST:
    output_gpkg = os.path.join(output_dir, f'{ref_source}_{product}-{prodtype}_{extract_type}-TEST.gpkg')

print(f"\nUsing input .gpkg {input_gpkg_path}")
print(f"\nOutput .gpkg: {output_gpkg}\n")

# print(input_gpkg_path)
# sys.exit()

##########################################
# Given an input gpkg schema, output gpkg/crs and some parameters, create an empty output gpkg
# *NOTE this will have problems with mistmatch if columns are different
def create_empty_output(ingpkg, gpkg, epsg, var_names, zonal_stats_list):

    in_gdf = gpd.read_file(input_gpkg_path, rows=1)

    # to hold input + output columns
    all_columns = {}

    # get input columns
    for col in in_gdf.columns:
        if col != 'geometry': 
            all_columns[col] = [np.nan]  # One NaN value
    
    # # same structure but no data, remove all rows, set to output epsg
    # empty_gdf = in_gdf.iloc[0:0].copy().to_crs(f"EPSG:{epsg}")
    # add the new columns from the vars and parameters
    if zonal_stats_list is not None: # extract type = zonal
        for var_name in var_names:
            for stat in zonal_stats_list:
                # new_cols.append(f'{var_name}_{stat}')
                all_columns[f'{var_name}_{stat}'] = [np.nan]#pd.Series(dtype='float')

    else:
        # new_cols = var_names
        for var_name in var_names:
            all_columns[var_name] = [np.nan]#pd.Series(dtype='float')

    # put empty_gdf in source crs, then convert to out epsg if its different 
    empty_gdf = gpd.GeoDataFrame(all_columns, geometry=in_gdf.geometry, 
                                                             # crs=f"EPSG:{epsg}")
                                                               crs = in_gdf.crs)
    if CRS.from_epsg(epsg) != empty_gdf.crs:
        empty_gdf = empty_gdf.to_crs(CRS.from_epsg(epsg))

    empty_gdf.to_file(gpkg, driver='GPKG')
        
    return

 
def julian_range_to_month(julian_range: str) -> str:
    """
    Converts a Julian range string like '2020001.2020031' to '2020.01',
    using the midpoint of the date range to determine the month.
    """
    start_str, end_str = julian_range.split(".")
    start_year, start_julian = int(start_str[:4]), int(start_str[4:])
    end_year, end_julian = int(end_str[:4]), int(end_str[4:])
    start_date = datetime(start_year, 1, 1) + timedelta(days=start_julian - 1)
    end_date = datetime(end_year, 1, 1) + timedelta(days=end_julian - 1)
    midpoint = start_date + (end_date - start_date) / 2
    return f"{midpoint.year}.{midpoint.month:02d}"

def transform_filename(filename: str) -> str:
    # Match Julian range
    match = re.search(r'\.(\d{7}\.\d{7})', filename)
    if match:
        julian_range = match.group(1)
        new_date_str = julian_range_to_month(julian_range)
        return filename.replace(julian_range, new_date_str)
    return filename

#*MW - adding extenson for vrt/tif
def filter_input_tifs(tif_paths: list, extension = '.tif') -> list:
    
    out_tifs = tif_paths.copy()
    for tif in tif_paths:
        # import pdb; pdb.set_trace()
        for band in exclude_bands:
            
            if tif.endswith(f'{band}{extension}'):
                out_tifs.remove(tif)
    return out_tifs

#* 
def fiona_crs_to_crs(fiona_crs):
    
    # for whatever reason, fiona or the gpkg may have been written with weird crs such that src.crs looks like e.g. {'init': 'epsg:4326'}
    # FutureWarning: '+init=<authority>:<code>' syntax is deprecated. 
    # '<authority>:<code>' is the preferred initialization method. When 
    # making the change, be mindful of axis order changes: https://pyproj4.github.io/pyproj/stable/gotchas.html#axis-order-changes-in-proj-6
    #* solution here to make function that will convert the fiona crs whether old or new, to CRS obj
    # gpkg_crs = fiona_crs_to_crs(src.crs)

    # try to grab epsg code from fiona_crs dict
    if not isinstance(fiona_crs, dict):
        print(f"ERROR with retrieving crs from fiona crs: {fiona_crs}. Expected dictionary")
        return

    try:
        fiona_epsg = fiona_crs['authority']
    except KeyError:
        fiona_epsg = fiona_crs['init']
    except KeyError:
        print(f"ERROR with retrieving crs from fiona crs: {fiona_crs}. Expected dictionary with 'authority' or 'init' key")
        return

    return CRS.from_string(fiona_epsg)
    

# Function to calculate zonal stats for a chunk of data
def calculate_zonal_stats(gpkg_path, layer_name, tif_paths, var_names, ndVal, 
                       bbox_geom, bbox_crs, out_gpkg, out_epsg=4326, lock=None):
    
    #################################
    print('in calculate_zonal_stats()')
    # Open the gpkg and filter features within the bounding box
    ## If no features are found in the bbox, return None
    ## The file is partially opened for the bbox_geom (defined grid)

    # import pdb; pdb. set_trace()
    # get output crs from epsg (or default to the one at top if None)
    if out_epsg is None:
        out_epsg = default_output_epsg
    out_crs = CRS.from_epsg(out_epsg)
    # print(f"Output epsg = {out_epsg}")
    # print(f"Output crs  = {out_crs}")
        
    # make sure they are in the same crs
    # try:# avoid having to dedent everything below try for now lol
    if True:
        # Get crs from fiona and convert to CRS object
        # Then read filtered features into gdf_chunk
        #################################
        with fiona.open(gpkg_path, layer=layer_name) as gpkg:
            gpkg_crs = fiona_crs_to_crs(gpkg.crs)
           
            filtered_features = [
                feature for feature in gpkg.filter(bbox=bbox_geom) 
                if shape(feature['geometry']).intersects(box(*bbox_geom))]

        if not filtered_features:
            # print("No features found in the bounding box.")
            return None
            
        # Convert filtered features into a GDF and apply CRS transformations 
        gdf_chunk = gpd.GeoDataFrame.from_features(filtered_features, crs=gpkg_crs)
        # gdf_chunk = gdf_chunk.head(100) #TEST
        #gdf_chunk.to_file('/panfs/ccds02/nobackup/projects/hls/HLS/mwooten3/hls_agb/data/vector/chunk.gpkg', f='GPKG')

        #################################    
        minx, miny, maxx, maxy = bbox_geom

        # Adding buffer to filter images (100 m)
        # Determine whether our bbox_geom is geographic or projected for buffering extent for raster clipping
        if bbox_crs.is_geographic:
            minx -= 0.0009
            maxx += 0.0009
            miny -= 0.0009
            maxy += 0.0009
        elif bbox_crs.is_projected: # assume meters
            minx -= 100
            maxx += 100
            miny -= 100
            maxy += 100

        # Now convert padded extents for chunk back to geom/gpd objs
        bbox_geom = box(minx, miny, maxx, maxy)
        bbox_gdf = gpd.GeoDataFrame(geometry=[bbox_geom], crs=bbox_crs) # necessary?

        #################################    
        #*TODO: find UTM zone (or general equal area projection for buffering - 
        # Save geometry and create buffer if runnign zonal stats
        # for now, move this to loop with .tifs)
        #* BUT, check to be sure this makes sense, because it might mess up since original_geometry part is in the input crs - but maybe it's ok, because we drop original geometry etc before converting to output proj? idk, just check
        # if zonal_stats_list is not None: # extract_type == 'zonal'
        #     # Save original (point) geometry
        #     gdf_chunk['original_geometry'] = gdf_chunk['geometry'].copy()
        #     # Create buffer with GEDI footprint size
        #     ### Approximate degrees equivalent to 25 meters (Earth ci ~111.32 km)
        #     if gdf_chunk.crs.is_geographic:
        #         gdf_chunk['geometry'] = gdf_chunk.geometry.buffer(0.000225)
        #     elif gdf_chunk.crs.is_projected:
        #         gdf_chunk['geometry'] = gdf_chunk.geometry.buffer(25)

        # instead of modifying gdf_chunk with extracted columns in iteration, 
        # create dict then concat when finished
        extracted_data = {}
        
        ################################# 
        # Loop through each raster file and calculate zonal statistics for the polygons
        for tif_path, var_name in zip(tif_paths, var_names):
            # print(f'{var_name} : {tif_path}')

            # Check if tif intersects chunk bbox
            with rasterio.open(tif_path) as src: 
                raster_bounds = src.bounds
                raster_geom = box(*raster_bounds)
                raster_crs = src.crs

                # reproject padded chunk bbox to projection of tif if necessary (also get new vars)
                if bbox_gdf.crs != raster_crs:
                    bbox_gdf = bbox_gdf.to_crs(raster_crs)
                    minx, miny, maxx, maxy = bbox_gdf.total_bounds
                    bbox_geom = box(minx, miny, maxx, maxy)

                #* NOW, regardless, everything (except maybe bbox_crs var, which will always be input bbox coords) is in the proj of .tif

                # Skip tif if extent doesn't intersect grid chunk
                if not bbox_geom.intersects(raster_geom):
                    continue 
                    
                # minx, miny, maxx, maxy = bbox_geom
                window = from_bounds(minx, miny, maxx, maxy, transform=src.transform)
                data = src.read(1, window=window)
                affine = src.window_transform(window)

            # Now, also reproject gdf_chunk if it differs from raster_crs
            if gdf_chunk.crs != raster_crs:
                gdf_chunk = gdf_chunk.to_crs(raster_crs)



            # Run zonal stats 
            if zonal_stats_list is not None: #  if extract_type == 'zonal'

                #* TODO: move this out of loop back up top, but with temp conversion to a projected - see notes above
                # for now, do it here so at least for HLS, chunk we will do projected buffer
                # Save original (point) geometry
                gdf_chunk['original_geometry'] = gdf_chunk['geometry'].copy()
                
                # Create buffer with GEDI footprint size
                ### Approximate degrees equivalent to 25 meters (Earth ci ~111.32 km)
                if gdf_chunk.crs.is_geographic:
                    gdf_chunk['geometry'] = gdf_chunk.geometry.buffer(0.000225)
                elif gdf_chunk.crs.is_projected:
                    gdf_chunk['geometry'] = gdf_chunk.geometry.buffer(25)
                    
                #########################
                # With buffer (zonal stats)
                zonal_stats_result = zonal_stats(gdf_chunk, 
                                                 data, 
                                                 affine=affine, 
                                                 nodata=ndVal,
                                                 stats=zonal_stats_list)

                # Save data extracted from above in the dictionary
                for stat in zonal_stats_list:
                    extracted_data[f'{var_name}_{stat}'] = [feature[stat] for feature in zonal_stats_result]

            # Run point query
            elif zonal_stats_list is None: # if extract_type == 'point'

                values = point_query(gdf_chunk, data, nodata=ndVal, affine=affine)

                # Save data extracted from above in the dictionary
                extracted_data[var_name] = values
            
        if TEST:
            import pdb; pdb.set_trace()
            
        #* Now add extracted values to gdf
        gdf_chunk = pd.concat([gdf_chunk, pd.DataFrame(extracted_data, index=gdf_chunk.index)], axis=1)

        #* IFF we want to drop rows whose new columsn (eg extracted data) are all NaNs
        if False:
            gdf_chunk = gdf_chunk.dropna(subset=extracted_data.keys(), how='all')

        # #################################
        # # Convert gdf_chunk back to points if zonal
        if zonal_stats_list is not None: # if extract_type == 'zonal':
            gdf_chunk['geometry'] = gdf_chunk['original_geometry']
            #Drop temporary column
            gdf_chunk.drop(columns=['original_geometry'], inplace=True)

        # make sure geometry column at the end of chunk_gdf (might not be necessary, might only matter for zonal)
        cols = [col for col in gdf_chunk.columns if col != 'geometry'] + ['geometry']
        gdf_chunk = gdf_chunk[cols]


        # Before writing, reproject output gdf_chunk to speciied output epsg
        if gdf_chunk.crs != out_crs:
            # import pdb; pdb.set_trace()
            gdf_chunk = gdf_chunk.to_crs(out_crs)

        #* TEMP
        if TEST:

            # write the columns for the output we are appending to input
            with open('_appendFrom_gpkg_cols.txt', 'w') as of:
                # of.write([f"{c}\n" for c in gdf_chunk.columns])
                of.write("".join([f"{c}\n" for c in gdf_chunk.columns]))

            # for c in gdf_chunk.columns: print(c)
            # write cols for output we're appending to
            with open('_appendTo_gpkg_cols.txt', 'w') as of:
                gdf2 = gpd.read_file(out_gpkg, rows=1)
                # of.write([f"{c}\n" for c in gdf_chunk.columns])
                of.write("".join([f"{c}\n" for c in gdf2.columns]))
                
            # return #* TEMP
        #* maybe comment out later
        print("GDF_CHUNK:")
        print(gdf_chunk)
        # print("OUT_GPKG:")
        # print(out_gpkg)

        # WORKAROUND for emtpy shp thing - if do it (eg ref_source is shm) check cols before appending
        # overwrite with cols that match but same data
        # THIS can cause lock problems when running in parallel - for now does not
        # matter because chm is only 3 grids and can run with 1 multiprocessor
        if True and ref_source.startswith('chm'):
            outg_gdf = gpd.read_file(out_gpkg)
            out_cols = [c for c in outg_gdf.columns if c in gdf_chunk.columns]

            if out_cols != outg_gdf.columns.to_list():
                outg_gdf = outg_gdf[out_cols]
                gdf_chunk.to_file(out_gpkg, driver='GPKG', mode='w')
                
        ################################# 
        # Write the results to the GeoPackage, using a lock to prevent concurrent writing issues
        output_layer_name = os.path.splitext(os.path.basename(out_gpkg))[0]

        if lock: # for concurrent writing
            with lock:
                gdf_chunk.to_file(out_gpkg, driver='GPKG', layer=output_layer_name, mode='a')

            del gdf_chunk
            return None#gdf_chunk
            
        else:
            gdf_chunk.to_file(out_gpkg, driver='GPKG', layer=output_layer_name, mode='a')
                
            del gdf_chunk
            return None#gdf_chunk

        del gdf_chunk
        return None#gdf_chunk
        
    # except Exception as e:
    #     print(f"Error in calculate_zonal_stats: {e}")
    #     del gdf_chunk
    #     return None

#################################
# Main function: Process input arguments and execute zonal statistics calculation
def main(args):
    
    #################################
    # Read the files and define variable names
    grid_path = args['grid']
    grid = gpd.read_file(grid_path)

    # tif_path = args['inputImage']
    # ndVal = -9999
    # var_name = os.path.splitext(os.path.basename(tif_path))[0]
    
    tif_paths = args['inputImages']
    ndVal = -9999
    if TEST: import pdb; pdb.set_trace()
    
    # remove and .tif/bands that we dont want to run with zonal stats - filter using exclude_bands
    ext = os.path.splitext(tif_paths[0])[1]
    tif_paths = filter_input_tifs(tif_paths, extension=ext)
    tif_paths = sorted(tif_paths) # sort in alpha order to ensue dates are together

    #################################
    # Get var_names (e.g. output column prefixes)
    var_names_original = [os.path.splitext(os.path.basename(tif_path))[0] for tif_path in tif_paths]
    var_names = [transform_filename(fn) for fn in var_names_original]
    # print(var_names[0])

    product = args['product']
    if product == 'S1-SAR':
        # .MonthlyMedian.vrt
        var_names = [os.path.basename(f).strip('.MonthlyMedian.vrt') for f in tif_paths]
        # import pdb; pdb.set_trace()
    
    gpkg_path = args['inputFc']
    layer_name = os.path.splitext(os.path.basename(gpkg_path))[0]

    ##################################
    # Create empty output gpkg to avoid race condition in parallel processing
    out_gpkg = args['outputFc']
    out_epsg = args['outputEpsg']

    if os.path.isfile(out_gpkg):
        # print(f"Output {out_gpkg} exists. Try again!")
        # return # dont run progam if exists
        print(f"Appending to {out_gpkg}, if possible")

    else: # create empty if it doesnt exist
        create_empty_output(gpkg_path, out_gpkg, out_epsg, var_names, zonal_stats_list)
        #*TEMP TEST
        gdf=gpd.read_file(out_gpkg, rows=5)
        with open('_empty_file_cols.txt', 'w') as of:
            of.write("".join([f"{c}\n" for c in gdf.columns]))
    ##################################
    
    # Convert crs: grid and gpkg should be in same crs
    with fiona.open(gpkg_path, layer=layer_name) as src:
        gpkg_crs = fiona_crs_to_crs(src.crs)
        # gpkg_schema = src.schema


    # Now convert grid to gpkg crs if they're different, pass grid.crs as bbox_crs
    if grid.crs != gpkg_crs:
        grid = grid.to_crs(gpkg_crs)

    bbox_crs = grid.crs

    ######################################
    ## Running sample chunk IF TEST
    # if False:
    # if True:  # if test
    if TEST:
        tif_path_0 = tif_paths#[0:2]
        var_name_0 = var_names#[0:2]
        
        # sample_tile = grid.iloc[[379]]  # keep it as GeoDataFrame
        # sample_tile = grid.iloc[[2109]]  # 2109 for 9 test .vrts near BH - keep it as GeoDataFrame
        sample_tile = grid.iloc[[1]] # rp lp
        # sample_tile = sample_tile.set_crs(epsg=epsg)  # ensure CRS is set if not already
        minx, miny, maxx, maxy = sample_tile.total_bounds
        bbox_geom = (minx, miny, maxx, maxy)

        # gdf=gpd.read_file(gpkg_path, rows=5)
        # with open('_empty_file_cols.txt', 'w') as of:
        #     of.write("".join([f"{c}\n" for c in gdf.columns]))
        
        # out_gdf = calculate_zonal_stats(gpkg_path, layer_name, tif_path_0, var_name_0, ndVal, bbox_geom, ref_crs)
        calculate_zonal_stats(gpkg_path, layer_name, tif_path_0, var_name_0, ndVal, 
                                  bbox_geom, bbox_crs, out_gpkg, out_epsg)
        
        return # TEMP

    # run in parallel
    manager = mp.Manager()
    lock = manager.Lock()
    
    if ref_source.startswith('chm'):
        n_proc = 1
    else:
        n_proc = 12
    # with mp.Pool(processes=12) as pool:
    with mp.Pool(processes=n_proc) as pool:
        # for index, row in grid.iloc[2109:2110].iterrows(): # TEST - alternate way to test with one chunk (?) via mp - if doing this, probably best to set processes to 1 (but may bnot matter)
        for index, row in grid.iterrows():

            minx, miny, maxx, maxy = row.geometry.bounds
            bbox_geom = (minx, miny, maxx, maxy)

            pool.apply_async(
                calculate_zonal_stats, 
                args=(gpkg_path, layer_name, tif_paths, var_names, ndVal, 
                           bbox_geom, bbox_crs, out_gpkg, out_epsg, lock))

        pool.close()
        pool.join()
        

if __name__ == "__main__":
       
    start_time = time.time()
    
    args = {
        'inputImages': tifs_list,
        'inputFc': input_gpkg_path,
        'grid': grid_path,
        'product': product,
        'outputFc': output_gpkg,
        'outputEpsg': default_output_epsg
        # 'grid': f'{dir_shp}/base/ls_unit50km.gpkg'
    }
    main(args)

    end_time = time.time()
    elapsed_time = end_time - start_time
     # Convert elapsed time to minutes and hours
    elapsed_minutes = elapsed_time / 60
    elapsed_hours = elapsed_minutes / 60

    print(f"Processing completed in {elapsed_minutes:.2f} minutes ({elapsed_hours:.2f} hours).")
    print(f"Output file: {output_gpkg}")


# if __name__ == "__main__":
#     parser.add_argument('--inputImage', required=True, help='Path to the input TIFF image.')
#     parser.add_argument('--inputFc', required=True, help='Path to the input GeoPackage file.')
#     parser.add_argument('--grid', required=True, help='Path to the grid GeoPackage file.')

#     args = parser.parse_args()
#     main(args)




################################
# ####### TEST for loop 

#     for index, row in grid.iterrows():
#         row_i = gpd.GeoDataFrame(pd.DataFrame(row).T,  
#                                  crs='EPSG:{}'.format(epsg))
#         # Get the bbox
#         bbox_df = row_i.geometry.bounds.iloc[0]
#         minx, miny, maxx, maxy = bbox_df['minx'], bbox_df['miny'], bbox_df['maxx'], bbox_df['maxy']
#         bbox_geom = (minx, miny, maxx, maxy)

#         ####################
#         ## Parallel process the bbox

#         gdf_concat = calculate_zonal_stats(gpkg_path,
#                                            layer_name,
#                                            tif_path,
#                                            var_name,
#                                            ndVal, 
#                                            bbox_geom)
""
