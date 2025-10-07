

# """
# Script to compute zonal statistics for GEDI L4A points against a multiband HLS raster.
# Steps:
# 1. Load GEDI L4A geospatial points (optionally buffer it into polygons - ## check errors on rastererization).
# 2. Load single band image - e.g. HLS composite raster for one band.
# 3. Reproject points to raster CRS and rasterize zones.
# 4. Compute per-zone statistics (e.g., mean) using xrspatial.zonal_stats.
# 5. Join results back to the input GeoDataFrame and export to GeoPackage.
# """

# In[1]:


import os
import glob
import numpy as np
import geopandas as gpd
import rioxarray
from rasterio import features
from shapely.geometry import Point
from xrspatial import zonal_stats
import pandas as pd
from pyproj import CRS



################################
# Functions
################################

def build_zones_xarr(gpkg_points_path,
                     img_path,
                     ZONE_ID_COL='zone',
                     BUFFER_METERS=0,
                     all_touched=False):
    """
    Create a zone-labeled xarray.DataArray aligned to `img_path` by rasterizing features
    from `gpkg_points_path`. Works with point or polygon geometries. If points are given,
    optional buffering creates circular polygons before rasterization.

    Parameters
    ----------
    gpkg_points_path : str
        Path to a GeoPackage with point or polygon geometries.
    img_path : str
        Path to the raster whose grid/CRS/transform will be used.
    ZONE_ID_COL : str, default 'zone'
        Column in the GPKG that holds the integer zone IDs. If missing, it will be created.
    BUFFER_METERS : float, default 0
        Buffer distance in meters to apply to geometries (useful for points).
    all_touched : bool, default False
        Pass-through to rasterio.features.rasterize (include all pixels touched by geometries).

    Returns
    -------
    zones_xarr : xarray.DataArray
        Integer-labeled zone raster aligned to `img_path` grid/CRS/transform.
        Background pixels are set to the value of the source raster's nodata.
    """

    # 1) Open reference raster and grab grid metadata
    da = rioxarray.open_rasterio(img_path).squeeze()
    if 'band' in da.dims:
        da = da.squeeze('band', drop=True)

    raster_crs = CRS.from_user_input(da.rio.crs)
    transform  = da.rio.transform()
    out_shape  = da.shape
    nodata     = da.rio.nodata

    # 2) Read vector data, ensure zone IDs exist, project to raster CRS
    gdf = gpd.read_file(gpkg_points_path)
    if ZONE_ID_COL not in gdf.columns:
        gdf = gdf.reset_index(drop=True)
        gdf[ZONE_ID_COL] = np.arange(1, len(gdf) + 1, dtype=np.int32)

    gdf_proj = gdf.to_crs(raster_crs)

    # Optional buffering (useful for points)
    if BUFFER_METERS and BUFFER_METERS > 0:
        gdf_proj["geometry"] = gdf_proj.geometry.buffer(BUFFER_METERS)

    # 3) Prepare (geometry, id) tuples and rasterize to the raster grid
    geom = list(zip(gdf_proj.geometry, gdf_proj[ZONE_ID_COL].astype(np.int32)))

    zones_arr = features.rasterize(
        geom,
        out_shape=out_shape,
        transform=transform,
        fill=nodata,           # keep background consistent with the raster's nodata
        nodata=nodata,         # used when masked=True
        masked=True,           # produce a masked array (background masked)
        dtype="int32",
        all_touched=all_touched
    )

    # 4) Wrap into xarray aligned to the raster (copies coords/attrs, not data)
    zones_xarr = da.copy(deep=False)
    zones_xarr.data = zones_arr

    return zones_xarr



def zonal_stats_raster(*,gpkg_points_path = None, zones_xarr = None, 
                       img_path = None, ZONE_ID_COL = 'zone', BUFFER_METERS =0, band_name = 'b'):
        
    ####################################################
    # 1) Open raster as xarray (single band)
    ####################################################
    
    da = rioxarray.open_rasterio(img_path).squeeze()  # [y, x]
    if 'band' in da.dims:
        da = da.squeeze('band', drop=True)
        
    raster_crs = CRS.from_user_input(da.rio.crs)
    transform  = da.rio.transform()
    out_shape  = da.shape  # (rows, cols)
    nodata     = da.rio.nodata
    
    da = da.where(da != da.rio.nodata)
    
    
    
    ####################################################
    # 2) Create zones array if not input
    ####################################################
    if zones_xarr is None:
        if gpkg_points_path is not None:
            zones_xarr = build_zones_xarr(
                gpkg_points_path=gpkg_points_path,
                img_path=img_path,
                ZONE_ID_COL=ZONE_ID_COL,
                BUFFER_METERS=BUFFER_METERS
            )
        else:
            raise ValueError("You must provide either zones_xarr or gpkg_points_path.")

    # Sanity check
    if zones_xarr.shape != da.shape:
        raise ValueError("zones_xarr shape does not match raster shape.")
    
    
    ####################################################
    # 7) Compute zonal stats using xrspatial (min/max/mean/etc.)
    ####################################################
   
    # If your raster has nodata, pass it so stats ignore it
    zs_df = zonal_stats(
        zones=zones_xarr,      # your integer-labeled zone raster
        values=da,             # the raster with values to summarize
        # stats_funcs=['mean', 'max', 'min', 'sum', 'std', 'var', 'count'],
        stats_funcs=['mean'],
        nodata_values=nodata,   # very important → ensures nodata is ignored
        return_type='pandas.DataFrame'
    )
    
    
    ####################################################
    # # Keep only "mean" and rename the band
    ####################################################
    # zs_df = pd.DataFrame({"mean": pd.Series(zs["mean"])}).rename_axis(ZONE_ID_COL).reset_index()
    
    # Rename column -> e.g. "evi2_mean"
    zs_df = zs_df.rename(columns={"mean": f"{band_name}_mean"})
    
    
    ####################################################
    # 9) create output with zone id and values
    ####################################################
    # out = (gdf[[ZONE_ID_COL]]
    #        .drop_duplicates()
    #        .merge(zs_df, on=ZONE_ID_COL, how="left"))
    
   
    
    return zs_df
    
    


################################################################
# Define input / output directories  
################################################################

###################
# GEDI points path
# Define reference year and tile of the GEDI points
ref_year = 2022
ref_tile = 89
gpkg_points_path = f'/projects/my-private-bucket/HLS-1DCNN-AGB/data/shp/gedi/test/l4a_t90km_t{ref_tile}_veg{ref_year}_outrm.gpkg'


###################
# Directory with your raster .vrt files
dir_img = '/projects/my-private-bucket/HLS-1DCNN-AGB/data/tif/HLS_composites/monthly/br_af_grid60km_prj_evi2_max/vrt/'
# Define excluded band suffixes
exclude_bands = ['ValidMask.vrt', 'count.vrt', 'yearDate.vrt', 'JulianDate.vrt']


###################
# Output gpkg with the zonal stats of the bands
output_gpkg_zonalstats_fn = f'/projects/my-private-bucket/HLS-1DCNN-AGB/data/shp/gedi/l4a_t90km_t{ref_tile}_veg{ref_year}_outrm_zonal_HLS.gpkg'

###################
# Zone id column name and point buffer
ZONE_ID_COL = 'zone'
BUFFER_METERS = 0


###################
# Get list of imgs
img_paths = glob.glob(os.path.join(dir_img, "*.vrt"))
img_paths = [
    f for f in img_paths
    if not any(f.endswith(ex) for ex in exclude_bands)
]



################################################################
# Read input gdf that where the results will be merged to  
################################################################
# Load the base gdf
gdf = gpd.read_file(gpkg_points_path)
# Create ID for the zonal stats
if ZONE_ID_COL not in gdf.columns:
    gdf = gdf.reset_index(drop=True)
    gdf[ZONE_ID_COL] = np.arange(1, len(gdf) + 1, dtype=np.int32)





################################################################
# Zonal stats  
################################################################


# Create gpkg array to do zonal_stats
zones_xarr = build_zones_xarr(gpkg_points_path,
                     img_paths[0],
                     ZONE_ID_COL='zone',
                     BUFFER_METERS=0,
                     all_touched=False)

error_band = None
error_idx = None

for i, img_path in enumerate(img_paths):
    try:
        os.chdir("/tmp")
        fname = os.path.basename(img_path)
        band_name = fname.replace(".vrt", "")
        print(f'[{i}] Zonal stats: {band_name}')

        out = zonal_stats_raster(zones_xarr = zones_xarr, 
                                 img_path = img_path, 
                                 ZONE_ID_COL = ZONE_ID_COL, 
                                 BUFFER_METERS = BUFFER_METERS, 
                                 band_name = band_name)
        gdf = gdf.merge(out, on=ZONE_ID_COL, how="left")

    except Exception as e:
        print("\n❌ Error at index:", i)
        print("Band name:", band_name)
        print("File:", img_path)
        print("Error:", e, "\n")

        error_band = band_name
        error_idx = i
        break


