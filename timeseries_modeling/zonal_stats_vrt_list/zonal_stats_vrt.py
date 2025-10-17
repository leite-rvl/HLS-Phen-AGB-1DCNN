#!/usr/bin/env python3
"""
write_gpkg_copy.py
Reads a GeoPackage (.gpkg) and writes it to an output directory.
Usage:
    python write_gpkg_copy.py <input_gpkg_path> <output_directory>
"""

import sys
import os
import geopandas as gpd

def main():
    if len(sys.argv) != 3:
        print("Usage: python write_gpkg_copy.py <input_gpkg_path> <output_directory>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_dir = sys.argv[2]

    if not os.path.exists(input_path):
        print(f"❌ Input file does not exist: {input_path}")
        sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read GeoPackage
    gdf = gpd.read_file(input_path)

    # Define output path
    output_path = os.path.join(output_dir, os.path.basename(input_path))

    # Write GeoPackage
    gdf.to_file(output_path, driver="GPKG")

    print(f"✅ GeoPackage saved to: {output_path}")

if __name__ == "__main__":
    main()

####################################################################################

# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-

# """
# Compute zonal statistics (mean) for GEDI L4A points/polygons against a folder of raster images.
# """

# import os
# import glob
# import argparse
# import numpy as np
# import geopandas as gpd
# import rioxarray
# from rasterio import features
# from xrspatial import zonal_stats
# import pandas as pd
# from pyproj import CRS


# def build_zones_xarr(gpkg_points_path, img_path, ZONE_ID_COL='zone', BUFFER_METERS=0, all_touched=False):
#     """Create a zone-labeled xarray aligned to the raster."""
#     da = rioxarray.open_rasterio(img_path).squeeze()
#     if 'band' in da.dims:
#         da = da.squeeze('band', drop=True)

#     raster_crs = CRS.from_user_input(da.rio.crs)
#     transform = da.rio.transform()
#     out_shape = da.shape
#     nodata = da.rio.nodata

#     gdf = gpd.read_file(gpkg_points_path)
#     if ZONE_ID_COL not in gdf.columns:
#         gdf = gdf.reset_index(drop=True)
#         gdf[ZONE_ID_COL] = np.arange(1, len(gdf) + 1, dtype=np.int32)

#     gdf = gdf.to_crs(raster_crs)

#     if BUFFER_METERS > 0:
#         gdf["geometry"] = gdf.geometry.buffer(BUFFER_METERS)

        
#      # (geom, value) pairs; ensure values are int32
#     shapes = list(zip(gdf.geometry, gdf[ZONE_ID_COL].astype(np.int32)))

#     zones_arr = features.rasterize(
#         shapes,
#         out_shape=out_shape,
#         transform=transform,
#         fill=nodata,
#         dtype="int32",
#         all_touched=all_touched
#     )

#     zones_xarr = da.copy(deep=False)
#     zones_xarr.data = zones_arr
#     return zones_xarr


# def zonal_stats_raster(zones_xarr, img_path, band_name='b', ZONE_ID_COL='zone'):
#     """Compute mean zonal stats for a single raster."""
#     da = rioxarray.open_rasterio(img_path).squeeze()
#     if 'band' in da.dims:
#         da = da.squeeze('band', drop=True)
#     nodata = da.rio.nodata
#     da = da.where(da != nodata)

#     zs_df = zonal_stats(
#         zones=zones_xarr,
#         values=da,
#         stats_funcs=['mean'],
#         nodata_values=nodata#,
#         # return_type='pandas.DataFrame'
#     )

#     zs_df = zs_df.rename(columns={'mean': f'{band_name}_mean'})
#     return zs_df


# def run_zonal_stats(gpkg_points_path, dir_img, output_path, ZONE_ID_COL='zone', BUFFER_METERS=0, exclude_bands=None):
#     """Run zonal stats for all rasters in directory and save merged output."""
#     if exclude_bands is None:
#         exclude_bands = []

#     img_paths = [
#         f for f in glob.glob(os.path.join(dir_img, "*.vrt"))
#         if not any(f.endswith(ex) for ex in exclude_bands)
#     ]

#     gdf = gpd.read_file(gpkg_points_path)
#     if ZONE_ID_COL not in gdf.columns:
#         gdf = gdf.reset_index(drop=True)
#         gdf[ZONE_ID_COL] = np.arange(1, len(gdf) + 1, dtype=np.int32)

#     zones_xarr = build_zones_xarr(gpkg_points_path, img_paths[0], ZONE_ID_COL, BUFFER_METERS)
    
#     for i, img_path in enumerate(img_paths, 1):
#         band_name = os.path.basename(img_path).replace('.vrt', '')
#         print(f"[{i}/{len(img_paths)}] Processing {band_name}...")
#         try:
#             out = zonal_stats_raster(zones_xarr, img_path, band_name, ZONE_ID_COL)
#             # import pdb; pdb.set_trace();
#             gdf = gdf.merge(out, on=ZONE_ID_COL, how='left')
#         except Exception as e:
#             print(f"❌ Error on {band_name}: {e}")

#     gdf.to_file(output_path)
#     print(f"✅ Saved: {output_path}")
#     return gdf

# def main():
#     parser = argparse.ArgumentParser(description="Compute zonal statistics for GEDI L4A points/polygons.")

#     parser.add_argument('--gpkg', required=True, help='Input GeoPackage with GEDI points or polygons.')
#     parser.add_argument('--dir-img', required=True, help='Directory with raster .vrt files.')
#     parser.add_argument('--output', required=True, help='Output GeoPackage to save results.')
#     parser.add_argument('--zone-id-col', default='zone', help='Zone ID column name.')
#     parser.add_argument('--buffer-m', type=float, default=0, help='Buffer distance in meters (for points).')
#     parser.add_argument('--exclude-bands', nargs='*', default=['ValidMask.vrt', 'count.vrt'],
#                         help='Raster suffixes to exclude (space-separated).')

#     args = parser.parse_args()

#     run_zonal_stats(
#         gpkg_points_path=args.gpkg,
#         dir_img=args.dir_img,
#         output_path=args.output,
#         ZONE_ID_COL=args.zone_id_col,
#         BUFFER_METERS=args.buffer_m,
#         exclude_bands=args.exclude_bands
#     )


# if __name__ == "__main__":
#     main()


