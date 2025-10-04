#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
import geopandas as gpd
import pandas as pd
import rioxarray
from rasterio import features
from pyproj import CRS
from xrspatial import zonal_stats


def _infer_band_label(raster_path: str) -> str:
    """
    Build a short, readable band label from the filename.
    Falls back to the filename stem if the known pattern isn't present.
    Example:
      HLS_89_01-01_12-30_2022_2022_percentile95.0evi2_Blue.vrt
      -> HLS.2022-01-01.2022-12-30.Blue
    """
    stem = os.path.basename(raster_path).replace(".vrt", "").replace(".tif", "")
    parts = stem.split("_")
    try:
        # ['HLS','89','01-01','12-30','2022','2022','percentile95.0evi2','Blue']
        label = f"{parts[0]}.{parts[4]}-{parts[2]}.{parts[5]}-{parts[3]}.{parts[-1]}"
        return label
    except Exception:
        return stem


def point_zonal_stats(
    gpkg_points_path: str,
    raster_path: str,
    layer: str | None = None,
    zone_id_col: str = "zone",
    buffer_m: float = 0.0,
    stats_funcs: list[str] = ["mean"],
    all_touched: bool = False,
    nodata: float | int | None = None,
    dropna_result: bool = False,
    output_gpkg: str | None = None,
) -> gpd.GeoDataFrame:
    """
    Compute zonal statistics for (optionally buffered) point features against a single-band raster
    using xrspatial.zonal_stats. Returns a GeoDataFrame with the stats columns joined.

    Parameters
    ----------
    gpkg_points_path : str
        Path to input GeoPackage with points (GEDI L4A).
    raster_path : str
        Path to single-band raster (VRT/TIF).
    layer : str | None
        Layer name in the GeoPackage, if needed.
    zone_id_col : str
        Name of unique integer ID column to attach (created if missing).
    buffer_m : float
        Buffer radius in meters (0 = use points as single-pixel zones).
        NOTE: If the raster CRS is geographic (degrees), this will buffer in degree units.
    stats_funcs : list[str]
        Stats to compute; e.g., ["mean","min","max","sum","std","var","count"].
    all_touched : bool
        Pass-through to rasterize(); if True, any pixel touched by polygon is burned.
    nodata : float|int|None
        NoData value for the raster. If None, tries da.rio.nodata.
    dropna_result : bool
        If True, drop rows where all computed stats are NaN.
    output_gpkg : str|None
        If provided, writes the resulting GeoDataFrame to this file.

    Returns
    -------
    gdf_stats : geopandas.GeoDataFrame
        Input features with new columns for each requested stat.
    """
    # ---- Load raster as xarray ----
    da = rioxarray.open_rasterio(raster_path).squeeze()  # -> [y, x] (drop band dim if present)
    if "band" in da.dims:
        da = da.squeeze("band", drop=True)

    raster_crs = CRS.from_user_input(da.rio.crs)
    transform = da.rio.transform()
    out_shape = da.shape
    r_nodata = da.rio.nodata if nodata is None else nodata

    # mask out nodata in the data array so stats ignore those cells
    if r_nodata is not None:
        da = da.where(da != r_nodata)

    # ---- Load points ----
    gdf = gpd.read_file(gpkg_points_path, layer=layer) if layer else gpd.read_file(gpkg_points_path)

    # Ensure a unique integer zone id
    if zone_id_col not in gdf.columns:
        gdf = gdf.reset_index(drop=True)
        gdf[zone_id_col] = np.arange(1, len(gdf) + 1, dtype=np.int32)

    # Reproject to raster CRS
    gdf_proj = gdf.to_crs(raster_crs)

    # Optional buffering (note on units if geographic CRS)
    if buffer_m and buffer_m > 0:
        # If CRS is projected, buffer in meters. If geographic, this is in degrees.
        # (For geographic buffering in meters, reproject to an appropriate projected CRS first.)
        gdf_proj["geometry"] = gdf_proj.geometry.buffer(buffer_m)

    # ---- Rasterize zones ----
    geom = list(zip(gdf_proj.geometry, gdf_proj[zone_id_col]))
    zones_arr = features.rasterize(
        geom,
        out_shape=out_shape,
        transform=transform,
        fill=r_nodata,
        nodata=r_nodata,
        all_touched=all_touched,
        dtype="int32",
    )

    # Align a zones xarray with the raster grid
    zones_xarr = da.copy(deep=False)
    zones_xarr.data = zones_arr

    # ---- Compute stats ----
    zs_df = zonal_stats(
        zones=zones_xarr,
        values=da,
        stats_funcs=stats_funcs,
        nodata_values=r_nodata,
        return_type="pandas.DataFrame",
    )

    # Name columns with band label prefix
    band_label = _infer_band_label(raster_path)
    rename_map = {s: f"{band_label}_{s}" for s in stats_funcs}
    zs_df = zs_df.rename(columns=rename_map)

    # ---- Join back to features ----
    out = (
        gdf[[zone_id_col]]
        .drop_duplicates()
        .merge(zs_df, left_on=zone_id_col, right_index=True, how="left")
    )
    gdf_stats = gdf.merge(out, on=zone_id_col, how="left")

    if dropna_result:
        # Drop rows where all stat columns are NaN
        stat_cols = list(rename_map.values())
        gdf_stats = gdf_stats.dropna(subset=stat_cols, how="all")

    if output_gpkg:
        gdf_stats.to_file(output_gpkg)

    return gdf_stats


# -------------------------
# Example usage:
# -------------------------
if __name__ == "__main__":
    l4a_path = "/projects/my-private-bucket/HLS-1DCNN-AGB/data/shp/gedi/test/l4a_t90km_t89_veg2020_outrm.gpkg"
    img_path = "/projects/my-private-bucket/HLS-1DCNN-AGB/data/tif/HLS_composites/yearly/br_af_grid90km_evi2_p95/tile_089/bands_vrt/HLS_89_01-01_12-30_2022_2022_percentile95.0evi2_Blue.vrt"

    out_gpkg = "/projects/my-private-bucket/HLS-1DCNN-AGB/data/shp/gedi/test/l4a_t90km_t89_veg2020_outrm_zonal.gpkg"

    gdf_result = point_zonal_stats(
        gpkg_points_path=l4a_path,
        raster_path=img_path,
        zone_id_col="zone",
        buffer_m=0.0,                  # set >0 to make per-point buffers (meters if projected)
        stats_funcs=["mean"],          # add more: ["mean","min","max","sum","std","var","count"]
        all_touched=False,
        nodata=None,                   # let it read from raster; or pass a value
        dropna_result=False,
        output_gpkg=out_gpkg
    )

    print(gdf_result.filter(regex="zone|_mean$").head())
