# TEMPLATE
# gdalbuildvrt -allow_projection_difference -vrtnodata -9999 /explore/nobackup/projects/hls/HLS/mwooten3/Composites/HLS.30m/_vrt/01/HLS.M30.2020.2022.01.v2.0.Fmask.vrt /explore/nobackup/projects/hls/HLS/mwooten3/Composites/HLS.30m/_vrt/2020001.2020031/HLS.M30.2020001.2020031.v2.0.Fmask.vrt /explore/nobackup/projects/hls/HLS/mwooten3/Composites/HLS.30m/_vrt/2021001.2021031/HLS.M30.2021001.2021031.v2.0.Fmask.vrt /explore/nobackup/projects/hls/HLS/mwooten3/Composites/HLS.30m/_vrt/2022001.2022031/HLS.M30.2022001.2022031.v2.0.Fmask.vrt


gdalbuildvrt \
  -srcnodata 0 -vrtnodata 0 \
  /projects/my-private-bucket/HLS-1DCNN-AGB/data/tif/VegetationMask/AtlanticForest/2024/vegmask_2024.vrt \
  /projects/my-private-bucket/HLS-1DCNN-AGB/data/tif/VegetationMask/AtlanticForest/2024/Tiles/*.tif



gdalbuildvrt \
  -srcnodata 0 -vrtnodata 0 \
  /projects/my-private-bucket/HLS-1DCNN-AGB/data/tif/VegetationMask/AtlanticForest/2022/vegmask_2022.vrt \
  /projects/my-private-bucket/HLS-1DCNN-AGB/data/tif/VegetationMask/AtlanticForest/2022/Tiles/*.tif

gdalbuildvrt \
  -srcnodata 0 -vrtnodata 0 \
  /projects/my-private-bucket/HLS-1DCNN-AGB/data/tif/VegetationMask/AtlanticForest/2020/vegmask_2020.vrt \
  /projects/my-private-bucket/HLS-1DCNN-AGB/data/tif/VegetationMask/AtlanticForest/2020/Tiles/*.tif


