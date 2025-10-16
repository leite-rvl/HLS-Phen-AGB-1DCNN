#!/usr/bin/env bash
set -euo pipefail

# Get current location of build script
# basedir=$(dirname "$(readlink -f "$0")")

# Make default dir the .sh command dir
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# # Create env once per machine; safe to call repeatedly
# ENV_PREFIX="/projects/env/maap-zs"
# if [ ! -d "$ENV_PREFIX" ]; then
#   conda env create -f environment.yml --prefix "$ENV_PREFIX"
# fi

# source activate "$ENV_PREFIX" || conda activate "$ENV_PREFIX"


# TESTING
GPKG="https://maap-ops-workspace.s3.amazonaws.com/rodrigo.leite/HLS-1DCNN-AGB/data/shp/gedi/test/l4a_t90km_t89_veg2022_outrm.gpkg" 
DIR_IMG="/projects/my-private-bucket/HLS-1DCNN-AGB/data/tif/HLS_composites/monthly/br_af_grid60km_prj_evi2_max/vrt_test/" 
OUTPUT="/projects/my-private-bucket/HLS-1DCNN-AGB/data/shp/gedi/l4a_t90km_t89_veg2022_outrm_zonal_HLS.gpkg" 
ZONE_ID_COL="zone" 
BUFFER_M=0 
EXCLUDE_BANDS="ValidMask.vrt count.vrt yearDate.vrt JulianDate.vrt" 

# With position arguments
# GPKG=$1
# DIR_IMG=$2
# OUTPUT=$3


conda run --live-stream --name python python zonal_stats_vrt.py \
  --gpkg       "${GPKG}" \
  --dir-img    "${DIR_IMG}" \
  --output     "${OUTPUT}" \
  --zone-id-col "${ZONE_ID_COL:-zone}" \
  --buffer-m    "${BUFFER_M:-0}" \
  --exclude-bands ${EXCLUDE_BANDS:-ValidMask.vrt count.vrt}


# python zonal_stats_vrt.py \
#   --gpkg       "${GPKG}" \
#   --dir-img    "${DIR_IMG}" \
#   --output     "${OUTPUT}" \
#   --zone-id-col "${ZONE_ID_COL:-zone}" \
#   --buffer-m    "${BUFFER_M:-0}" \
#   --exclude-bands ${EXCLUDE_BANDS:-ValidMask.vrt count.vrt}
