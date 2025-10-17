#!/usr/bin/env bash
set -euo pipefail

# --- USER CONFIG ---
VRT_URL="https://maap-ops-workspace.s3.amazonaws.com/rodrigo.leite/HLS-Phen-AGB-1DCNN/data/tif/HLS_composites/monthly/br_af_grid60km_prj_evi2_max/vrt_test/HLS.2018.01.maxevi2.Blue.vrt"


# --- DEFINE OUTPUT FILE ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_FILE="${SCRIPT_DIR}/vrt_test_result.txt"

# --- TEST GDAL VRT READ ---
echo "🔍 Testing if GDAL can open the VRT file:" | tee "$OUTPUT_FILE"
echo "   $VRT_URL" | tee -a "$OUTPUT_FILE"
echo "-------------------------------------------" | tee -a "$OUTPUT_FILE"

if gdalinfo "$VRT_URL" > "$OUTPUT_FILE.tmp" 2>&1; then
    echo "✅ SUCCESS: GDAL successfully read the VRT." | tee -a "$OUTPUT_FILE"
    echo "Output summary:" | tee -a "$OUTPUT_FILE"
    grep -E "Size is|Driver:|Coordinate System" "$OUTPUT_FILE.tmp" | tee -a "$OUTPUT_FILE" || head -n 10 "$OUTPUT_FILE.tmp" | tee -a "$OUTPUT_FILE"
else
    echo "❌ ERROR: GDAL could not open the VRT." | tee -a "$OUTPUT_FILE"
    echo "Detailed log:" | tee -a "$OUTPUT_FILE"
    cat "$OUTPUT_FILE.tmp" | tee -a "$OUTPUT_FILE"
fi

# # Cleanup temp file
# rm -f "$OUTPUT_FILE.tmp"

echo
echo "📄 Results saved to: $OUTPUT_FILE"


# #!/usr/bin/env bash
# set -euo pipefail

# # Get current location of build script
# # basedir=$(dirname "$(readlink -f "$0")")

# # Make default dir the .sh command dir
# SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# cd "$SCRIPT_DIR"

# # # Create env once per machine; safe to call repeatedly
# # ENV_PREFIX="/projects/env/maap-zs"
# # if [ ! -d "$ENV_PREFIX" ]; then
# #   conda env create -f environment.yml --prefix "$ENV_PREFIX"
# # fi

# # source activate "$ENV_PREFIX" || conda activate "$ENV_PREFIX"


# # TESTING
# GPKG="https://maap-ops-workspace.s3.amazonaws.com/rodrigo.leite/HLS-1DCNN-AGB/data/shp/gedi/test/l4a_t90km_t89_veg2022_outrm.gpkg" 
# DIR_IMG="https://maap-ops-workspace.s3.amazonaws.com/rodrigo.leite/HLS-1DCNN-AGB/data/tif/HLS_composites/monthly/br_af_grid60km_prj_evi2_max/vrt_test/" 
# OUTPUT="/projects/my-private-bucket/HLS-1DCNN-AGB/data/shp/gedi/l4a_t90km_t89_veg2022_outrm_zonal_HLS.gpkg" 
# ZONE_ID_COL="zone" 
# BUFFER_M=0 
# EXCLUDE_BANDS="ValidMask.vrt count.vrt yearDate.vrt JulianDate.vrt" 

# # With position arguments
# # GPKG=$1
# # DIR_IMG=$2
# # OUTPUT=$3


# conda run --live-stream --name python python zonal_stats_vrt.py \
#   --gpkg       "${GPKG}" \
#   --dir-img    "${DIR_IMG}" \
#   --output     "${OUTPUT}" #\
#   # --zone-id-col "${ZONE_ID_COL:-zone}" \
#   # --buffer-m    "${BUFFER_M:-0}" \
#   # --exclude-bands ${EXCLUDE_BANDS:-ValidMask.vrt count.vrt}


# # python zonal_stats_vrt.py \
# #   --gpkg       "${GPKG}" \
# #   --dir-img    "${DIR_IMG}" \
# #   --output     "${OUTPUT}" \
# #   --zone-id-col "${ZONE_ID_COL:-zone}" \
# #   --buffer-m    "${BUFFER_M:-0}" \
# #   --exclude-bands ${EXCLUDE_BANDS:-ValidMask.vrt count.vrt}
