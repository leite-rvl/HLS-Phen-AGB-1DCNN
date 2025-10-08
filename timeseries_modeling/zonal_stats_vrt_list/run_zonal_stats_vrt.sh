#!/usr/bin/env bash
set -euo pipefail

# Create env once per machine; safe to call repeatedly
ENV_PREFIX="/projects/env/maap-zs"
if [ ! -d "$ENV_PREFIX" ]; then
  conda env create -f environment.yml --prefix "$ENV_PREFIX"
fi

source activate "$ENV_PREFIX" || conda activate "$ENV_PREFIX"

python zonal_stats_vrt_list.py \
  --gpkg       "${GPKG}" \
  --dir-img    "${DIR_IMG}" \
  --output     "${OUTPUT}" \
  --zone-id-col "${ZONE_ID_COL:-zone}" \
  --buffer-m    "${BUFFER_M:-0}" \
  --exclude-bands ${EXCLUDE_BANDS:-ValidMask.vrt count.vrt}
