import os
import time

import fiona

from pathlib import Path
from datetime import datetime
from calendar import monthrange
import multiprocessing as mp

start_time = time.time()

# Input tiles
inputFc = '/projects/my-private-bucket/HLS/data/shp/atlantic_forest/tiles/br_af_grid90km_prj.gpkg'
# inputFc = '/projects/HLS/data/shp/atlantic_forest/tiles/ls_unit50km_utm23s.gpkg'

# Output directory
outdir = '/projects/my-private-bucket/HLS/data/tif/monthly_composites_ndvi_test_terminal/'


# Functions 
def wrapper_composite(params):

    FOCAL_TILE = params.get('FOCAL_TILE')
    SAT_API = params.get('SAT_API')
    MS_COMP_TYPE = params.get('MS_COMP_TYPE')
    YEAR = params.get('YEAR')
    MIN_N_FILT_RESULTS = params.get('MIN_N_FILT_RESULTS')
    SEASON_START = params.get('SEASON_START')
    SEASON_STOP = params.get('SEASON_STOP')
    # INDEX_FN = params.get('INDEX_FN')
    # INDEX_LYR = params.get('INDEX_LYR')
    STAT = params.get('STAT')
    STAT_PCT = params.get('STAT_PCT')
    TARGET_SPECTRAL = params.get('TARGET_SPECTRAL')

    INDEX_FN =  params.get('INDEX_FN') #'https://maap-ops-workspace.s3.amazonaws.com/shared/montesano/databank/boreal_tiles_v004.gpkg'
    INDEX_LYR = params.get('INDEX_LYR') # 'boreal_tiles_v004'
    YEAR_START, YEAR_STOP = (YEAR, YEAR)
    HLS_PRODUCT = params.get('HLS_PRODUCT') #HLS_PRODUCT = 'H30'
    MAX_CLOUDS = params.get('MAX_CLOUDS') #MAX_CLOUDS = 0

    OUTDIR = params.get('OUTDIR') #'/projects/my-private-bucket/tmp/mask_test_keep_snow'
    
    args = f"--in_tile_fn {INDEX_FN} \
        --in_tile_layer {INDEX_LYR} \
        --sat_api {SAT_API} \
        --tile_buffer_m 0 \
        --in_tile_num {FOCAL_TILE} \
        --output_dir {OUTDIR} \
        -sy {YEAR_START} -ey {YEAR_STOP} -smd {SEASON_START} -emd {SEASON_STOP} -mc {MAX_CLOUDS} \
        --composite_type {MS_COMP_TYPE} \
        --hls_product {HLS_PRODUCT} \
        --thresh_min_ndvi -1 \
        --min_n_filt_results {MIN_N_FILT_RESULTS} \
        --stat {STAT} \
        --stat_pct {STAT_PCT} \
        --target_spectral_index {TARGET_SPECTRAL}"
    args += " --do_indices"
    #args += " --search_only"
    #args += " --rangelims_red 0.01 0.1" # the default now effectively no limit [-1e9, 1e9]
    
    cmd = f'python /projects/my-private-bucket/code/icesat2_boreal/lib/build_ms_composite.py {args}'
    #!echo $cmd
    # get_ipython().system('eval $cmd')
    os.system(cmd)

    fn = f'{OUTDIR}/{MS_COMP_TYPE}_{FOCAL_TILE}_{SEASON_START}_{SEASON_STOP}_{YEAR_START}_{YEAR_STOP}_{STAT}{TARGET_SPECTRAL}.tif'
    if STAT == 'percentile':
        fn = f'{OUTDIR}/{MS_COMP_TYPE}_{FOCAL_TILE}_{SEASON_START}_{SEASON_STOP}_{YEAR_START}_{YEAR_STOP}_{STAT}{STAT_PCT}{TARGET_SPECTRAL}.tif'
    #rescaled_multiband_fn = os.path.join(os.path.dirname(fn), os.path.basename(fn).replace('.tif','_rescaled_3band_temp.tif'))
    # plotlib.rescale_multiband_for_plot(fn, rescaled_multiband_fn, bandlist = [5,7,3], pct=[20,90], nodata=-9999.0) 

    return fn

def main(args):

    # ---------------------------------------------
    # Set default parameters 
    # ---------------------------------------------
    SAT_API = 'https://cmr.earthdata.nasa.gov/stac/LPCLOUD'
    MS_COMP_TYPE = 'HLS'
    HLS_PRODUCT = 'H30'
    
    STAT = 'max'
    # TARGET_SPECTRAL = 'ndvi'
    # STAT_PCT = 50.0

    # STAT = 'percentile'
    TARGET_SPECTRAL = 'ndvi'
    STAT_PCT = 50.0

    
    MIN_N_FILT_RESULTS = 5
    MAX_CLOUDS = 0
    
    # ---------------------------------------------
    # Define input fc
    # ---------------------------------------------
    INDEX_FN = args['inputFc'] #'/projects/HLS/data/shp/atlantic_forest/tiles/ls_unit50km_utm23s.gpkg'

    # Get first layer name
    layer_names = fiona.listlayers(INDEX_FN)
    
    # Get the first layer name
    first_layer_name = layer_names[0]
    
    INDEX_LYR = first_layer_name# 'ls_unit50km_utm23s'
    BASE_OUTDIR = args['outdir'] #'/projects/HLS/data/tif/monthly_composites'
    
    
    # In[15]:
    
    
    # Define list of tiles and years
    #TODO define argument in a way that if it is defined as all it process all the tiles in the gpkg
    with fiona.open(INDEX_FN, layer=0) as src:
        tiles = [feature["properties"]["tile_num"] for feature in src]


    # tiles = [5]
    years = [2022]
    
    
    
    # Create dictionary of params
    params = {
        'SAT_API': SAT_API,
        'HLS_PRODUCT': HLS_PRODUCT,
        'MS_COMP_TYPE': MS_COMP_TYPE,
        'MAX_CLOUDS': MAX_CLOUDS,
        'MIN_N_FILT_RESULTS': MIN_N_FILT_RESULTS,
        'STAT': STAT,
        'STAT_PCT': STAT_PCT,
        'TARGET_SPECTRAL': TARGET_SPECTRAL,
        'INDEX_FN': INDEX_FN,
        'INDEX_LYR': INDEX_LYR,
    }
    
    
    
    manager = mp.Manager()
    lock = manager.Lock()
    # ---------------------------------------------
    # Create parameter list
    # ---------------------------------------------

    # Testing
    tiles_run = tiles[1:5]
    params_list = []
    for tile in tiles_run: 
        for year in years:
            for month in range(2, 4):
            # for month in range(1, 13):
                start_day = f"{month:02d}-01"
                end_day = f"{month:02d}-{monthrange(year, month)[1]:02d}"
        
                # Compose output directory for this specific run
                outdir = f'{BASE_OUTDIR}/tile_{tile:03d}/{year}/{month:02d}'
                os.makedirs(outdir, exist_ok=True)
        
                run_params = params.copy()
                run_params.update({
                    'FOCAL_TILE': tile,
                    'YEAR': year,
                    'SEASON_START': start_day,
                    'SEASON_STOP': end_day,
                    'OUTDIR': str(outdir)
                })
        
                params_list.append(run_params)

    
    # ---------------------------------------------
    # Run wrapper composite
    # ---------------------------------------------
   
    
    # import pdb; pdb.set_trace()
    
    with mp.Pool(processes= mp.cpu_count() - 1) as pool:
        fn_list = pool.map(wrapper_composite, params_list)


if __name__ == "__main__":
    args = {
        'inputFc': inputFc,
        'outdir': outdir
    }
    
    main(args)

    end_time = time.time()
    elapsed_time = end_time - start_time
     # Convert elapsed time to minutes and hours
    elapsed_minutes = elapsed_time / 60
    elapsed_hours = elapsed_minutes / 60

    # print(f"Processing completed in {elapsed_time:.2f} seconds.")
    print(f"Processing completed in {elapsed_minutes:.2f} minutes.")
    print(f"Processing completed in {elapsed_hours:.2f} hours.")


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Run HLS monthly composites")
#     parser.add_argument('--inputFc', type=str, required=True, help="Path to input GeoPackage with tiles")
#     parser.add_argument('--outdir', type=str, required=True, help="Directory to save composite outputs")

#     args = parser.parse_args()

#     main({'inputFc': args.inputFc, 'outdir': args.outdir})

    
# nohup python /projects/my-private-bucket/HLS/IGARSS25/ms_composite_month.py > /projects/my-private-bucket/HLS/ms_composite_2022test_log.txt 2>&1 &
