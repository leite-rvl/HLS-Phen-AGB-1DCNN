import os
import time
import fiona
from pathlib import Path
from datetime import datetime
from calendar import monthrange
import multiprocessing as mp

import copy
import json

# Input tiles
inputFc = '/projects/my-private-bucket/HLS/data/shp/atlantic_forest/tiles/br_af_grid30km_prj.gpkg'
# inputFc = '/projects/my-private-bucket/HLS/data/shp/atlantic_forest/tiles/br_af_grid90km_prj.gpkg'
# inputFc = '/projects/HLS/data/shp/atlantic_forest/tiles/ls_unit50km_utm23s.gpkg'

# Output directory
# outdir = '/projects/my-private-bucket/HLS/data/tif/monthly_composites_ndvi_test_terminal/'
outdir = '/projects/data/tmp/monthly_composites_ndvi_test_terminal2'

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

     # Define output filename
    fn = f'{OUTDIR}/{MS_COMP_TYPE}_{FOCAL_TILE}_{SEASON_START}_{SEASON_STOP}_{YEAR_START}_{YEAR_STOP}_{STAT}{TARGET_SPECTRAL}.tif'
    if STAT == 'percentile':
        fn = f'{OUTDIR}/{MS_COMP_TYPE}_{FOCAL_TILE}_{SEASON_START}_{SEASON_STOP}_{YEAR_START}_{YEAR_STOP}_{STAT}{STAT_PCT}{TARGET_SPECTRAL}.tif'

    # Skip if output already exists
    if os.path.exists(fn):
        print(f"[{FOCAL_TILE}] Output exists, skipping: {fn}")
        return fn
        
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

    
    # cmd = f'python /projects/my-private-bucket/code/icesat2_boreal/lib/build_ms_composite.py {args}' 
    cmd = f'python /projects/my-private-bucket/code/icesat2_boreal/lib/build_ms_composite_multip.py {args}'
    #!echo $cmd
    # get_ipython().system('eval $cmd')
    os.system(cmd)

    # fn = f'{OUTDIR}/{MS_COMP_TYPE}_{FOCAL_TILE}_{SEASON_START}_{SEASON_STOP}_{YEAR_START}_{YEAR_STOP}_{STAT}{TARGET_SPECTRAL}.tif'
    # if STAT == 'percentile':
    #     fn = f'{OUTDIR}/{MS_COMP_TYPE}_{FOCAL_TILE}_{SEASON_START}_{SEASON_STOP}_{YEAR_START}_{YEAR_STOP}_{STAT}{STAT_PCT}{TARGET_SPECTRAL}.tif'
    #rescaled_multiband_fn = os.path.join(os.path.dirname(fn), os.path.basename(fn).replace('.tif','_rescaled_3band_temp.tif'))
    # plotlib.rescale_multiband_for_plot(fn, rescaled_multiband_fn, bandlist = [5,7,3], pct=[20,90], nodata=-9999.0) 
    if not os.path.exists(fn):
        raise RuntimeError(f"Expected output not found: {fn}")
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
    
    
    
    # Define list of tiles and years
    #TODO define argument in a way that if it is defined as all it process all the tiles in the gpkg
    with fiona.open(INDEX_FN, layer=0) as src:
        tiles = [feature["properties"]["tile_num"] for feature in src]

    
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
    
    
    
    
    # ---------------------------------------------
    # Create parameter list
    # ---------------------------------------------

    # Testing
    # tiles_run = tiles[0:1]
    # Multiprocess 2022, months range(2,5)
    ## Processing completed in 5.07 minutes.
    ## Processing completed in 0.08 hours.
    ## Error in month 02 (output not exported)
    
    tiles_run = tiles[1:2]
    
   
    
    params_list = []
    for tile in tiles_run: 
        for year in years:
            for month in range(2, 5):
            # for month in range(1, 13):
                start_day = f"{month:02d}-01"
                end_day = f"{month:02d}-{monthrange(year, month)[1]:02d}"
        
                # Compose output directory for this specific run
                outdir = f'{BASE_OUTDIR}/tile_{tile:03d}/{year}/{month:02d}'
                os.makedirs(outdir, exist_ok=True)
        
                # run_params = params.copy()
                # run_params = copy.deepcopy(params)
                # run_params.update({
                    # 'FOCAL_TILE': tile,
                    # 'YEAR': year,
                    # 'SEASON_START': start_day,
                    # 'SEASON_STOP': end_day,
                    # 'OUTDIR': str(outdir)
                # })
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

                    
                    'FOCAL_TILE': tile,
                    'YEAR': year,
                    'SEASON_START': start_day,
                    'SEASON_STOP': end_day,
                    'OUTDIR': str(outdir)
                    }
        
                params_list.append(params)


    #TEST
    # import pdb; pdb.set_trace()
    # for params in params_list:
    #     output_file = wrapper_composite(params)
    params = params_list[0]
    fn = wrapper_composite(params)
    # ---------------------------------------------
    # Run wrapper composite
    # ---------------------------------------------

    def safe_wrapper(params):
        try:
            return wrapper_composite(params)
        except Exception as e:
            print(f"[ERROR] Failed to process tile {params.get('FOCAL_TILE')} "
                  f"season {params.get('SEASON_START')}–{params.get('SEASON_STOP')} "
                  f"year {params.get('YEAR')}\n        Reason: {e}")
            return {'failed': True, 'params': params, 'error': str(e)}

    # You can adjust the number of processes as needed
    num_procs = max(mp.cpu_count() - 5, 1)
    
    results = []
    with Pool(processes=num_procs) as pool:
        results = pool.map(safe_wrapper, params_list)

# Save only failed cases
failed_params = [r['params'] for r in results if isinstance(r, dict) and r.get('failed')]

if failed_params:
    failed_file = f'{BASE_OUTDIR}/failed_params.json'
    with open(failed_file, "w") as f:
        json.dump(failed_params, f, indent=4)
    print(f"\n{len(failed_params)} composites failed. Saved details to: {failed_file}")
    
    failed_params = []

    for params in params_list:
        try:
            output_file = wrapper_composite(params)
        except Exception as e:
            print(f"[ERROR] Failed to process tile {params.get('FOCAL_TILE')} season {params.get('SEASON_START')}–{params.get('SEASON_STOP')} year {params.get('YEAR')}")
            print(f"        Reason: {str(e)}")
            failed_params.append(params)

    # Save failed params to file if any
    if failed_params:
        failed_file = f'{BASE_OUTDIR}/{tile:03d}/failed_params_{tile:03d}.json'
        with open(failed_file, "w") as f:
            json.dump(failed_params, f, indent=4)
        print(f"\n {len(failed_params)} composites failed. Saved details to: {failed_file}")
    
    # manager = mp.Manager()
    # lock = manager.Lock()
    # fn_outputs = []
    # with mp.Pool(processes = 2) as pool:
    # # with mp.Pool(processes = mp.cpu_count() - 10) as pool:
    #     for params in params_list:
    #         fn = pool.apply_async(wrapper_composite,args=(params,))
    #         fn_outputs.append(fn)
    #     pool.close()
    #     pool.join()
         
    # with mp.Pool(processes=mp.cpu_count() - 5) as pool:
    #     fn_list = pool.map(wrapper_composite, params_list)

if __name__ == "__main__":
    
    start_time = time.time()

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

    
# nohup python /projects/my-private-bucket/HLS/IGARSS25/ms_composite_month_multip.py > /projects/my-private-bucket/HLS/ms_composite_2022test_log.txt 2>&1 &

# rerun
# with open("failed_params.json") as f:
#     failed_params = json.load(f)
#     for p in failed_params:
#         wrapper_composite(p)
