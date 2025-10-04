#!/usr/bin/env python
# coding: utf-8

# # Search for HLS data using `pystac_client`

# ### For debugging HLS querying from LPDAAC STAC

# In[1]:


pip install cogeo_mosaic


# In[10]:


import os
import sys
sys.path.append('/projects/code/icesat2_boreal/lib')
import plotlib
from mosaiclib import *


# In[11]:


import rasterio
from rasterio.plot import show_hist, show
import numpy as np
import matplotlib.pyplot as plt


# ## Tests
# 
# ### NOTE:
# if `cogeo_mosaic` has been installed, it bring with it the newest `rio_tiler` which will break this script.  
# You want `rio_tiler==6.6.1`

# In[12]:


import CovariateUtils


# In[13]:


import geopandas as gpd
from mosaiclib import *


# In[14]:


import rio_tiler


# In[15]:


rio_tiler.__version__


# In[16]:


boreal_tiles = gpd.read_file('https://maap-ops-workspace.s3.amazonaws.com/shared/montesano/databank/boreal_tiles_v004.gpkg')
#boreal_tiles.explore()


# ## Wrap `build_ms_composite.py` for testing

# In[17]:


# Landsat Collection 2 SR
SAT_API = 'https://landsatlook.usgs.gov/stac-server'
MS_COMP_TYPE = 'LC2SR' 

# HLS
SAT_API = 'https://cmr.earthdata.nasa.gov/stac/LPCLOUD'
MS_COMP_TYPE = 'HLS'

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
    #args += " --do_indices"
    #args += " --search_only"
    #args += " --rangelims_red 0.01 0.1" # the default now effectively no limit [-1e9, 1e9]
    
    cmd = f'python /projects/code/icesat2_boreal/lib/build_ms_composite.py {args}'
    #!echo $cmd
    get_ipython().system('eval $cmd')

    fn = f'{OUTDIR}/{MS_COMP_TYPE}_{FOCAL_TILE}_{SEASON_START}_{SEASON_STOP}_{YEAR_START}_{YEAR_STOP}_{STAT}{TARGET_SPECTRAL}.tif'
    if STAT == 'percentile':
        fn = f'{OUTDIR}/{MS_COMP_TYPE}_{FOCAL_TILE}_{SEASON_START}_{SEASON_STOP}_{YEAR_START}_{YEAR_STOP}_{STAT}{STAT_PCT}{TARGET_SPECTRAL}.tif'
    #rescaled_multiband_fn = os.path.join(os.path.dirname(fn), os.path.basename(fn).replace('.tif','_rescaled_3band_temp.tif'))
    # plotlib.rescale_multiband_for_plot(fn, rescaled_multiband_fn, bandlist = [5,7,3], pct=[20,90], nodata=-9999.0) 

    return fn


# ## Test Build MS Composite (new version of 3.1.2)
# 
# TODO:  
# + Fix 0 values coming out of comp when a pixel slice is all NaN. Need this to be NaN.
# + adjust red_rangelims : change for winter comps; are they even needed in any comps?
# + adjust min ndvi thresh
# + test for tiles with ocean
# + test for Landsat 5/7 comps
# 
# ### Learned details:  
# + put flag that would turn on snow/ice masking for summer comps - this way, that nodata could be masked all pixels down through the stack across all seasons. But dont do this for winter, fall seasons.

# ## Current testing: `build_ms_composite.py`
# 

# ### Edited the wrapper function and the arguments list to edit everything in the dictionary

# In[36]:


STAT = 'max'
TARGET_SPECTRAL = 'ndvi'
MIN_N = 10


params_test1 = {
    
    'SAT_API': 'https://cmr.earthdata.nasa.gov/stac/LPCLOUD',
    'HLS_PRODUCT': 'H30',
    'MS_COMP_TYPE': 'HLS', # Composite type?
    
    'MAX_CLOUDS': 0,
    'MIN_N_FILT_RESULTS': 10, # Minimum number of filtered scenes required to build the composite
    
    'YEAR': 2024,
    'SEASON_START': '07-01', # Start date for the time window for the composite
    'SEASON_STOP': '08-31',
    
    'STAT': 'max', 
    'STAT_PCT': 50.0,
    'TARGET_SPECTRAL': 'ndvi',
     
    'INDEX_FN': 'https://maap-ops-workspace.s3.amazonaws.com/shared/montesano/databank/boreal_tiles_v004.gpkg', # Path to the GeoPackage tile grid
    'INDEX_LYR': 'boreal_tiles_v004', # Layer name inside the gpkg that defines tile geometry
    'FOCAL_TILE': 3366, # ID of the tile to process (index),
    
    'OUTDIR': f'/projects/HLS/IGARSS25/data/tif/tmp/test3_{STAT}{TARGET_SPECTRAL}_anyseason_n{MIN_N:02}'
     
}



# In[37]:


fn = wrapper_composite(params_test1)


# In[34]:


fn


# ## Plot

# In[35]:


rescaled_multiband_fn = fn.replace('.tif', '_rescaled_3band_temp.tif')

plotlib.rescale_multiband_for_plot(fn, rescaled_multiband_fn, bandlist = [6,4,3], pct=[20,98], nodata=-9999.0) 



fig,ax=plt.subplots(figsize=(10,10))
with rasterio.open(rescaled_multiband_fn) as src:
    #print(src.profile)
    show(src.read(),transform=src.transform, ax=ax, title=os.path.basename(rescaled_multiband_fn))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ### Test with tiles in Atlantic Forest Brazil

# In[40]:


STAT = 'max'
TARGET_SPECTRAL = 'ndvi'
MIN_N = 10


params_test_af = {
    
    'SAT_API': 'https://cmr.earthdata.nasa.gov/stac/LPCLOUD',
    'HLS_PRODUCT': 'H30',
    'MS_COMP_TYPE': 'HLS', # Composite type?
    
    'MAX_CLOUDS': 20,
    'MIN_N_FILT_RESULTS': MIN_N, # Minimum number of filtered scenes required to build the composite
    
    'YEAR': 2024,
    'SEASON_START': '05-01', # Start date for the time window for the composite
    'SEASON_STOP': '08-31',
    
    'STAT': STAT, 
    'STAT_PCT': 50.0,
    'TARGET_SPECTRAL': TARGET_SPECTRAL,
     
    'INDEX_FN': '/projects/HLS/data/shp/atlantic_forest/tiles/ls_unit50km_utm23s.gpkg', # Path to the GeoPackage tile grid
    'INDEX_LYR': 'ls_unit50km_utm23s', # Layer name inside the gpkg that defines tile geometry
    'FOCAL_TILE': 5, # ID of the tile to process (index),
    
    'OUTDIR': f'/projects/HLS/IGARSS25/data/tif/tmp/test1_af_{STAT}{TARGET_SPECTRAL}_anyseason_n{MIN_N:02}'
     
}


# In[41]:


fn = wrapper_composite(params_test_af)


# In[42]:


fn


# In[43]:


rescaled_multiband_fn = fn.replace('.tif', '_rescaled_3band_temp.tif')

plotlib.rescale_multiband_for_plot(fn, rescaled_multiband_fn, bandlist = [6,4,3], pct=[20,98], nodata=-9999.0) 



fig,ax=plt.subplots(figsize=(10,10))
with rasterio.open(rescaled_multiband_fn) as src:
    #print(src.profile)
    show(src.read(),transform=src.transform, ax=ax, title=os.path.basename(rescaled_multiband_fn))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ## TO EDIT

# In[ ]:





# In[55]:


MIN_N = 10 # minimum number of images

#OUTDIR = f'/projects/my-private-bucket/tmp/test_percentile_n{MIN_N:02}'
#OUTDIR = f'/projects/my-private-bucket/tmp/test_max_anyseason_n{MIN_N:02}'

args_dict_septoct = {
    'SAT_API': 'https://cmr.earthdata.nasa.gov/stac/LPCLOUD', 
    'MS_COMP_TYPE': 'HLS',
    'YEAR': 2024, 
    'MIN_N_FILT_RESULTS': MIN_N,
    'SEASON_START': '09-01',
    'SEASON_STOP': '10-31',
    'OUTDIR': '',
    'TARGET_SPECTRAL': ''
}
args_dict_julyaug = {
    'SAT_API': 'https://cmr.earthdata.nasa.gov/stac/LPCLOUD', 
    'MS_COMP_TYPE': 'HLS',
    'YEAR': 2024, 
    'MIN_N_FILT_RESULTS': MIN_N, ###< ----------------------
    'SEASON_START': '07-01',
    'SEASON_STOP': '08-31',
    'OUTDIR': '',
    'TARGET_SPECTRAL': ''
}
args_dict_marchapril = {
    'SAT_API': 'https://cmr.earthdata.nasa.gov/stac/LPCLOUD', 
    'MS_COMP_TYPE': 'HLS',
    'YEAR': 2024, 
    'MIN_N_FILT_RESULTS': MIN_N, ###< ----------------------
    'SEASON_START': '03-01',
    'SEASON_STOP': '04-30',
    'OUTDIR': '',
    'TARGET_SPECTRAL': ''
}
#print(OUTDIR)


# In[57]:


INDEX_FN =  'https://maap-ops-workspace.s3.amazonaws.com/shared/montesano/databank/boreal_tiles_v004.gpkg'
INDEX_LYR = 'boreal_tiles_v004'


YEAR_LIST = [2019] #[2016,2022,2025]
FOCAL_TILE_LIST = [3366] # [3273] #[41248] #[936] #,3266,,3458,3457,37743]
STAT = 'max'
TARGET_SPECTRAL = 'ndvi'
# OUTDIR = f'/projects/my-private-bucket/tmp/test_{STAT}{TARGET_SPECTRAL}_anyseason_n{MIN_N:02}'
OUTDIR = f'/projects/HLS/IGARSS25/data/tif/tmp/test4_{STAT}{TARGET_SPECTRAL}_anyseason_n{MIN_N:02}'


# In[66]:


# Build a list of parameter dictionaries (`params_list`) for all combinations of:
# - years in YEAR_LIST
# - focal tiles in FOCAL_TILE_LIST
# - predefined seasonal periods (March–April, July–August, September–October)

params_list = []
for YEAR in YEAR_LIST:
    for FOCAL_TILE in FOCAL_TILE_LIST:
        params_d_list = []
        # params_d_list.append(args_dict_marchapril.copy())
        params_d_list.append(args_dict_julyaug.copy())
        # params_d_list.append(args_dict_septoct.copy())

        for params in params_d_list:
            
            params['OUTDIR'] = OUTDIR
            params['TARGET_SPECTRAL'] = TARGET_SPECTRAL
            params['FOCAL_TILE'] = FOCAL_TILE
            params['YEAR'] = YEAR
            params['STAT'] = STAT
            params['STAT_PCT'] = 50.0
            
            params_list += [params]   


# In[67]:


len(params_list)


# In[68]:


params_list


# In[61]:


from multiprocessing import Pool
from functools import partial


# #### Multiprocess...

# In[62]:


get_ipython().run_cell_magic('time', '', 'with Pool(processes=3) as pool:\n    fn_list = pool.map(wrapper_composite, params_list)\n')


# #### Single process....

# In[18]:


# params = args_dict_marchapril.copy()

# params['FOCAL_TILE'] = 3266 #41248 # in 2021 this comp is bad - see if still bad #3458 #3457 #37743 #936
# params['YEAR'] = 2022
# params['STAT'] = 'max'
# params['STAT_PCT'] = 50.0


# In[19]:


# %%time

# fn = wrapper_composite(params)


# In[64]:


fn = fn_list[0] #= '/projects/my-private-bucket/tmp/mask_test_keep_snow/HLS_3458_03-01_04-30_2024_2024_percentile50.0.tif'
fn


# In[72]:


rescaled_multiband_fn = fn.replace('.tif', '_rescaled_3band_temp.tif')

# plotlib.rescale_multiband_for_plot(fn, rescaled_multiband_fn, bandlist = [6,4,3], pct=[20,98], nodata=-9999.0) 
plotlib.rescale_multiband_for_plot(fn, rescaled_multiband_fn, bandlist = [6,4,3], pct=[20,98], nodata=-9999.0) 



fig,ax=plt.subplots(figsize=(10,10))
with rasterio.open(rescaled_multiband_fn) as src:
    #print(src.profile)
    show(src.read(),transform=src.transform, ax=ax, title=os.path.basename(rescaled_multiband_fn))


# In[29]:


import rasterio

with rasterio.open(fn) as src:
    print("number of bands:", src.count)
    print("Band descriptions:", src.descriptions)

