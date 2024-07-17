#!/bin/python

import numpy as np 
import pandas as pd
import xarray as xr
import glob
from parallax_correction_funcs import haversine, parallax_correction, era_Tlapse_height, call_parallax_era
import sys

####### satellite position
lat_sat = 0
lon_sat = 41.5

####### prescribe time controls
startyear = int(sys.argv[1])
endyear = int(sys.argv[2])
#mon_arr = np.array([int(sys.argv[3])])
mon_arr = np.array([6,7,8,9])
#mon_arr = np.array([1,2,3,4,5,6,7,8,9,10,11,12])
mon_dict = {1:'jan',2:'feb',3:'mar',4:'apr',5:'may',6:'jun',7:'jul',8:'aug',9:'sep',10:'oct',11:'nov',12:'dec'}

####### relevant directories
df_path = '/data/keeling/a/dchug2/f/projects_final/nwindia_SM_convection/tracking/monthly_clouds/'
save_path = '/data/keeling/a/dchug2/f/projects_final/nwindia_SM_convection/tracking/monthly_clouds/llax_correction/'

####### loop over years
year = startyear
while year<=endyear:
    
    # loop over array of desired months
    for mon in mon_arr:
        
        # load relevant csv
        df = pd.read_csv(df_path + str(year) + '_' + mon_dict[mon] + '_core241_thres253_coolrate8.csv')
        
        # extract reqiured data for parallax correction
        t_cloud = df['Tb'].values
        lon_cloud = df['lon'].values
        lat_cloud = df['lat'].values
        
        # calculate correction
        km_coords = np.array([call_parallax_era(mon, year, row.Tb, 
                                                row.lon, row.lat, 
                                                lon_sat, lat_sat) for row in df.itertuples()])
        # update dataframe
        df['lat_plax'] = df['lat'] - km_coords[:,3]
        df['lon_plax'] = df['lon'] - km_coords[:,2]
        
        # save dataframe
        df = df.set_index('time')
        df.to_csv(save_path + str(year) + '_' + mon_dict[mon] + '_core241_thres253_coolrate8.csv')
        del df
        
    year = year + 1
