import numpy as np
from netCDF4 import Dataset
import xarray as xr
import pandas as pd
import datetime as DT
import time
from calendar import monthrange
import os
import sys
import logging

# Setup logging to display information about the processing steps
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def filenamegen(year, mon, day, utc):
    """
    Generates filenames for MERG IR data files based on the given date and time.

    Parameters:
    year (int): Year of the data file.
    mon (int): Month of the data file.
    day (int): Day of the data file.
    utc (int): UTC hour of the data file.

    Returns:
    tuple: A tuple containing two filenames corresponding to the current and previous hours.
    """
    # Create datetime objects for current and previous hours
    date = DT.datetime(year, mon, day, utc)
    dateminus = date - DT.timedelta(hours=1)

    # Format datetime objects into string representations for filenames
    date_str = date.strftime('%Y%m%d%H')
    dateminus_str = dateminus.strftime('%Y%m%d%H')

    # Define the directory where the data files are stored
    datadir = '/data/keeling/a/dchug2/a/GPM_MERGIR/nwindia_clipped/'

    # Construct filenames using the formatted datetime strings
    fname = os.path.join(datadir, f'merg_{date_str}_4km-pixel.nc4')
    fname_minus = os.path.join(datadir, f'merg_{dateminus_str}_4km-pixel.nc4')

    return fname, fname_minus

def initialize_variables():
    """
    Initialize variables used for processing.
    Returns a dictionary containing empty arrays and counters.
    """
    return {
        'Tb_field': np.empty((1,11,11)),
        'lat_field': np.empty((1,11)),
        'lon_field': np.empty((1,11)),
        'cenlat': np.array([], dtype='float64'),
        'cenlon': np.array([], dtype='float64'),
        'savetime': np.array([], dtype='datetime64[ns]'),
        'cenTb': np.array([]),
        'counter': 0,
        'fname_minus_exists': 0,
        'fname_exists': 0,
        'isnan_center': 0,
        'valid_retrievals': 0,
        'local_minima': 0,
        'cold_core': 0,
        'vicinity': 0,
        'cooler_30min': 0,
        'cool_rate': 0
    }

def process_file(year, mon, day, hour, constants, variables):
    """
    Processes a single MERG IR data file for a given hour.

    Parameters:
    year (int): Year of the data file.
    mon (int): Month of the data file.
    day (int): Day of the data file.
    hour (int): UTC hour of the data file.
    constants (dict): Dictionary containing constants used for processing.
    variables (tuple): Tuple containing arrays and dictionaries for storing data and counts.
    """
    # Generate filenames for current and previous hours
    fname, fname_minus = filenamegen(year, mon, day, hour)

    try:
        # Load datasets using xarray
        xr_T = xr.load_dataset(fname)
        xr_Tminus = xr.load_dataset(fname_minus)
    except FileNotFoundError:
        logging.warning(f'File not found for year: {year}, month: {mon}, day: {day}, hour: {hour}. Skipping...')
        variables['fname_exists'] += 1
        return

    # Extract latitude and longitude arrays from the datasets
    lat_arr = xr_T.lat.values
    lon_arr = xr_T.lon.values

    # Process two time steps: current and previous hour
    for step in [0, 1]:
        Tb, Tb_minus30, Tb_minus60, timeobj = get_time_step_data(xr_T, xr_Tminus, step)

        # Iterate over latitude and longitude indices within the valid range
        for latidx in range(constants['vicsize'] - 1, len(lat_arr) - (constants['vicsize'] - 1)):
            for lonidx in range(constants['vicsize'] - 1, len(lon_arr) - (constants['vicsize'] - 1)):
                process_grid_point(latidx, lonidx, Tb, Tb_minus30, Tb_minus60, timeobj, lat_arr, lon_arr, constants, variables)


def get_time_step_data(xr_T, xr_Tminus, step):
    """
    Extracts temperature data arrays for a specific time step.

    Parameters:
    xr_T (xarray.Dataset): Dataset containing current time step data.
    xr_Tminus (xarray.Dataset): Dataset containing previous time step data.
    step (int): Step index (0 for current hour, 1 for previous hour).

    Returns:
    tuple: Arrays containing temperature data and time objects for the specified time step.
    """
    if step == 0:
        Tb = xr_T.Tb.isel(time=0).values
        Tb_minus30 = xr_Tminus.Tb.isel(time=1).values
        Tb_minus60 = xr_Tminus.Tb.isel(time=0).values
        timeobj = xr_T.time.isel(time=0).values
    else:
        Tb = xr_T.Tb.isel(time=1).values
        Tb_minus30 = xr_T.Tb.isel(time=0).values
        Tb_minus60 = xr_Tminus.Tb.isel(time=1).values
        timeobj = xr_T.time.isel(time=1).values
    return Tb, Tb_minus30, Tb_minus60, timeobj

def process_grid_point(latidx, lonidx, Tb, Tb_minus30, Tb_minus60, timeobj, lat_arr, lon_arr, constants, variables):
    """
    Processes a single grid point within the MERG IR data.

    Parameters:
    latidx (int): Index of latitude array for the grid point.
    lonidx (int): Index of longitude array for the grid point.
    Tb (np.ndarray): Temperature data array for the current hour.
    Tb_minus30 (np.ndarray): Temperature data array for the previous hour (30 minutes ago).
    Tb_minus60 (np.ndarray): Temperature data array for the previous hour (60 minutes ago).
    timeobj (np.datetime64): Time object corresponding to the current hour.
    lat_arr (np.ndarray): Array of latitude values.
    lon_arr (np.ndarray): Array of longitude values.
    constants (dict): Dictionary containing constants used for processing.
    variables (tuple): Tuple containing arrays and dictionaries for storing data and counts.
    """
    # Extract vicinity and domain arrays centered at the grid point
    vicinity = Tb[latidx - constants['vicrad']:latidx + constants['vicrad'] + 1, 
                  lonidx - constants['vicrad']:lonidx + constants['vicrad'] + 1]
    domain = Tb[latidx - constants['domrad']:latidx + constants['domrad'] + 1, 
                lonidx - constants['domrad']:lonidx + constants['domrad'] + 1]
    center = vicinity[constants['vicrad'], constants['vicrad']]

    # Check conditions for data validity and features of interest
    if np.isnan(center):
        variables['isnan_center'] += 1
        return

    if np.sum(np.isnan(domain)) > np.floor((constants['domsize'] ** 2)/2):
        variables['valid_retrievals'] += 1
        return

    if not is_local_minima(center, vicinity):
        variables['local_minima'] += 1
        return

    if center > constants['coldcore']:
        variables['cold_core'] += 1
        return

    if not is_vicinity_warmer_than_threshold(vicinity, constants['cloud_thres'], constants['vicsize'], constants['domsize']):
        variables['vicinity'] += 1
        return

    # Calculate temperature thresholds for comparison
    ts3, ts2, ts1 = center, np.nanmin(Tb_minus30[latidx - constants['vicrad']:latidx + constants['vicrad'] + 1, lonidx - constants['vicrad']:lonidx + constants['vicrad'] + 1]), np.nanmin(Tb_minus60[latidx - constants['vicrad']:latidx + constants['vicrad'] + 1, lonidx - constants['vicrad']:lonidx + constants['vicrad'] + 1])

    # Check conditions based on temperature thresholds
    if ts3 > ts2:
        variables['cooler_30min'] += 1
        return

    if (ts3 - ts1) > constants['coolrate']:
        variables['cool_rate'] += 1
        return

    # Store data for the grid point if all conditions are met
    append_data(latidx, lonidx, timeobj, lat_arr, lon_arr, ts3, vicinity, constants, variables)

def is_local_minima(center, vicinity):
    """
    Checks if the center temperature is a local minima within its vicinity.

    Parameters:
    center (float): Temperature at the center grid point.
    vicinity (np.ndarray): Temperature values in the vicinity of the center grid point.

    Returns:
    bool: True if the center temperature is a local minima, False otherwise.
    """
    test1 = vicinity.copy() - center
    test1[test1 < 0] = 0
    return np.sum(test1 == 0) == 1

def is_vicinity_warmer_than_threshold(vicinity, cloud_thres, vicsize, domsize):
    """
    Checks if the vicinity temperatures are warmer than a given threshold.

    Parameters:
    vicinity (np.ndarray): Temperature values in the vicinity of the center grid point.
    cloud_thres (float): Temperature threshold for clouds.

    Returns:
    bool: True if the vicinity is warmer than the threshold, False otherwise.
    """
    peri = int((vicsize - 1) / 2 - (domsize - 1) / 2)
    return not (np.any(vicinity[:peri, :] < cloud_thres) or np.any(vicinity[-peri:, :] < cloud_thres) or np.any(vicinity[peri:-peri, :peri] < cloud_thres) or np.any(vicinity[peri:-peri, -peri:] < cloud_thres))

def append_data(latidx, lonidx, timeobj, lat_arr, lon_arr, ts3, vicinity, constants, variables):
    """
    Appends data for a valid grid point to the storage arrays.

    Parameters:
    latidx (int): Index of latitude array for the grid point.
    lonidx (int): Index of longitude array for the grid point.
    timeobj (np.datetime64): Time object corresponding to the current hour.
    lat_arr (np.ndarray): Array of latitude values.
    lon_arr (np.ndarray): Array of longitude values.
    ts3 (float): Temperature at the center grid point for the current hour.
    vicinity (np.ndarray): Temperature values in the vicinity of the center grid point.
    variables (tuple): Tuple containing arrays and dictionaries for storing data and counts.
    """
    if variables['counter'] == 0:
        variables['Tb_field'] = np.expand_dims(vicinity, 
                                               axis=0)
        variables['lat_field'] = np.expand_dims(lat_arr[latidx - constants['vicrad']:latidx + constants['vicrad'] + 1], 
                                                axis=0)
        variables['lon_field'] = np.expand_dims(lon_arr[lonidx - constants['vicrad']:lonidx + constants['vicrad'] + 1], 
                                                axis=0)
        variables['cenlat'] = lat_arr[latidx]
        variables['cenlon'] = lon_arr[lonidx]
        variables['savetime'] = timeobj
        variables['cenTb'] = ts3
    else:
        variables['Tb_field'] = np.append(variables['Tb_field'], 
                                          np.expand_dims(vicinity, axis=0), 
                                          axis=0)
        variables['lat_field'] = np.append(variables['lat_field'], 
                                           np.expand_dims(lat_arr[latidx - constants['vicrad']:latidx + constants['vicrad'] + 1], axis=0), 
                                           axis=0)
        variables['lon_field'] = np.append(variables['lon_field'], 
                                           np.expand_dims(lon_arr[lonidx - constants['vicrad']:lonidx + constants['vicrad'] + 1], axis=0), 
                                           axis=0)
        variables['cenlat'] = np.append(variables['cenlat'], lat_arr[latidx])
        variables['cenlon'] = np.append(variables['cenlon'], lon_arr[lonidx])
        variables['savetime'] = np.append(variables['savetime'], timeobj)
        variables['cenTb'] = np.append(variables['cenTb'], ts3)
    variables['counter'] += 1

def save_csv(year, mon, coldcore, cloud_thres, coolrate, variables, savepath):
    """
    Saves the processed data to a CSV file.

    Parameters:
    year (int): Year of the data.
    mon (int): Month of the data.
    coldcore (float): Temperature threshold for identifying cold cores.
    cloud_thres (float): Temperature threshold for identifying cloud cover.
    coolrate (float): Cooling rate threshold for identifying cooling regions.
    variables (tuple): Tuple containing arrays and dictionaries for storing data and counts.
    savepath (str): Path where the CSV file will be saved.

    Returns:
    str: Filename of the saved CSV file.
    """
    mon_dict = {1:'jan', 2:'feb', 3:'mar', 4:'apr', 5:'may', 6:'jun', 7:'jul', 8:'aug', 9:'sep', 10:'oct', 11:'nov', 12:'dec'}
    result = pd.DataFrame({'lat': np.array(variables['cenlat']), 'lon': np.array(variables['cenlon']), 
                           'Tb': np.array(variables['cenTb']), 'time': np.array(variables['savetime'])})
    result = result.set_index('time')
    save_fname = f"{year}_{mon_dict[mon]}_core{int(coldcore)}_thres{int(cloud_thres)}_coolrate{int(abs(coolrate))}.csv"
    result.to_csv(os.path.join(savepath, save_fname))
    return save_fname

def save_to_netcdf(year, mon, coldcore, cloud_thres, coolrate, variables, savepath):
    """
    Saves the processed data to a NetCDF file.

    Parameters:
    year (int): Year of the data.
    mon (int): Month of the data.
    coldcore (float): Temperature threshold for identifying cold cores.
    cloud_thres (float): Temperature threshold for identifying cloud cover.
    coolrate (float): Cooling rate threshold for identifying cooling regions.
    variables (tuple): Tuple containing arrays and dictionaries for storing data and counts.
    savepath (str): Path where the NetCDF file will be saved.
    """
    mon_dict = {1:'jan', 2:'feb', 3:'mar', 4:'apr', 5:'may', 6:'jun', 7:'jul', 8:'aug', 9:'sep', 10:'oct', 11:'nov', 12:'dec'}
    if variables['cenlat'].size > 0:
        s_n = np.arange(-5, 6)
        w_e = np.arange(-5, 6)
        N = np.arange(1, len(variables['cenlat'])+1)
        save_fname = f"{year}_{mon_dict[mon]}_core{int(coldcore)}_thres{int(cloud_thres)}_coolrate{int(abs(coolrate))}.nc"
        new_ds = xr.Dataset({'Tb': (('N', 's_n', 'w_e'), variables['Tb_field']),
                             'lat': (('N', 's_n'), variables['lat_field']),
                             'lon': (('N', 'w_e'), variables['lon_field']),
                             'cenlat': (('N'), variables['cenlat']),
                             'cenlon': (('N'), variables['cenlon']),
                             'time': (('N'), variables['savetime'])},
                            coords={'N': N, 's_n': s_n, 'w_e': w_e})
        new_ds.to_netcdf(os.path.join(savepath, save_fname))

def main():
    """
    Main function to process MERG IR data files over a range of years and months.
    """
    if len(sys.argv) != 6:
        logging.error('Usage: python3 code.py <start year> <end year> <cold core temp> <cloud thres> <cooling rate>')
        return
    
    # Parse command-line arguments
    # Extract arguments and convert to appropriate types
    year, mon = map(int, sys.argv[1:3])
    coldcore, cloud_thres, coolrate = map(float, sys.argv[3:])

    # declare constants
    constants = {'vicsize': 11, 'domsize': 5, 'valid': 15, 'coldcore': coldcore, 
                 'cloud_thres': cloud_thres, 'coolrate': coolrate}
    constants['vicrad'] = int((constants['vicsize']-1)/2)
    constants['domrad'] = int((constants['domsize']-1)/2)
    savepath = "/data/keeling/a/dchug2/f/projects_final/nwindia_SM_convection/tracking/monthly_clouds/test_github_07162024/"
    
    # Iterate over the specified years and months
    mon_day = monthrange(year, mon)[1]
    variables = initialize_variables()
    for day in range(1, mon_day + 1):
        for hour in range(24):
            process_file(year, mon, day, hour, constants, variables)
    
    # Save the processed data to CSV and NetCDF files
    save_fname = save_csv(year, mon, coldcore, cloud_thres, coolrate, variables, savepath)
    logging.info(f'Saved data to {save_fname}')
    save_to_netcdf(year, mon, coldcore, cloud_thres, coolrate, variables, savepath)
    logging.info(f'Saved NetCDF to {save_fname.replace(".csv", ".nc")}')
    
    logging.info('Finished processing all files.')

if __name__ == "__main__":
    main()