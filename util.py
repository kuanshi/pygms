import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

def plot_filter_effect(data_array, filter_min, filter_max, parameter_name, target_value=None):
    """
    Plots a histogram of a given data array with vertical lines for filter bounds and optional target value.

    Parameters:
    -----------
    data_array : array-like
        The dataset to analyze (e.g., magnitudes, Vs30 values, distances).
    filter_min : float
        Lower bound of the filter.
    filter_max : float
        Upper bound of the filter.
    parameter_name : str
        Name of the parameter (used for axis labeling and title).
    target_value : float, optional
        The target value to be shown as a blue vertical line.
    """
    data_array = np.asarray(data_array)
    total_count = len(data_array)
    passed_count = np.sum((data_array >= filter_min) & (data_array <= filter_max))

    plt.figure(figsize=(8, 5))
    plt.hist(data_array, bins=50, color='skyblue', edgecolor='black', alpha=0.7)

    plt.axvline(filter_min, color='red', linestyle='--', label=f'Min {parameter_name} = {filter_min:.2f}')
    plt.axvline(filter_max, color='green', linestyle='--', label=f'Max {parameter_name} = {filter_max:.2f}')
    
    if target_value is not None:
        plt.axvline(target_value, color='blue', linestyle='-', linewidth=2, label=f'Target {parameter_name} = {target_value:.2f}')

    plt.xlabel(parameter_name)
    plt.ylabel('Number of Ground Motions')
    plt.title(f'{parameter_name} Distribution\nTotal: {total_count}, Passing filter: {passed_count}')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


# kz-250722: csv input option for periods
def load_user_defined_periods(infile):
    # initialize
    peridos = []
    # load the periods in the input file path
    infilename, infileext = os.path.splitext(infile)
    file_formats = ['.txt','.csv']
    if infileext not in file_formats:
        print('load_user_defined_periods: WARNING - the suggested period file format is {}.'.format(file_formats))
    # try to load the file anyway
    try:
        df_periods = pd.read_csv(infile,header=None)
    except:
        print('load_user_defined_periods: ERROR - could not load the period file {}.'.format(infile))
        return peridos
    # convert to list
    periods = df_periods[0].tolist()
    # return
    return periods

# kz-250724: check file location and format
def file_validation(infile,file_format=[]):
    # initialize
    outfile = None
    # infile is not defined
    if infile is None:
        return outfile
    # infile format check
    if len(file_format)>0:
        infilename, infileext = os.path.splitext(infile)
        if infileext not in file_format:
            print('file_validation: ERROR - the file format {} is not valid in {}'.format(infileext,file_format))
            return outfile
    # infile location check
    if os.path.isfile(infile):
        outfile = infile
    else:
        print('file_validation: ERROR - the file {} not found'.format(infile))
    # return
    return outfile


# kz-250724: parse user-flatfile
def parse_user_flatfile(infile,user_filters,data_style='PEER_NGA'):
    # initialize (default parameters)
    dict_gmdb = {
        'Periods': [],
        'SA': [],
        'DS575': [],
        'DS595': [],
        'Filename': [],
        'LowestUsableFrequency': [],
        'Magnitude': [],
        'Vs30': [],
        'Distance': [],
    }
    # load the flatfile first
    infilename, infileext = os.path.splitext(infile)
    if infileext in ['.xlsx','.xls']:
        df_flatfile = pd.read_excel(infile,header=0)
    elif infileext in ['.txt','.csv']:
        df_flatfile = pd.read_csv(infile,header=0)
    else:
        print('parse_user_flatfile: ERROR - only supporting txt, csv, xlsx, and xls formats now')
        return dict_gmdb, user_filters
    # parse default parameters
    if data_style == 'PEER_NGA':
        # PEER NGA flatfile format: https://peer.berkeley.edu/research/databases/databases
        dict_gmdb.update({'Magnitude': df_flatfile['Earthquake Magnitude']})
        dict_gmdb.update({'Vs30': df_flatfile['Vs30 (m/s) selected for analysis']})
        # kz: confirming with md that the distance filter is based on epicenter distance as implemented in 2c84f0c
        dict_gmdb.update({'Distance': df_flatfile['EpiD (km)']})
        # get column names with "TxxxS"
        sa_cols = [x for x in df_flatfile.columns if x.startswith('T') and x.endswith('S')]
        if len(sa_cols) < 1:
            print('parse_user_flatfile: ERROR - please check the column names for Sa(T), e.g., T0.010S.')
            return dict_gmdb, user_filters
        dict_gmdb.update({'SA': df_flatfile[sa_cols]})
        dict_gmdb.update({'Periods': pd.DataFrame.from_dict({'Periods': [float(x[1:-1]) for x in sa_cols]})})
        dict_gmdb.update({'LowestUsableFrequency': df_flatfile['Lowest Usable Freq - H1 (Hz)']})
        # file name
        file_cols = [x for x in df_flatfile.columns if x.startswith('File Name')]
        filenames = []
        for i in range(len(df_flatfile.index)):
            filenames.append(df_flatfile[file_cols].iloc[i].tolist())
        dict_gmdb.update({'Filename': pd.DataFrame.from_dict({'Filename': filenames})})
        # DS575 and DS959 are not available in NGA-West2 official flatfile...
    else:
        print('parse_user_flatfile: please revise util.parse_user_flatfile to add the new data_style.')
        return dict_gmdb, user_filters
    # parse parameters required by user-filters
    if user_filters is None:
        # no user-filter defined
        pass
    else:
        user_filtering_vars = list(user_filters.keys())
        # check user-filters are available in the input flatfile
        for cur_var in user_filtering_vars:
            if cur_var not in df_flatfile.columns:
                print('parse_user_flatfile: WARNING - user-filter {} not found in the input flatfile - will be ignored.'.format(cur_var))
                del user_filters[cur_var]
        # add the available filtering parameters to the dictionary
        dict_gmdb.update({cur_var: df_flatfile[cur_var].tolist()})
    # return
    return dict_gmdb, user_filters


        