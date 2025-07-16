import matplotlib.pyplot as plt
import numpy as np

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
