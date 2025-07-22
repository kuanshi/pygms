# PyGMS: Python Ground Motion Selection

## Overview

PyGMS is a Python package for selecting and scaling ground motion records to match target response spectra. The package implements various target spectra including Conditional Spectrum (CS), Conditional Mean Spectrum (CMS), and CS with duration constraints (CS+duration). 

The implementation provides a user-friendly interface through JSON configuration files.

## Features

- Support for multiple target intensity measure types:
  - Conditional Spectrum (CS)
  - Conditional Mean Spectrum (CMS)
  - Conditional Spectrum plus duration constraints (CS+duration)
- Integration with USGS hazard API for ground motion prediction equations
- Built-in significant duration models
- Includes dataset of NGA-West ground motions (Sa from 0.01s to 10s, Ds5-75, and Ds5-95)

## Structure

The package consists of:

- **pygms.py**: Main script that interfaces with user-defined JSON configuration files
- **Classes**:
  - `TargetIntensityMeasure`: Processes target spectra configuration
  - `GroundMotionSelection`: Handles the ground motion selection algorithm
- **Modules**:
  - `USGSHazardGMM`: Fetches USGS API for different ground motion prediction equations
  - `SignificantDurationModel`: Implements significant duration models
  - `CorrelationModel`: Implements intensity-measure correlation models

## Usage

### Basic Usage

1. Create a JSON configuration file (see examples below)
2. Run the selection algorithm using the provided Jupyter notebook [examples.ipynb](https://github.com/bakerjw/pygms/blob/main/examples.ipynb)  or directly with the pygms.py Python script:

```python
import pygms

# Load configuration from JSON file
pygms.run_selection("example1.json")
```

## Configuration Options

Configuration is handled through a JSON file with the following main sections:

1. **Target Intensity Measure**: Define the target spectrum type, periods, and conditional intensity measure
2. **Ground Motion Models**: Select models and parameters (magnitude, distance, Vs30)
3. **Selection Parameters**: Define scaling approach, number of ground motions, and error weights

## Examples

The package includes three example configurations:

1. **example1.json**: CS+duration target
2. **example2.json**: CS only target
3. **example3.json**: CMS only target

## Selection Algorithm

The ground motion selection algorithm follows these steps:

1. Sample a set of pseudo spectra from the target distribution
2. Perform initial selection to find the best match to individual pseudo spectra (minimizing weighted total error)
3. Apply a greedy algorithm to optimize the weighted total error for:
   - Mean and standard deviation of spectral acceleration (Sa)
   - Mean and standard deviation of significant duration (for CS+duration targets)

## Output

The default output is a CSV file containing:
1. Ground motion filenames (RSN identifiers)
2. Unscaled spectra values
3. Significant duration values
4. Scaling factors

## References

Baker, J. W., and Lee, C. (2018). “[An Improved Algorithm for Selecting Ground Motions to Match a Conditional Spectrum](https://www.jackwbaker.com/Publications/Baker_Lee_(2018)_GM_Selection,_JEE.pdf).” Journal of Earthquake Engineering, 22(4), 708–723. This algorithm is also implemented in [Matlab](https://github.com/bakerjw/CS_Selection).


## License

[Add license information]

## Contributors

- Kuanshi Zhong (Implementation in Python)
- Jack Baker (Development of a previous implementation in Matlab)
- Marguerite Dubertret (Performance and functional improvements)

## Contact

- Kuanshi Zhong [zhongki@ucmail.uc.edu]
- Jack Baker [bakerjw@stanford.edu]

## TODO
- Integrating site-specific seismic disaggregation info
- Multi-threading/processing function
