# Fresh Convective Cloud Growth Detection from Geostationary Satellite Data

This repository contains a Python script designed to detect fresh convective cloud growth using geostationary satellite data. The script processes data files, applies various checks and criteria to identify convective clouds, and saves the results in both CSV and NetCDF formats.

![Demo](assets/clouds_2010_01_08.gif)

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Prerequisites

Ensure you have Python 3.6+ installed.

### Dependencies

Install the required Python packages using `pip`:

```
pip install numpy netCDF4 xarray pandas
```

### Required Dataset

This repository is set up to work with the NCEP/CPC Level 3 Merged Infrared Brightness Temperatures dataset. See [Dataset Page](https://gpm.nasa.gov/data/directory/ncepcpc-level-3-merged-infrared-brightness-temperatures-0) for more details. The preferred method to download global dataset files is through the [NASA Earthdata portal](https://search.earthdata.nasa.gov/search).

As a teaser, this repo provides two days of hourly input for two different regions - [Subtropical South America](https://github.com/divyanshchg007/fresh-convective-cloud-detection/tree/b6890f7800a02d563c5a59efcde11adfa15264a5/data/SouthAmerica_Subtrop_clipped) and [Northwest India](https://github.com/divyanshchg007/fresh-convective-cloud-detection/tree/b6890f7800a02d563c5a59efcde11adfa15264a5/data/NWIndia_clipped). These regional files were clipped from the global dataset using [Climate Data Operators (CDO) library](https://code.mpimet.mpg.de/projects/cdo). 

### Cloning the Repository

Clone this repository to your local machine using:

```
git clone https://github.com/yourusername/fresh-convective-cloud-detection.git
cd fresh-convective-cloud-detection
```

## Usage

### Running the Script

The core script 'extract_fresh_dcc.py' requires five command-line arguments: 'year', 'mon', 'cold core temperature', 'cloud threshold temperature', and 'cooling rate'. It is currently setup to run for an entire month and save the corresponding monthly list of the time and location of where convective clouds initiate.

1. Set up the necessary parameters and paths:
    - Ensure the 'data_dir' variable in 'extract_fresh_dcc.py' points to your data directory containing geostationary satellite data.
    - Set the 'savepath' variable in to specify where outputs (CSV files and NetCDF files) should be saved.

2. Input the desired 'year', 'mon', 'cold core temp', 'cloud thres', and 'cooling rate' in 'extract_fresh_dcc.py' with appropriate values.
    - cold core temp: Specifies the temperature threshold for identifying cold cores in deep convective clouds. Units: K.
    - cloud thres: Defines the temperature threshold to distinguish convective clouds from other atmospheric features. Units: K.
    - cooling rate: Sets the rate of temperature change used to identify fresh cloud growth. Units: K/hour.

3. Optional: Adjust the cloud definition parameters in the constants dictionary, such as the size of the cloud detection window defined by 'domsize'. 

4. Run the main script with the following arguments:

```
python extract_fresh_dcc.py <year> <month> <cold core temp> <cloud thres> <cooling rate>
```

5. Make sure to apply parallax correction to the detected cloud locations using the llax_correction.py script. Modify the location of the satellite (latitude, longitude) based on the region. Currently set to Meteosat for analysis over NW India.

### Example

```
python extract_fresh_dcc.py 2020 7 235 253 -10
```

This example processes data from July of 2020 with a cold core temperature of 235K, a cloud threshold temperature of 253K, and a cooling rate of 10K/hour.

### Directory Structure

Ensure the following directories exists:
    - data_dir
    - savepath

Adjust the paths in the 'filenamegen' function and the 'savepath' variable.

## Contributing

Contributions are welcome! Please fork this repository and submit pull requests to contribute. Ensure your code follows best practices and is well-documented.

### Steps to Contribute:
    - Fork this repository.
    - Create a new branch (git checkout -b feature-branch).
    - Make your changes.
    - Commit your changes (git commit -am 'Add new feature').
    - Push to the branch (git push origin feature-branch).
    - Create a new Pull Request.

## License

This project is licensed under the MIT License. See the LICENSE file for details.
