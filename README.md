# Fresh Convective Cloud Growth Detection from Geostationary Satellite Data

This repository contains a Python script designed to detect fresh convective cloud growth using geostationary satellite data. The script processes data files, applies various checks and criteria to identify convective clouds, and saves the results in both CSV and NetCDF formats.

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

### Cloning the Repository

Clone this repository to your local machine using:

```
git clone https://github.com/yourusername/fresh-convective-cloud-detection.git
cd fresh-convective-cloud-detection
```

## Usage

### Running the Script

The script requires five command-line arguments: 'start year', 'end year', 'cold core temperature', 'cloud threshold temperature', and 'cooling rate'. It is currently setup to run for an entire month and save the corresponding monthly list of the time and location of where convective clouds initiate.

1. Set up the necessary parameters and paths:
    - Ensure the 'data_dir' variable in 'extract_fresh_dcc.py' points to your data directory containing geostationary satellite data.
    - Set the 'savepath' variable to specify where outputs (CSV files and NetCDF files) should be saved.

2. Replace 'start year', 'end year', 'cold core temp', 'cloud thres', and 'cooling rate' in 'extract_fresh_dcc.py' with appropriate values.
    - cold core temp: Specifies the temperature threshold for identifying cold cores in deep convective clouds. Units: K.
    - cloud thres: Defines the temperature threshold to distinguish convective clouds from other atmospheric features. Units: K.
    - cooling rate: Sets the rate of temperature change used to identify fresh cloud growth. Units: K/hour.

3. Run the main script with the following arguments:

```
python extract_fresh_dcc.py <start year> <end year> <cold core temp> <cloud thres> <cooling rate>
```

### Example

```
python extract_fresh_dcc.py 2020 2021 235 253 -10
```

This example processes data from 2020 to 2021 with a cold core temperature of 235K, a cloud threshold temperature of 253K, and a cooling rate of 10K/hour.

### Directory Structure

Ensure the following directories exists:
    - data_dir
    - savepath

Adjust the paths in the 'filenamegen' function and 'savepath' variable.

## Contributing

Contributions are welcome! Please fork this repository and submit pull requests to contribute. Ensure your code follows best practices and is well-documented.

Steps to Contribute:
    - Fork this repository.
    - Create a new branch (git checkout -b feature-branch).
    - Make your changes.
    - Commit your changes (git commit -am 'Add new feature').
    - Push to the branch (git push origin feature-branch).
    - Create a new Pull Request.

## License

This project is licensed under the MIT License. See the LICENSE file for details.
