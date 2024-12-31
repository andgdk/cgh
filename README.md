# Computer Generated Hologram Generation with Scilab

## Requirements

-   Tested with Scilab 2025.0.0 (version 6.1.0 or higher)
-   Image Processing and Computer Vision Toolbox (IPCV)

## Installation

1. Download and install Scilab from [Scilab Official Website](https://www.scilab.org/download).
2. Install the Image Processing and Computer Vision Toolbox (IPCV) dependency via the Scilab console:
    ```scilab
    atomsInstall("IPCV")
    ```
3. Clone this repository:
    ```sh
    git clone https://github.com/andgdk/cgh.git
    cd cgh
    ```

## Usage

1. Open Scilab.
2. Navigate to the project directory.
3. Run the main script:
    ```scilab
    exec('index.sci');
    ```

## Project Structure

-   `index.sci`: Main script file.
-   `src/`: Source code directory.
-   `data/`: Data files directory with examples.
-   `results/`: Directory for output results.

## Examples

TODO: Include examples of how to use the project here.

## TODO

-   Consider creating Jupyter notebook using the Scilab kernel to provide examples
-   Explain manufacturing constraints (binary DOE element)
-   Add user option to compute CGH with 8-step quantized sampling/continuous phase in addition to binary phase
-   Design further abort criteria for iteration
-   Plot CGH graphs as subgraphs
-   Investigate how signatures of functions can be improved
