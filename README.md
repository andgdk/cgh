> This repository is archived

# Computer Generated Hologram Generation with Scilab (CGH)

## Generate CGH and Signal Reconstructions

Generate binary phase masks from black and white images and compute a
reconstruction from the phase mask into the signal plane.

<img alt="Resulting computation of the OSS example" src="assets/plot-example.png" width="700px" />

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
    exec('index.sce');
    ```
4. Observe the resulting graphs.
5. Add your own black and white signal images.
6. Update the signal image name variable in `index.sce` and execute the program.

## Project Structure

-   `index.sce`: Main script file.
-   `src/`: Source code directory.
-   `data/`: Data files directory with examples.
-   `results/`: Directory for output results.

## Examples

TODO: Include examples of how to use the project here.
