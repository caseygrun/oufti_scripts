
This is the repository for MATLAB scripts adapted from Christina Lin and Daniel Lee's research on the dynamics of T3SS expression in _Pseudomonas aeruginosa_. Scripts have been adapted by Casey Grun for general-purpose use. MIT License. 

The scripts manipulate output files from Oufti, created by the Jacobs-Wagner lab and restructured by Casey Grun (github.com/caseygrun/oufti).
The scripts require MATLAB R2018b to run. Later versions may or may not work. The following toolboxes are also required; you can install them with "Home" > "Add-Ons" > "Get Add-Ons" from within MATLAB. 

- Curve fitting toolbox
- Mapping toolbox

Included in this repository are several `functions` and usage examples. To use the usage examples, download the contents of [this entire folder](https://yale.box.com/s/0k148uadhiav9ky5xmi39pqwosichru3) to a subdirectory called `example_data` in the same folder as this document. See the script `examples.m` for examples of how to use these functions (`examples.m` relies on `example_data`, so download it first). Make sure the `functions` directory is in your MATLAB Path (Right click > "Add folder to MATLAB Path"). Then you should be able to run `examples.m`, or ask for help on any function in `functions` with `help FUNCTION_NAME` on the MATLAB command line. 

In general, you will first use Oufti meshes to calculate MFIs per cell with the `parse_MFI_lineage` function. This requires an Oufti mesh (`.mat`) file and image (`.tif`) files for each channel (e.g. GFP, RFP, etc.). This script works for movies or for single images, and generates two kinds of files: a .CSV file (one per colony) with MFI and linage information, and .DAT files, one per colony per channel per frame with just the MFIs. From there, you can use `plot_MFIs` to plot a histogram of MFIs in each frame of each colony in a dataset. `plot_gaussian` can fit a Gaussian mixture model for a single frame and tell you what % of cells in an image are ON. If you have movies (or single images) of constitutively ON and OFF controls, you can use `calculate_threshold` to calculate, at each time point, a threshold MFI whereupon 50% of cells at that MFI will be ON or OFF. Then use `classify_colonies` to classify each cell as ON or OFF, compared to the controls. All of this is demonstrated in `examples.m`
