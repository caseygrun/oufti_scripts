
This is the repository for MATLAB scripts adapted from Christina Lin and Daniel Lee's research on studying the dynamics of T3SS expression in Pseudomonas aeruginosa.
The scripts manipulate output files from Oufti, created by the Jacobs-Wagner lab and restructured by Casey Grun (github.com/caseygrun/oufti).
The scripts require MATLAB R2018b to run. Later versions may or may not work. The following toolboxes are also required: 

- Curve fitting toolbox
- Mapping toolbox

Included in this repository are several `functions` and some `example_data`. See the script `examples.m` for examples of how to use these functions. Make sure the `functions` directory is in your MATLAB Path ("Right click > Add folder to MATLAB Path"). Then you should be able to run `examples.m`, or ask for help on any function with `help FUNCTION_NAME`. 

In general, you will first use Oufti meshes to calculate MFIs per cell with the `parse_MFI_lineage` function. This requires an Oufti mesh (`.mat`) file and image (`.tif`) files for each channel (e.g. GFP, RFP, etc.). From there, `plot_gaussian` can fit a Gaussian mixture model for a single frame and tell you what % of cells in an image are ON. If you have movies (or single images) of constitutively ON and OFF controls, you can use `calculate_threshold` to calculate, at each time point, a threshold MFI whereupon 50% of cells at that MFI will be ON or OFF. 

Three simple analysis scripts are also provided; one, for training of control data to determine a threshold for T3SS-ON and T3SS-OFF cells as described in the full manuscript ('threshold_calculate.m'). Once this model is "trained", another script, 'classify_colonies.m', will calculate a percentage of T3SS-ON cells based on their fluoresence intensities, given a list of mean fluoresence intensities (e.g., as output by 'GFP_MFI_parse.m'). Finally, 'tree_analyze.m' performs change point analysis along lineage trajectories as also described in the accompanying manuscript, parsing .csv data files in the format of '/example_data_and_outputs/outputs/lineage/gfpmfi_lineage.csv'. 
