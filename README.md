# U-Net D2min Analysis

AI-based nucleus segmentation and tracking tool with downstream analysis of D2min

This tool was developed and modified to measure and analysis the collective migration of two-dimensional cell sheets that are imaged over time using phase contrast microscopy. The whole workflow has three main components: 

1. U-Net nucleus segmentation (trained_model folder) 
2. Nucleus tracking (UNetCrocker_Tracker_v2.m)
3. Analyze nucleus motion compared to their close neighbors (their neighbors are founding using and calculate how willing they are to follow their neighbors’ motion (D2min). 

The current version of the code was adapted for analysis of the data in:

Shuyao Gu, Rachel M. Lee, Zackery Benson, Chenyi Ling, Michele I. Vitolo, Stuart S. Martin, Joe Chalfoun, Wolfgang Losert, "Label-free Cell Tracking Enables Collective Motion Phenotyping in Epithelial Monolayers"

https://doi.org/10.1101/2021.12.14.472148

## Dependencies and Installation
All scripted were tested in MATLAB 2020a with the Image Processing Toolbox, Deep Learning Toolbox, and Computer Vision Toolbox.

## Workflow

### `/trained_model/Inferencing_2D_with_overlap.m`
The `trained_model` folder contains the trained network in MATLAB ((`UNet_single_cell_seg_norm_256_weights_1_3_Workspace.mat`) to segment the nucleus of single stem cells in phase contrast, as well as the scripts in MATLAB to run the inferencing code (`Inferencing_2D_with_overlap.m`). This is the raw output from the neural network, one still need to launch some morphological cleaning on the output and some object separation methods (like FogBank or Watershed in Fiji) to get the final result. Here we included the Fiji macro `Watershed_macro.ijm` that has been used to segment objects.

This code should work without GPU support, but it will be much faster if you have a GPU on your machine (strongly recommend).

User Inputs
- `read_imDir`: path to the raw images to segment
- `read_img_common_name`: common name of the images
- `read_file_ext`: extension of the name of the images
- `save_pxDir`: path to save the segmented masks

### `UNetCrocker_Tracker_v3.m`
This code uses `track.m` function to track the motion of individual nucleus over time. Segmented binary mask from the previous step is needed as input, and this code extract the centroid of each region of interest and performs tracking. The output tracks are stored in variable `tr`.

User Inputs
- `number_Frames`: number of frames in the input image sequence
- `pathname`: path to the segmented binary mask to perform tracking
- `savefolder`: path to save tracks
- `file_short`: common name of the images
- `bwareaopen_size`: any region of interested smaller than this size is removed as noise, this value should be smaller than the average nucleus size of your cells
- `xmin`, `xmax`: min and max of pixels in the x direction of the input images
- `ymin`, `ymax`: min and max of pixels in the y direction of the input images

### `find_n_neighbors2D.m`
Finds the nearest n neighbors of each cell nuclei. 

User Inputs
- `r`: how many nearest neighbors to use
- `loadname`: name of the tracks to load 
- `rootdir`: path to tracks
- `pathname`: path to save neighbors

### `D2min2D_v5.m`
Uses output from `find_n_neighbors2D.m` to calculate the D2min for each nucleus at each time t

User Inputs
- `findneighborname`: name of the neighbors to input
- `findneighbordir`: path to the neighbors to input
- `savedfilename`: name of the D2min file to save
- `save_path`: path to save D2min file
- `total_frame`: number of frames to calculate D2min
