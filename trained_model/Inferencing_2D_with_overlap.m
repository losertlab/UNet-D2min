

%% Set Paths
read_imDir = 'example/sample_image/'; % path to the raw images to segment
read_img_common_name = 'sample'; %common name of the images
read_file_ext = '.tif';
save_pxDir = 'example/segmented/'; % path to save the segmented masks
add_read_image_to_output_name = 1; % Set this variable to 0 if output images are to have only common_name_indx.tif, if set to 1 then output image name will be read_image_name_common_name_indx.tif
save_img_common_name = 'inferenced'; % if not empty, the program automatically adds an _ between common name and image id 
zero_pad = 4; % zero padding added to the image name
output_class = 'uint8';
save_file_ext = '.tif';

%% Load Model
load("trained_model/UNet_single_cell_seg_norm_256_weights_1_3_Workspace.mat",'Unet_seg') % path to the saved model and its name
tile_size = Unet_seg.Layers(1,1).InputSize; % Get the input size image used during training the network
neighbor_size = 128; % pixels to add to the tiled image from surrounding neighboring neighbors before inferencing
normalization_choice = 'zscore'; % 'zerocenter', 'zscore', 'none'

%% Set Output Directory
if ~exist(save_pxDir,'dir')
    mkdir(save_pxDir);
end

%% Read Images
input_images = dir(fullfile(read_imDir,['*' read_img_common_name '*' read_file_ext]));
input_images = {input_images.name}';
nb_images = length(input_images);

%% Inference
for im = 1:nb_images
    % Read image
    display("Inferencing image " + num2str(im) + input_images{im})
    I = imread(fullfile(read_imDir,input_images{im}));
    I = normalize_image(I,normalization_choice); % Normalize image
    I = im2uint8(I);
    [nb_rows,nb_cols] = size(I);
    
    % Pad image by tile_size to make it easier to navigate
    Ipad = padarray(I,[neighbor_size neighbor_size]);
    
    % Compute the number of tiles per dimension
    ni = nb_rows / tile_size(1);
    nj = nb_cols / tile_size(2);
    % Check if we need to handle the last column and the last row, when the number of tiles with overlap does not equal the width or hight of the image
    handle_row = 0;
    handle_col = 0;
    if floor(ni) ~= ni, handle_row = 1; end
    if floor(nj) ~= nj, handle_col = 1; end
    ni = floor(ni); % make sure number is int
    nj = floor(nj); % make sure number is int
    
    % Initialize the segmented image
    S = zeros(nb_rows,nb_cols,'uint8');
    for i = 1:ni
        for j = 1:nj
            % Extract the tile with extra neighboring pixels
            i_min = (i-1)*tile_size(1)+1;
            i_max = i*tile_size(1);
            j_min = (j-1)*tile_size(2)+1;
            j_max = j*tile_size(2);
            I1 = Ipad(i_min:i_max+2*neighbor_size, j_min:j_max+2*neighbor_size);
            % Segment using UNet model
            S1 = semanticseg(I1,Unet_seg,'OutputType','uint8');
            % Save in correct place in the whole image
            S(i_min:i_max, j_min:j_max) = S1(neighbor_size+1:neighbor_size+tile_size(1),neighbor_size+1:neighbor_size+tile_size(2));
        end
        % Handle last column
        if handle_col
            i_min = (i-1)*tile_size(1)+1;
            i_max = i*tile_size(1);
            j_min = nb_cols-tile_size(2)+1;
            j_max = nb_cols;
            I1 = Ipad(i_min:i_max+2*neighbor_size, j_min:j_max+2*neighbor_size);
            % Segment using UNet model
            S1 = semanticseg(I1,Unet_seg,'OutputType','uint8');
            % Save in correct place in the whole image
            S(i_min:i_max, j_min:j_max) = S1(neighbor_size+1:neighbor_size+tile_size(1),neighbor_size+1:neighbor_size+tile_size(2));
        end
    end
    
    % Handle last row
    if handle_row
        for j = 1:nj
            % Extract the tile with extra neighboring pixels
            i_min = nb_rows-tile_size(1)+1;
            i_max = nb_rows;
            j_min = (j-1)*tile_size(2)+1;
            j_max = j*tile_size(2);
            I1 = Ipad(i_min:i_max+2*neighbor_size, j_min:j_max+2*neighbor_size);
            % Segment using UNet model
            S1 = semanticseg(I1,Unet_seg,'OutputType','uint8');
            % Save in correct place in the whole image
            S(i_min:i_max, j_min:j_max) = S1(neighbor_size+1:neighbor_size+tile_size(1),neighbor_size+1:neighbor_size+tile_size(2));
        end
        % Handle last column and last row, that lower corner tile!
        if handle_col
            i_min = nb_rows-tile_size(1)+1;
            i_max = nb_rows;
            j_min = nb_cols-tile_size(2)+1;
            j_max = nb_cols;
            I1 = Ipad(i_min:i_max+2*neighbor_size, j_min:j_max+2*neighbor_size);
            % Segment using UNet model
            S1 = semanticseg(I1,Unet_seg,'OutputType','uint8');
            % Save in correct place in the whole image
            S(i_min:i_max, j_min:j_max) = S1(neighbor_size+1:neighbor_size+tile_size(1),neighbor_size+1:neighbor_size+tile_size(2));
        end
    end
    % Save image to disk
    S = S-1; % put background to 0
    input_image_name = [];
    if add_read_image_to_output_name, input_image_name = input_images{im}; end
    save_image_to_disk(S,save_pxDir,save_img_common_name,im,zero_pad,input_image_name,save_file_ext,output_class)
end
