

function save_image_to_disk(I,file_path,common_name,img_id,zero_pad,input_image_name,file_ext,output_class)

% Itinialize allowable skipped inputs
if nargin < 4, img_id = 1; end
if nargin < 5, zero_pad = 3; end
if nargin < 6, input_image_name = []; end
if nargin < 8 , output_class = 'uint8'; end

% Get input_file name and extension if specified
if ~isempty(input_image_name), [~,init_name,input_file_ext] = fileparts(input_image_name); end

% Use same extension as input if extension is skipped and input image name is used
if ~isempty(input_image_name) && (nargin < 7 || isempty(file_ext)), file_ext = input_file_ext; end

% Use default extension if not specified at this stage
if nargin < 7 && isempty(file_ext), file_ext = '.tif'; end

% Use default common name "img_" if nothing is pecified
if isempty(common_name) && isempty(input_image_name), save_img_name = ['img_' sprintf(['%0' num2str(zero_pad) 'd'],img_id) file_ext]; end

% Use same input name if no common name is specified
if isempty(common_name) && ~isempty(input_image_name), save_img_name = [init_name file_ext]; end

% Use same input name with common name is specified
if ~isempty(common_name) && ~isempty(input_image_name), save_img_name = [common_name '_' init_name file_ext]; end

% Use common name with zero padding if only common is used
if ~isempty(common_name) && isempty(input_image_name), save_img_name = [common_name '_' sprintf(['%0' num2str(zero_pad) 'd'],img_id) file_ext]; end

% Cast to appropriate class
if string(output_class) == "uint8", I = im2uint8(I); end
if string(output_class) == "uint16", I = im2uint16(I); end
if string(output_class) == "single", I = im2single(I); end
if string(output_class) == "double", I = im2double(I); end

% Write image to disk
imwrite(I,fullfile(file_path, save_img_name))





