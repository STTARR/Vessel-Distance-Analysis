%% image_alignment_FRONTIERS
% This script is used to perform an intensity-based, semi-automated image
% registration between two histological images of similar morphology
% (either as a restained sample or serial section).

% Prior to running the script, root_static and root_moving need to be
% changed to point to the folders containing the static (reference image to
% align to) and moving (images that need to be aligned) images,
% respectively. scale_factor specifies the factor by which the image is
% downsampled, and should not be less than 1. Ideally, with a powerful
% computer, the transformation matrix should be computed on an unresized
% image (scale_factor=1), however we found a scale factor of 9 is sufficient
% to get a decent registration on most modern computers and laptops. Note,
% this does not affect the resolution of the output image. Downsampling is
% only performed on an internal intermediate image to calculate the
% transformation matrix, which is then applied on the full size moving
% images. In the event that registration fails, a manual alignment can be
% done.

% There are 4 inputs for image registration. file_static is the static
% image and file_moving is the moving image, which are used to calculate
% the transformation matrix. Note, after calculating this matrix, the
% images have not yet been registered. files_to_transform is a list of one
% or more images that the transformation matrix must be applied to.
% Typically, these are the separated channels of an IF image, or the
% separated R-G-B components of a brightfield image. static_to_rewrite is
% one or more images that no registration needs to be applied to. It is
% simply read and rewritten to have consistent write parameters and
% metadata with the rest of the images, as some image analysis software
% such as Definiens may have trouble loading improperly configured images.
% All input images must be single-channel grayscale images. If a scanner
% reconstructs a color image, the channels must be split into individual
% stain components.

% The outputs are transformedImage and static_to_rewrite. transformedImage
% are the transformed images, aligned to static_to_rewrite's spatial
% coordinates. static_to_rewrite is the rewritten static image. All images
% are written as .tifs in the Aligned_Images subbfolder


% Defining global variables
clear;
root_static='K:\000000\000000_Definiens\images\Dan Cojacari\re_aligned_201plus_images\temp_out';
root_moving='K:\000000\000000_Definiens\images\Dan Cojacari\re_aligned_201plus_images\temp_out';
scale_factor=9;

%Select static and moving image(s)
file_static=uigetfile([root_static '\*.tif'],'Select static image');
file_moving=uigetfile([root_moving '\*.tif'],'Select moving image');
files_to_transform=uigetfile([root_moving '\*.tif'],'Select images to apply transformation to','multiselect','on');
static_to_rewrite=uigetfile([root_static '\*.tif'],'Select static images to read and rewrite','multiselect','on');

%Creation of an output folder if it doesn't already exist
aligned_folder=fullfile(root_moving,'Aligned_Images\');
if ~isdir(aligned_folder)
    mkdir(aligned_folder)
end

%Read static and moving images
static=imread(fullfile(root_static,file_static));
moving=imread(fullfile(root_moving,file_moving));

%Calculate transformation matrix using RegisterImages_FRONTIERS (see top
%comments in RegisterImages_FRONTIERS for details)
[saved_tform,needSave] = RegisterImages_FRONTIERS(imadjust(static), imadjust(moving),'ScalingFactor',1/scale_factor);

%In the event that only 1 image needs to be transformed, it's name will
%need to be converted to a cell array
if ~iscell(files_to_transform);
    files_to_transform=cellstr(files_to_transform);
end

%apply and write transformation to moving images
for i=1:length(files_to_transform);

    curr=imread(fullfile(root_moving,files_to_transform{i}));
    transformedImage=imwarp(curr,affine2d(saved_tform),'OutputView',imref2d(size(static)),'FillValues',[0]);
    imwrite(transformedImage,fullfile(aligned_folder,files_to_transform{i}),'compression','jpeg','rowsperstrip',8);
end
clear i

%In the event that only 1 static image is specified to be rewritten, it's name will
%need to be converted to a cell array
if ~iscell(static_to_rewrite);
    static_to_rewrite=cellstr(static_to_rewrite);
end

%rewrite the static image (no transformation applied) to have consistent
%metadata with the other .tifs (e.g. compression type, number of rows to
%include in each image strip)
for i=1:length(static_to_rewrite);

    imwrite(imread(fullfile(root_static,static_to_rewrite{i})),fullfile(aligned_folder,static_to_rewrite{i}),'compression','jpeg','rowsperstrip',8);
    %Definiens can't read images generated from czi-to-tiff
    %conversion tool; need to read in and rewrite
end
