% Reset the orientation for the images

% set the center to the center of each images
% batch_3Dto4D;
batch_reset_orientation_pCASL;

%realign the functional images to the first functional image of eachsubject
batch_realign_pCASL;

%coreg the functional images to the anatomical image
batch_coreg_pCASL;

batch_filtering_pCASL;

%smooth the coreged functional images
batch_smooth_pCASL;

%create perfusion mask
batch_create_mask_pCASL;