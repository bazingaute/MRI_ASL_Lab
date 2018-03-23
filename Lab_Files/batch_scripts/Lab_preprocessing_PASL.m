% set the center to the center of each images
batch_reset_orientation_PASL;

%set the center to AC
batch_setorigin_PASL;

%realign the functional images to the first functional image of eachsubject
batch_realign_PASL;

%coreg the functional images to the anatomical image
batch_coreg_PASL;

%smooth the coreged functional images
batch_smooth_PASL;

%create perfusion mask
batch_create_mask_PASL;