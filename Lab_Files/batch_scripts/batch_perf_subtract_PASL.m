% Toolbox for batch processing ASL perfusion based fMRI data.
% All rights reserved.
% Ze Wang @ TRC, CFN, Upenn 2004
%
% Batch calculation for the perfusion signals.

for sb = 1:PAR.nsubs % for each subject

    sprintf('Calculate perfusion and CBF signals for subject #%g ... %g subjects left...\n',sb,length(PAR.subjects)-sb)
    for c=1:PAR.ncond

        P=[];
        %ptmp=spm_get('files', PAR.condirs{s,1}, [PAR.subjects{s} '*' PAR.confilters{1} '*img']);
        %P=ptmp(1,:);
        % creating a mask image for removing background
        maskimg = spm_select('FPList', PAR.condirs{sb,c},    ['^mask_perf_cbf\w*\.img$']);
        Filename=spm_select('ExtFPList', PAR.condirs{sb,c}, ['^sr' PAR.confilters{c} '\w*\.img$']);


        Ptmp=spm_select('FPList', char(PAR.M0dirs{sb,c}), ['^sr' PAR.M0filters{c} '\w*\.img$']);
        if  length(deblank( Ptmp (1,:)))~=1
            M0img= Ptmp (1,:);

        end

        %asl_perf_subtract(Filename,FirstimageType, SubtractionType,...
        %   SubtractionOrder,Flag,
        %Timeshift,AslType,labeff,MagType,
        %   Labeltime,Delaytime,Slicetime,TE,M0img,M0seg,maskimg)
        
        % These are the values that impact the actual calculation of the
        % perfusion 
        FirstImage=0; % What is the first image? Control or label
        SubtractionType = 0; 
        SubtractionOrder=0;
        % Flag is more complex, these settings ask what outputs you would
        % like to save and create
        Flag = [1 1 1 0 0 1 0]; % For PASL
        % Flag = [1 1 1 0 0 1 0 1 1]; % For pCASL
        Timeshift = 0.5;
        AslType = 0; % Is the ASL PASL, CASL, or pCASL? 0=PASL, 1=CASL, 2=pCASL 
        labeff=0.95; % What is the label efficiency? 
        MagType = 1; % What field strength? 0 = 1.5T, 1=3T
        Labeltime= 2; % What is the labeling time?
        DelayTime = 1.5; % What is the post labeling delay? 
        SliceTime = 45; % To correct for slice timing delays, what is the slice time (only for 2D acquisitions
        TE=17; % What is the TE for the acquisition?
   
        asl_perf_subtract(Filename, FirstImage, SubtractionType, ...
            SubtractionOrder,     Flag,...
            Timeshift,     AslType,      labeff, MagType,...
            Labeltime, DelayTime, SliceTime, TE, M0img, [], maskimg);
        fprintf('\n%40s%30s','',' ');
    end
end

