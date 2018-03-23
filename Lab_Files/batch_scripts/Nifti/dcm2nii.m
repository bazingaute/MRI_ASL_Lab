function varargout = dcm2nii(src, dataFolder, varargin)
% DCM2NII converts Siemens dicom files into NIfTI or img/hdr files. 
% 
% DCM2NII(dcmSource, niiFolder, outFormat, MoCoOption, subjName)
% 
% The input arguments are all optional:
%  1. dicom source. It can be a zip file or a folder (including subfolders)
%     containing dicom files, or dicom file(s). It can also contain
%     wildcards like 'run1_*' for all files start with 'run1_'.
%  2. data folder to save result files.
%  3. output file format:
%      0 or 'nii'    for uncompressed single file, run1.nii.
%      1 or 'nii.gz' for single file, run1.nii.gz (default).
%      2 or 'img'    for triplet files (img/hdr/mat, no compression).
%      3 or 'img.gz' for compressed triplet files (img/hdr/mat with gz).
%  4. MoCo series options:
%      0 create files for both original and MoCo series.
%      1 ignore MoCo series if both present (default).
%      2 ignore original series if both present.
%     Note that if only one of the two series is present, it will be
%     created always.
%  5. subject name to create nii. The code will do only one subject each
%     time. If you have a folder or zip file for multiple subjects (not
%     recommended), you can specifiy which subject to create with the 5th
%     input. The name should be the entry for the 'Last name' at scanner
%     console (dicominfo field 'PatientName.FamilyName').
% 
% Typical examples to use:
%  dcm2nii('D:/myProj/zip/subj1.zip', 'D:/myProj/subj1/data'); % zip file
%  dcm2nii('D:/myProj/subj1/dicom/', 'D:/myProj/subj1/data'); % folder
% 
% Less useful examples:
%  dcm2nii('D:/myProj/dicom/', 'D:/myProj/subj2/data', [], [], 'subj2');
%  dcm2nii('D:/myProj/dicom/run2*', 'D:/myProj/subj/data');
% 
% File dialog will appear if any of first two inputs is left out or empty.
% 
% If the first input is a zip file, such as those downloaded from dicom
% server, DCM2NII will extract files into a temp folder, create NIfTI files
% into the data folder, and then delete the temp folder. For this reason,
% it is better keep the zip file as a backup.
% 
% If a folder is the data source, DCM2NII will convert all files in the
% folder and subfolders.
% 
% If there is DTI series, bval and bvec files will be generated for FSL. A
% Matlab data file, dcmHeaders.mat, will be saved into the data folder too.
% This file contains dicom header from the first file for created series
% and some information from last file in field LastFile. For DTI series,
% B_value and DiffusionGradientDirection for all directions are saved into
% the mat file. For MoCo series, motion parameters, RBMoCoTrans and
% RBMoCoRot, are also saved.
% 
% Some information, such as TE, phase encoding direction and effective
% dwell time are stored in descrip of nii header. These are useful for
% fieldmap B0 unwarp correction. Acquisition start time and date are also
% stored, and this may be useful if one wants to align the functional data
% to some physiological recording, like pulse, respiration or ECG.
% 
% The output file names adopt ProtocolName of each series used on scanner
% console. If both original and MoCo series are created, '_MoCo' will be
% appended for MoCo series. For phase image, such as those from field map,
% '_phase' will be appended to the name. In case of name conflict,
% SeriesNumber, such as '_s005', will be appended to make file names
% unique. You are suggested to use short and descriptive ProtocolName on
% the Siemens scanner console. The name better contains only letters,
% numbers and underscores.
% 
% Besides Image Processing Toolbox from Matlab, this needs Jimmy Shen's
% NIfTI toolbox at http://research.baycrest.org/~jimmy/NIfTI. For
% convenience, the related code are included in this file.
% 
% Tested for dicom data from syngo MR versions 'B17' and '2004A 4VA25A'.

% Thanks to:
% Jimmy Shen's Tools for NIfTI and ANALYZE image,
% Chris Rorden's dcm2nii pascal source code,
% Przemyslaw Baranski for direction cosine matrix to quaternions. 

% History (yymmdd):
% 130512 Publish to CCBBI users (Xiangrui Li).
% 130513 Convert img from uint16 to int16 if range allows;
%        Expand output file format to img/hdr/mat.
% 130515 Change creation order to acquisition order (more natural).
%        If MoCo series is included, append _MoCo in file names.
% 130516 Use SpacingBetweenSlices, if exists, for SliceThickness. 
% 130518 Use NumberOfImagesInMosaic in CSA header (work for older data).
% 130604 Add scl_inter/scl_slope and special naming for fieldmap.
% 130614 Work out the way to get EffectiveEchoSpacing for B0 unwarp.
% 130616 Add needed dicom field check, so it won't err later.
% 130618 Reorient if non-mosaic or slice_dim is still 3 and no slice flip.
% 130619 Simplify DERIVED series detection. No '_mag' in fieldmap name.
% 130629 Improve the method to get phase direction;
%        Permute img dim1&2 (no -90 rotation) & simplify xform accordingly.
% 130711 Make MoCoOption smarter: create nii if only 1 of 2 series exists.
% 130712 Remove 5th input (allHeader). Save memory by using partial header.
% 130712 Bug fix: dim_info with reorient. No problem since no EPI reorient.
% 130715 Use 2 slices for xform. No slice flip needed except revNum mosaic.
% 130716 Take care of lower/upper cases for output file names;
%        Apply scl_slope and inter to img if range allows and no rounding;
%        Save motion parameters, if any, into dcmHeader.mat.
% 130722 Ugly fix for isMos, so it works for '2004A 4VA25A' phase data;
%        Store dTE instead of TE if two TE are used, such as fieldmap.
% 130724 Add two more ways for dwell time, useful for '2004A 4VA25A' dicom.
% 130801 Can't use DERIVED since MoCo series can be labeled as DERIVED.
% 130807 Check PixelSpacing consistency for a series;
%        Prepare to publish to Matlab Central.
% 130809 Add 5th input for subjName, so one can choose a subject.
% 130813 Store ImageComments, if exists and meaningful, into aux_file.
% 130818 Expand source to dicom file(s) and wildcards like run1*.dcm.
%        Update fields in dcmHeader.mat, rather than overwriting the file.
%        Include save_nii etc in the code for easy distribution.

%% Deal with output format first, so error out if invalid
if nargin<3 || isempty(varargin{1}), fmt = 1; % default .nii.gz
else fmt = varargin{1};
end

if (isnumeric(fmt) && (fmt==0 || fmt==1)) || ...
      (ischar(fmt) && ~isempty(regexpi(fmt, 'nii')))
    ext = '.nii';
elseif (isnumeric(fmt) && (fmt==2 || fmt==3)) || (ischar(fmt) && ...
        (~isempty(regexpi(fmt, 'hdr')) || ~isempty(regexpi(fmt, 'img'))))
    ext = '.img';
else
    error(' Invalid output file format (the 3rd input).');
end

if (isnumeric(fmt) && mod(fmt,2)) || ...
        (ischar(fmt) && ~isempty(regexpi(fmt, '.gz')))
    ext = [ext '.gz']; % gzip file
end

%% Deal with MoCo option
if nargin<4 || isempty(varargin{2})
    MoCo = 1; % only original series by default
else
    MoCo = varargin{2};
    if ~any(MoCo==0:2)
        error(' Invalid MoCoOption. The 4th input must be 0, 1 or 2.');
    end
end

%% Deal with 5th input: we do one subject once
if nargin<5 || isempty(varargin{3})
    subjProvided = false; subj = '';
else 
    subjProvided = true; subj = varargin{3};
    if ~ischar(subj), error(' Invalid subject name.');end
end

%% Deal with data source
varargout = {};
srcIsZip = false;
if nargin<1 || isempty(src)
    src = questdlg('Select data source:', 'Dicom to NIfTI', ...
                   'Zip File', 'Folder', 'Dicom Files', 'Zip File');
    if isempty(src), return; end
    if strcmp(src, 'Folder')
        src = uigetdir(cd, 'Select root folder containing DICOM files');
        if isnumeric(src), return; end
        dcmFolder = src;
    elseif strcmp(src, 'Zip File')
        [src, dcmFolder] = uigetfile('*.zip', ...
            'Select a zip file containing DICOM files');
        if isnumeric(src), return; end
        srcIsZip = true;
        src = fullfile(dcmFolder, src);
    else % one or more dicom files
        [src, dcmFolder] = uigetfile('*.dcm', ...
            'Select one or more DICOM files', 'MultiSelect', 'on');
        if isnumeric(src), return; end
        src = cellstr(src); % in case only 1 file selected
        fnames = fullfile(dcmFolder, src);
    end
elseif isnumeric(src)
    error('Invalid dicom source.');    
elseif ~exist(src, 'file') % like input: run1*.dcm
    fnames = dir(src);
    if isempty(fnames), error('%s does not exist.', src); end
    dcmFolder = fullfile(getfield(what(fileparts(src)), 'path'));
    fnames = fullfile(dcmFolder, {fnames.name});    
elseif isdir(src) % folder
    dcmFolder = src;
elseif ischar(src) % 1 dicom or zip file
    dcmFolder = fullfile(getfield(what(fileparts(src)), 'path'));
    if isdicom(src)
        fnames = dir(src);
        fnames = fullfile(dcmFolder, {fnames.name});
    else % likely zip file
        srcIsZip = true;
    end
elseif iscellstr(src) % multiple files
    dcmFolder = fullfile(getfield(what(fileparts(src)), 'path'));
    n = length(src);
    fnames = src;
    for i = 1:n
        foo = dir(src{i});
        if isempty(foo), error('%s does not exist.', src{i}); end
        fnames{i} = fullfile(dcmFolder, foo.name); 
    end
else 
    error('Unknown dicom source.');
end
dcmFolder = fullfile(getfield(what(dcmFolder), 'path'));

%% Deal with dataFolder
if nargin<2 || isempty(dataFolder)
    dataFolder = uigetdir(dcmFolder, 'Select a folder to save data files');
    if dataFolder==0, return; end
end
if ~isdir(dataFolder), mkdir(dataFolder); end
dataFolder = fullfile([getfield(what(dataFolder), 'path'), filesep]);

disp('Xiangrui.Li@gmail.com ''s dcm2nii.m for Siemens MR');

%% Unzip if zip file is the source
tic;
if srcIsZip
    if ~iszip(src), error('Unknown source file.'); end
    [~, fname] = fileparts(src);
    dcmFolder = sprintf('%stmpDcm%s/', dataFolder, fname);
    if ~isdir(dcmFolder), mkdir(dcmFolder); end
    disp(['Extracting files from ' fname '.zip ...']);

    cmd = sprintf('unzip -qq -o %s -d %s', src, dcmFolder);
    err = system(cmd); % first try system unzip
    if err, unzip(src, dcmFolder); end % Matlab's unzip is too slow
    drawnow;
end 

%% Get all file names including those in subfolders, if not specified
if ~exist('fnames', 'var')
    dirs = genpath(dcmFolder);
    dirs = textscan(dirs, '%s', 'Delimiter', pathsep);
    dirs = dirs{1}; % cell str
    fnames = {};
    for i = 1:length(dirs)
        curFolder = [dirs{i} filesep];
        foo = dir(curFolder); % all files and folders
        foo([foo.isdir]) = []; % remove folders
        foo = strcat(curFolder, {foo.name});
        fnames = [fnames foo]; %#ok<*AGROW>
    end
end
nFile = length(fnames);
if nFile<1, error(' No files found in the folder.'); end

%% Check each file,  store header in cell array h
% cell is more flexible than struct for this purpose. 
% Empty cells take only a little memory.
% these are to check valid dicom. EchoTime can exclude PhoenixZIPReport.
flds1 = {'Filename' 'SeriesNumber' 'InstanceNumber' 'ProtocolName' ...
    'PatientName' 'ImageType' 'Columns' 'Rows' 'EchoTime' 'RepetitionTime' ...
    'ImageOrientationPatient' 'ImagePositionPatient' 'PixelSpacing' ...
    'SliceThickness' 'InPlanePhaseEncodingDirection' 'SoftwareVersion' ...
    'Private_0029_1010' 'Private_0029_1020'};

% these are for files other than 1st. If extra needed, only add it here
flds2 = {'Filename' 'AcquisitionTime' 'Columns' 'Rows' 'ProtocolName' ...
    'ImagePositionPatient' 'ImageOrientationPatient' 'PixelSpacing' ...
    'Private_0019_100c' 'Private_0019_100e' ... % these two for DTI
    'Private_0019_1025' 'Private_0019_1026'}; % these two for MOCO

junk = {'\MEAN' '\DUMMY IMAGE' '\FMRI\TTEST' '\FMRI\DESIGN' ... % GLM
        '\DIFFUSION\ADC\' '\DIFFUSION\FA\' '\DIFFUSION\TRACEW\'}; % DTI

h = {}; % in case of no dicom files at all
fprintf('Validating %g files ...\n', nFile);
for k = nFile:-1:1 % reversed order may help pre-allocate cell arrays
    fname = fnames{k};
    if ~isdicom(fname), continue; end
    s = dicominfo(fname); % takes ~20 ms, bottle neck for a lot of files
    if any(~isfield(s, flds1)), continue; end % also exclude some DTI junk
    if isType(s, junk),  continue; end
    i = s.SeriesNumber; j = s.InstanceNumber; 
    subj1 = s.PatientName.FamilyName;
    
%     if s.SamplesPerPixel>1
%         fprintf(2, ' Not implemented for SamplesPerPixel>1: %s skipped.', ...
%             s.ProtocolName);
%         continue;
%     end
%     
%     if ~strcmp(s.PatientPosition, 'HFS')
%         fprintf(2, 'Not Head first supine: %s skipped.', s.ProtocolName);
%         continue;
%     end
        
    % if not the same subject, do the first only
    if isempty(subj)
        subj = subj1; % store it for later check
    elseif ~strcmpi(subj, subj1)
        if ~subjProvided
            fprintf(2, ' %s is for a different subject %s. Skipped.\n', ...
                fname, subj1);
        end
        continue;
    end

    % For fieldmap mag image, we use the one with short TE, which has
    % better quality. This also skips repeated copy of a file.
    try %#ok ignore the error if the cell in h hasn't been filled.
       if s.EchoTime >= h{i}{j}.EchoTime, continue; end
    end
    
    if j>1 % only several fields for rest: to save memory if lots of files
        a = struct;
        for m = 1:length(flds2)
            if isfield(s, flds2{m}), a.(flds2{m}) = s.(flds2{m}); end
        end
        h{i}{j} = a;
    else % store full header only for first file
        h{i}{j} = s;
    end    
end
if nargout, varargout = {subj}; end % return subject ID

%% Check headers: remove file-missing and dim-inconsistent runs etc
nRun = length(h);
if nRun<1, error('No dicom files found for %s.', subj); end
keep = true(1, nRun); % true for useful runs
isMoCo = false(1, nRun); % deal moco together later
for i = 1:nRun
    if isempty(h{i}), keep(i) = 0; continue; end % no file for this run
    ind = cellfun(@isempty, h{i});
    if any(ind) % missing file(s)
        k = find(~ind, 1);
        fprintf(2, ['%s, Series %g, Instance %g file missing. ' ...
            'Run skipped.\n'], h{i}{k}.ProtocolName, i, find(ind,1));
        keep(i) = 0; continue; % skip
    end
    
    s = h{i}{1}; % header for 1st file
    isMoCo(i) = isType(s, '\MOCO\');
    
    % check dimension, orientation, pixel size consistency
    for j = 2:length(h{i})
        s1 = h{i}{j};
        err = ~isequal([s.Columns s.Rows], [s1.Columns s1.Rows]);
        err = err || sum(abs(s.PixelSpacing-s1.PixelSpacing)) > 0.01;
        if err || (sum(abs(s1.ImageOrientationPatient - ...
               s.ImageOrientationPatient)) > 0.01); % 0.01 is arbituary     
            fprintf(2, ['Inconsistent pixel size, image orientation ' ...
             'and/or dimension for subject %s, %s, Series %g. ' ...
             'Run skipped.\n'], subj, s.ProtocolName, s.SeriesNumber);
            keep(i) = 0; break; % skip
        end
    end   
end

ind = find(isMoCo); % decide MoCo after checking all series
for i = 1:length(ind)
    if MoCo==1 && keep(ind(i)-1) % in case original skipped, keep MOCO
        keep(ind(i)) = 0; continue; % skip MOCO
    elseif MoCo==2 && keep(ind(i)) % in case MOCO skipped, keep original
        keep(ind(i)-1) = 0; % skip previous series (original)
    end
end
h = h(keep); % remove all unwanted series once

%% Generate unique file names
% Unique names are in format of ProtocolName_s007. Special cases are: 
%  for phase image, such as field_map phase, append '_phase' to the name;
%  for MoCo series, append '_MoCo' to the name if both series are needed.
nRun = length(h); % update it, since we may have removed some
rNames = cell(1, nRun);
for i = 1:nRun
    s = h{i}{1};
    a = s.ProtocolName;
    a(~isstrprop(a, 'alphanum')) = '_'; % make str valid for field name
    if isType(s, '\P\'), a = [a '_phase']; end % phase image
    if MoCo==0 && isType(s, '\MOCO\'), a = [a '_MoCo']; end
    a = sprintf('%s_s%03g', a, s.SeriesNumber);
    rNames{i} = a;
end
rNames = genvarname(rNames); % add 'x' if started with a digit, and more

% After following sort, we need compare only neighboring names. Remove
% _s007 if there is no conflict. Have to ignore case for Windows & MAC
fnames = rNames; % copy it, keep letter cases
[rNames, iRuns] = sort(lower(rNames)); 
for i = 1:nRun
    a = rNames{i}(1:end-5); % remove _s003
    % no conflict with both previous and next name
    if nRun==1 || ... % only one run
         (i==1    && ~strcmpi(a, rNames{2}(1:end-5))) || ... % first
         (i==nRun && ~strcmpi(a, rNames{i-1}(1:end-5))) || ... % last
         (i>1 && i<nRun && ~strcmpi(a, rNames{i-1}(1:end-5)) ...
         && ~strcmpi(a, rNames{i+1}(1:end-5))); % middle ones
        fnames{iRuns(i)}(end+(-4:0)) = [];
    end
end

%% Now ready to convert nii run by run
fprintf('Converting %g series into %s: subject %s\n', nRun, ext, subj);
for i = 1:nRun
    nFile = length(h{i});
    fprintf(' %s %4g\n', fnames{i}, nFile); % show info and progress
    
    h{i}{1}.LastFile = h{i}{nFile}; % store last header into 1st
    s = h{i}{1}; % 1st header (with last) is needed for nii info
    img = dicomread(s); % initialize with proper data type and img size
    img(:, :, 2:nFile) = 0; % pre-allocate
    for j = 2:nFile, img(:,:,j) = dicomread(h{i}{j}); end % read rest files
 
    % Save FSL bval and bvec files for DTI data
    if all(isfield(s.LastFile, {'Private_0019_100c' 'Private_0019_100e'}))
        bval = zeros(nFile, 1);
        bvec = zeros(nFile, 3);
        for j = 2:nFile % no these values for 1st file
            bval(j)   = h{i}{j}.Private_0019_100c;
            bvec(j,:) = h{i}{j}.Private_0019_100e;
        end
        h{i}{1}.B_value = bval; % store all into header of 1st file
        h{i}{1}.DiffusionGradientDirection = bvec;
        
        m = reshape(s.ImageOrientationPatient, 3, 2); % row & col dir cos
        m(:,3) = null(m'); % slice direction
        bvec = bvec * m; % diffusionGradient * dir cos matrix
        bvec(:,2) = -bvec(:,2); % change y sign

        fid = fopen([dataFolder fnames{i} '.bval'], 'w');
        fprintf(fid, '%g ', bval); % one row
        fclose(fid);

        str = repmat('%.15g ', 1, nFile);
        fid = fopen([dataFolder fnames{i} '.bvec'], 'w');
        fprintf(fid, [str '\n'], bvec); % 3 rows by # direction cols
        fclose(fid);
    end
    
    % Store motion parameters for MoCo series
    if all(isfield(s, {'Private_0019_1025' 'Private_0019_1026'}))
        trans = zeros(nFile, 3);
        rotat = zeros(nFile, 3);
        for j = 1:nFile
            trans(j,:) = h{i}{j}.Private_0019_1025;
            rotat(j,:) = h{i}{j}.Private_0019_1026;
        end
        h{i}{1}.RBMoCoTrans = trans;
        h{i}{1}.RBMoCoRot = rotat;
    end
    h{i} = h{i}{1}; % keep 1st dcm header only
    
    if isMosaic(s), img = mos2vol(img, s); end % convert mosaic to volume
    img = permute(img, [2 1 3 4]); % dicom img is row-major
    if s.BitsStored<16 && isa(img, 'uint16')
        img = int16(img); % use int16 if lossless. seems always true
    end
    
    nii = make_nii(img); % need NIfTI toolbox
    nii = set_nii_header(nii, s); % a little complicate here
    save_nii(nii, [dataFolder fnames{i} ext]); % need NIfTI toolbox
end

h = cell2struct(h, fnames, 2); % convert into struct
fname = [dataFolder 'dcmHeaders.mat'];
if exist(fname, 'file') % only update fields if file exists
    S = load(fname);
    for i = 1:length(fnames), S.h.(fnames{i}) = h.(fnames{i}); end
    h = S.h; %#ok
end
save(fname, 'h', '-v7'); % -v7 better compatibility
fprintf('Elapsed time by dcm2nii.m is %.1f seconds\n\n', toc);
if srcIsZip, rmdir(dcmFolder, 's'); end % delete tmp dicom folder
return;

%% Subfunction
% Return true if any of keywords is in s.ImageType
function tf = isType(s, keywords)
keywords = cellstr(keywords);
for i = 1:length(keywords)
    key = strrep(keywords{i}, '\', '\\'); % for regexp
    tf = ~isempty(regexp(s.ImageType, key, 'once'));
    if tf, return; end
end

% Set most nii header. Also reorient img is needed.  s is dicom header for
% the first file, with partial header of last file in field 'LastFile'.
function nii = set_nii_header(nii, s)
dim = nii.hdr.dime.dim(2:4); % image dim, set by make_nii
pixdim = [s.PixelSpacing(:)' s.SliceThickness];
if isfield(s, 'SpacingBetweenSlices')
    pixdim(3) = s.SpacingBetweenSlices; % in case gap used
end
R44 = reshape(s.ImageOrientationPatient, 3, 2);
R44(:,3) = null(R44'); % this 3x3 matrix contains a lot of info!
[~, ixyz] = max(abs(R44)); % orientation info: perm of 1:3
R44 = R44 * diag(pixdim); % apply vox size
R44(4,4) = 1; % now it is 4x4
fps_dim = 1:3; % freq_dim, phase_dim, slice_dim
if strcmpi(s.InPlanePhaseEncodingDirection, 'ROW'), fps_dim = [2 1 3]; end
isMos = isMosaic(s);

% dim_info byte: freq_dim, phase_dim, slice_dim low to high, each 2 bits
nii.hdr.hk.dim_info = fps_dim * [1 4 16]'; % update later if reorient
nii.hdr.dime.xyzt_units = 10; % mm (2) and seconds (8)
nii.hdr.dime.pixdim(2:5) = [pixdim s.RepetitionTime/1000]; % voxSize and TR
if isfield(s, 'ImageComments') && ~isType(s, '\MOCO\')
    nii.hdr.hist.aux_file = s.ImageComments; % char[24], info only
elseif isfield(s, 'SeriesDescription')
    nii.hdr.hist.aux_file = s.SeriesDescription;
end
seq = asc_header(s, 'tSequenceFileName'); % like '%SiemensSeq%\ep2d_bold'
[~, seq] = strtok(seq, '\'); seq = strtok(seq, '\'); % like 'ep2d_bold'
id = ''; if isfield(s, 'PatientID'), id = s.PatientID; end
nii.hdr.hk.db_name = [seq ';' id]; % char[18], optional

% save some useful info in descrip: optional
descrip = '';
if isfield(s, 'AcquisitionTime')
    descrip = sprintf('time=%s;', s.AcquisitionTime(1:11));
end
if isfield(s.LastFile, 'AcquisitionTime')
    descrip = sprintf('%stimeN=%s;', descrip, s.LastFile.AcquisitionTime(1:11)); 
end
if isfield(s, 'AcquisitionDate')
    descrip = sprintf('%sdate=%s;', descrip, s.AcquisitionDate);
end
TE0 = asc_header(s, 'alTE[0]')/1000; % s.EchoTime stores only 1 TE
dTE = asc_header(s, 'alTE[1]')/1000 - TE0; % TE difference for fieldmap
if ~isempty(dTE)
    descrip = sprintf('dTE=%g;%s', abs(dTE), descrip);
elseif ~isempty(TE0)
    descrip = sprintf('TE=%g;%s', TE0, descrip);
end

% This relies on the fact that EPI is always stored in mosaic. Maybe it is
% better to use MRAcquisitionType (2D or 3D) for this?
if isMos % Get dwell time & slice timing info, reverse slices if needed
    hz = csa_header(s, 'BandwidthPerPixelPhaseEncode');
    dwell = 1000 ./ hz / dim(fps_dim(2)); % in ms
    if isempty(dwell) % true for syngo MR 2004A
        % ppf = [1 2 4 8 16] represent [4 5 6 7 8] 8ths PartialFourier
        % ppf = asc_header(s, 'sKSpace.ucPhasePartialFourier');
        lns = asc_header(s, 'sKSpace.lPhaseEncodingLines');
        dur = csa_header(s, 'SliceMeasurementDuration');
        dwell = dur ./ lns; % ./ (log2(ppf)+4) * 8;
    end
    if isempty(dwell) % next is not accurate, so as last resort
        dur = csa_header(s, 'RealDwellTime') * 1e-6; % ns to ms
        dwell = dur * asc_header(s, 'sKSpace.lBaseResolution');
    end
    if ~isempty(dwell)
        descrip = sprintf('dwell=%.3g;%s', dwell, descrip);
    end
    
    t = sort(csa_header(s, 'MosaicRefAcqTimes'));
    dur = min(diff(t)) / 1000; % 2.5ms diff, in seconds
    if isempty(dur), dur = s.RepetitionTime/dim(3)/1000; end % for old data
    nii.hdr.dime.slice_duration = dur; 
    sc = asc_header(s, 'sSliceArray.ucMode'); % 1 2 or 4
    if isempty(sc), sc = 0; % unknown
    elseif sc==4, sc = 5 - mod(dim(3),2)*2; % 3/5 interleaved odd/even
    end
    nii.hdr.dime.slice_code = sc;
    nii.hdr.dime.slice_end = dim(3)-1; % 0-based, slice_start default to 0

    % Following should be the only situation to flip slices. The method 
    % using SliceNormalVector is not right: it will flip normal Sag slices.
    % The sPosition in asc_header is for slices after following flip.
    % The actual keyword is sSliceArray.ucImageNumbSag, ...Cor or ...Tra
    if asc_header(s, 'sSliceArray.ucImageNumb') % reversed numbering
        nii.img = flipdim(nii.img, 3); % flip slices before reorient
    end
end

% Phase encoding direction
iPhase = ixyz(fps_dim(2)); % phase axis index
phPos = csa_header(s, 'PhaseEncodingDirectionPositive'); % 0 or 1
if isempty(phPos), phPos=1; fprintf(2, ' Guessed phase direction.\n'); end
if iPhase ~= 3, phPos = ~phPos; end % dicom LPS vs nii RAS
phDir = 'xyz'; phDir = phDir(iPhase);
if ~phPos, phDir = ['-' phDir]; end % phase dir such as 'y' or '-y'
descrip = sprintf('phase=%s;%s', phDir, descrip);
nii.hdr.hist.descrip = descrip; % char[80], drop from end if exceed

% data slope and intercept: apply to img if no rounding error 
nii.hdr.dime.scl_slope = 1; % default scl_inter is 0
scl = isfield(s, {'RescaleSlope' 'RescaleIntercept'});
if any(scl)
    slope = 1; inter = 0;
    if scl(1), slope = s.RescaleSlope; end
    if scl(2), inter = s.RescaleIntercept; end
    val = sort([nii.hdr.dime.glmax nii.hdr.dime.glmin] * slope + inter);
    dClass = class(nii.img);
    if isa(nii.img, 'float') || (val(1)>=intmin(dClass) && ...
            val(2)<=intmax(dClass) && mod(slope,1)==0 && mod(inter,1)==0)
        nii.img = nii.img * slope + inter; % apply to img if no rounding
    else
        nii.hdr.dime.scl_slope = slope;
        nii.hdr.dime.scl_inter = inter;
    end
end

% Transformation matrix: important feature for nii
if isMos % mosaic EPI or DTI
    x = [dim(1:2)'*[1 1]/2; 0 dim(3)-1; 1 1]; % ijk 1st&last slice center
    y = [mosaic_sPosition(s, [0 dim(3)-1]); 1 1]; % pos for slice center
else % non-mosaic image
    x = [zeros(2); 0 dim(3)-1; 1 1]; % upper-left corner of 1st&last slices
    y = [s.ImagePositionPatient s.LastFile.ImagePositionPatient; 1 1];
end
if isequal(y(:,1), y(:,2)) % fake 2nd slice in case of only 1 slice
    x(:,2) = [0 0 1 0]'; % 2nd slice top-left corner without translation
    y(:,2) = R44(:,3); % pos for 2nd slice without translation
end
R = [R44(:,1:2) y] / [eye(4,2) x]; % also take care of slice direction
R(1:2,:) = -R(1:2,:); % dicom LPS vs nifti RAS

if any(isnan(R(:))) % unlikely: skip reorient and won't set s/qform_code
	fprintf(2, ' Failed to compute transformation matrix.\n');
    return;
end

% Reorient for mosaic with nVols<3 & (no slice_dim change & no slice flip).
% If FSL can read related info from nii header, we can always reorient.
[~, perm] = sort(ixyz); % may permute 3 dimensions in this order
ind4 = ixyz + [0 4 8]; % index in 4xn matrix
flip = R(ind4)<0; % may flip an axis if true
if (nii.hdr.dime.dim(5)<3 || (perm(3)==3 && ~flip(3))) && dim(3)>1 && ...
        (any(flip) || ~isequal(perm, 1:3)) % skip if already standard view
    ijk = (dim-1) .* flip; % ijk for original [0 0 0] after possible flip
    rotM = zeros(4); % only 0 1 or -1 for 3x3 part
    rotM(ind4) = sign(R(ind4)); % assign 1 or -1 to 3 max
    rotM(:,4) = [ijk(perm) 1]; % ijk for original [0 0 0] after flip&perm
    R = R / rotM; % this is xform matrix after flip & perm
    
    % update header accordingly
    nii.hdr.hk.dim_info = ixyz(fps_dim) * [1 4 16]'; % useful for EPI only
    nii.hdr.dime.dim(2:4) = dim(perm);
    nii.hdr.dime.pixdim(2:4) = pixdim(perm);
    if flip(3) && isMos && sc>0 % for future if we reorient mosaic
        nii.hdr.dime.slice_code = sc+mod(sc,2)*2-1; % 1<->2, 3<->4, 5<->6
    end
    
    for k = 1:3, if flip(k), nii.img = flipdim(nii.img, k); end; end
    nii.img = permute(nii.img, [perm 4]); % permute img after flip
end

nii.hdr.hist.sform_code = 1; % SCANNER_ANAT
nii.hdr.hist.srow_x = R(1,:);
nii.hdr.hist.srow_y = R(2,:);
nii.hdr.hist.srow_z = R(3,:);

nii.hdr.hist.qform_code = 1; % SCANNER_ANAT
nii.hdr.hist.qoffset_x = R(1,4);
nii.hdr.hist.qoffset_y = R(2,4);
nii.hdr.hist.qoffset_z = R(3,4);

R = normc(R(1:3, 1:3)); % for quaternion
proper = round(det(R)); % always 1 if reorient, otherwise can be -1
if proper<0, R(:,3) = -R(:,3); end
nii.hdr.dime.pixdim(1) = proper; % -1 or 1 

q = dcm2quat(R); % 3x3 dir cos matrix to quaternion
if q(1)<0, q = -q; end % as MRICron
nii.hdr.hist.quatern_b = q(2);
nii.hdr.hist.quatern_c = q(3);
nii.hdr.hist.quatern_d = q(4);
return;
% hdr.hist.magic, glmax, glmin will be taken care of by save_nii.
% magic: 'ni1', hdr/img pair; 'n+1', single nii file, empty for ANALYZE. 
% Not used: char data_type[10], char regular


%% Subfunction, reshape mosaic into volume, remove padded zeros
function img = mos2vol(img, s)
nSlice = csa_header(s, 'NumberOfImagesInMosaic'); % not work for some data
if isempty(nSlice), nSlice = asc_header(s, 'sSliceArray.lSize'); end
nMos = ceil(sqrt(nSlice)); % always NxN tiles
[nc, nr, nv] = size(img); % number of col, row and vol
sz = [nc nr] / nMos; % slice size

% Get index in vol for one mosaic: not elegant, but brief
[rm, cm] = ind2sub([nc nr], 1:nc*nr); % row, col sub in mosaic
rv = mod(rm-1, sz(1)) + 1; % row index in vol
cv = mod(cm-1, sz(2)) + 1; % col index in vol
sv = floor((cm-1)/sz(2))+1 + floor((rm-1)/sz(1))*nMos; % slice index
iv = sub2ind([sz nMos^2], rv, cv, sv); % singlar index in vol

img = reshape(img, [nc*nr nv]); % one col per mosaic
img(iv, :) = img; % change to vol order
img = reshape(img, [sz nMos^2 nv]); % vol now
img(:, :, nSlice+1:end, :) = []; % remove padded slices        

%% Subfunction, return a parameter from CSA header
% Return nx1 array if the val is numeric, a string if the value is str, or
% a cell str if multiple str. The last one may not happen at all.
function val = csa_header(s, key)
% This table is unnecessary. But it may be safer and faster, and also makes
% main code more readable by avoiding those private tags
keys = {'Private_0019_100a' 'NumberOfImagesInMosaic';
        'Private_0019_100b' 'SliceMeasurementDuration';
        'Private_0019_100c' 'B_value';
        'Private_0019_100e' 'DiffusionGradientDirection';
        'Private_0019_1018' 'RealDwellTime';
        'Private_0019_1025' 'RBMoCoTrans';
        'Private_0019_1026' 'RBMoCoRot';
        'Private_0019_1028' 'BandwidthPerPixelPhaseEncode';
        'Private_0019_1029' 'MosaicRefAcqTimes';
        'Private_0051_100b' 'AcquisitionMatrixText';
        'Private_0051_1011' 'ImaPATModeText';
        'Private_0051_1013' 'PositivePCSDirections'};
ind = find(strcmp(keys(:,2), key));
if ~isempty(ind) && isfield(s, keys{ind,1}) % use main header instead
    val = s.(keys{ind,1});
    if isnumeric(val), val = double(val); end
    return;
end

csa = [s.Private_0029_1010' s.Private_0029_1020']; % image and series info
ind = strfind(csa, [0 key 0]); % only whole word, offset by -1
val = {};
for j = 1:length(ind) % if keyword repeated, try each till we get val
    i0 = ind(j) + 76; % skip name vm vr syngodt (64+4+4+4)
    n = typecast(csa(i0+(1:4)), 'int32'); % number of items
    if n<=0 || n>512, continue; end % give up if weird, 512 arbituary
    i0 = i0 + 8; % skip nitems xx (4+4)
    for i = 1:n % often times, n=6, but only 1st has len>0
        len = typecast(csa(i0+(1:4)), 'int32'); % # of bytes for item
        if len<=0 || len>512, break; end
        i0 = i0 + 16; % skip len len & int32 * 2 junk
        val{i,1} = char(csa(i0+(1:len-1))); % drop null
        i0 = i0 + ceil(double(len)/4)*4; % multiple 4-byte
    end
    if ~isempty(val), break; end
end

foo = str2double(val); % try to convert into numeric
if isempty(foo) || all(~isnan(foo)) % return numeric if not NaN
    val = foo; 
elseif length(val) == 1 % return char if single cell. Seems always
    val = val{1};
end

%% Subfunction, Convert 3x3 direction cosine matrix to quaternion
% Simplied from Quaternions by Przemyslaw Baranski 
function q = dcm2quat(R)
if ~isequal(size(R), [3 3]), error('R must be a 3x3 matrix.'); end

q = zeros(1, 4);
q(1) = sqrt(max(trace(R)+1, 0)) / 2; % if negative, zero it
q(2) = sqrt(1 + R(1,1) - R(2,2) - R(3,3)) / 2;
q(3) = sqrt(1 + R(2,2) - R(1,1) - R(3,3)) / 2;
q(4) = sqrt(1 + R(3,3) - R(2,2) - R(1,1)) / 2;
[m, ind] = max(q);

switch ind
    case 1
        q(2) = (R(3,2) - R(2,3)) /m/4;
        q(3) = (R(1,3) - R(3,1)) /m/4;
        q(4) = (R(2,1) - R(1,2)) /m/4;
    case 2
        q(1) = (R(3,2) - R(2,3)) /m/4;
        q(3) = (R(1,2) + R(2,1)) /m/4;
        q(4) = (R(3,1) + R(1,3)) /m/4;
    case 3
        q(1) = (R(1,3) - R(3,1)) /m/4;
        q(2) = (R(1,2) + R(2,1)) /m/4;
        q(4) = (R(2,3) + R(3,2)) /m/4;
    case 4
        q(1) = (R(2,1) - R(1,2)) /m/4;
        q(2) = (R(3,1) + R(1,3)) /m/4;
        q(3) = (R(2,3) + R(3,2)) /m/4;
end

%% Subfunction: return slice position in CSA series ASC header.
% Position points to slice center, not upper-left corner as dicom.
function sPos = mosaic_sPosition(s, iSlice)
nSlice = asc_header(s, 'sSliceArray.lSize');
if any(iSlice>nSlice-1), error('Index exceeds number of slices.'); end
ori = ['Sag'; 'Cor'; 'Tra']; % xyz
n = length(iSlice);
sPos = zeros(3, n); % zero values omitted in ASC header
for i = 1:n
    for j = 1:3
        key = sprintf('[%g].sPosition.d%s', iSlice(i), ori(j,:));
        foo = asc_header(s, ['sSliceArray.asSlice' key]);
        if isempty(foo), continue; end
        sPos(j,i) = foo;
    end
end

%% Subfunction: return parameter in CSA series ASC header.
function val = asc_header(s, key)
str = char(s.Private_0029_1020');
k0 = strfind(str, '### ASCCONV BEGIN ###');
k  = strfind(str, '### ASCCONV END ###');
str = str(k0:k); % avoid key before BEGIN and after END
k = strfind(str, [char(10) key]); % start with new line: safer
if isempty(k), val = []; return; end
str = strtok(str(k(1):end), char(10)); % the line
[~, str] = strtok(str, '='); % '=' and the vaule
str = strtrim(strtok(str, '=')); % remvoe '=' and space 

if strncmp(str, '""', 2) % str parameter
    val = str(3:end-2);
elseif strncmp(str, '"', 1) % str parameter for version like 2004A
    val = str(2:end-1);
elseif strncmp(str, '0x', 2) % hex parameter, convert to decimal
    val = sscanf(str(3:end), '%x', 1);
else % decimal
    val = str2double(str);
end

%% Subfunction: needed only for old version data.
function tf = isMosaic(s)
tf = isType(s, '\MOSAIC');
if tf || ~isType(s, '\P\'), return; end % to be safer
% the ugly fix below is only for syngo MR 2004A 4VA25A phase image, which
% is not labeled as MOSAIC in ImageType. Siemens bug I believe
v = strtok(s.SoftwareVersion, 'syngo MR');
if strncmp(v, '2004', 4) % I know only 'syngo MR 2004A 4VA25A' for now
    % for 'syngo MR B17' fieldmap img, lSize>1 even it is not masaic
    tf = asc_header(s, 'sSliceArray.lSize')>1;
end

%% Subfunction: return true if the file is a zip file.
function tf = iszip(fname)
fid = fopen(fname);
if fid<0, tf = false; return; end
sig = fread(fid, 4);
fclose(fid);
tf = isequal(sig, [80 75 3 4]'); % zip file signature

%% Subfunction: simplified from Jimmy Shen's NIfTI toolbox
function save_nii(nii, fileprefix)
%  Check file extension. If .gz, unpack it into temp folder
if length(fileprefix) > 2 && strcmp(fileprefix(end-2:end), '.gz')
    if ~strcmp(fileprefix(end-6:end), '.img.gz') && ...
            ~strcmp(fileprefix(end-6:end), '.hdr.gz') && ...
            ~strcmp(fileprefix(end-6:end), '.nii.gz')
        error('Please check filename.');
    end
    
    v = version;
    if str2double(v(1:3)) < 7.1 || ~usejava('jvm')
        error(['Please use MATLAB 7.1 (with java) and above, ' ...
            'or run gunzip outside MATLAB.']);
    else
        gzFile = 1;
        fileprefix = fileprefix(1:end-3);
    end
end

filetype = 1;

%  Note: fileprefix is actually the filename you want to save
if findstr('.nii',fileprefix) & strcmp(fileprefix(end-3:end), '.nii') %#ok<*AND2,*FSTR>
    filetype = 2;
    fileprefix(end-3:end)='';
end

if findstr('.hdr',fileprefix) & strcmp(fileprefix(end-3:end), '.hdr')
    fileprefix(end-3:end)='';
end

if findstr('.img',fileprefix) & strcmp(fileprefix(end-3:end), '.img')
    fileprefix(end-3:end)='';
end

write_nii(nii, filetype, fileprefix);

%  gzip output file if requested
if exist('gzFile', 'var') && gzFile
    if filetype == 1
        gzipOS([fileprefix, '.img']);
        gzipOS([fileprefix, '.hdr']);
    elseif filetype == 2
        gzipOS([fileprefix, '.nii']);
    end
end

if filetype == 1
    %  So earlier versions of SPM can also open it with correct originator
    M=[[diag(nii.hdr.dime.pixdim(2:4)) -[nii.hdr.hist.originator(1:3).* ...
        nii.hdr.dime.pixdim(2:4)]'];[0 0 0 1]]; %#ok
    save([fileprefix '.mat'], 'M');
end
return					% save_nii

%% Subfunction: use system gzip if available (faster)
function gzipOS(fname)
persistent useMatlabgzip
if isempty(useMatlabgzip), [useMatlabgzip, ~] = system('gzip -h'); end
if useMatlabgzip
    gzip(fname);
    delete(fname);
else
    system(['gzip -f ' fname]); % overwrite
end
return % gzipOS

%% Subfunction: simplified from Jimmy Shen's NIfTI toolbox
function write_nii(nii, filetype, fileprefix)
hdr = nii.hdr;
switch double(hdr.dime.datatype),
    case   1,    hdr.dime.bitpix = int16(1 );    precision = 'ubit1';
    case   2,    hdr.dime.bitpix = int16(8 );    precision = 'uint8';
    case   4,    hdr.dime.bitpix = int16(16);    precision = 'int16';
    case   8,    hdr.dime.bitpix = int16(32);    precision = 'int32';
    case  16,    hdr.dime.bitpix = int16(32);    precision = 'float32';
    case  32,    hdr.dime.bitpix = int16(64);    precision = 'float32';
    case  64,    hdr.dime.bitpix = int16(64);    precision = 'float64';
    case 128,    hdr.dime.bitpix = int16(24);    precision = 'uint8';
    case 256,    hdr.dime.bitpix = int16(8 );    precision = 'int8';
    case 511,    hdr.dime.bitpix = int16(96);    precision = 'float32';
    case 512,    hdr.dime.bitpix = int16(16);    precision = 'uint16';
    case 768,    hdr.dime.bitpix = int16(32);    precision = 'uint32';
    case 1024,   hdr.dime.bitpix = int16(64);    precision = 'int64';
    case 1280,   hdr.dime.bitpix = int16(64);    precision = 'uint64';
    case 1792,   hdr.dime.bitpix = int16(128);   precision = 'float64';
    otherwise
        error('This datatype is not supported');
end

hdr.dime.glmax = round(double(max(nii.img(:))));
hdr.dime.glmin = round(double(min(nii.img(:))));

if filetype == 2
    fid = fopen([fileprefix '.nii'], 'w');
    if fid < 0, error('Cannot open file %s.nii.',fileprefix); end
    hdr.dime.vox_offset = 352;
        
    hdr.hist.magic = 'n+1';
    save_nii_hdr(hdr, fid);
    skip_bytes = double(hdr.dime.vox_offset) - 348;
    fwrite(fid, zeros(1, skip_bytes), 'uint8');
else
    fid = fopen([fileprefix '.hdr'], 'w');
    if fid < 0, error('Cannot open file %s.hdr.', fileprefix); end
    hdr.dime.vox_offset = 0;
    hdr.hist.magic = 'ni1';
    save_nii_hdr(hdr, fid);
        
    fclose(fid);
    fid = fopen([fileprefix '.img'], 'w');
end

fwrite(fid, nii.img, precision);
fclose(fid);
return;					% write_nii

%% Subfunction: simplified from Jimmy Shen's NIfTI toolbox
function save_nii_hdr(hdr, fid)
fseek(fid, 0 ,'bof');
fwrite(fid, hdr.hk.sizeof_hdr(1),    'int32');	% must be 348.
fwrite(fid, padChar(hdr.hk.data_type, 10), 'uchar');
fwrite(fid, padChar(hdr.hk.db_name, 18), 'uchar');
fwrite(fid, hdr.hk.extents(1),       'int32');
fwrite(fid, hdr.hk.session_error(1), 'int16');
fwrite(fid, hdr.hk.regular(1),       'uchar');	% might be uint8
fwrite(fid, hdr.hk.dim_info(1),      'uchar');

fwrite(fid, hdr.dime.dim(1:8),        'int16');
fwrite(fid, hdr.dime.intent_p1(1),  'float32');
fwrite(fid, hdr.dime.intent_p2(1),  'float32');
fwrite(fid, hdr.dime.intent_p3(1),  'float32');
fwrite(fid, hdr.dime.intent_code(1),  'int16');
fwrite(fid, hdr.dime.datatype(1),     'int16');
fwrite(fid, hdr.dime.bitpix(1),       'int16');
fwrite(fid, hdr.dime.slice_start(1),  'int16');
fwrite(fid, hdr.dime.pixdim(1:8),   'float32');
fwrite(fid, hdr.dime.vox_offset(1), 'float32');
fwrite(fid, hdr.dime.scl_slope(1),  'float32');
fwrite(fid, hdr.dime.scl_inter(1),  'float32');
fwrite(fid, hdr.dime.slice_end(1),    'int16');
fwrite(fid, hdr.dime.slice_code(1),   'uchar');
fwrite(fid, hdr.dime.xyzt_units(1),   'uchar');
fwrite(fid, hdr.dime.cal_max(1),    'float32');
fwrite(fid, hdr.dime.cal_min(1),    'float32');
fwrite(fid, hdr.dime.slice_duration(1), 'float32');
fwrite(fid, hdr.dime.toffset(1),    'float32');
fwrite(fid, hdr.dime.glmax(1),        'int32');
fwrite(fid, hdr.dime.glmin(1),        'int32');

fwrite(fid, padChar(hdr.hist.descrip, 80), 'uchar');
fwrite(fid, padChar(hdr.hist.aux_file, 24), 'uchar');
fwrite(fid, hdr.hist.qform_code,    'int16');
fwrite(fid, hdr.hist.sform_code,    'int16');
fwrite(fid, hdr.hist.quatern_b,   'float32');
fwrite(fid, hdr.hist.quatern_c,   'float32');
fwrite(fid, hdr.hist.quatern_d,   'float32');
fwrite(fid, hdr.hist.qoffset_x,   'float32');
fwrite(fid, hdr.hist.qoffset_y,   'float32');
fwrite(fid, hdr.hist.qoffset_z,   'float32');
fwrite(fid, hdr.hist.srow_x(1:4), 'float32');
fwrite(fid, hdr.hist.srow_y(1:4), 'float32');
fwrite(fid, hdr.hist.srow_z(1:4), 'float32');
fwrite(fid, padChar(hdr.hist.intent_name, 16), 'uchar');
fwrite(fid, padChar(hdr.hist.magic, 4), 'uchar');

if ~isequal(ftell(fid), 348), warning('Header size is not 348 bytes.'); end
return;  % save_nii_hdr

%% Subfunction: pad or chop char to correct length. Called by save_nii_hdr
function buf = padChar(ch, len)
len1 = length(ch);
if len1 >= len,  buf = ch(1:len);
else buf = [ch zeros(1, len-len1, 'uint8')];
end

%% Subfunction: simplified from Jimmy Shen's NIfTI toolbox
function nii = make_nii(img)
switch class(img)
    case 'uint8',       datatype = 2;
    case 'int16',       datatype = 4;
    case 'int32',       datatype = 8;
    case 'single',
        if isreal(img), datatype = 16;
        else            datatype = 32;
        end
    case 'double',
        if isreal(img), datatype = 64;
        else            datatype = 1792;
        end
    case 'int8',        datatype = 256;
    case 'uint16',      datatype = 512;
    case 'uint32',      datatype = 768;
    otherwise
        error('Datatype is not supported by make_nii.');
end

dims = size(img);
dims = [length(dims) dims ones(1,8)];
dims = dims(1:8);

hdr.hk.sizeof_hdr       = 348;			% must be 348!
hdr.hk.data_type        = '';
hdr.hk.db_name          = '';
hdr.hk.extents          = 0;
hdr.hk.session_error    = 0;
hdr.hk.regular          = 'r';
hdr.hk.dim_info         = 0;

hdr.dime.dim = dims;
hdr.dime.intent_p1 = 0;
hdr.dime.intent_p2 = 0;
hdr.dime.intent_p3 = 0;
hdr.dime.intent_code = 0;
hdr.dime.datatype = datatype;

mx = round(double(max(img(:))));
mn = round(double(min(img(:))));
switch datatype
    case 2,     img =  uint8(img);  bitpix = 8;
    case 4,     img =  int16(img);  bitpix = 16;
    case 8,     img =  int32(img);  bitpix = 32;
    case 16,    img = single(img);  bitpix = 32;
    case 32,    img = single(img);  bitpix = 64;
    case 64,    img = double(img);  bitpix = 64;
    case 128,   img =  uint8(img);  bitpix = 24;
    case 256,   img =   int8(img);  bitpix = 8;
    case 511,
        img = double(img);
        img = (img - mn)/(mx - mn);
        img = single(img);          bitpix = 96;
        mx = 1;
        mn = 0;
    case 512,   img = uint16(img);  bitpix = 16;
    case 768,   img = uint32(img);  bitpix = 32;
    case 1792,  img = double(img);  bitpix = 128;
    otherwise
        error('Datatype is not supported by make_nii.');
end

hdr.dime.bitpix = bitpix;
hdr.dime.slice_start = 0;
hdr.dime.pixdim = [0 ones(1,7)];
hdr.dime.vox_offset = 0;
hdr.dime.scl_slope = 0;
hdr.dime.scl_inter = 0;
hdr.dime.slice_end = 0;
hdr.dime.slice_code = 0;
hdr.dime.xyzt_units = 0;
hdr.dime.cal_max = 0;
hdr.dime.cal_min = 0;
hdr.dime.slice_duration = 0;
hdr.dime.toffset = 0;
hdr.dime.glmax = mx;
hdr.dime.glmin = mn;

hdr.hist.descrip = '';
hdr.hist.aux_file = 'none';
hdr.hist.qform_code = 0;
hdr.hist.sform_code = 0;
hdr.hist.quatern_b = 0;
hdr.hist.quatern_c = 0;
hdr.hist.quatern_d = 0;
hdr.hist.qoffset_x = 0;
hdr.hist.qoffset_y = 0;
hdr.hist.qoffset_z = 0;
hdr.hist.srow_x = zeros(1,4);
hdr.hist.srow_y = zeros(1,4);
hdr.hist.srow_z = zeros(1,4);
hdr.hist.intent_name = '';
hdr.hist.magic = '';
hdr.hist.originator = zeros(1,5);

nii.hdr = hdr;
nii.img = img;
