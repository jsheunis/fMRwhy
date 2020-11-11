function options = fmrwhy_bids_setupQcSubDirs(bids_dir, sub, options)

if isempty(options)
    options = struct;
end

% Derivatives directories
options.deriv_dir = fullfile(bids_dir, 'derivatives');
options.preproc_dir = fullfile(options.deriv_dir, 'fmrwhy-preproc');
options.qc_dir = fullfile(options.deriv_dir, 'fmrwhy-qc');

% Subject directories
options.sub_dir_preproc = fullfile(options.preproc_dir, ['sub-' sub]);
options.sub_dir_qc = fullfile(options.qc_dir, ['sub-' sub]);
options.sub_dir_BIDS = fullfile(bids_dir, ['sub-' sub]);

% Create and copy content, if necessary
if ~exist(options.sub_dir_preproc, 'dir')
    % Create new preprod deriv dir
    mkdir(options.sub_dir_preproc);
    % grab BIDS subject dir contents
    dir_contents = dir(options.sub_dir_BIDS);
    % Loop through files/folders
    % If session, loop through session folder contents, and:
    %   - if anat or func, loop through content and copy or gunzip (nii.gz)
    %   - if other, don't copy
    % If anat or func (i.e. no session), loop through folder contents, and:
    %   - copy or gunzip (nii.gz) all
    for i = 1:numel(dir_contents)
        % If sessions exist
        if (dir_contents(i).isdir == 1) && contains(dir_contents(i).name, 'ses')
            sesdir_contents = dir(fullfile(dir_contents(i).folder, dir_contents(i).name));
            % Loop through session dir
            for j = 1:numel(sesdir_contents)
                % if there are anat or func directories, copy/gunzip files accordingly
                if (sesdir_contents(j).isdir == 1) && (contains(sesdir_contents(j).name, 'anat') || contains(sesdir_contents(j).name, 'func'))
                    anatfunc_dir = fullfile(sesdir_contents(j).folder, sesdir_contents(j).name);
                    deriv_subdir = fullfile(options.sub_dir_preproc, dir_contents(i).name, sesdir_contents(j).name);
                    copyGunzipDirContents(anatfunc_dir, deriv_subdir);
                else
                    % Don't copy/gunzip other directories/files; only anat and func used for fmrwhy_bids_workflowQC;
                end
            end
        % If sessions DO NOT exist
        elseif (dir_contents(i).isdir == 1) && (contains(dir_contents(i).name, 'anat') || contains(dir_contents(i).name, 'func'))
            anatfunc_dir = fullfile(dir_contents(i).folder, dir_contents(i).name);
            deriv_subdir = fullfile(options.sub_dir_preproc, dir_contents(i).name);
            copyGunzipDirContents(anatfunc_dir, deriv_subdir);
        else
            % Ignore other file types
        end
    end
end


if ~exist(options.sub_dir_qc, 'dir')
    % create sub dir in qc derivatives dir
    mkdir(options.sub_dir_qc);
    dir_contents = dir(options.sub_dir_BIDS);
    
    for i = 1:numel(dir_contents)
        % If sessions exist
        if (dir_contents(i).isdir == 1) && contains(dir_contents(i).name, 'ses')
            % create session dir
            mkdir(fullfile(options.sub_dir_qc, dir_contents(i).name));
            % Loop through session dir
            for j = 1:numel(sesdir_contents)
                % if there are anat or func directories, copy/gunzip files accordingly
                if (sesdir_contents(j).isdir == 1) && (contains(sesdir_contents(j).name, 'anat') || contains(sesdir_contents(j).name, 'func'))
                    mkdir(fullfile(options.sub_dir_qc, dir_contents(i).name, sesdir_contents(j).name));
                end
            end
        % If sessions DO NOT exist
        elseif (dir_contents(i).isdir == 1) && (contains(dir_contents(i).name, 'anat') || contains(dir_contents(i).name, 'func'))
            mkdir(fullfile(options.sub_dir_qc, dir_contents(i).name));
        else
            % Ignore other file types
        end
    end
end

function copyGunzipDirContents(srcDir, destDir)
    % disp(['source: ' srcDir])
    % disp(['destination: ' destDir])
    srcDir_contents = dir(srcDir);
    % assignin('base', 'srcDir_contents', srcDir_contents);
    for i = 1:numel(srcDir_contents)
        srcDir_file = fullfile(srcDir_contents(i).folder, srcDir_contents(i).name);
        if contains(srcDir_contents(i).name, 'sub')
            if contains(srcDir_contents(i).name, 'nii.gz')
                gunzip(srcDir_file, destDir);
            else
                if ~exist(destDir, 'dir')
                    mkdir(destDir)
                end
                copyfile(srcDir_file, destDir);
            end
        end
    end
