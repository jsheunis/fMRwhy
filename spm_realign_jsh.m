function output = spm_realign_jsh(ref_fn, functional4D_fn)

func_spm = spm_vol(functional4D_fn);
Ndyn = numel(func_spm);
% Filenames to realign
fns={};

if ref_fn == 0
    for i = 1:Ndyn
        fns{i} = [functional4D_fn ',' num2str(i) ];
    end
else
    fns{1} = [ref_fn ',1'];
    for i = 1:Ndyn
        fns{i+1} = [functional4D_fn ',' num2str(i) ];
    end
end

% Data
realign_estimate_reslice = struct;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.data={fns'};
% Eoptions
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % register to first
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
% Roptions
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [1 0]; % images [2..n]
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
% Run
cfg_util('run',realign_estimate_reslice.matlabbatch);
[d, f, e] = fileparts(functional4D_fn);
[dref, fref, eref] = fileparts(ref_fn);
% Output
output = struct;
output.rfunctional_fn = [d filesep 'r' f e];
if ref_fn == 0
    output.mp_fn = [d filesep 'rp_' f '.txt'];
else
    output.mp_fn = [dref filesep 'rp_' fref '.txt'];
end
output.MP = load(output.mp_fn);