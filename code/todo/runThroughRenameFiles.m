files_dir = '/Users/jheunis/Documents/PYTHON/drin/sample_data';


for i = 5:140
    for j = 1:3
        source = [files_dir filesep  'fMRI_tapping_d' num2str(i) '_e' num2str(j) '.nii'];
        destination = [files_dir filesep  'fMRI_tapping_d' sprintf('%03d', i) '_e' num2str(j) '.nii'];
        movefile(source, destination)
    end
    
end

