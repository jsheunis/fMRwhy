function spm_createNII_jsh(template_spm, img, fn, descrip)
new_spm = template_spm;
new_spm.fname = fn;
new_spm.private.dat.fname = fn;
new_spm.descrip = descrip;
new_spm.private.dat.dim = new_spm(1).private.dat.dim(1:3);
new_spm.n = [1 1];
new_spm.pinfo(1) = 1; % https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;12fa60a.1205
spm_write_vol(new_spm,img);