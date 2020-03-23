function fmrwhy_util_gridPlot(img_fn, k)

% Simple function to

img_spm = spm_vol(img_fn);
img = spm_read_vols(img_spm);



tmp={};
for i=1:3
    tmp{i} = zeros(3,64);
end
j = 0;
for i = 1:64
    j=j+1;
    VoxelSubscriptsX = [i 1 1];
    VoxelSubscriptsY = [1 i 1];
    VoxelSubscriptsZ = [1 1 i];
    tempX = [VoxelSubscriptsX' ; 1];
    tempY = [VoxelSubscriptsY' ; 1];
    tempZ = [VoxelSubscriptsZ' ; 1];
    WorldSpaceCoordinatesX = img_spm.mat * tempX;
    WorldSpaceCoordinatesY = img_spm.mat * tempY;
    WorldSpaceCoordinatesZ = img_spm.mat * tempZ;
    WorldSpaceCoordinatesX = WorldSpaceCoordinatesX(1:3);
    WorldSpaceCoordinatesY = WorldSpaceCoordinatesY(1:3);
    WorldSpaceCoordinatesZ = WorldSpaceCoordinatesZ(1:3);
    tmp{1}(:,j) = WorldSpaceCoordinatesX;
    tmp{2}(:,j) = WorldSpaceCoordinatesY;
    tmp{3}(:,j) = WorldSpaceCoordinatesZ;
end


colordata = zeros(64,64);
for i = 1:64
    for j = 1:64
        colordata(i,j) = img(i,j,k);
    end
end
xx = tmp{1}(1,:);
yy = tmp{2}(2,:);

figure; surf(xx, yy, colordata)
view(2);
title('surf')
xlabel('X-coordinates')
ylabel('Y-coordinates')

figure; pcolor(xx, yy, colordata)
view(2);
title('pcolor')
xlabel('X-coordinates')
ylabel('Y-coordinates')






