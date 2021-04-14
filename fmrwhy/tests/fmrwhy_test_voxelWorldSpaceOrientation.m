tmp = {};
for i = 1:3
    tmp{i} = zeros(3, 64);
end
j = 0;
for i = 1:64
    j = j + 1;
    VoxelSubscriptsX = [i 1 1];
    VoxelSubscriptsY = [1 i 1];
    VoxelSubscriptsZ = [1 1 i];
    tempX = [VoxelSubscriptsX'; 1];
    tempY = [VoxelSubscriptsY'; 1];
    tempZ = [VoxelSubscriptsZ'; 1];
    WorldSpaceCoordinatesX = template_spm.mat * tempX;
    WorldSpaceCoordinatesY = template_spm.mat * tempY;
    WorldSpaceCoordinatesZ = template_spm.mat * tempZ;
    WorldSpaceCoordinatesX = WorldSpaceCoordinatesX(1:3);
    WorldSpaceCoordinatesY = WorldSpaceCoordinatesY(1:3);
    WorldSpaceCoordinatesZ = WorldSpaceCoordinatesZ(1:3);
    tmp{1}(:, j) = WorldSpaceCoordinatesX;
    tmp{2}(:, j) = WorldSpaceCoordinatesY;
    tmp{3}(:, j) = WorldSpaceCoordinatesZ;
end

k = 1;
colordata = zeros(64, 64);
for i = 1:64
    for j = 1:64
        colordata(i, j) = img(i, j, k);
    end
end

xx = tmp{1}(1, j);
yy = tmp{2}(2, j);
