fprintf('\nRecPC: Reconstruction from Partial Circulant Data, ver 1.0\n');
fprintf('Wotao Yin (Rice) and Junfeng Yang (Nanjing)\n');
fprintf('Compile some MEX files ...');
    cd utilities;
    mex -O Compute_Ux_Uy.c;
    mex -O Compute_Wx_Wy.c;
    mex -O Compute_rhs_DxtU_DytU.c;
    cd ..;
fprintf('Done!\n');

fprintf('Add subfolders to path ...');
    pathstr=fileparts(mfilename('fullpath'));
    addpath(genpath(pathstr));
fprintf('Done!\n');

