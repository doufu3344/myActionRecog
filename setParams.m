% data information
info.nact = 3;
info.nsbj = 1;
info.ntms = 1;
info.train = [1, 3, 5, 7, 9];
info.test = [2, 4, 6, 8, 10];
info.vidpath = 'E:\datasets\02 MSR Action 3D Dataset\';
info.dstippath = 'DSTIP\';
info.normalpath = 'Normals\';
info.featurespath = 'Features\';
info.coefficients = 'Coefficients\';

% dstip information
stip.npoints = 500; % number of intrest points
stip.gauss_size = 5; % size of Gaussian lowpass filter
stip.gauss_variance = 1; % Gauss variance, sigma
stip.gabor_scale = 10; % temproal scale of Gabor filter, tao
stip.gabor_Ormig = 0.6/stip.gabor_scale; % w
stip.noise_thresh = 1;  % larger value removes more noise (too large value may hurt fast movement)
stip.Noisesup_all = 1;  % noise suppressor, 1-all frames

%stip.boarder=15;  % leave enough space to extract a cuboid around the DSTIPs
%stip.noiseSuppressor=1; %if apply the noise suppressor: set to 0 when there is no backgroud for faster process.

cuboid.spatioscale = 600;
cuboid.xpoints = 3;
cuboid.ypoints = 3;
cuboid.temporalscale = 3;

sparse.maxiteration = 100;
sparse.threshold = 500000;
sparse.dictdimen = [108 100];
sparse.lamda = 0.01;
sparse.yita = 0.01;
sparse.miu = 0.01;
sparse.omiga = [7 15 20];