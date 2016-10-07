global constK Box datafile detRatio_dense sideLength max_dist_sqr max_Acut
global indicator dim reinitialize datafile1 mat_L mat_U
global array_m array_u array_detRatio globalCount

%Initialization
globalCount = 0; love = 0;
dim = 3; N = 2*8^3; indicator = ones(1, N); totalSteps = 2*N;
array_m = [];  array_u = [];  array_detRatio = [];
sigma = 0.135*sideLength;
% ParticleConfig = Pconfig; OrbitalConfig = Oconfig;

[ParticleConfig, OrbitalConfig, Box, sigma, constK, sideLength, max_dist_sqr, max_Acut] = read_configuration (N);
A = compute_matrix (OrbitalConfig, OrbitalConfig);
Ainitial = A;