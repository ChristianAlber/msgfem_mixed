addpath('2D_Q_RT0')
addpath('Functions')
addpath('plotting')
addpath('Simulation')
addpath('nulls')
addpath('SparseNullOrth')

warning('on','all')

% Define all relevant parameters for which we want to run the Simulation
nx = 60;

% Define parameters for the mesh
ny = 180;
hx = 1/nx;
hy = hx;

% build the fine mesh
tic
msh = buildMesh(hx,hy,nx,ny);
toc
disp(['running time (building the mesh): ',num2str(toc)])

% coeffs = ones(msh.dim,msh.nelem);
coeffs = getCoeffExample3();
% rhs = getSourceNeu(msh);
rhs = getSourceExample3(msh);
bvals = zeros(1,length(msh.bfaces));
Nc_xs = 4;
Nc_ys = 3 * Nc_xs;
extensions = 5;
nlocs = 5:5:15;
% gammas = 1/max(coeffs,[],'all');
gammas = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3];

% Run the simulation
filename = 'data/Example3_gamma_test.mat';
[v_err,p_err] = runSimulation(msh, coeffs, rhs, bvals, Nc_xs, Nc_ys, extensions, nlocs, gammas, filename);

