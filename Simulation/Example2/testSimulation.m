addpath('../../2D_Q_RT0')
addpath('../../Functions')
addpath('../../plotting')
addpath('../../nulls')
addpath('../../SparseNullOrth')
addpath('../../data')

warning('on','all')

% Define all relevant parameters for which we want to run the Simulation
nx = 240;

% Define parameters for the mesh
ny = nx;
hx = 1/nx;
hy = hx;

% build the fine mesh
tic
msh = buildMesh(hx,hy,nx,ny);
toc
disp(['running time (building the mesh): ',num2str(toc)])

coeffs = getCoeffExample2(msh);
rhs = getSourceExample2(msh);
bvals = zeros(1,length(msh.bfaces));
Nc_xs = 6;
Nc_ys = Nc_xs;
extensions = 4:4:16;
nlocs = 2:2:20;
gammas = 1/max(coeffs,[],'all');
% Run the simulation
filename = '../../data/Example2_Ven_l_nloc.mat';
[v_err,p_err] = runSimulation(msh, coeffs, rhs, bvals, Nc_xs, Nc_ys, extensions, nlocs, gammas, filename);

