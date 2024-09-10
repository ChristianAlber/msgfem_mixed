% This code produces the data in Figure 1.

addpath('../../2D_Q_RT0')
addpath('../../Functions')
addpath('../../plotting')
addpath('../../nulls')
addpath('../../SparseNullOrth')
addpath('../../data')

warning('on','all')

% Define all relevant parameters for which we want to run the Simulation
nx = 100;

% Define parameters for the mesh
ny = nx;
hx = 1/nx;
hy = hx;

% build the fine mesh
tic
msh = buildMesh(hx,hy,nx,ny);
toc
disp(['running time (building the mesh): ',num2str(toc)])

coeffs = ones(msh.dim,msh.nelem);
rhs = getSourceNeu(msh,1);
bvals = zeros(1,length(msh.bfaces));
Nc_xs = 2;
Nc_ys = Nc_xs;
extensions = 2:2:8;
nlocs = 4:4:20;
gammas = 1;

% Run the simulation
filename = '../../data/testrun.mat';
[v_err,p_err] = runSimulation(msh, coeffs, rhs, bvals, Nc_xs, Nc_ys, extensions, nlocs, gammas, filename);

