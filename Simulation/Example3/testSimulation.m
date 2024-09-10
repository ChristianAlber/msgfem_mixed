addpath('../../2D_Q_RT0')
addpath('../../Functions')
addpath('../../plotting')
addpath('../../nulls')
addpath('../../SparseNullOrth')
addpath('../../data')

warning('on','all')

% Define all relevant parameters for which we want to run the Simulation
fine_fac = 4;
nx = fine_fac * 60;

% Define parameters for the mesh
ny = nx;
hx = 1/nx;
hy = hx;

% build the fine mesh
tic
msh = buildMesh(hx,hy,nx,ny);
toc
disp(['running time (building the mesh): ',num2str(toc)])

coeffs = getCoeffExample3();

coeffs = coeffs(1,:);

coeffs_tmp = reshape(coeffs,[nx/fine_fac,ny/fine_fac]);
coeffs_inter = repelem(coeffs_tmp, fine_fac,fine_fac);
coeffs_inter = reshape(coeffs_inter,[nx*ny,1]);
coeffs = [coeffs_inter, coeffs_inter]';

% rhs = getSourceNeu(msh);
rhs = getSourceExample3(msh);
bvals = zeros(1,length(msh.bfaces));
Nc_xs = 6;
Nc_ys = Nc_xs;
extensions = 5:5:15;
nlocs = 2:4:18;
gammas = 1/max(coeffs,[],'all');
%gammas = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3];

% Run the simulation
filename = '../../data/Example3_nloc_l_venr_test.mat';
[v_err,p_err] = runSimulation(msh, coeffs, rhs, bvals, Nc_xs, Nc_ys, extensions, nlocs, gammas, filename);

