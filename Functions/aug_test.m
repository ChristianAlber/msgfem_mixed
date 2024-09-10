addpath('2D_Q_RT0')
addpath('Functions')
addpath('plotting')



% Note that we use the functions 1/h as basis functions for the pressure on
% cells with sidelength h. Hence on a cell Ki, the function value f and the
% coefficient fi in the basis are related by fi = f * h. Hence, we multiply
% the source term by the mesh size for computations. Converting back from
% the computed pressure basis coefficients pi to the real pressure values,
% we have: p = pi / h. Converting back from the basis coefficients of the
% velocity ui, we have to divide by h again, because the dof is given by an
% integral over the face: 
% ui = int_{E_i}u*n ~~ u * n * h  -> u * n ~~ ui / h


% Define parameters for the mesh
hx = 0.1;
hy = 0.1;
nx = 10;
ny = 10;

% Define parameter for augmented form
gamma = 0;

% build the fine mesh
tic
msh = buildMesh(hx,hy,nx,ny);
toc
disp(['running time (building the mesh): ',num2str(toc)])

% Define coeffs of the matrix A: Seems like code only supports diagonal
% matrix A, hence two entries.
coeffs = ones(msh.dim,msh.nelem);

% Build fine matrix for mixed formulation
tic
[K, Kloc, ~, Aug] = buildFineMatrix_aug(msh, coeffs, gamma);
toc
disp(['running time (building the fine matrix): ',num2str(toc)])

% Define the right hand side
f = zeros(msh.elengdof,1);
% f(msh.nfaces+1) = 1;
% f(msh.nfaces+2) = 1;
% f(msh.nfaces+3) = -1;
% f(msh.elengdof) = -1;

% We scale the source function by h to obtain the values in the basis
f((msh.nfaces+1):msh.elengdof) = - getSourceNeu(msh) * msh.hx;

p = getSourceNeu(msh) /(2 * pi * pi);
% Define boundary conditions for normal flux
bc_dofs = msh.bfaces;

% Need to satisfy the compatibility condition
bc_vals = zeros(1,length(bc_dofs));

% Add augmented form
f(1:msh.nfaces) = f(1:msh.nfaces) + gamma * K(1:msh.nfaces,(msh.nfaces+1):msh.elengdof) * f((msh.nfaces+1):msh.elengdof);

% apply the boundary conditions
[K, f] = apply_bc(K,f, [bc_dofs; bc_vals]);
% change RHS due to augmented term


% solve the system
sol = K \ f;

% scale by h to go from basis values to true values
sol = sol / msh.hx;

% Split into velocity part and pressure part
uh = sol(1:msh.nfaces);
ph = sol((msh.nfaces+1):msh.elengdof);

% normalize for neumann problems
p = p - mean(p);
ph = ph -mean(ph);


% visualize pressure
pressure = reshape(ph, [msh.nx,msh.ny])';

%ph = ph - mean(ph);
plotField(ph, msh, 'test', 'pressure.vtk'); 
plotField(f((msh.nfaces+1):msh.elengdof), msh, 'test', 'f.vtk'); 

% Sanity check
err = sum(abs(p-ph))/sum(abs(p));
sprintf('relative error = %e', err)

p_tmp = p./ph;


% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Initialize everything for the multiscale method % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Initialize parameters
Nc_x = 2;
Nc_y = Nc_x;
overlap = 2;
extension = 1;

% build the overlapping subdomains
tic
    msh = overlapDomains(msh, Nc_x, Nc_y, overlap, extension, 1);
toc
disp(['running time (Constructing overlapping Subdomains): ',num2str(toc)])

% build the partition of unity
tic
    msh =  partitionOfUnity(msh);
toc
disp(['running time (Constructing partition of unity): ',num2str(toc)])


% Testing partition of unity for pressure:
u = ones(length(msh.Emg{1}.elements),1);
x = msh.Emg{1}.Xp * u;
% x = reshape(x,[28,28]);


% build coarse space
tic
    basis = buildCoarseSpace(msh, coeffs, f);
toc
disp(['running time (Building coarse space): ',num2str(toc)])

% test
x = basis{1}.par_p;
% x = reshape(x,[8,8]);
