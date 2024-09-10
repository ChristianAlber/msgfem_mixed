addpath('2D_Q_RT0')
addpath('Functions')
addpath('plotting')
addpath('nulls')
addpath('SparseNullOrth')

% Note that we use the functions 1/h as basis functions for the pressure on
% cells with sidelength h. Hence on a cell Ki, the function value f and the
% coefficient fi in the basis are related by fi = f * h. Hence, we multiply
% the source term by the mesh size for computations. Converting back from
% the computed pressure basis coefficients pi to the real pressure values,
% we have: p = pi / h. Converting back from the basis coefficients of the
% velocity ui, we have to divide by h again, because the dof is given by an
% integral over the face: 
% ui = int_{E_i}u*n ~~ u * n * h  -> u * n ~~ ui / h


warning('off','all')


% !!!!Attention!!!!: currently only works if coarse domains all have the same size,
% would have to adjust computation of the coarse VRT spaces in
% buildCoarseSpace!

% Define parameter for rhs
alpha = 1;

% Define parameters for the mesh
nx = 150;
ny = nx;
hx = 1/nx;
hy = hx;

% build the fine mesh
tic
msh = buildMesh(hx,hy,nx,ny);
toc
disp(['running time (building the mesh): ',num2str(toc)])

% test
% signs = signs_edges_Q(msh.elems2nodes');

% Define coeffs of the matrix A: Seems like code only supports diagonal
% matrix A, hence two entries.
coeffs = ones(msh.dim,msh.nelem);

% Build fine matrix for mixed formulation
tic
[K, ~, ~] = buildFineMatrix(msh, coeffs);
toc
disp(['running time (building the fine matrix): ',num2str(toc)])

% % % Tests for entries of matrix (condition)
% % disp(max(abs(K(1:msh.nfaces,1:msh.nfaces)), [], 'all'))
% % disp(max(abs(K(msh.nfaces+1:msh.elengdof, 1:msh.nfaces)), [], 'all'))
% % disp(1/condest(K))

% Define the right hand side
f = zeros(msh.elengdof,1);
ftmp = zeros(msh.elengdof,1);
% f(msh.nfaces+1) = 1;
% f(msh.nfaces+2) = 1;
% f(msh.nfaces+3) = -1;
% f(msh.elengdof) = -1;

% We scale the source function by h to obtain the values in the basis
f((msh.nfaces+1):msh.elengdof) = - getSourceTest(msh) * msh.hx;
ftmp((msh.nfaces+1):msh.elengdof) = - getSourceTest(msh) * msh.hx;

p = getSourceTest(msh) /(pi * pi);
% p = getSourceTest(msh);

% Define boundary conditions for normal flux
% Need to satisfy the compatibility condition

% bc_dofs = msh.bfaces;
% msh.bvals = zeros(1,length(msh.bfaces));

right_bdofs = find(msh.face_midpoints(1,:)>1-1e-6);
left_bdofs = find(msh.face_midpoints(1,:)<1e-6);
top_bdofs = find(msh.face_midpoints(2,:)>1-1e-6);
bottom_bdofs = find(msh.face_midpoints(2,:)<1e-6);

right_vals =  pi *ones(1,length(right_bdofs)) *msh.hx *0.5;
left_vals =  - pi *ones(1,length(left_bdofs)) *msh.hx * 0.5;
% right_vals = - ones(1,length(right_bdofs)) *msh.hx *0.5;
% left_vals = - ones(1,length(left_bdofs)) *msh.hx * 0.5;
top_vals = zeros(1,length(top_bdofs));
bottom_vals = zeros(1,length(bottom_bdofs));

bc_dofs = [right_bdofs left_bdofs top_bdofs bottom_bdofs];
msh.bvals = [right_vals left_vals top_vals bottom_vals];



% FREE_tmp = setdiff(msh.elengdof, msh.bfaces);
% disp(1/condest(K(FREE_tmp, FREE_tmp)));

% apply the boundary conditions
[Kb, f] = apply_bc(K,f, [bc_dofs; msh.bvals]);
condest(K)
condest(Kb)
% solve the system
sol = Kb \ f;


% Split into velocity part and pressure part
uh = sol(1:msh.nfaces);
ph = sol((msh.nfaces+1):msh.elengdof) / msh.hx;

% normalize for neumann problems
p = p - mean(p);
ph = ph -mean(ph);

% visualize pressure
pressure = reshape(ph, [msh.nx,msh.ny])';

plotField(p, msh, 'test', 'p.vtk'); 
plotField(ph, msh, 'test', 'ph.vtk'); 
% plotField(f((msh.nfaces+1):msh.elengdof), msh, 'test', 'f.vtk'); 

% Sanity check
err = norm(p-ph)/norm(p);
sprintf('relative error true fine = %e', err)

% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Initialize everything for the multiscale method % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % Initialize parameters
% Nc_x = 2;           % Number of subdomains in x direction
% Nc_y = Nc_x;        % Number of subdomains in y direction
% overlap = 2;        % Number of layers used to construct overlapping domains    
% extension = 5;      % Number of layers used to construct oversampling domains
% nloc = 10;          % Number of eigenmodes used in each subdomain
% gamma = 0.0001;          % Parameter for augmented formulation
% bool_RTv = 1;       % 1 if VMS contains coarse RT0, otherwise 0
% 
% 
% % build the overlapping subdomains
% tic
%     msh = overlapDomains(msh, Nc_x, Nc_y, overlap, extension, 1);
% toc
% disp(['running time (Constructing overlapping Subdomains): ',num2str(toc)])
% 
% % build the partition of unity
% tic
%     msh =  partitionOfUnity(msh);
% toc
% disp(['running time (Constructing partition of unity): ',num2str(toc)])
% 
% % Build coarse space
% tic
%     [basis] = buildCoarseSpace(msh, coeffs, ftmp, nloc);
% toc
% disp(['running time (Building coarse space): ',num2str(toc)])
% 
% % Assemble coarse space
% [K_aug, ~, ~] = buildFineMatrix_aug(msh, coeffs, gamma);
% 
% tic
% [KH, fH, psol_v, psol_p, RHv, RHp] = assembleCoarseSpace(msh, basis, K_aug, ftmp, gamma, bool_RTv);
% toc
% disp(['running time (Assembling coarse space): ',num2str(toc)])
% 
% % Solve the coarse problem
% solH = KH \ fH;
% 
% % Split into pressure and velocity part, add particular function and vector
% % field to get multiscale solution.
% 
% if not(bool_RTv)
%     basis{1}.totalSizeRTv = 0;            % Turn of coarse RT velocity
% end
% sizev = size(RHv,2);
% sizep = size(RHp,2);
% PH = psol_p + RHp * solH((sizev+1):(sizev+sizep));
% PH = PH - mean(PH);
% UH = psol_v + RHv * solH(1:(sizev));
% 
% % scale by h to go from basis values to true values
% PH = PH / msh.hx;
% 
% % Print relative error
% p_err = norm(PH-ph)/norm(ph);
% sprintf('relative pressure error multiscale fine = %e', p_err)
% 
% diff_u = uh - UH;
% v_err = sqrt((diff_u' * K(1:msh.nfaces, 1:msh.nfaces) * diff_u) / (uh' * K(1:msh.nfaces, 1:msh.nfaces) * uh) );
% sprintf('relative velocity error multiscale fine = %e', v_err)
% 
% % Plot the pressure solution and error
% plotField(PH, msh, 'test', 'PH.vtk');
% plotField(abs(PH-ph), msh, 'test', 'abs(PH-ph).vtk');
% 
% % Calculate velocity x component and error in x component
uh_X = (- 1/2 * uh(msh.elems2faces(1,:)) - 1/2 * uh(msh.elems2faces(3,:))) /msh.hx;
% UH_X = (- 1/2 * UH(msh.elems2faces(1,:)) - 1/2 * UH(msh.elems2faces(3,:))) / msh.hx;
% 
% % Plot velocity x component and error in x component
plotField(uh_X, msh, 'uh_X', 'uh_X.vtk');
plotVectorField(uh, msh, 'vf_fine.vtk');
% plotField(UH_X, msh, 'UH_X', 'UH_X.vtk');
% plotField(abs(UH_X-uh_X), msh, 'UH_X_err', 'abs(UH_X-uh_X).vtk');
% 
% % Calculate velocity y component and error in x component
% uh_Y = ( 1/2 * uh(msh.elems2faces(2,:)) + 1/2 * uh(msh.elems2faces(4,:))) /msh.hx;
% UH_Y = ( 1/2 * UH(msh.elems2faces(2,:)) + 1/2 * UH(msh.elems2faces(4,:))) / msh.hx;
% 
% % Plot velocity y component and error in x component
% plotField(UH_Y, msh, 'UH_Y', 'UH_Y.vtk');
% plotField(abs(UH_Y-uh_Y), msh, 'UH_Y_err', 'abs(UH_Y-uh_Y).vtk');
% 
