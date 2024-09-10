% This code runs Example 2 from the paper

addpath('../../2D_Q_RT0')
addpath('../../Functions')
addpath('../../plotting')
addpath('../../nulls')
addpath('../../SparseNullOrth')

warning('on','all')


% !!!!Attention!!!!: currently only works if coarse domains all have the same size,
% would have to adjust computation of the coarse VRT spaces in
% buildCoarseSpace!


% Define parameters for the mesh
nx = 240;
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
%coeffs = ones(msh.dim,msh.nelem);
coeffs = getCoeffExample2(msh);

% Plott permeability coefficients
plotField(coeffs(1,:), msh, 'test', 'coeffs.vtk');

% Build fine matrix for mixed formulation
tic
[K, ~, ~] = buildFineMatrix(msh, coeffs);
toc
disp(['running time (building the fine matrix): ',num2str(toc)])

% Define the right hand side
f = zeros(msh.elengdof,1);
ftmp = zeros(msh.elengdof,1);

% We scale the source function by h to obtain the values in the basis
f((msh.nfaces+1):msh.elengdof) = - getSourceExample2(msh) * msh.hx;
ftmp((msh.nfaces+1):msh.elengdof) = - getSourceExample2(msh) * msh.hx;

% p = getSourceNeu(msh, alpha) /(2 * pi * pi * alpha * alpha);

% Define boundary conditions for normal flux
% Need to satisfy the compatibility condition

bc_dofs = msh.bfaces;
msh.bvals = zeros(1,length(msh.bfaces));

% apply the boundary conditions
[Kb, f] = apply_bc(K,f, [msh.bfaces; msh.bvals]);

% solve the system
% sol = Kb \ f;
sol = zeros(size(Kb,1), 1);
sol(1:end-1) = Kb(1:end-1,1:end-1) \ f(1:end-1);

% Split into velocity part and pressure part
uh = sol(1:msh.nfaces);
ph = sol((msh.nfaces+1):msh.elengdof) / msh.hx;

% normalize for neumann problems
ph = ph -mean(ph);



plotField(ph, msh, 'test', 'ph_fine.vtk'); 
plotField(f((msh.nfaces+1):msh.elengdof), msh, 'test', 'f.vtk'); 

% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Initialize everything for the multiscale method % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Initialize parameters
Nc_x = 6;           % Number of subdomains in x direction
Nc_y = Nc_x;        % Number of subdomains in y direction
overlap = 2;        % Number of layers used to construct overlapping domains    
extension = 8;      % Number of layers used to construct oversampling domains
nloc = 20;          % Number of eigenmodes used in each subdomain
gamma = 0.001;          % Parameter for augmented formulation
bool_RTv = 1;       % 1 if VMS contains coarse RT0, otherwise 0


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

% Build coarse space
tic
    [basis] = buildCoarseSpace(msh, coeffs, ftmp, nloc);
toc
disp(['running time (Building coarse space): ',num2str(toc)])

% Assemble coarse space
[K_aug, ~, ~] = buildFineMatrix_aug(msh, coeffs, gamma);

tic
[KH, fH, psol_v, psol_p, RHv, RHp] = assembleCoarseSpace(msh, basis, K_aug, ftmp, gamma, bool_RTv);
toc
disp(['running time (Assembling coarse space): ',num2str(toc)])

% Solve the coarse problem
solH = KH \ fH;

% Split into pressure and velocity part, add particular function and vector
% field to get multiscale solution.

if not(bool_RTv)
    basis{1}.totalSizeRTv = 0;            % Turn of coarse RT velocity
end
sizev = size(RHv,2);
sizep = size(RHp,2);
PH = psol_p + RHp * solH((sizev+1):(sizev+sizep));
PH = PH - mean(PH);
UH = psol_v + RHv * solH(1:(sizev));

% scale by h to go from basis values to true values
PH = PH / msh.hx;

% Print relative error
err = norm(PH-ph)/norm(ph);
sprintf('relative pressure error multiscale fine = %e', err)

diff_u = uh - UH;
err = sqrt((diff_u' * K(1:msh.nfaces, 1:msh.nfaces) * diff_u) / (uh' * K(1:msh.nfaces, 1:msh.nfaces) * uh) );
sprintf('relative velocity error multiscale fine = %e', err)



% Plot the pressure solution and error
plotField(PH, msh, 'test', 'PH.vtk');
plotField(abs(PH-ph), msh, 'test', 'abs(PH-ph).vtk');

% Calculate velocity x component and error in x component
uh_X = (- 1/2 * uh(msh.elems2faces(1,:)) - 1/2 * uh(msh.elems2faces(3,:))) /msh.hx;
UH_X = (- 1/2 * UH(msh.elems2faces(1,:)) - 1/2 * UH(msh.elems2faces(3,:))) / msh.hx;

% Plot velocity x component and error in x component
plotField(UH_X, msh, 'UH_X', 'UH_X.vtk');
plotField(uh_X, msh, 'uh_X', 'uh_fine_X.vtk');
plotField(abs(UH_X-uh_X), msh, 'UH_X_err', 'abs(UH_X-uh_X).vtk');

% Calculate velocity y component and error in x component
uh_Y = ( 1/2 * uh(msh.elems2faces(2,:)) + 1/2 * uh(msh.elems2faces(4,:))) /msh.hx;
UH_Y = ( 1/2 * UH(msh.elems2faces(2,:)) + 1/2 * UH(msh.elems2faces(4,:))) / msh.hx;

% Plot velocity y component and error in x component
plotField(UH_Y, msh, 'UH_Y', 'UH_Y.vtk');
plotField(uh_Y, msh, 'uh_Y', 'uh_fine_Y.vtk');
plotField(abs(UH_Y-uh_Y), msh, 'UH_Y_err', 'abs(UH_Y-uh_Y).vtk');

% Save setup
% save('Visualization/setup.mat', 'nx' , 'ny', 'coeffs', 'rhs', 'bvals', 'Nc_x', 'Nc_y', 'extension', 'nloc', 'gamma');

% for i = 1:20
% u = msh.Emg{1}.Rv' * basis{1}.Ven(:,i);
% u_X = (- 1/2 * u(msh.elems2faces(1,:)) - 1/2 * u(msh.elems2faces(3,:))) /msh.hx;
% plotField(u_X, msh, 'u_x', strcat('u_X_',num2str(i),'.vtk'));
% end

