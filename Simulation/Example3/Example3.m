addpath('../../2D_Q_RT0')
addpath('../../Functions')
addpath('../../plotting')
addpath('../../nulls')
addpath('../../SparseNullOrth')

warning('on','all')


% !!!!Attention!!!!: currently only works if coarse domains all have the same size,
% would have to adjust computation of the coarse VRT spaces in
% buildCoarseSpace!

% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Computation of a reference solution on finer mesh %
% % % % % % % % % % % % % % % % % % % % % % % % % % % 

finer_fac = 8;                              % We use finer_fac*60 fine element in x direction for the overkill solution
fine_fac = finer_fac / 2;
% Define parameters for the mesh
nx = finer_fac * 60;
ny = nx;
hx = 1/nx;
hy = hx;

% build the fine mesh
tic
msh_finer = buildMesh(hx,hy,nx,ny);
toc
disp(['running time (building the mesh): ',num2str(toc)])

% test
% signs = signs_edges_Q(msh.elems2nodes');

% Define coeffs of the matrix A: Seems like code only supports diagonal
% matrix A, hence two entries.
%coeffs = ones(msh.dim,msh.nelem);
coeffs = getCoeffExample3();
coeffs = coeffs(1,:);

coeffs_tmp = reshape(coeffs,[nx/finer_fac,ny/finer_fac]);
coeffs_inter = repelem(coeffs_tmp, finer_fac,finer_fac);
coeffs_inter = reshape(coeffs_inter,[nx*ny,1]);
coeffs = [coeffs_inter, coeffs_inter]';
coeffs_finer = coeffs;

% Build fine matrix for mixed formulation
tic
[K_finer, ~, ~] = buildFineMatrix(msh_finer, coeffs);
toc
disp(['running time (building the fine matrix): ',num2str(toc)])

% Define the right hand side
f = zeros(msh_finer.elengdof,1);
ftmp = zeros(msh_finer.elengdof,1);

% We scale the source function by h to obtain the values in the basis
f((msh_finer.nfaces+1):msh_finer.elengdof) = - getSourceExample3(msh_finer) * msh_finer.hx;
ftmp((msh_finer.nfaces+1):msh_finer.elengdof) = - getSourceExample3(msh_finer) * msh_finer.hx;

% Define boundary conditions for normal flux
% Need to satisfy the compatibility condition

bc_dofs = msh_finer.bfaces;
msh_finer.bvals = zeros(1,length(msh_finer.bfaces));

% apply the boundary conditions
[Kb, f] = apply_bc(K_finer,f, [msh_finer.bfaces; msh_finer.bvals]);

% solve the system
% sol = Kb \ f;
sol = zeros(size(Kb,1), 1);
sol(1:end-1) = Kb(1:end-1,1:end-1) \ f(1:end-1);

% Split into velocity part and pressure part
u = sol(1:msh_finer.nfaces);
p = sol((msh_finer.nfaces+1):msh_finer.elengdof) / msh_finer.hx;

u_X = (- 1/2 * u(msh_finer.elems2faces(1,:)) - 1/2 * u(msh_finer.elems2faces(3,:))) /msh_finer.hx;
u_Y = ( 1/2 * u(msh_finer.elems2faces(2,:)) + 1/2 * u(msh_finer.elems2faces(4,:))) /msh_finer.hx;
% normalize for neumann problems
p = p -mean(p);

% M = reshape(p, [nx,ny]);
% A = conv2(M,ones(2,2),'valid')/4;
% A = A(1:2:end,1:2:end);
% p = reshape(A, [nx*ny /4,1]);



% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % End of reference sol computation on finer mesh  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Define parameters for the mesh
nx = fine_fac * 60;
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
coeffs = getCoeffExample3();

coeffs = coeffs(1,:);

coeffs_tmp = reshape(coeffs,[nx/fine_fac,ny/fine_fac]);
coeffs_inter = repelem(coeffs_tmp, fine_fac,fine_fac);
coeffs_inter = reshape(coeffs_inter,[nx*ny,1]);
coeffs = [coeffs_inter, coeffs_inter]';

% Plott permeability coefficients
plotField(log10(coeffs(1,:)), msh, 'test', 'log10(coeffs).vtk');

% Build fine matrix for mixed formulation
tic
[K, ~, ~] = buildFineMatrix(msh, coeffs);
toc
disp(['running time (building the fine matrix): ',num2str(toc)])

% Define the right hand side
f = zeros(msh.elengdof,1);
ftmp = zeros(msh.elengdof,1);

% We scale the source function by h to obtain the values in the basis
f((msh.nfaces+1):msh.elengdof) = - getSourceExample3(msh) * msh.hx;
ftmp((msh.nfaces+1):msh.elengdof) = - getSourceExample3(msh) * msh.hx;

% p = getSourceNeu(msh, alpha) /(2 * pi * pi * alpha * alpha);

% Define boundary conditions for normal flux
% Need to satisfy the compatibility condition

bc_dofs = msh.bfaces;
msh.bvals = zeros(1,length(msh.bfaces));

% apply the boundary conditions
[Kb, f] = apply_bc(K,f, [msh.bfaces; msh.bvals]);

% solve the system
sol = Kb \ f;

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
extension = 15;      % Number of layers used to construct oversampling domains
nloc = 6;          % Number of eigenmodes used in each subdomain
gamma = 1/max(coeffs,[],'all');          % Parameter for augmented formulation
%  gamma = 0;
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
plotField(abs(UH_Y-uh_Y), msh, 'UH_Y_err', 'abs(UH_Y-uh_Y).vtk');

% plot velocity as vector field
plotVectorField(uh, msh, 'vf_fine.vtk');
plotVectorField(UH, msh, 'vf_coarse.vtk');


% Interpolate pressure solutions to the finer mesh
p_tmp = reshape(ph,[nx,ny]);
ph_inter = repelem(p_tmp, 2,2);
ph_inter = reshape(ph_inter,[4 * nx*ny,1]);
p_tmp = reshape(PH,[nx,ny]);
PH_inter = repelem(p_tmp, 2,2);
PH_inter = reshape(PH_inter,[4 * nx*ny,1]);
p_err_finer_fine = norm(p-ph_inter)/norm(p);
sprintf('relative pressure error finer fine = %e', norm(p-ph_inter)/norm(p))
p_err_finer_multi = norm(p-PH_inter)/norm(p);
sprintf('relative pressure error finer multiscale = %e', norm(p-PH_inter)/norm(p))
p_err_fine_multi = norm(ph_inter-PH_inter)/norm(ph_inter);
sprintf('relative pressure error fine multiscale = %e', norm(ph_inter-PH_inter)/norm(ph_inter))

% Interpolate velocity solutions to the finer mesh
uh_X_inter = reshape(repelem(reshape(uh_X,[nx,ny]), 2,2),[4 * nx*ny,1]);
uh_Y_inter = reshape(repelem(reshape(uh_Y,[nx,ny]), 2,2),[4 * nx*ny,1]);
UH_X_inter = reshape(repelem(reshape(UH_X,[nx,ny]), 2,2),[4 * nx*ny,1]);
UH_Y_inter = reshape(repelem(reshape(UH_Y,[nx,ny]), 2,2),[4 * nx*ny,1]);

% Rough estimates for the velocity error. 
diff = ((uh_X_inter - u_X) ./ coeffs_finer(1,:)')' * (uh_X_inter - u_X) + ((uh_Y_inter - u_Y) ./ coeffs_finer(2,:)')' * (uh_Y_inter - u_Y);
normalize = (u_X ./ coeffs_finer(1,:)')' * u_X + (u_Y ./ coeffs_finer(2,:)')' * u_Y;
err = sqrt(diff/normalize);
sprintf('estimate relative velocity error finer fine = %e', err)
diff = ((UH_X_inter - u_X) ./ coeffs_finer(1,:)')' * (UH_X_inter - u_X) + ((UH_Y_inter - u_Y) ./ coeffs_finer(2,:)')' * (UH_Y_inter - u_Y);
normalize = (u_X ./ coeffs_finer(1,:)')' * u_X + (u_Y ./ coeffs_finer(2,:)')' * u_Y;
err = sqrt(diff/normalize);
sprintf('estimate relative velocity error finer multiscale = %e', err)
% diff = ((uh_X_inter - UH_X_inter) ./ coeffs_finer(1,:)')' * (uh_X_inter - UH_X_inter) + ((uh_Y_inter - UH_Y_inter) ./ coeffs_finer(2,:)')' * (uh_Y_inter - UH_Y_inter);
% normalize = (uh_X_inter ./ coeffs_finer(1,:)')' * uh_X_inter + (uh_Y_inter ./ coeffs_finer(2,:)')' * uh_Y_inter;
% err = sqrt(diff/normalize);
% sprintf('estimate relative velocity error fine multiscale = %e', err)


% Velocity errors
% Interpolate onto finer mesh
uh_inter = interpolation(uh, msh, msh_finer);
UH_inter = interpolation(UH, msh, msh_finer);

diff_u = u - uh_inter;
err = sqrt((diff_u' * K_finer(1:msh_finer.nfaces, 1:msh_finer.nfaces) * diff_u) / (u' * K_finer(1:msh_finer.nfaces, 1:msh_finer.nfaces) * u) );
sprintf('relative velocity error finer fine = %e', err)

diff_u = u - UH_inter;
err = sqrt((diff_u' * K_finer(1:msh_finer.nfaces, 1:msh_finer.nfaces) * diff_u) / (u' * K_finer(1:msh_finer.nfaces, 1:msh_finer.nfaces) * u) );
sprintf('relative velocity error finer multiscale = %e', err)


% Save setup
% save('LaTeX/Visualization/Example3/setup.mat', 'nx' , 'ny', 'coeffs', 'rhs', 'bvals', 'Nc_x', 'Nc_y', 'extension', 'nloc', 'gamma', 'p_err_finer_fine', 'p_err_finer_multi', 'p_err_fine_multi');

