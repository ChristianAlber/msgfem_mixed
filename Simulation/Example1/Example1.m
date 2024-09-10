% This code runs Example 1 from the paper 
% @article{alber2024mixed,
%   title={A Mixed Multiscale Spectral Generalized Finite Element Method},
%   author={Alber, Christian and Ma, Chupeng and Scheichl, Robert},
%   journal={arXiv preprint arXiv:2403.16714},
%   year={2024}
% }

addpath('../../2D_Q_RT0')
addpath('../../Functions')
addpath('../../plotting')
addpath('../../nulls')
addpath('../../SparseNullOrth')

warning('on','all')


% !!!!Attention!!!!: currently only works if coarse domains all have the same size,
% would have to adjust computation of the coarse VRT spaces in
% buildCoarseSpace!

% Define parameter for rhs
alpha = 1;

% Define parameters for the mesh
nx = 100;
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
f((msh.nfaces+1):msh.elengdof) = - getSourceNeu(msh, alpha) * msh.hx;
ftmp((msh.nfaces+1):msh.elengdof) = - getSourceNeu(msh, alpha) * msh.hx;

p = getSourceNeu(msh, alpha) /(2 * pi * pi * alpha * alpha);

% Define boundary conditions for normal flux
% Need to satisfy the compatibility condition

bc_dofs = msh.bfaces;
msh.bvals = zeros(1,length(msh.bfaces));

% FREE_tmp = setdiff(msh.elengdof, msh.bfaces);
% disp(1/condest(K(FREE_tmp, FREE_tmp)));

% apply the boundary conditions
[Kb, f] = apply_bc(K,f, [msh.bfaces; msh.bvals]);

% solve the system
% sol = Kb \ f;
sol = zeros(size(Kb,1), 1);
sol(1:end-1) = Kb(1:end-1,1:end-1) \ f(1:end-1);

% scale by h to go from basis values to true values
sol = sol;

% Split into velocity part and pressure part
uh = sol(1:msh.nfaces);
ph = sol((msh.nfaces+1):msh.elengdof) / msh.hx;

% normalize for neumann problems
p = p - mean(p);
ph = ph -mean(ph);

% visualize pressure
pressure = reshape(ph, [msh.nx,msh.ny])';

%ph = ph - mean(ph);
% plotField(ph, msh, 'test', 'pressure.vtk'); 
% plotField(f((msh.nfaces+1):msh.elengdof), msh, 'test', 'f.vtk'); 

% Sanity check
err = norm(p-ph)/norm(p);
sprintf('relative error true fine = %e', err)

% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Initialize everything for the multiscale method % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Initialize parameters
Nc_x = 4;           % Number of subdomains in x direction
Nc_y = Nc_x;        % Number of subdomains in y direction
overlap = 2;        % Number of layers used to construct overlapping domains    
extension = 4;      % Number of layers used to construct oversampling domains
nloc = 16;          % Number of eigenmodes used in each subdomain
gamma = 0;          % Parameter for augmented formulation
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
p_err = norm(PH-ph)/norm(ph);
sprintf('relative pressure error multiscale fine = %e', p_err)

diff_u = uh - UH;
v_err = sqrt((diff_u' * K(1:msh.nfaces, 1:msh.nfaces) * diff_u) / (uh' * K(1:msh.nfaces, 1:msh.nfaces) * uh) );
sprintf('relative velocity error multiscale fine = %e', v_err)

% % % Plot the pressure solution and error
% plotField(psol_p /msh.hx, msh, 'test', 'psol_p.vtk');
% plotField(abs(psol_p /msh.hx - ph), msh, 'test', 'abs(psol_p-ph).vtk');
% psol_v_X = (- 1/2 * psol_v(msh.elems2faces(1,:)) - 1/2 * psol_v(msh.elems2faces(3,:))) /msh.hx;
% psol_v_Y = ( 1/2 * psol_v(msh.elems2faces(2,:)) + 1/2 * psol_v(msh.elems2faces(4,:))) /msh.hx;
% plotField(psol_v_X, msh, 'X', 'psol_v_X.vtk');
% plotField(psol_v_Y, msh, 'Y', 'psol_v_Y.vtk');
% 
% % plotField(PH, msh, 'test', 'PH.vtk');
% plotField(abs(PH-ph), msh, 'test', "abs(PH-ph)_nloc="+num2str(nloc)+".vtk");
% % % 
% % % % Calculate velocity x component and error in x component
% uh_X = (- 1/2 * uh(msh.elems2faces(1,:)) - 1/2 * uh(msh.elems2faces(3,:))) /msh.hx;
% UH_X = (- 1/2 * UH(msh.elems2faces(1,:)) - 1/2 * UH(msh.elems2faces(3,:))) / msh.hx;
% % % 
% % % % Plot velocity x component and error in x component
% % plotField(UH_X, msh, 'UH_X', 'UH_X.vtk');
% plotField(abs(UH_X-uh_X), msh, 'UH_X_err', "abs(UH_X-uh_X)_nloc="+num2str(nloc)+".vtk");
% % % 
% % % % Calculate velocity y component and error in x component
% uh_Y = ( 1/2 * uh(msh.elems2faces(2,:)) + 1/2 * uh(msh.elems2faces(4,:))) /msh.hx;
% UH_Y = ( 1/2 * UH(msh.elems2faces(2,:)) + 1/2 * UH(msh.elems2faces(4,:))) / msh.hx;
% % % 
% % % % Plot velocity y component and error in x component
% % plotField(UH_Y, msh, 'UH_Y', 'UH_Y.vtk');
% plotField(abs(UH_Y-uh_Y), msh, 'UH_Y_err', "abs(UH_Y-uh_Y)_nloc="+num2str(nloc)+".vtk");
% 
% plotField(abs(psol_v_X-uh_X), msh, 'X', 'abs(psol_v_X-uh_X).vtk');
% plotField(abs(psol_v_Y-uh_Y), msh, 'Y', 'abs(psol_v_Y-uh_Y).vtk');

% Save setup
% save('Visualization/Example1/setup.mat', 'nx', 'coeffs', 'Nc_x', 'extension', 'nloc', 'gamma');

% Test for simulations

% cS = cell(msh.numDomain,1);          
% cntv = 1;                               % global counter for velocity
% cntp = 1;                               % global counter for pressure
% 
%     % Put multiscale space together.
%     % First for velocity, care about indexing
% nloc = 20;
% for j = 1:msh.numDomain
% 
%     cS{j}.par_vtilde = basis{j}.par_vtilde;
%     
%     cS{j}.Vtilde = basis{j}.Vtilde(:,1:nloc);
%     
%     cS{j}.Ven = basis{j}.Ven(:,1:nloc);
% 
%     cS{j}.VMS = orth([basis{j}.Vtilde(:,1:nloc) basis{j}.Ven(:,1:nloc)]);
% 
%     cS{j}.numModesv = length(cS{j}.VMS(1,:));
% 
%     idv = cntv : (cntv + cS{j}.numModesv - 1);
% 
%     cS{j}.globalIdsv = idv;
% 
%     cntv = idv(end) + 1;
% 
%     cS{j}.par_ptilde = basis{j}.par_ptilde;
% 
%     cS{j}.PMS = orth([basis{j}.Ptilde(:, 1:nloc) basis{j}.PRT]);  
% 
%     cS{j}.numModesp = length(cS{j}.PMS(1,:));
% 
%     idp = cntp : (cntp + cS{j}.numModesp - 1);
% 
%     cS{j}.globalIdsp = idp;
% 
%     cntp = idp(end) + 1;
% end
% 
% for j = 1:msh.numDomain   
%     cS{j}.totalSizev = cntv - 1;
%     cS{j}.totalSizep = cntp - 1;
%     cS{j}.totalSizeRTv = basis{j}.totalSizeRTv;
% end
% 
% cS{1}.VRT_all = basis{1}.VRT_all;
