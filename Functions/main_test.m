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


warning('on','all')


% !!!!Attention!!!!: currently only works if coarse domains all have the same size,
% would have to adjust computation of the coarse VRT spaces in
% buildCoarseSpace!

% Define parameter for rhs
alpha = 1;

% Define parameters for the mesh
nx = 128;
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
ftmp3 = zeros(msh.elengdof,1);
% f(msh.nfaces+1) = 1;
% f(msh.nfaces+2) = 1;
% f(msh.nfaces+3) = -1;
% f(msh.elengdof) = -1;

% We scale the source function by h to obtain the values in the basis
f((msh.nfaces+1):msh.elengdof) = - getSourceNeu(msh, alpha) * msh.hx;
ftmp3((msh.nfaces+1):msh.elengdof) = - getSourceNeu(msh, alpha) * msh.hx;

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
sol = Kb \ f;

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
Nc_x = 8;           % Number of subdomains in x direction
Nc_y = Nc_x;        % Number of subdomains in y direction
overlap = 2;        % Number of layers used to construct overlapping domains    
extension = 8;      % Number of layers used to construct oversampling domains
nloc = 10;          % Number of eigenmodes used in each subdomain
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
%   Computation of the Coarse Space
    %
    % IN:   msh     The mesh
    %       coeffs  coefficients of permeability field
    %       f       Right hand side of the problem (velocity+pressure part)
    %       nloc    Number of eigenmodes per subdomain
    %
    % OUT:  basis   contains particular function and multiscale spaces     

    % Build fine matrix for mixed formulation

    msh_loc = buildMesh(hx,hy,nx/Nc_x,ny/Nc_y);
    FREE = setdiff(1:msh_loc.elengdof, msh_loc.bfaces);
    [K_en, ~, ~] = buildFineMatrix(msh_loc, ones(msh_loc.dim,msh_loc.nelem));
    K_en = K_en(FREE,FREE);
    
    Klocal = cell(msh.numDomain,1);         % local matrices on oversampling domains
    Klocal_Omg = cell(msh.numDomain,1);     % local matrix on oversampling domain (only assemble elements of overlapping domain)
    Klocal_Par = cell(msh.numDomain,1);     % local matrix on oversampling domain (only assemble elements of partition domain)
    basis = cell(msh.numDomain,1);          

    Vtmp = cell(msh.numDomain,1);           % tmp storage for velocity basis
    Ptmp = cell(msh.numDomain,1);           % tmp storage for pressure basis

    cntv = 1;                               % global counter for velocity
    cntp = 1;                               % global counter for pressure
    cntRTv = 1;                             % Counter for coarse RT0 velocities
    
    % Waitbar to see how far we are
    h = waitbar(0,'Building the Coarse Space');
    f = ftmp3;
    for j = 1 : msh.numDomain

        % Assemble the matrix for the local problem 
        [jac,detj] = getDeterminants(msh.Emg{j}.elems2localnodes,msh.Emg{j}.nodes2coord,msh.dim,msh.Emg{j}.nel);
        signs = signs_edges_Q(msh.Emg{j}.elems2localnodes');
        MPV = mpv_matrix_RT0_Darcy(msh.Emg{j}.elems2localfaces',jac,detj,signs,coeffs(:,msh.Emg{j}.elements), msh.hx);
        MPV0 = mat2cell(MPV,size(MPV,1),size(MPV,2),ones(size(MPV,3),1));
        Klocal{j} = assemble(MPV0,msh.Emg{j}.elems2localfaces2,1:(length(msh.Emg{j}.dofp)+length(msh.Emg{j}.dofv)),1:msh.Emg{j}.nel,msh.Emg{j}.totdof);
        
%       % Find out which elements to assemble for right hand matrix of the
%       % Eigenvalue problem (only assemble elements of the overlapping domains)
        el2assemble = find(ismember(msh.Emg{j}.elements, msh.Omg{j}.elements) == 1); 
        Klocal_Omg{j} = assemble(MPV0,msh.Emg{j}.elems2localfaces2,1:(length(msh.Emg{j}.dofp)+length(msh.Emg{j}.dofv)),el2assemble,msh.Emg{j}.totdof);
            
        % Assemble matrix for the partition domain (only assemble elements of the partition domain)
        el2assemble = find(ismember(msh.Emg{j}.elements, msh.Par{j}.elements) == 1); 
        Klocal_Par{j} = assemble(MPV0,msh.Emg{j}.elems2localfaces2,1:(length(msh.Emg{j}.dofp)+length(msh.Emg{j}.dofv)),el2assemble,msh.Emg{j}.totdof);
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % % % Compute the local particular functions  % % % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % 

        % Allocate particular solution
        par_sol = zeros(msh.Emg{j}.totdof, 1);
        
        % Define right hand side
        flocal = f([msh.Emg{j}.dofv', msh.Emg{j}.dofp']);    
        
        %  Take into account boundary conditions: Zero Dirichlet in
        %  interior, prescribed normal flux on global boundary
        if not(isempty(msh.Emg{j}.gBlfaces))
            [K_tmp, flocal] = apply_bc(Klocal{j},flocal, [msh.Emg{j}.gBlfaces; msh.Emg{j}.gBvals]);
        else
            K_tmp = Klocal{j};
        end
        
        % Solve for local particular solution
        par_sol = K_tmp \ flocal;  

        % extract the velocity part and the pressure part.
        basis{j}.par_v = par_sol(1:length(msh.Emg{j}.dofv));
        basis{j}.par_p = par_sol(length(msh.Emg{j}.dofv) + (1:length(msh.Emg{j}.dofp)));
          
        % multiply the computed velocity by the partition of unity
        basis{j}.par_vtilde = msh.Emg{j}.Xv * basis{j}.par_v;    
        
        % Normalize particular pressure (subtract mean on partition domain)
%         basis{j}.par_ptilde = basis{j}.par_p - mean(basis{j}.par_p(ismember(msh.Emg{j}.elements, msh.Par{j}.elements)));
        basis{j}.par_ptilde = basis{j}.par_p - mean(basis{j}.par_p(msh.Par{j}.lelements));
        
        % restrict the particular pressure to the partition domain (=0 outside of partition domain)
        basis{j}.par_ptilde = msh.Emg{j}.Res_p * basis{j}.par_ptilde;
        
    
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % Solve the local eigenvalues problems for the velocities %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        % Build matrices for the Eigenvalue problem (MA * V = lambda * MB * V)       
        % get number of velocity dofs and pressure dofs in oversampling
        % domain
        ndofv = length(msh.Emg{j}.dofv);
        ndofp = length(msh.Emg{j}.dofp);
        
        % Build Matrix for the A-harmonic constraint
        % Find basis for space Z, consisting of divergence free vector fields
        % with everywhere vanishing normal flux. 
        % C_flux contains constraint for normal flux
        tmp = speye(ndofv);
        C_flux = tmp(union(msh.Emg{j}.gBlfaces,msh.Emg{j}.iBlfaces),1:ndofv);
        
        % Columns of Z are the basis for the above mentioned space
        Z = nulls([Klocal{j}(ndofv+1:ndofv+ndofp, 1:ndofv); C_flux]);
        [~, dimZ] = size(Z);
        
        
        % Constraints are enforced by matrix C
        C = Z' * Klocal{j}(1:ndofv,1:ndofv);
        
        % Build nullmatrices
        Zpp = sparse(ndofp,ndofp);
        Zzz = sparse(dimZ,dimZ);
        Zzv = sparse(dimZ,ndofv);
        Zpz = sparse(ndofp,dimZ);
        Zpv = sparse(ndofp,ndofv);
        
        % Build matrices for EVP
        MA = [Klocal{j}(1:ndofv,1:ndofv) Klocal{j}(1:ndofv,(ndofv+1):(ndofv+ndofp)) C'; Klocal{j}((ndofv+1):(ndofv+ndofp),1:ndofv) Zpp Zpz; C Zpz' Zzz];
        MB = [Klocal_Omg{j}(1:ndofv,1:ndofv) Zpv' Zzv';Zpv Zpp Zpz; Zzv Zpz' Zzz];
        
        % Have to constrain considered free dofs
        
        % A harmonic vector fields have vanishing normal trace on global
        % Neumann boundary
        FREE = setdiff(1:msh.Emg{j}.totdof, msh.Emg{j}.gBlfaces );
        
        % Test functions in the a harmonic constraint have vanishing normal
        % trace on the whole boundary of the oversampling domain
        FREE = union(FREE, msh.Emg{j}.totdof + (1:dimZ));
        
        % Allocate space for Eigenfunctions
        Vtmp{j} = zeros(msh.Emg{j}.totdof + ndofv,nloc);
        
        % Solve generalized Eigenproblem using eigs
        [Vtmp{j}(FREE,:), D] = eigs(MA(FREE, FREE), MB(FREE, FREE), nloc, 'sm', 'Tolerance', 1e-8);
        
        % Forget about the Lagrange multipliers, only care about the
        % eigenvectors
        basis{j}.V = Vtmp{j}(1:ndofv, 1:nloc);
        
        % Store the eigenvalues
        D = 1 ./ D;
        basis{j}.D = diag(D);
        
        % Multiply the velocites by partition of unity
        basis{j}.Vtilde = msh.Emg{j}.Xv * basis{j}.V;


        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % % Obtain pressure approximation space from velocities % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        
        
        % Allocate space for pressure basis
        Ptmp{j} = zeros(ndofp,nloc);
        
        % Constrained to zero
        con_v = setdiff(1:length(msh.Emg{j}.dofv), msh.Par{j}.lfaces);                              % constrained velocity dofs
        con_p = length(msh.Emg{j}.dofv) + setdiff(1:length(msh.Emg{j}.dofp), msh.Par{j}.lelements); % constrained pressure dofs
        con_z = [con_v con_p];                                                                      % all constrained dofs
        FREE = setdiff(1:msh.Emg{j}.totdof, con_z);                                                 % free dofs
        
        % For every velocity basis function, find the corresponding
        % pressure basis function
        for i = 1 : nloc

            tmp_sol = zeros(msh.Emg{j}.totdof, 1);
        
            % Define right hand side
            ftmp = zeros(msh.Emg{j}.totdof, 1); 

            % Remove dofs which are on oversampling domain, but not on
            % partition domain, by constraining dof to 0 for all these
            % dofs. For dofs on boundary of partition domain, prescribe
            % boundary condition given from the eigenvector of the above
            % eigenproblem. 
            % Solve for local particular solution
            % Second variant 
            FREEv = setdiff(msh.Par{j}.lfaces, msh.Par{j}.Blfaces);
            rhs_tmp = - Klocal_Par{j}(FREEv,msh.Par{j}.lfaces) * basis{j}.V(msh.Par{j}.lfaces,i);
            B_par = Klocal_Par{j}(FREEv,ndofv + msh.Par{j}.lelements);
            Ptmp{j}(msh.Par{j}.lelements,i) = B_par \ rhs_tmp;
            
            % normalize by subtracting mean of partition domain
            Ptmp{j}(msh.Par{j}.lelements,i) = Ptmp{j}(msh.Par{j}.lelements,i) - mean(Ptmp{j}(msh.Par{j}.lelements,i));
        end
      
        % Add obtained pressure to basis
        basis{j}.Ptilde = orth(Ptmp{j});
        
        
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % Build enrichment space from pressures % % % % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        
        
        % Allocate space for enrichment basis
        Vtmp{j} = zeros(ndofv,nloc);
        
        % Constrained to zero
        con_v = setdiff(1:length(msh.Emg{j}.dofv), msh.Par{j}.lfaces);
        con_p = length(msh.Emg{j}.dofv) + setdiff(1:length(msh.Emg{j}.dofp), msh.Par{j}.lelements);
        con_z = [con_v con_p];
        FREE = setdiff(1:msh.Emg{j}.totdof, [con_z msh.Par{j}.Blfaces']);
        
        % For every pressure basis function, find the corresponding
        % enrichment velocity
        for i = 1 : nloc
            tmp_sol = zeros(msh.Emg{j}.totdof, 1);
        
            % Define right hand side
            ftmp = zeros(msh.Emg{j}.totdof, 1); 
            ftmp((ndofv+1) : (ndofv+ndofp)) = basis{j}.Ptilde(:, i);

            % Remove dofs which are on oversampling domain, but not on
            % partition domain, by constraining dof to 0 for all these
            % dofs. For dofs on boundary of partition domain, prescribe
            % boundary condition given from the eigenvector of the above
            % eigenproblem. 
            
            %tests 
%             FREE_tmp = setdiff(FREE, msh.Par{j}.Blfaces);
%             K_tmp = Klocal_Par{j};
%             K_tmp(1:ndofv, 1:ndofv) = K_tmp(1:ndofv, 1:ndofv);
%             K_tmp = K_tmp(FREE,FREE);
            % Solve for local particular solution
            tmp_sol(FREE) = K_en \ ftmp(FREE);  
            
%             % TODO: Might be better to change to matrix without coefficients 
%             [K_tmp, ftmp] = apply_bc(Klocal_Par{j},ftmp, [msh.Par{j}.Blfaces'; zeros(1,length(msh.Par{j}.Blfaces))]);
% 
%             % Solve for local particular solution
%             tmp_sol(FREE) = K_tmp(FREE,FREE) \ ftmp(FREE);  
            
%             % Sanity check: Prints something if the computed velocity
%             % doesn't have the right divergence
%             if abs(mean(basis{j}.Ptilde(:, i) - Klocal{j}((ndofv+1):(ndofv+ndofp), 1:ndofv) * tmp_sol(1:ndofv))) > 1e-6
%                 disp('Inf-Sup')
%                 disp(abs(mean(basis{j}.Ptilde(:, i) - Klocal{j}((ndofv+1):(ndofv+ndofp), 1:ndofv) * tmp_sol(1:ndofv))))
%             end
            


            % extract the  velocity part.
            Vtmp{j}(1:ndofv,i) = tmp_sol(1:ndofv);  
        end
        
        % Add obtained velocities to basis
        basis{j}.Ven = Vtmp{j};
        
        
       
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % % Add coarse grid Raviart Thomas space of 0th order % % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        % Pressure basis function:
        % We build all coarse grid RT0 pressures 
        basis{j}.PRT = msh.Emg{j}.Res_p * ones(msh.Emg{j}.nel,1);
        
        % Velocity basis function: 
        % We build all coarse grid RT0 velocities
        basis{j}.VRT = [];
        
      % First, the velocity associated to the right boundary edge
        % If current partition domain is not on the right boundary  
        if not (mod(j, msh.Nc_x) == 0)
            cntRTv = cntRTv +1;
            tmp = zeros(msh.nfaces,1);
            % Find all vertical faces in partition domain j 
            faces1 = intersect(msh.vfaces, msh.Par{j}.gfaces);
            % Interpolate the coarse RT0 velocity at faces of partition
            % domain j
            tmp(faces1) = msh.face_midpoints(1, faces1) - msh.Par{j}.xmin;
            
            % Find all vertical faces in partition domain j+1 
            faces2 = intersect(msh.vfaces, msh.Par{j+1}.gfaces);
            % Interpolate the coarse RT0 velocity at faces of partition
            % domain j+1
            tmp(faces2) = msh.Par{j+1}.xmax - msh.face_midpoints(1, faces2);
            
            % Add obtained velocity to the struct basis
            basis{j}.VRT = [basis{j}.VRT tmp];
        end
        
        % Then, the velocity associated to the upper boundary edge
        % If current partition domain is not on the upper boundary  
        if j <= msh.Nc_y * (msh.Nc_x - 1)
            cntRTv = cntRTv +1;
            tmp = zeros(msh.nfaces,1);
            % Find all horizontal faces in partition domain j 
            faces1 = intersect(msh.hfaces, msh.Par{j}.gfaces);
            % Interpolate the coarse RT0 velocity at faces of partition
            % domain j
            tmp(faces1) = msh.face_midpoints(2, faces1) - msh.Par{j}.ymin;
            
            % Find all horizontal faces in partition domain j+Nc_x 
            faces2 = intersect(msh.hfaces, msh.Par{j+ msh.Nc_x}.gfaces);
            % Interpolate the coarse RT0 velocity at faces of partition
            % domain j+Nc_x
            tmp(faces2) = msh.Par{j+msh.Nc_x}.ymax - msh.face_midpoints(2, faces2);
            
            % Add obtained velocity to the struct basis
            basis{j}.VRT = [basis{j}.VRT tmp];
        end
        
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % % % Build the Multiscale Spaces % % % % % % % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        % Put multiscale space together.
        % First for velocity, care about indexing
        disp(j)
        basis{j}.VMS = orth([basis{j}.Vtilde basis{j}.Ven]);
        
        basis{j}.numModesv = length(basis{j}.VMS(1,:));

        idv = cntv : (cntv + basis{j}.numModesv - 1);

        basis{j}.globalIdsv = idv;

        cntv = idv(end) + 1;
        
        % For pressure
%         With coarse RT0
        basis{j}.PMS = orth([basis{j}.Ptilde basis{j}.PRT]); 

%         Without coarse RT0
%         basis{j}.PMS = orth(basis{j}.Ptilde); 

        
        basis{j}.numModesp = length(basis{j}.PMS(1,:));

        idp = cntp : (cntp + basis{j}.numModesp - 1);

        basis{j}.globalIdsp = idp;

        cntp = idp(end) + 1;
        
        % update waitbar 
        waitbar(j/msh.numDomain);
        
    end
    
        % close the waitbar;
        close(h);
    
    
    basis{1}.VRT_all = [];
    % Store size of the multiscale spaces
    for j = 1 : msh.numDomain
        basis{j}.totalSizev = cntv - 1;
        basis{j}.totalSizep = cntp - 1;
        basis{j}.totalSizeRTv = cntRTv - 1;
        basis{1}.VRT_all = [basis{1}.VRT_all basis{j}.VRT] ;
    end
%     basis{1}.VRT_all = orth(basis{1}.VRT_all);

% Assemble coarse space
[K_aug, ~, ~] = buildFineMatrix_aug(msh, coeffs, gamma);

tic
[KH, fH, psol_v, psol_p, RHv, RHp] = assembleCoarseSpace(msh, basis, K_aug, ftmp3, gamma, bool_RTv);
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

% % Plot the pressure solution and error
% plotField(PH, msh, 'test', 'PH.vtk');
% plotField(abs(PH-ph), msh, 'test', 'abs(PH-ph).vtk');
% 
% % Calculate velocity x component and error in x component
% uh_X = (- 1/2 * uh(msh.elems2faces(1,:)) - 1/2 * uh(msh.elems2faces(3,:))) /msh.hx;
% UH_X = (- 1/2 * UH(msh.elems2faces(1,:)) - 1/2 * UH(msh.elems2faces(3,:))) / msh.hx;
% 
% % Plot velocity x component and error in x component
% plotField(UH_X, msh, 'UH_X', 'UH_X.vtk');
% plotField(abs(UH_X-uh_X), msh, 'UH_X_err', 'abs(UH_X-uh_X).vtk');
% 
% % Calculate velocity y component and error in x component
% uh_Y = ( 1/2 * uh(msh.elems2faces(2,:)) + 1/2 * uh(msh.elems2faces(4,:))) /msh.hx;
% UH_Y = ( 1/2 * UH(msh.elems2faces(2,:)) + 1/2 * UH(msh.elems2faces(4,:))) / msh.hx;

% Plot velocity y component and error in x component
% plotField(UH_Y, msh, 'UH_Y', 'UH_Y.vtk');
% plotField(abs(UH_Y-uh_Y), msh, 'UH_Y_err', 'abs(UH_Y-uh_Y).vtk');

% Save setup
% save('Visualization/setup.mat', 'nx', 'coeffs', 'rhs', 'bvals', 'Nc_x', 'extension', 'nloc', 'gamma');

