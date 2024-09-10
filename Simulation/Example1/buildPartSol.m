% Here we build the coarse space
function [basis] = buildCoarseSpace(msh, coeffs, f, nloc)
    %   Computation of the Coarse Space
    %
    % IN:   msh     The mesh
    %       coeffs  coefficients of permeability field
    %       f       Right hand side of the problem (velocity+pressure part)
    %       nloc    Number of eigenmodes per subdomain
    %
    % OUT:  basis   contains particular function and multiscale spaces     
    
    nlocp = nloc;                           % Number of used pressure basis functions per subdomain
    
    % Prebuild the matrix for velocity reconstruction. Only need to do this
    % once 
    msh_loc = buildMesh(msh.hx,msh.hy,msh.nx/msh.Nc_x,msh.ny/msh.Nc_y);
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
    
    time_nullspace = 0;
    time_evp = 0;
    for j = 1 : msh.numDomain
        
        % get number of velocity dofs and pressure dofs in oversampling
        % domain
        ndofv = length(msh.Emg{j}.dofv);
        ndofp = length(msh.Emg{j}.dofp);

        % Assemble the matrix for the local problem 
        [jac,detj] = getDeterminants(msh.Emg{j}.elems2localnodes,msh.Emg{j}.nodes2coord,msh.dim,msh.Emg{j}.nel);
        signs = signs_edges_Q(msh.Emg{j}.elems2localnodes');
        MPV = mpv_matrix_RT0_Darcy(msh.Emg{j}.elems2localfaces',jac,detj,signs,coeffs(:,msh.Emg{j}.elements), msh.hx);
        MPV0 = mat2cell(MPV,size(MPV,1),size(MPV,2),ones(size(MPV,3),1));
        Klocal{j} = assemble(MPV0,msh.Emg{j}.elems2localfaces2,1:(length(msh.Emg{j}.dofp)+length(msh.Emg{j}.dofv)),1:msh.Emg{j}.nel,msh.Emg{j}.totdof);
        
%       % Find out which elements to assemble for right hand matrix of the
%       % Eigenvalue problem (only assemble elements of the overlapping domains)
        el2assemble = find(ismember(msh.Emg{j}.elements, msh.Omg{j}.elements) == 1); 
        
        % !!!!!!!!!!!!!! Comment out one of the following two blocks !!!!!!!!!!!!!!
        % No partition of unity in Eigenvalue problem     
        Klocal_Omg{j} = assemble(MPV0,msh.Emg{j}.elems2localfaces2,1:(length(msh.Emg{j}.dofp)+length(msh.Emg{j}.dofv)),el2assemble,msh.Emg{j}.totdof);
        
        % Parition of unity in Eigenvalue problem
%         Aug = augmentedForm(msh.Emg{j}.elems2localfaces',signs, msh.hx, 1.0);   % Copute part involving scalar product of divergences
%         Aug = augmentedFormWeighted(msh.Emg{j}.elems2localfaces',signs, msh.hx, coeffs(1,msh.Emg{j}.elements));     % Use weighted div norm
%         MPV_aug = MPV(1:4,1:4,:) + Aug;                                         % Add the local augmented matrices to the usual local matrices
%         MPV0_aug = mat2cell(MPV_aug,size(MPV_aug,1),size(MPV_aug,2),ones(size(MPV_aug,3),1));
%         Klocal_Omg{j} = assemble(MPV0_aug,msh.Emg{j}.elems2localfaces2,1:(length(msh.Emg{j}.dofp)+length(msh.Emg{j}.dofv)),el2assemble,msh.Emg{j}.totdof);
%         Klocal_Omg{j}(1:ndofv,1:ndofv) =  msh.Emg{j}.Xv' * Klocal_Omg{j}(1:ndofv,1:ndofv) * msh.Emg{j}.Xv;       %involve partition of unity
        
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
        
        %  Need to comment out one of the following two options !!!!!!!!!!!
        %  Option 1:
        %  Take into account boundary conditions: Zero Dirichlet in
        %  interior, prescribed normal flux on global boundary
        
        %/ 
%         if not(isempty(msh.Emg{j}.gBlfaces))
%             [K_tmp, flocal] = apply_bc(Klocal{j},flocal, [msh.Emg{j}.gBlfaces; msh.Emg{j}.gBvals]);
%         else
%             K_tmp = Klocal{j};
%         end
        K_tmp = Klocal{j};
        %/
        
        %  Option 2: Working now
        %  Neumann boundary conditions everywhere. Constant Neumann in
        %  interior
        
        %/
%         % Find out right orientation for outward unit normal
%         orient_temp = ones(1,length(msh.Emg{j}.iBgfaces));
%         orient_temp(find(msh.face_midpoints(1,msh.Emg{j}.iBgfaces)> msh.Emg{j}.xmax -1e-6)) = -1;
%         orient_temp(find(msh.face_midpoints(2,msh.Emg{j}.iBgfaces)< msh.Emg{j}.ymin +1e-6)) = -1;
% 
%         % Compute constant for interior boundary conditions  
%         const_tmp =  0.5 * msh.hx / (length(msh.Emg{j}.iBgfaces)) * sum(flocal) * orient_temp;
%         
%         % Apply boundary conditions
%         [K_tmp, flocal] = apply_bc(Klocal{j},flocal, [msh.Emg{j}.gBlfaces msh.Emg{j}.iBlfaces'; msh.Emg{j}.gBvals const_tmp]);
        %/
        
        % Solve for local particular solution
        par_sol = K_tmp \ flocal;  
    
        % extract the velocity part and the pressure part.
        basis{j}.par_v = par_sol(1:length(msh.Emg{j}.dofv));
        basis{j}.par_p = par_sol(length(msh.Emg{j}.dofv) + (1:length(msh.Emg{j}.dofp)));
        
%         disp("f")
%         ftmmmp =  f([msh.Emg{j}.dofv', msh.Emg{j}.dofp']);
%         disp(sum(ftmmmp));
%         disp("normal trace")
%         disp(orient_temp * basis{j}.par_v(msh.Emg{j}.iBlfaces) )
        % multiply the computed velocity by the partition of unity
        basis{j}.par_vtilde = msh.Emg{j}.Xv * basis{j}.par_v;    
        
        % Normalize particular pressure (subtract mean on partition domain)
%         basis{j}.par_ptilde = basis{j}.par_p - mean(basis{j}.par_p(ismember(msh.Emg{j}.elements, msh.Par{j}.elements)));
        basis{j}.par_ptilde = basis{j}.par_p;
        
        % restrict the particular pressure to the partition domain (=0 outside of partition domain)
        basis{j}.par_ptilde = msh.Emg{j}.Xp * basis{j}.par_ptilde;
        
    
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % Solve the local eigenvalues problems for the velocities (old) %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        %\
        % Build matrices for the Eigenvalue problem (MA * V = lambda * MB * V)       
        
        
        % Build Matrix for the A-harmonic constraint
        % Find basis for space Z, consisting of divergence free vector fields
        % with everywhere vanishing normal flux. 
        % C_flux contains constraint for normal flux
        tmp = speye(ndofv);
        C_flux = tmp(union(msh.Emg{j}.gBlfaces,msh.Emg{j}.iBlfaces),1:ndofv);
        
        % Columns of Z are the basis for the above mentioned space
        Z = nulls([Klocal{j}(ndofv+1:ndofv+ndofp, 1:ndofv); C_flux]);
        
        % Save nullspace for tests
        basis{j}.Z=Z;
        
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
        
%         disp(size(MA))
%         disp(condest(MA))
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
        
        %\
        
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % Solve the local eigenvalues problems for the velocities (new) %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        % \
%         tmp = speye(ndofv);
%         C_in = tmp(msh.Emg{j}.iBlfaces,1:ndofv);
%         C_gl = tmp(msh.Emg{j}.gBlfaces,1:ndofv);
%         
%         tic;
%         Z = nulls([Klocal{j}(ndofv+1:ndofv+ndofp, 1:ndofv);C_gl]);      % Divergence free and zero flux on global boundary
%         Z_flux = Z * nulls(C_in * Z);                                   % Divergence free and zero flux everywhere
%         t_tmp = toc;
%         time_nullspace = time_nullspace + t_tmp;
%         
%         % Build matrices for EVP
%         A = Klocal{j}(1:ndofv,1:ndofv);
%         A_test = Z' * A * Z_flux;
%         disp(size(A_test));
%         MA = [Z' * A * Z A_test; A_test' sparse(size(Z_flux,2),size(Z_flux,2))];
% %         MA = [Z' * A * Z Z' * A * Z_flux; Z_flux' * A * Z sparse(size(Z_flux,2),size(Z_flux,2))];
%         MB = [Z' * Klocal_Omg{j}(1:ndofv,1:ndofv) * Z sparse(size(Z,2),size(Z_flux,2));sparse(size(Z_flux,2),size(Z,2))  sparse(size(Z_flux,2),size(Z_flux,2))];
%         
%         tic;
%         [Z_tmp, D] = eigs(MA, MB, nloc, 'sm', 'Tolerance', 1e-8);
%         t_tmp = toc;
%         time_evp = time_evp + t_tmp;
%         
%         Vtmp{j} = Z * Z_tmp(1:size(Z,2),:);
        
        %\
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % Solve the local eigenvalues problems for the velocities (test) %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        % \
%         tmp = speye(ndofv);
%         C_gl = tmp(msh.Emg{j}.gBlfaces,1:ndofv);
%         tic;
%         Z_gl = nulls([Klocal{j}(ndofv+1:ndofv+ndofp, 1:ndofv);C_gl]);      % Divergence free and zero flux on global boundary
%         tmp = speye(ndofv);
%         C_flux = tmp(union(msh.Emg{j}.gBlfaces,msh.Emg{j}.iBlfaces),1:ndofv);
%         Z_flux = nulls([Klocal{j}(ndofv+1:ndofv+ndofp, 1:ndofv);C_flux]);   % Divergence free and zero flux everywhere
%         t_tmp = toc;
%         time_nullspace = time_nullspace + t_tmp;
%         
%         % Build matrices for EVP
%         A = Klocal{j}(1:ndofv,1:ndofv);
%         A_test = Z_gl' * A * Z_flux;
%         MA = [Z_gl' * A * Z_gl A_test; A_test' sparse(size(Z_flux,2),size(Z_flux,2))];
%         MB = [Z_gl' * Klocal_Omg{j}(1:ndofv,1:ndofv) * Z_gl sparse(size(Z_gl,2),size(Z_flux,2));sparse(size(Z_flux,2),size(Z_gl,2))  sparse(size(Z_flux,2),size(Z_flux,2))];
%         
%         issymmetric(MA)
%         disp(size(MA));
%         disp(condest(MA));
%         
%         tic;
%         [Z_tmp, D] = eigs(MA, MB, nloc, 'sm', 'Tolerance', 1e-8);
%         t_tmp = toc;
%         time_evp = time_evp + t_tmp;
%         
%         Vtmp{j} = Z_gl * Z_tmp(1:size(Z_gl,2),:);
%         
        %\
        
        
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
        Ptmp{j} = zeros(ndofp,nlocp);
%         Ptmp2{j} = zeros(ndofp,nloc);
        % Constrained to zero
        con_v = setdiff(1:length(msh.Emg{j}.dofv), msh.Par{j}.lfaces);                              % constrained velocity dofs
        con_p = length(msh.Emg{j}.dofv) + setdiff(1:length(msh.Emg{j}.dofp), msh.Par{j}.lelements); % constrained pressure dofs
        con_z = [con_v con_p];                                                                      % all constrained dofs
        FREE = setdiff(1:msh.Emg{j}.totdof, con_z);                                                 % free dofs
        
        % For every velocity basis function, find the corresponding
        % pressure basis function
        for i = 1 : nlocp
            tmp_sol = zeros(msh.Emg{j}.totdof, 1);
        
            % Define right hand side
            ftmp = zeros(msh.Emg{j}.totdof, 1); 

            % Remove dofs which are on oversampling domain, but not on
            % partition domain, by constraining dof to 0 for all these
            % dofs. For dofs on boundary of partition domain, prescribe
            % boundary condition given from the eigenvector of the above
            % eigenproblem. 
            
            
            % Two versions: Comment out one
            
            %/
            % First variant
%             [K_tmp, ftmp] = apply_bc(Klocal_Par{j},ftmp, [msh.Par{j}.Blfaces'; basis{j}.V(msh.Par{j}.Blfaces,i)']);
% 
%             % Solve for local particular solution
%             tmp_sol(FREE) = K_tmp(FREE,FREE) \ ftmp(FREE);  
%             
%             % extract the  pressure part.
%             Ptmp{j}(1:ndofp,i) = tmp_sol((ndofv+1) : (ndofv+ndofp));
            
%             Sanity check: Prints something if the computed velocity is
%             different from the initial eigenvelocity
%             if abs(mean(tmp_sol(msh.Par{j}.lfaces) - basis{j}.V(msh.Par{j}.lfaces,i))) > 1e-6
%                 disp('Pressure from Velocity')
%                 disp(abs(mean(tmp_sol(msh.Par{j}.lfaces) - basis{j}.V(msh.Par{j}.lfaces,i))))
%             end
            %/
            
            
            %/
            % Second variant: Just solve the first equation (on partition domain).
            FREEv = setdiff(msh.Par{j}.lfaces, msh.Par{j}.Blfaces);
            rhs_tmp = - Klocal_Par{j}(FREEv,msh.Par{j}.lfaces) * basis{j}.V(msh.Par{j}.lfaces,i);
            B_par = Klocal_Par{j}(FREEv,ndofv + msh.Par{j}.lelements);
            Ptmp{j}(msh.Par{j}.lelements,i) = B_par \ rhs_tmp;
            %/
            
            %/
%             % Third variant: Just solve the first equation (on oversampling domain).
%             FREEv = setdiff(msh.Emg{j}.lfaces, msh.Emg{j}.Blfaces);
%             rhs_tmp = - Klocal{j}(FREEv,msh.Emg{j}.lfaces) * basis{j}.V(msh.Emg{j}.lfaces,i);
%             B_par = Klocal{j}(FREEv,ndofv+1:ndofv+ndofp);
%             sol_tmp = B_par \ rhs_tmp;
%             % restrict to partition domain, extend by zero
%             Ptmp{j}(msh.Par{j}.lelements,i) = sol_tmp(msh.Par{j}.lelements);     
            %/
            
            % normalize by subtracting mean of partition domain
            Ptmp{j}(msh.Par{j}.lelements,i) = Ptmp{j}(msh.Par{j}.lelements,i) - mean(Ptmp{j}(msh.Par{j}.lelements,i));
        end
        
        
        % Add obtained pressure to basis
        basis{j}.Ptilde = Ptmp{j};
        
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % Build enrichment space from pressures % % % % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        
        
        % Allocate space for enrichment basis
        Vtmp{j} = zeros(ndofv,nlocp);
        
        % Constrained to zero
        con_v = setdiff(1:length(msh.Emg{j}.dofv), msh.Par{j}.lfaces);
        con_p = length(msh.Emg{j}.dofv) + setdiff(1:length(msh.Emg{j}.dofp), msh.Par{j}.lelements);
        con_z = [con_v con_p];
        FREE = setdiff(1:msh.Emg{j}.totdof, con_z);
        
        % For every pressure basis function, find the corresponding
        % enrichment velocity
        for i = 1 : nlocp
            tmp_sol = zeros(msh.Emg{j}.totdof, 1);
        
            % Define right hand side
            ftmp = zeros(msh.Emg{j}.totdof, 1); 
            ftmp((ndofv+1) : (ndofv+ndofp)) = basis{j}.Ptilde(:, i);

            % Remove dofs which are on oversampling domain, but not on
            % partition domain, by constraining dof to 0 for all these
            % dofs. For dofs on boundary of partition domain, prescribe
            % boundary condition given from the eigenvector of the above
            % eigenproblem. 
            
            % First variant:
%             [K_tmp, ftmp] = apply_bc(Klocal_Par{j},ftmp, [msh.Par{j}.Blfaces'; zeros(1,length(msh.Par{j}.Blfaces))]);
% 
%             % Solve for local particular solution
%             tmp_sol(FREE) = K_tmp(FREE,FREE) \ ftmp(FREE);  

            % Second variant:
            FREE = setdiff(1:msh.Emg{j}.totdof, [con_z msh.Par{j}.Blfaces']);
            tmp_sol(FREE) = K_en \ ftmp(FREE);
            
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

        basis{j}.VMS = orth([basis{j}.Vtilde basis{j}.Ven]);
%         basis{j}.VMS = orth([basis{j}.Vtilde]);
         
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
    disp(['running time (Building Nullspaces): ',num2str(time_nullspace)])
    disp(['running time (solving EVP): ',num2str(time_evp)])
end

