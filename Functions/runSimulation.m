function [v_err,p_err] = runSimulation(msh, coeffs, rhs, bvals, Nc_xs, Nc_ys, extensions, nlocs, gammas, filename )
% This function runs test with the given parameters
    nx = msh.nx;
    
    

    % Build fine matrix for mixed formulation
    tic
    [K, ~, ~] = buildFineMatrix(msh, coeffs);
    toc
    disp(['running time (building the fine matrix): ',num2str(toc)])
    
    % Build matrix corresponding to scalar product of divergences, need for
    % error calculations
    [K_aug, ~, ~] = buildFineMatrix_aug(msh, coeffs, 1);
    AUG = K_aug-K;

    % Define the right hand side
    f = zeros(msh.elengdof,1);
    ftmp = zeros(msh.elengdof,1);

    % We scale the source function by h to obtain the values in the basis
    f((msh.nfaces+1):msh.elengdof) = - rhs * msh.hx;
    ftmp((msh.nfaces+1):msh.elengdof) = - rhs * msh.hx;

    % Define boundary conditions for normal flux
    % Need to satisfy the compatibility condition

    bc_dofs = msh.bfaces;
    msh.bvals = bvals;
    
    p = getSourceNeu(msh, 1) /(2 * pi * pi);
    p = p - mean(p);

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
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % % % % % % % % % multiscale method % % % % % % % % % 
    % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    v_err = zeros(length(Nc_xs), length(extensions), length(nlocs), length(gammas));
    div_err = zeros(length(Nc_xs), length(extensions), length(nlocs), length(gammas));
    p_err = zeros(length(Nc_xs), length(extensions), length(nlocs), length(gammas));
    
    for i_Nc_x = 1:length(Nc_xs)
        for i_extension = 1:length(extensions)
            
            Nc_x = Nc_xs(i_Nc_x);               % Number of subdomains in x direction
            Nc_y = Nc_ys(i_Nc_x);                        % Number of subdomains in y direction
            overlap = 2;                        % Number of layers used to construct overlapping domains    
            extension = extensions(i_extension); % Number of layers used to construct oversampling domains
            
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
            
%             % Build coarse space
%             tic
%                 [basis] = buildCoarseSpace(msh, coeffs, ftmp, max(nlocs));
%             toc
%             disp(['running time (Building coarse space): ',num2str(toc)])
%             
%             for i_nloc = 1:length(nlocs)
%                 
%                 nloc = nlocs(i_nloc);                   % Number of eigenmodes used in each subdomain
%                 
%                 cS = cell(msh.numDomain,1);          
%                 cntv = 1;                               % global counter for velocity
%                 cntp = 1;                               % global counter for pressure
%                 
%                     % Put multiscale space together.
%                     % First for velocity, care about indexing
%                 for j = 1:msh.numDomain
%                    
%                     cS{j}.par_vtilde = basis{j}.par_vtilde;
%                     
%                     cS{j}.VMS = orth([basis{j}.Vtilde(:,1:nloc) basis{j}.Ven(:,1:nloc)]);
% 
%                     cS{j}.numModesv = length(cS{j}.VMS(1,:));
% 
%                     idv = cntv : (cntv + cS{j}.numModesv - 1);
% 
%                     cS{j}.globalIdsv = idv;
% 
%                     cntv = idv(end) + 1;
%                     
%                     cS{j}.par_ptilde = basis{j}.par_ptilde;
% 
%                     cS{j}.PMS = orth([basis{j}.Ptilde(:, 1:nloc) basis{j}.PRT]);  
% 
%                     cS{j}.numModesp = length(cS{j}.PMS(1,:));
% 
%                     idp = cntp : (cntp + cS{j}.numModesp - 1);
% 
%                     cS{j}.globalIdsp = idp;
% 
%                     cntp = idp(end) + 1;
%                 end
%                 
%                 for j = 1:msh.numDomain   
%                     cS{j}.totalSizev = cntv - 1;
%                     cS{j}.totalSizep = cntp - 1;
%                     cS{j}.totalSizeRTv = basis{j}.totalSizeRTv;
%                 end
%                 
%                 cS{1}.VRT_all = basis{1}.VRT_all;
                    
            for i_nloc = 1:length(nlocs)
                
                nloc = nlocs(i_nloc);                   % Number of eigenmodes used in each subdomain
                
                % Build coarse space
                tic
                    [cS] = buildCoarseSpace(msh, coeffs, ftmp, nloc);
                toc
                disp(['running time (Building coarse space): ',num2str(toc)])
                
                for i_gamma = 1:length(gammas) 

                    gamma = gammas(i_gamma);            % Parameter for augmented formulation
                    bool_RTv = 1;                       % 1 if VMS contains coarse RT0, otherwise 0

                    % Assemble coarse space
                    [K_aug, ~, ~] = buildFineMatrix_aug(msh, coeffs, gamma);

                    tic
                    [KH, fH, psol_v, psol_p, RHv, RHp] = assembleCoarseSpace(msh, cS, K_aug, ftmp, gamma, bool_RTv);
                    toc
                    disp(['running time (Assembling coarse space): ',num2str(toc)])

                    % Solve the coarse problem
                    solH = KH \ fH;

                    % Split into pressure and velocity part, add particular function and vector
                    % field to get multiscale solution.

                    if not(bool_RTv)
                        cS{1}.totalSizeRTv = 0;            % Turn of coarse RT velocity
                    end
                    sizev = size(RHv,2);
                    sizep = size(RHp,2);
                    PH = psol_p + RHp * solH((sizev+1):(sizev+sizep));
                    PH = PH - mean(PH);
                    UH = psol_v + RHv * solH(1:(sizev));

                    % scale by h to go from basis values to true values
                    PH = PH / msh.hx;
                    
                    disp(Nc_x)
                    disp(extension)
                    disp(nloc)
                    
                    % Print relative error
                    err = norm(PH-ph)/norm(ph);
                    sprintf('relative pressure error multiscale fine = %e', err)
                    p_err(i_Nc_x, i_extension, i_nloc, i_gamma) = err;

                    diff_u = uh - UH;
                    err = sqrt((diff_u' * K(1:msh.nfaces, 1:msh.nfaces) * diff_u) / (uh' * K(1:msh.nfaces, 1:msh.nfaces) * uh) );
                    sprintf('relative velocity error multiscale fine = %e', err)
                    v_err(i_Nc_x, i_extension, i_nloc, i_gamma) = err;
                    
           
                    err = sqrt((diff_u' * AUG(1:msh.nfaces, 1:msh.nfaces) * diff_u) / (uh' * AUG(1:msh.nfaces, 1:msh.nfaces) * uh) );
                    sprintf('relative divergence error multiscale fine = %e', err)
                    div_err(i_Nc_x, i_extension, i_nloc, i_gamma) = err;
                    
                    K_all = K + AUG;
                    err = sqrt((diff_u' * K_all(1:msh.nfaces, 1:msh.nfaces) * diff_u) / (uh' * K_all(1:msh.nfaces, 1:msh.nfaces) * uh) );
                    sprintf('relative full velocity error multiscale fine = %e', err)
                    %full_err(i_Nc_x, i_extension, i_nloc, i_gamma) = err;
                    
                    sprintf('true pressure error = %e', norm(PH-p)/norm(p))
                    
                    save(filename, 'nx', 'coeffs', 'rhs', 'bvals', 'Nc_xs', 'extensions', 'nlocs', 'gammas', 'v_err', 'p_err', 'div_err' );
                end
            end
        end
    end
    v_err = squeeze(v_err);
    div_err = squeeze(div_err);
    p_err = squeeze(p_err);
end

