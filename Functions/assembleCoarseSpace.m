function [KH, fH, psol_v, psol_p, RHv, RHp] = assembleCoarseSpace(msh, cS, K, f, gamma, bool_RTv)
    %   Assembly of the Coarse Space
    %
    % IN:   msh     The mesh
    %       cS      The coarse space
    %       K       Global (augmented!) problem matrix
    %       f       Right hand side of the problem (velocity+pressure part)
    %       gamma   Augmentation parameter
    %
    % OUT:  KH      Multiscale problem matrix
    %       fH      Multiscale right hand side
    %       psol_v  Particular velocity
    %       psol_p  Particular pressure
    %       RHv     Restriction matrix from Vh to VMS
    %       RHp     Restriction matrix from Ph to PMS

    if not(bool_RTv)
        cS{1}.totalSizeRTv =0;              % Turn off coarse RT velocities
    end


    RHv = sparse(msh.nfaces, cS{1}.totalSizev + cS{1}.totalSizeRTv);
    RHp = sparse(msh.nelem, cS{1}.totalSizep);
    psol_v = zeros(msh.nfaces, 1);
    psol_p = zeros(msh.nelem, 1);
    
    % Since we solve the augmented form, we change the right hand side
    % according to the problem formulation. 
    f(1:msh.nfaces) = f(1:msh.nfaces) + gamma * K(1:msh.nfaces,(msh.nfaces+1):msh.elengdof) * f((msh.nfaces+1):msh.elengdof);

    for j = 1 : msh.numDomain

        % Get global indices for velocity and pressure basis
        idv = cS{j}.globalIdsv;
        idp = cS{j}.globalIdsp;

        % Build restriction matrices to multiscale spaces
        RHv(:,idv) = msh.Emg{j}.Rv' * cS{j}.VMS;
        if bool_RTv
            RHv(:, (cS{1}.totalSizev+1): (cS{1}.totalSizev + cS{1}.totalSizeRTv)) = cS{1}.VRT_all;
        end
        RHp(:,idp) = msh.Emg{j}.Rp' * cS{j}.PMS;

        % Compute global particular function and vector field
        psol_v = psol_v + msh.Emg{j}.Rv' * cS{j}.par_vtilde;
        psol_p = psol_p + msh.Emg{j}.Rp' * cS{j}.par_ptilde; 
        
    end
    
    % Test whether orth basis different result
%     RHv = sporth(RHv);
%     RHp = sporth(RHp);
    
    % Get new sizes of multiscale spaces
    sizev = size(RHv,2);
    sizep = size(RHp,2);
    
    % Define coarse RHS
    fH = zeros(sizev + sizep,1);
    

    % Restrict the FEM matrices to the multiscale basis
    % AH corresponds to the form a(.,.)
    % BH corresponds to the form b(.,.)
    AH = RHv' * K(1:msh.nfaces, 1:msh.nfaces) * RHv;
    BH = RHp' * K(msh.nfaces+1:msh.elengdof, 1:msh.nfaces) * RHv;
    Z0 = sparse(sizep, sizep);
    KH = [AH BH'; BH Z0];
    
    % Compute the right hand side for the multiscale problem, i.e. restrict
    % to the multiscale basis and subtract the particular stuff. 
    fH(1:sizev) = RHv' * f(1:msh.nfaces) - RHv' * K(1:msh.nfaces,1:msh.nfaces) * psol_v - RHv' * K(1:msh.nfaces,(msh.nfaces+1):msh.elengdof) * psol_p;
    fH((sizev+1):(sizev + sizep)) = RHp' * f((msh.nfaces+1):msh.elengdof) - RHp' * K((1+msh.nfaces):msh.elengdof,1:msh.nfaces) * psol_v;


end