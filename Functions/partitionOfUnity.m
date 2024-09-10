function msh = partitionOfUnity(msh)

      xi_all = zeros(msh.nnode,1);

        for ip = 1 :  msh.numDomain
            xi_all(msh.Omg{ip}.dof) = xi_all(msh.Omg{ip}.dof) + 1;
        end

        xi_all = 1 ./ xi_all;

        for ip = 1 :  msh.numDomain
            msh.Emg{ip}.xi = xi_all(msh.Emg{ip}.dofb);
            [msh.Emg{ip}.procbnd,IA] = setdiff(msh.Emg{ip}.dofb,msh.Emg{ip}.dof);
            msh.Emg{ip}.procbnd_local = IA;
            
            [IC,IB] = setdiff(msh.Emg{ip}.dofb,msh.Omg{ip}.dof);            
            msh.Emg{ip}.xi(IB) = zeros(length(IC),1);
            msh.Emg{ip}.X = sparse(1:length(msh.Emg{ip}.xi),1:length(msh.Emg{ip}.xi), msh.Emg{ip}.xi);
            
            % build the partition of unity matrices for velocity and
            % pressure separately
            msh.Emg{ip}.xip = mean(msh.Emg{ip}.xi(msh.Emg{ip}.elems2localnodes));
            msh.Emg{ip}.xiv = mean(msh.Emg{ip}.xi(msh.Emg{ip}.faces2localnodes));
            msh.Emg{ip}.Xp = sparse(1:length(msh.Emg{ip}.xip),1:length(msh.Emg{ip}.xip), msh.Emg{ip}.xip);
            msh.Emg{ip}.Xv = sparse(1:length(msh.Emg{ip}.xiv),1:length(msh.Emg{ip}.xiv), msh.Emg{ip}.xiv);
%             msh.Emg{ip}.Xp = sparse
        end
        
%         for ip = 1 :  msh.numDomain
%             msh.Omg{ip}.xi = xi_all(msh.Omg{ip}.dofb);
%             [msh.Omg{ip}.procbnd,ia] = setdiff(msh.Omg{ip}.dofb,msh.Omg{ip}.dof);
%             msh.Omg{ip}.procbnd_local = ia;
%             msh.Omg{ip}.xi(ia) = zeros(length(msh.Omg{ip}.procbnd),1);
%             msh.Omg{ip}.X = sparse(diag(msh.Omg{ip}.xi));
%         end
    
    
end
