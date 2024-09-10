function [MPVA, MPV, MPV0] = buildFineMatrix(msh, coeffs)
%BUILDFINEMATRIX builds the fine matrix for the Darcy problem
%   Detailed explanation goes here

[jac,detj] = getDeterminants(msh.elems2nodes,msh.nodes2coord,msh.dim,msh.nelem);
signs = signs_edges_Q(msh.elems2nodes');
% MPV = mpv_matrix_RT0(msh.elems2faces',jac,detj,signs,coeffs);
MPV = mpv_matrix_RT0_Darcy(msh.elems2faces',jac,detj,signs,coeffs, msh.hx);
%test
% for i = 1:msh.nelem
%     for j = 1:4
%         MPV(5,j,i) = sign(MPV(5,j,i)) * 2/msh.hx;
%         MPV(j,5,i) = sign(MPV(j,5,i)) * 2/msh.hx;
%     end
% end
MPV0 = mat2cell(MPV,size(MPV,1),size(MPV,2),ones(size(MPV,3),1));
MPVA = assemble(MPV0,msh.elems2faces2,1:msh.elengdof,1:msh.nelem,msh.elengdof);
end

