m = 1;
dom = 1;

% x = basis{dom}.Ptilde(:,m);
% x = flip(reshape(x,[sqrt(msh.Emg{dom}.nel),sqrt(msh.Emg{dom}.nel)]));
% imagesc(x);
for m = 63:63
    x = zeros(msh.nelem,1);
    x(msh.Emg{dom}.elements) = basis{dom}.Ptilde(:,m);
    plotField(x, msh, 'test', "pressure_basis_m="+num2str(m)+".vtk");
end


% for m = 1:10
%     tmp = zeros(msh.nfaces,1);
%     tmp(msh.Emg{dom}.dofv) = basis{dom}.Vtilde(:,m);
%     x = zeros(msh.nelem,1);
%     x(msh.Emg{dom}.elements) = (- 1/2 * tmp(msh.elems2faces(1,msh.Emg{dom}.elements)) - 1/2 * tmp(msh.elems2faces(3,msh.Emg{dom}.elements))) / msh.hx;
%     plotField(x, msh, 'test', "vel_basis_x_m="+num2str(m)+".vtk");
% end
% 
% 
% for m = 1:10
%     tmp = zeros(msh.nfaces,1);
%     tmp(msh.Emg{dom}.dofv) = basis{dom}.Vtilde(:,m);
%     x = zeros(msh.nelem,1);
%     x(msh.Emg{dom}.elements) = ( 1/2 * tmp(msh.elems2faces(2,msh.Emg{dom}.elements)) + 1/2 * tmp(msh.elems2faces(4,msh.Emg{dom}.elements))) / msh.hx;
%     plotField(x, msh, 'test', "vel_basis_y_m="+num2str(m)+".vtk");
% end