
for j = 1:nloc
    ven = msh.Emg{1}.Rv' * basis{1}.Ven(:,j);
    ven_X = (- 1/2 * ven(msh.elems2faces(1,:)) - 1/2 * ven(msh.elems2faces(3,:))) /msh.hx;
    plotField(ven_X, msh, 'X', "Ven"+num2str(j)+".vtk");
    
    vtilde = msh.Emg{1}.Rv' * basis{1}.Vtilde(:,j);
    vtilde_X = (- 1/2 * vtilde(msh.elems2faces(1,:)) - 1/2 * vtilde(msh.elems2faces(3,:))) /msh.hx;
    plotField(vtilde_X, msh, 'X', "Vtilde"+num2str(j)+".vtk");
    
    vdiff_X = abs(ven_X - vtilde_X);
    plotField(vdiff_X, msh, 'X', "Vdiff"+num2str(j)+".vtk");
    
end