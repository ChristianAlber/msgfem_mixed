function plotVectorField(vec, msh, filename)
    % Creates the vtk file for the given vector field vec
    
    % x components
    u = - 1/2 * vec(msh.elems2faces(1,:)) - 1/2 * vec(msh.elems2faces(3,:));
    
    % y components
    v = 1/2 * vec(msh.elems2faces(2,:)) + 1/2 * vec(msh.elems2faces(4,:));
    
    
    vtkwrite(filename,'unstructured_grid', msh.midpoints(:,1)',msh.midpoints(:,2)',zeros(1,msh.nelem),'vectors','v',u,v,zeros(1,msh.nelem),'scalars', 'p', zeros(1,msh.nelem));

end

