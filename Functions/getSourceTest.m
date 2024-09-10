function rhs = getSourceTest(msh)
%Computes the source term
rhs = zeros(msh.nelem,1);
for i = 1:msh.nelem
    x = msh.midpoints(i,1);
    rhs(i) = pi * pi * sin(pi *x);
%     rhs(i)= x -1/2;
end

end


