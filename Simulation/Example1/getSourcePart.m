function rhs = getSourcePart(msh, alpha)
%Computes the source term
rhs = zeros(msh.nelem,1);
for i = 1:msh.nelem
    x = msh.midpoints(i,1);
    y = msh.midpoints(i,2);
%     rhs(i) = 1 - 2 *x;
    rhs(i) = exp(- 10 * ((x-0.25)^2 +(y-0.5)^2));
% if x < 0.2 && y<0.2
%     rhs(i) =1;
% end
% if x > 0.8 && y>0.8
%     rhs(i) =-1;
% end
%     rhs(i) = 2* (1-x-y);
end

end


