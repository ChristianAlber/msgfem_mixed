function rhs = getSourceNeu(msh, alpha)
%Computes the source term
rhs = zeros(msh.nelem,1);
for i = 1:msh.nelem
    x = msh.midpoints(i,1);
    y = msh.midpoints(i,2);
%     rhs(i) = 1 - 2 *x;
    rhs(i) = 2 * alpha * alpha *pi * pi * cos(alpha * pi *x) * cos(alpha * pi *y);
% if x < 0.2 && y<0.2
%     rhs(i) =1;
% end
% if x > 0.8 && y>0.8
%     rhs(i) =-1;
% end
%     rhs(i) = 2* (1-x-y);
end

end


