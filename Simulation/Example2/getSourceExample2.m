function rhs = getSourceExample2(msh)
% rhs = zeros(msh.nelem,1);
% for i = 1:msh.nelem
%     x = msh.midpoints(i,1);
%     y = msh.midpoints(i,2);
%     if x < 0.1 && y < 0.1
%         rhs(i) = 100;
%     end
%     
%     if x > 0.9 && y > 0.9
%         rhs(i) = 100;
%     end
%     
%     if x < 0.1 && y > 0.9
%         rhs(i) = -100;
%     end
%     
%     if x > 0.9 && y < 0.1
%         rhs(i) = -100;
%     end
rhs = getSourceNeu(msh, 1);

end

