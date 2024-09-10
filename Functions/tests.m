% L = msh.Lx;
% F = reshape(f(msh.nfaces+1:msh.elengdof),[100,100])';
% 
% N = length(F);
% dx = L / (N - 1);
% p = 0:N-1;
% q = 0:N-1;
% [p,q] = meshgrid(p,q);
% B = dct2(F);
% A = dx^2 * B ./ (2 * cos(pi * p / N) + 2 * cos(pi * q / N) - 4); A(1,1) = 0;
% U = idct2(A);
% 
% plotField(reshape(U',[100*100,1]), msh, 'test', 'U.vtk'); 

dim=2

[nelem,nnode,elems2nodes,nodes2coord] = create_mesh_Q_2D_P1(hx,hy,nx,ny);
% adjust: facewise data for RT0 element
% ToDo: elems2faces = get_edges(elems2nodes')
% etc.
%[elems2faces,edges2nodes] = get_edges3(elems2nodes');
[elems2faces,edges2nodes] = get_edges2(elems2nodes');
elems2faces = elems2faces';
nfaces = max(max(elems2faces));
elengdof = nfaces + nelem;
figure(1); show_mesh2(elems2nodes',nodes2coord'); title('mesh');%...

coeffs = rand(dim,nelem) * 10;
faces2nodes = edges2nodes;

[Mloc,signs] = setup_stiffness_matrix_Darcy_2D_Q_RT0(...
            coeffs,nelem,elems2nodes,elems2faces,faces2nodes,nodes2coord);
elems2nodes2 = [elems2faces; nfaces+(1:nelem)];
Mloc0 = mat2cell(Mloc,size(Mloc,1),size(Mloc,2),ones(size(Mloc,3),1));
M = assemble(Mloc0,elems2nodes2,1:elengdof,1:nelem,elengdof);

% Define the right hand side
f = zeros(msh.elengdof,1);
%f(msh.nfaces+1) = 1;
%f(msh.elengdof) = -1;
f((msh.nfaces+1):msh.elengdof) = getSource(msh);

% apply the boundary conditions
[M, f] = apply_bc(M,f, [bc_dofs; bc_vals]);

% solve the system
sol = M \ f;

% Split into velocity part and pressure part
uhtest = sol(1:msh.nfaces);
phtest = sol((msh.nfaces+1):msh.elengdof);

phtest = phtest - mean(phtest);
plotField(phtest, msh, 'test', 'pressure_test.vtk'); 