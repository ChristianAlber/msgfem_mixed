
% nx=100;
% ny=100;
% dim = 2;           % physical space dimension (only ndim = 2 supported)
% 
% % (mean) value of the coefficient
% 
% %Create mesh and show our figure, setup input data, currently random coeffs
% coeffs_mean = 1;
% [nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_test_Q_P1(nx,ny);
% [elems2faces,faces2nodes] = get_edges2(elems2nodes');
% elems2faces = elems2faces';
% nfaces = max(max(elems2faces));
% elengdof = nfaces + nelem;
% elems2faces2 = [elems2faces; nfaces+(1:nelem)];
% figure(1); show_mesh2(elems2nodes',nodes2coord'); title('mesh');%...
% coeffs = ones(dim,nelem);
% 
% %Vectorized
% tic         
% [jac,detj] = getDeterminants(elems2nodes,nodes2coord,dim,nelem);
% %TED NOTE:  Modify signs by standard rules
% %signs = ones(nelem,4);
% signs = signs_edges_Q(elems2nodes');
% %signs = [1 1 -1 -1; 1 1 -1 -1; 1 1 -1 -1; 1 1 -1 -1];
% %test = mass_matrix_RT0(elems2faces,jac,detj,signs,coeffs)
% [MPV] = mpv_matrix_RT0(elems2faces',jac,detj,signs,coeffs);
% time(3) = toc;
% tic
% MPV0 = mat2cell(MPV,size(MPV,1),size(MPV,2),ones(size(MPV,3),1));
% MPVA = assemble(MPV0,elems2faces2,1:elengdof,1:nelem,elengdof);
% time(4) = toc;
% 
% f = zeros(1,nelem);
% f(1:2 * nx) = 1;
% f(nx/2 *ny:(nx/2+1)*ny) = -1;
% g = zeros(1,elengdof - nelem);
% rhs = [g,f];
% 
% bc_dof = union(1:nx, nx * (ny-1) +1:nx*ny);
% bc_dof = union(bc_dof,1:nx:nx*ny);
% bc_dof = union(bc_dof,nx:nx:nx*ny);
% bc_dof = bc_dof + (elengdof - nelem);
% 
% sol = zeros(elengdof,1);
% free = setdiff(1:elengdof, bc_dof);
% sol(free) = MPVA(free,free) \ rhs(free)';
% p = sol(elengdof - nelem +1:elengdof);
% p = reshape(p,[nx,ny]); 

hx = 0.01;
hy = 0.01;
nx=100;
ny =100;

[nelem,ngnodes,elems2nodes,nodes2coord] = create_mesh_Q_2D_P1(hx,hy,nx,ny);
[elems2edges, edges2nodes]=get_edges2(elems2nodes');
