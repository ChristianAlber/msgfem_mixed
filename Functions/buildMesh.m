function msh = buildMesh(hx, hy, nx, ny)
    %BUILDMESH: This functions builds the fine mesh
    %  Input:   nx, ny: Number of elements in x and y direction
    %           hx, hy: Meshsize in x and y direction
    %  Output:  msh: the fine mesh
    msh.nx = nx;
    msh.ny = ny;
    msh.Nx = nx+1;
    msh.Ny = nx+1;
    
    msh.hx = hx;
    msh.hy = hy;

    % Compute sidelengths
    msh.Lx = msh.nx * msh.hx;
    msh.Ly = msh.ny * msh.hy;

    % We only consider two dimensional problems
    msh.dim = 2;

    % Create the mesh using the fast implementation of the package
    [msh.nelem,msh.nnode,msh.elems2nodes,msh.nodes2coord] = create_mesh_Q_2D_P1(hx,hy,nx,ny);

    % Construct arrays containing information which faces belong to elements;
    % which nodes belong to faces
    [msh.elems2faces, msh.faces2nodes] = get_edges2(msh.elems2nodes');
    msh.elems2faces = msh.elems2faces';

    %  Compute number of faces in the mesh
    msh.nfaces = max(max(msh.elems2faces));

    % Compute number of total dofs (both velocity and pressure dofs)
    msh.elengdof = msh.nfaces + msh.nelem;

    % Contains information, which dofs are related to elem
    msh.elems2faces2 = [msh.elems2faces; msh.nfaces+(1:msh.nelem)];

    % Determine boundary faces
    msh.bfaces = 1:2:(2 * msh.nx-1);                                                                    % add bottom boundary faces
    msh.bfaces = union(msh.bfaces, 2:(2 * msh.nx + 1): 2 + (2 * msh.nx + 1) * (msh.ny-1))     ;         % add left boundary faces
    msh.bfaces = union(msh.bfaces, (msh.nfaces-msh.nx+1):msh.nfaces);                                   % add upper boundary faces
    msh.bfaces = union(msh.bfaces, (2 * msh.nx + 1):(2 * msh.nx + 1):((2 * msh.nx + 1)* (msh.ny))); 	% add right boundary faces

    
    % Determine vertical and horizontal faces 
    msh.vfaces = find(abs(msh.nodes2coord(1,msh.faces2nodes(:,1))-msh.nodes2coord(1,msh.faces2nodes(:,2)))< 1e-7);
    msh.hfaces = setdiff(1:msh.nfaces, msh.vfaces);
    
    % Determine face midpoints
    msh.face_midpoints = 1/2 * msh.nodes2coord(:,msh.faces2nodes(:, 1)) + 1/2 *msh.nodes2coord(:,msh.faces2nodes(:, 2));

    
    % Determine midpoints of elements
    % old non vectorized version
%     for i = 1:msh.nelem
%         msh.midpoints(i,1) = mean(msh.nodes2coord(1, msh.elems2nodes(:,i)));
%         msh.midpoints(i,2) = mean(msh.nodes2coord(2, msh.elems2nodes(:,i)));
%     end

    % vectorized version of Chupeng
    msh.midpoints(:,1) = mean(reshape(msh.nodes2coord(1, msh.elems2nodes), 4, msh.nelem));
    msh.midpoints(:,2) = mean(reshape(msh.nodes2coord(2, msh.elems2nodes), 4, msh.nelem));
    
    
    % some renaming that the function overlapDomains works
end

