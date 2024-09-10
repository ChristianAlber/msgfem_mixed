function u_finer= interpolation(u, msh, msh_finer)
% Interpolates the vector field u, which lives on the space V_h
% corresponding to msh to the finer space V_{h/2} corresponding
% to msh_finer

u_finer = zeros(msh_finer.nfaces, 1);

row_ids = [];
col_ids = [];
vals = [];
for j = 1:msh.nfaces
    % Find coordinates of coarse edge j's endpoints
    node1_id = msh.faces2nodes(j,1);
    node2_id = msh.faces2nodes(j,2);
    node1_x = msh.nodes2coord(1,node1_id);
    node1_y = msh.nodes2coord(2,node1_id);
    node2_x = msh.nodes2coord(1,node2_id);
    node2_y = msh.nodes2coord(2,node2_id);

    % Distinguish cases of vertical and horizontal coarse edge

    if abs(node1_x - node2_x) < 1e-6
        % vertical edge

        % Find all relevant fine edges that are affected by the 
        % considered coarse basis function. There are at most 6.
        rel_edge_ids = find(abs(msh_finer.face_midpoints(1,msh_finer.vfaces) - node1_x) < msh_finer.hx + 1e-6  & abs(msh_finer.face_midpoints(2,msh_finer.vfaces) - (node1_y + node2_y)/2 ) < msh_finer.hy/2 + 1e-6);
        rel_edge_ids = msh_finer.vfaces(rel_edge_ids);
        %disp(rel_edge_ids);
        for i = 1:size(rel_edge_ids,2)
            row_ids(end+1) = rel_edge_ids(i);
            col_ids(end+1) = j; 
            if abs(msh_finer.face_midpoints(1,rel_edge_ids(i)) - node1_x) < 1e-6
                % If fine edge is on coarse edge
                vals(end+1) = 0.5;
            else
                vals(end+1) = 0.25;
            end        
        end
    else
        % horizontal edge

        % Find all relevant fine edges that are affected by the 
        % considered coarse basis function. There are at most 6.
        rel_edge_ids = find(abs(msh_finer.face_midpoints(2,msh_finer.hfaces) - node1_y) < msh_finer.hy + 1e-6  & abs(msh_finer.face_midpoints(1,msh_finer.hfaces) - (node1_x + node2_x)/2 ) < msh_finer.hx/2 +1e-6);
        rel_edge_ids = msh_finer.hfaces(rel_edge_ids);
        %disp(rel_edge_ids);
        for i = 1:size(rel_edge_ids,2)
            row_ids(end+1) = rel_edge_ids(i);
            col_ids(end+1) = j; 
            if abs(msh_finer.face_midpoints(2,rel_edge_ids(i)) - node1_y) < 1e-6
                % If fine edge is on coarse edge
                vals(end+1) = 0.5;
            else
                vals(end+1) = 0.25;
            end        
        end
    end
end

P = sparse(row_ids, col_ids, vals, msh_finer.nfaces, msh.nfaces);
u_finer = P * u;

end