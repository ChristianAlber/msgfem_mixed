function msh = overlapDomains(msh, Nc_x, Nc_y, overlap, extension, plotPartition)
    
    if mod(msh.nx, Nc_x) ~= 0
        error('The number of fine elements %d is not divisible by the number of coarse elements %d.', msh.nx, Nc_x);
    end

    if mod(msh.ny, Nc_y) ~= 0
        error('The number of fine elements %d is not divisible by the number of coarse elements %d.', msh.ny, Nc_y);
    end
       
    % Build Partitioned Mesh

    msh.numDomain = Nc_x * Nc_y;
    msh.Nc_x = Nc_x;
    msh.Nc_y = Nc_y;
    
    msh.Par = cell(msh.numDomain,1);                % coarse Grid partition
    msh.Omg = cell(msh.numDomain,1);                % overlapping subdomains
    msh.Emg = cell(msh.numDomain,1);                % oversampling subdomains

    part = zeros(msh.nelem,1);

    for ie = 1 : msh.nelem

        ind_x = floor(Nc_x* msh.midpoints(ie,1)/msh.Lx)  +  1;
        ind_y = floor(Nc_y* msh.midpoints(ie,2)/msh.Ly)  +  1;

        part(ie) = Nc_x * (ind_y-1) + ind_x;

    end
    
    

    for ip = 1 : msh.numDomain
        msh.Omg{ip}.elements = find(part == ip);
        msh.Par{ip}.elements = msh.Omg{ip}.elements;
        msh.Omg{ip}.overlap = [];
    end
    
    % Find all boundary faces of the partition domains
    
     for ip = 1 : msh.numDomain
        % First, we find the maximum and minimum of coordinates of
        % midpoints of elements in the partition
        min_midpoint_el = min(msh.midpoints(msh.Par{ip}.elements,:)); 
        max_midpoint_el = max(msh.midpoints(msh.Par{ip}.elements,:));
        
        % Next, find all faces in the partition
        tmp_faces = unique(msh.elems2faces(:,msh.Par{ip}.elements));
        msh.Par{ip}.gfaces = tmp_faces;                                 % Contains all faces, indexed globally
        
        % Compute all midpoints of faces in the partition
        tmp_faces_midpoints = 1/2 * msh.nodes2coord(:,msh.faces2nodes(tmp_faces, 1)) + 1/2 *msh.nodes2coord(:,msh.faces2nodes(tmp_faces, 2));
        
        % Compute maximum and minimum coordinates of edges
        tmp_max = max(tmp_faces_midpoints');
        msh.Par{ip}.xmax = tmp_max(1);
        msh.Par{ip}.ymax = tmp_max(2);
        tmp_min= min(tmp_faces_midpoints');
        msh.Par{ip}.xmin = tmp_min(1);
        msh.Par{ip}.ymin = tmp_min(2);
        
        
        % Boundary faces are faces whose midpoint is outside of the region
        % determined by the maximum midpoints of elements
        tmp_bd_left = find(tmp_faces_midpoints(1,:) < (min_midpoint_el(1) -1e-5));
        tmp_bd_right = find(tmp_faces_midpoints(1,:) > (max_midpoint_el(1) + 1e-5));
        tmp_bd_bottom = find(tmp_faces_midpoints(2,:) < (min_midpoint_el(2) - 1e-5));
        tmp_bd_top = find(tmp_faces_midpoints(2,:) > (max_midpoint_el(2) + 1e-5));
        
        % save the obtained faces of the partition, they are global.
        msh.Par{ip}.Bgfaces = tmp_faces(union(union(union(tmp_bd_left, tmp_bd_right), tmp_bd_top), tmp_bd_bottom));
        %msh.Par{ip}.face_midpoints = tmp_faces_midpoints(:,union(union(union(tmp_bd_left, tmp_bd_right), tmp_bd_top), tmp_bd_bottom));
     end
     
     
    
    if(plotPartition == true)
    
        partition.name = 'Partition_Data';

        partition.data = zeros(msh.nelem,1);

        for i = 1:msh.numDomain
            ie = msh.Omg{i}.elements;
            partition.data(ie) = i*ones(length(ie),1);
        end

        matlab2vtk('partition_test2D_initalpartition.vtk','Partion_Test2D Initial Partition', msh,'quad',[],[],partition);
        
    end
    
    

% Construct Overlapping Domain

msh.sizeOverlap = overlap;
msh.sizeExtension = extension;

Ke = (1/4)*ones(4,1); Ke_all = repmat(Ke,1,msh.nelem);

indx_i = msh.elems2nodes;
indx_j = repmat(1:msh.nelem,4,1);

M  = sparse(indx_i,indx_j,Ke_all);
% disp(M)
msh.nodeCnt = zeros(msh.nnode,1);

% constructing the overlapping domains from the non-overlapping domains

    for ip = 1 : msh.numDomain

        for i = 1 : msh.sizeOverlap

            e = zeros(msh.nelem,1); % Initial element marker to zero

            e(msh.Omg{ip}.elements) = ones(length(msh.Omg{ip}.elements),1); % Mark Element in current partition with 1

            tmp_elements = msh.Omg{ip}.elements;

            n = M * e;
            
            if(i == msh.sizeOverlap)
                tmpNodes = msh.Omg{ip}.nodes;
            end

            msh.Omg{ip}.nodes = find(n > 0);
            

            if (i == msh.sizeOverlap)
                msh.Omg{ip}.dof = msh.Omg{ip}.nodes;

                %msh.Omg{ip}.R = sparse(1:length(msh.Omg{ip}.dof),msh.Omg{ip}.dof,ones(1,length(msh.Omg{ip}.dof)),length(msh.Omg{ip}.dof),msh.nnode);
            end

            n(msh.Omg{ip}.nodes) = ones(length(msh.Omg{ip}.nodes),1);

            e = M'*n;

            msh.Omg{ip}.elements = find(e > 0);


        end

        e = zeros(msh.nelem,1);

        e(msh.Omg{ip}.elements) = ones(length(msh.Omg{ip}.elements),1);

        n = M * e;

        msh.Omg{ip}.nodes = find(n > 0);
        
        
        [msh.Omg{ip}.iBglobal,msh.Omg{ip}.iBlocal] = setdiff(msh.Omg{ip}.nodes,tmpNodes);
        
        msh.Omg{ip}.dofb = msh.Omg{ip}.nodes;

        msh.Omg{ip}.R = sparse(1:length(msh.Omg{ip}.dofb),msh.Omg{ip}.dofb,ones(1,length(msh.Omg{ip}.dofb)),length(msh.Omg{ip}.dofb),msh.nnode);

  %     msh.nodeCnt(msh.Omg{ip}.dof) =  msh.nodeCnt(msh.Omg{ip}.dof) +  ones(length(msh.nodeCnt(msh.Omg{ip}.dof)),1);

        msh.Omg{ip}.nodes2coord = msh.nodes2coord(:,msh.Omg{ip}.nodes);

    end
    
    
    
    % constructing the extension domains from the overlapping domains
    
     for ip = 1 : msh.numDomain
         
         msh.Emg{ip}.elements = msh.Omg{ip}.elements;      
         
       if(msh.sizeExtension > 0.5)

        for i = 1 : msh.sizeExtension

            e = zeros(msh.nelem,1); % Initial element marker to zero

            e(msh.Emg{ip}.elements) = ones(length(msh.Emg{ip}.elements),1); % Mark Element in current partition with 1

            tmp_elements = msh.Emg{ip}.elements;

            n = M * e;
            
           

            msh.Emg{ip}.nodes = find(n > 0);
            
            if(i == msh.sizeExtension)
                tmpNodes = msh.Emg{ip}.nodes;
            end
            

            if (i == msh.sizeExtension)
                msh.Emg{ip}.dof = msh.Emg{ip}.nodes;

                %msh.Omg{ip}.R = sparse(1:length(msh.Omg{ip}.dof),msh.Omg{ip}.dof,ones(1,length(msh.Omg{ip}.dof)),length(msh.Omg{ip}.dof),msh.nnode);
            end

            n(msh.Emg{ip}.nodes) = ones(length(msh.Emg{ip}.nodes),1);

            e = M'*n;

            msh.Emg{ip}.elements = find(e > 0);

        end     
        
 
        

        e = zeros(msh.nelem,1);

        e(msh.Emg{ip}.elements) = ones(length(msh.Emg{ip}.elements),1);

        n = M * e;

        msh.Emg{ip}.nodes = find(n > 0);
        
        
        [msh.Emg{ip}.iBglobal,msh.Emg{ip}.iBlocal] = setdiff(msh.Emg{ip}.nodes,tmpNodes);
        

        msh.Emg{ip}.dofb = msh.Emg{ip}.nodes;

        % Build Restriction operators
        msh.Emg{ip}.R = sparse(1:length(msh.Emg{ip}.dofb),msh.Emg{ip}.dofb,ones(1,length(msh.Emg{ip}.dofb)),length(msh.Emg{ip}.dofb),msh.nnode);
        

        msh.nodeCnt(msh.Emg{ip}.dof) =  msh.nodeCnt(msh.Emg{ip}.dof) + ones(length(msh.nodeCnt(msh.Emg{ip}.dof)),1);

        msh.Emg{ip}.nodes2coord = msh.nodes2coord(:,msh.Emg{ip}.nodes);
        
        % Find out which are the interior boundary elements
        msh.Emg{ip}.iBgelements = setdiff(msh.Emg{ip}.elements, tmp_elements);
        % Find out which are the interior boundary faces
        % Obtain all faces corresponding to interior boundary elements
        unique_faces = unique(msh.elems2faces(:,msh.Emg{ip}.iBgelements ));
        % Check for every such face, whether two corresponding nodes are
        % boundary nodes
        valid_nodes = ismember(msh.faces2nodes(unique_faces, :), msh.Emg{ip}.iBglobal);
        valid_faces = and(valid_nodes(:,1), valid_nodes(:,2));
        msh.Emg{ip}.iBgfaces = unique_faces(valid_faces);
        
       else
           
        msh.Emg{ip}.nodes = msh.Omg{ip}.nodes;
               
        msh.Emg{ip}.iBglobal = msh.Omg{ip}.iBglobal;
            
        msh.Emg{ip}.iBlocal = msh.Omg{ip}.iBlocal;
        
        msh.Emg{ip}.dofb = msh.Omg{ip}.dofb;
        
        msh.Emg{ip}.dof = msh.Omg{ip}.dof;

        msh.Emg{ip}.R = msh.Omg{ip}.R;

        msh.nodeCnt(msh.Emg{ip}.dof) =  msh.nodeCnt(msh.Emg{ip}.dof) +  ones(length(msh.nodeCnt(msh.Emg{ip}.dof)),1);

        msh.Emg{ip}.nodes2coord = msh.nodes2coord(:,msh.Emg{ip}.nodes);
       end

     end
     

    


     
  
     for i = 1 :msh.numDomain
    
        msh.Emg{i}.elems2nodes = msh.elems2nodes(:,msh.Emg{i}.elements);
        msh.Emg{i}.nel = length(msh.Emg{i}.elements);
        msh.Emg{i}.nnode = length(msh.Emg{i}.nodes);
        
        msh.Emg{i}.lelem = msh.Emg{i}.elems2nodes;

% Old version        
%         for j  = 1 : length(msh.Emg{i}.nodes)
% 
% %             counter =  counter  + 1;
% % 
% %             id =  find(msh.Emg{i}.elem == msh.Emg{i}.nodes(j));
% % 
% %             msh.Emg{i}.lelem(id) =  counter*ones(length(id),1);
% 
%          msh.Emg{i}.lelem(msh.Emg{i}.lelem == msh.Emg{i}.nodes(j)) = j;
%     
%         end 
       

% Updated by Chupeng recently

        lIndex = 1 : length(msh.Emg{i}.nodes);
        g2l_map = [msh.Emg{i}.nodes lIndex'];
        [~, loc] = ismember(msh.Emg{i}.lelem, g2l_map(:,1));
        msh.Emg{i}.lelem(:) = g2l_map(loc,2);
   
     end

     
     
%      for ip = 1 : msh.numDomain
%          msh.Omg{ip}.dofv = int16.empty;
%          msh.Omg{ip}.dofp =  msh.Omg{ip}.elements + msh.nfaces;
%          for j = 1: length(msh.Omg{ip}.elements)
%             msh.Omg{ip}.dofv = union(msh.Omg{ip}.dofv, msh.elems2faces(:,msh.Omg{ip}.elements(j)));
%          end
%          
%          msh.Emg{ip}.dofv = int16.empty;
%          msh.Emg{ip}.dofp = msh.Emg{ip}.elements + msh.nfaces;
%          for j = 1:msh.Emg{ip}.nel
%             msh.Emg{ip}.dofv = union(msh.Emg{ip}.dofv, msh.elems2faces(:,msh.Emg{ip}.elements(j)));
%          end
%      end
     
     % Modification of Christian for mixed problem. Here we determine the
     % velocity dofs (dofv) and the pressure dofs (dofp) of the subdomains
     
     for ip = 1 : msh.numDomain
         msh.Omg{ip}.dofp = msh.Omg{ip}.elements + msh.nfaces;
         msh.Omg{ip}.dofv = sort(unique(msh.elems2faces(:,msh.Omg{ip}.elements)));
         msh.Emg{ip}.dofp = msh.Emg{ip}.elements + msh.nfaces;
         msh.Emg{ip}.dofv = sort(unique(msh.elems2faces(:,msh.Emg{ip}.elements)));
     end
     
     %create key from global nodes to local nodes
     % the following sections are slow and need to be vectorized.
     
     for ip = 1 : msh.numDomain
         msh.Emg{ip}.nodeKey = containers.Map(msh.Emg{ip}.nodes,1:length(msh.Emg{ip}.nodes));
     end
     
     % create elems2localnodes 
     for ip = 1 :  msh.numDomain
        msh.Emg{ip}.elems2localnodes = zeros('like',msh.Emg{ip}.elems2nodes);
        for j = 1 : msh.Emg{ip}.nel
            for k = 1 : 4
            	msh.Emg{ip}.elems2localnodes(k,j) = msh.Emg{ip}.nodeKey(msh.Emg{ip}.elems2nodes(k,j));
            end
        end
     end
     
    
     % create faces2localnodes
     for ip = 1 :  msh.numDomain
        msh.Emg{ip}.faces2localnodes = zeros(2,length(msh.Emg{ip}.dofv));
        for j = 1 : length(msh.Emg{ip}.dofv)
            for k = 1 : 2
                msh.Emg{ip}.faces2localnodes(k,j) = msh.Emg{ip}.nodeKey(msh.faces2nodes(msh.Emg{ip}.dofv(j),k));
            end
        end
     end
     
     % create key from global faces to local faces
     % the following sections are slow and need to be vectorized.
     
     for ip = 1 : msh.numDomain
         msh.Emg{ip}.faceKey = containers.Map(msh.Emg{ip}.dofv,1:length(msh.Emg{ip}.dofv));
     end
     
     % create elems2localfaces

     for ip = 1 :  msh.numDomain
        msh.Emg{ip}.elems2faces = msh.elems2faces(:,msh.Emg{ip}.elements);
        msh.Emg{ip}.elems2localfaces = zeros('like',msh.Emg{ip}.elems2faces);
        for j = 1 : msh.Emg{ip}.nel
            for k = 1 : 4
            	msh.Emg{ip}.elems2localfaces(k,j) = msh.Emg{ip}.faceKey(msh.Emg{ip}.elems2faces(k,j));
            end
        end
        % create elems2localfaces2
        msh.Emg{ip}.elems2localfaces2 = [msh.Emg{ip}.elems2localfaces ;length(msh.Emg{ip}.dofv) + (1:msh.Emg{ip}.nel)];
        
        % total dof number
        msh.Emg{ip}.totdof = length(msh.Emg{ip}.dofp)+length(msh.Emg{ip}.dofv);
     
     end
     
     % create interior boundary faces local
     for ip = 1 : msh.numDomain
         msh.Emg{ip}.iBlfaces = arrayfun(@(x) msh.Emg{ip}.faceKey(x),msh.Emg{ip}.iBgfaces);
         % Also create locally indexed boundary faces for the partition
         % domains. Locally here means with respect to the indexing in the
         % oversampling domain Emg. 
         msh.Par{ip}.Blfaces = arrayfun(@(x) msh.Emg{ip}.faceKey(x),msh.Par{ip}.Bgfaces);
         
         
         % Also create locally indexed faces for the partition
         % domains. Locally here means with respect to the indexing in the
         % oversampling domain Emg. 
         msh.Par{ip}.lfaces = arrayfun(@(x) msh.Emg{ip}.faceKey(x),msh.Par{ip}.gfaces);
         msh.Emg{ip}.lfaces = arrayfun(@(x) msh.Emg{ip}.faceKey(x),msh.Emg{ip}.dofv);  
     end
     
     % create locally (wrt oversampling domain) index elements of partition
     % domain
     for ip = 1 : msh.numDomain
         msh.Par{ip}.lelements =  find(ismember(msh.Emg{ip}.elements, msh.Par{ip}.elements));  
     end
     
     % create boundary faces of the global problem
     for ip = 1 : msh.numDomain
         tmp_bool = ismember(msh.bfaces, msh.Emg{ip}.elems2faces);
         msh.Emg{ip}.gBgfaces = msh.bfaces(tmp_bool);
         msh.Emg{ip}.gBvals = msh.bvals(tmp_bool);
         msh.Emg{ip}.gBlfaces = arrayfun(@(x) msh.Emg{ip}.faceKey(x),msh.Emg{ip}.gBgfaces);
         msh.Emg{ip}.Blfaces = union(msh.Emg{ip}.gBlfaces,msh.Emg{ip}.iBlfaces);
     end
     
%      for ip = 1 : msh.numDomain
%          tmp_bool = ismember(msh.bfaces, msh.Emg{ip}.elems2faces(:,ismember(msh.Emg{ip}.elements, msh.Par{ip}.elements)));
%          msh.Par{ip}.gBgfaces = msh.bfaces(tmp_bool);
%          msh.Par{ip}.gBvals = msh.bvals(tmp_bool);
%          msh.Par{ip}.gBlfaces = arrayfun(@(x) msh.Emg{ip}.faceKey(x),msh.Par{ip}.gBgfaces);         % Local in the sense of Emg
%      end
     
     % create restriction matrix for pressure
     for ip = 1 : msh.numDomain
        msh.Emg{ip}.Res_p = sparse(1:length(msh.Emg{ip}.elements),1:length(msh.Emg{ip}.elements), ismember(msh.Emg{ip}.elements, msh.Par{ip}.elements));
         
        msh.Emg{ip}.Rv = sparse(1:length(msh.Emg{ip}.dofv),msh.Emg{ip}.dofv,ones(1,length(msh.Emg{ip}.dofv)),length(msh.Emg{ip}.dofv),msh.nfaces);
        msh.Emg{ip}.Rp = sparse(1:length(msh.Emg{ip}.elements),msh.Emg{ip}.elements,ones(1,length(msh.Emg{ip}.elements)),length(msh.Emg{ip}.elements),msh.nelem);
     end
     
     % Compute maximum and minimum coordinates of oversampling domains
      % Find all boundary faces of the partition domains
    
     for ip = 1 : msh.numDomain
        % First, we find the maximum and minimum of coordinates of
        % midpoints of elements in the partition
        min_midpoint_el = min(msh.midpoints(msh.Par{ip}.elements,:)); 
        max_midpoint_el = max(msh.midpoints(msh.Par{ip}.elements,:));
   
        
        % Compute maximum and minimum coordinates of edges
        tmp_max = max(msh.nodes2coord(:,msh.Emg{ip}.elems2nodes),[],2);
        msh.Emg{ip}.xmax = tmp_max(1);
        msh.Emg{ip}.ymax = tmp_max(2);
        tmp_min= min(msh.nodes2coord(:,msh.Emg{ip}.elems2nodes),[],2);
        msh.Emg{ip}.xmin = tmp_min(1);
        msh.Emg{ip}.ymin = tmp_min(2);
     end
     
end