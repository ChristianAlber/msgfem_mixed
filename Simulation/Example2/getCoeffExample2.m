function coeffs = getCoeffExample2(msh)
% This function returns a high contrast coefficient
    alpha = 0.02; 
    contrast = 1000;
    coeffs = ones(msh.dim,msh.nelem);
    for i = 1:msh.nelem
        x = msh.midpoints(i,1);
        y = msh.midpoints(i,2);
        
        % Two big quadrilaterals
        if x < 0.7 && x > 0.3
            if y < 0.75 && y > 0.65
                coeffs(:,i) = contrast;
            end
            if y < 0.45 && y > 0.35
                coeffs(:,i) = contrast;
            end
        end  
        
        % Small rectangles 
        rect_midpoints = union([ 0.15:0.1:0.85 ; 0.15 *ones(1,8)]', [ 0.15:0.1:0.85 ; 0.85 *ones(1,8)]', 'rows'); 
        rect_midpoints = union(rect_midpoints, [ 0.15 *ones(1,8); 0.15:0.1:0.85 ]', 'rows');
        rect_midpoints = union(rect_midpoints, [ 0.85 *ones(1,8); 0.15:0.1:0.85 ]', 'rows');
        
        for j = 1:28
            if abs(x - rect_midpoints(j,1)) < alpha && abs(y - rect_midpoints(j,2)) < alpha
                coeffs(:,i) = contrast;
            end
        end
    end
end

