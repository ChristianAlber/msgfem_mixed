function [Aug] = augmentedFormWeighted(elems, signs, h, coeffs)
% This function builds the matrix corresponding to the augmented form, ie.e
% gamma (divu,divv)
%
% IN:   elems    elements by: edges in 2D / faces in 3D
%       B_K_det  affine map determinants
%       signs    RT basis function signs per element
%       coeffs   coefficient A  
%
% OUT:  STIFF    the stiffness matrix
%
% Modified for Quadrilaterals December 2016


 nelems   = size(elems,1);        % number of elements
 [ip,w,nip] = intquadrilateral(2);
[~,~,nbasis] = basis_RT0(ip);

Aug = zeros(nbasis,nbasis,nelems);
    for m=1:nbasis
        for k=m:nbasis
            Aug(m,k,:) = 1/(h * h) * 4 * signs(:,m) .* signs(:,k) ./ coeffs'; 
        end
    end
% copy symmetric entries of the local matrices
Aug = copy_triu(Aug);
%system('free');
%whos
end