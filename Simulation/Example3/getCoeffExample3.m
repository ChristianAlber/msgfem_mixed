function coeffs = getCoeffExample4()
% Returns horizontal layer j of the SPE10 benchmark. In y-direction, we
% truncate from 220 elements to 180.

j = 1;
coeffs = zeros(2, 60 * 60);

% Load the data
SPE10_perm = load('SPE10_perm.mat', 'SPE10_perm').SPE10_perm;

%Extract layer j
coeff_layer = SPE10_perm(1:60, 1:60, j);
coeffs(1,:) = reshape(coeff_layer,[1,60 * 60]);
coeffs(2,:) = coeffs(1,:);
end

