% subdom = 21;
% tmp = zeros(msh.nelem,1);
% diff = abs(ph-PH);
% tmp(msh.Par{subdom}.elements) = diff(msh.Par{subdom}.elements);
% plotField(tmp, msh, 'test', 'abs(ph-PH).vtk');


% subdom = 1;
% tmp = zeros(msh.nelem,1);
% tmp(msh.Emg{subdom}.elements) = basis{subdom}.Ptilde(:,1);
% plotField(tmp, msh, 'test', 'ptilde.vtk');

cnt = 0;
for i = 1:msh.numDomain
    cnt = cnt + size(basis{i}.PMS,2); 
end