function plotPOU(msh)

    for i = 1 : msh.numDomain

        scalar.name = 'xi';
        scalar.data = msh.Omg{i}.R' * msh.Omg{i}.xi;
        matlab2vtk(strcat('POU_',int2str(i),'.vtk'),'POU', msh,'quad',scalar,[], []);
    end

end
