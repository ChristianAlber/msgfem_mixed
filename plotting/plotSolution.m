function plotSolution(val,msh, name,fileName)

    field.name = name;
    field.data = val;
    matlab2vtk (fileName,name, msh,'quad',field,[], []);

end