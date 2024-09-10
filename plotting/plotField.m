
function plotField(val,msh, name, filename)

    field.name = 'Perm';
    field.data = val;
    matlab2vtk (filename, name, msh,'quad',[],[], field);

end