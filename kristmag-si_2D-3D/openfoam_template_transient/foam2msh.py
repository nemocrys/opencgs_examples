import meshio

# Use:
# foamToVTK -latestTime -legacy -ascii
input = "VTK/sim-OF_2D_SIMPLE_fixedCrysInter_294.vtk"
output = "VTK/of-result.msh"

mesh = meshio.read(input)

print("\n\nMSH file:\n")
print("\npoints = \n", mesh.points)
print("\ncells = \n", mesh.cells)
print("\npoint_data = \n", mesh.point_data)
print("\ncell_data= \n", mesh.cell_data)
print("\nfield_data = \n", mesh.field_data)
print("\ncell_sets = \n", mesh.cell_sets)
print("\ncells_dict = \n", mesh.cells_dict)


mesh.write(output, file_format="gmsh22", binary=False)
