import gmsh
from mpi4py       import MPI
import numpy as np

# /home/thomas/Documents/PhD/Mouse_experiment_dolfinx/

# Code made for parallel running
mesh_comm = MPI.COMM_WORLD
model_rank = 0


gmsh.initialize()


#Characteristic length
lc = 5e-6


[Length, Width, Height] = [6e-4, 2.5e-4, 4e-5]
xc,yc,zc,dx,dy,dz, r = 6e-4/2, 0, 0, 0, 0, 4e-5, 1.5e-4

print([Length, Width, Height])
gdim = 3

if mesh_comm.rank == model_rank:
	box      = gmsh.model.occ.addBox(0, 0, 0, Length, Width, Height)
	cylinder = gmsh.model.occ.addCylinder(xc,yc,zc,dx,dy,dz, r,tag=1000,angle=np.pi)
	gmsh.model.occ.synchronize()

	# Remove duplicate entities and synchronize
	gmsh.model.occ.removeAllDuplicates()
	gmsh.model.occ.synchronize()

	surfaces, volumes = [gmsh.model.getEntities(d) for d in [ gdim-1, gdim]]
	print(volumes)

	# Volumes
	gmsh.model.addPhysicalGroup(volumes[0][0], [volumes[0][1]], -1)
	gmsh.model.setPhysicalName(volumes[0][0], -1, 'Half_Cylinder')
	gmsh.model.addPhysicalGroup(volumes[1][0], [volumes[1][1]], -1)
	gmsh.model.setPhysicalName(volumes[1][0], -1, 'Box')

	# # 1 = loading, 2 = top minus loading, 3 = bottom, 4 = left, 5 = right, 6 = Front, 7 = back
	bottom_marker, front_marker, back_marker, left_marker, right_marker, top_marker, indenter_marker = 3, 6, 7, 4, 5, 2, 1
	bottom, front, back, left, right, top, indenter = [],[],[],[],[],[],[]
	boundaries = gmsh.model.getBoundary(volumes, oriented=False)
	for boundary in boundaries:
		center_of_mass = gmsh.model.occ.getCenterOfMass(boundary[0], boundary[1])
		if np.isclose(center_of_mass[1], Width):
			back.append(boundary[1])
		elif np.isclose(center_of_mass[1], 0):
			front.append(boundary[1])
		elif np.isclose(center_of_mass[0], 0):
			left.append(boundary[1])
		elif np.isclose(center_of_mass[0], Length):
			right.append(boundary[1])
		elif np.isclose(center_of_mass[2], 0):
			bottom.append(boundary[1])
		elif np.isclose(center_of_mass[2], Height) and center_of_mass[1]>Width/3: 
			top.append(boundary[1])
		else:
			indenter.append(boundary[1])
	# mark the surfaces
	gmsh.model.addPhysicalGroup(boundaries[0][0], bottom, bottom_marker)
	gmsh.model.setPhysicalName(boundaries[0][0], bottom_marker, 'bottom')
	gmsh.model.addPhysicalGroup(boundaries[0][0], front, front_marker)
	gmsh.model.setPhysicalName(boundaries[0][0], front_marker, 'front')
	gmsh.model.addPhysicalGroup(boundaries[0][0], back, back_marker)
	gmsh.model.setPhysicalName(boundaries[0][0], back_marker, 'back')
	gmsh.model.addPhysicalGroup(boundaries[0][0], left, left_marker)
	gmsh.model.setPhysicalName(boundaries[0][0], left_marker, 'left')
	gmsh.model.addPhysicalGroup(boundaries[0][0], right, right_marker)
	gmsh.model.setPhysicalName(boundaries[0][0], right_marker, 'right')
	gmsh.model.addPhysicalGroup(boundaries[0][0], top, top_marker)
	gmsh.model.setPhysicalName(boundaries[0][0], top_marker, 'top')
	gmsh.model.addPhysicalGroup(boundaries[0][0], indenter, indenter_marker)
	gmsh.model.setPhysicalName(boundaries[0][0], indenter_marker, 'indenter')
	gmsh.model.occ.synchronize()
	gmsh.write('Geom_2reelle_8EP.geo_unrolled')


	# Local refinement
	# Create the adaptive mesh around the indenter: https://jsdokken.com/src/tutorial_gmsh.html
	# # The final step is to set mesh resolutions. We will use GMSH Fields to do this. One can alternatively set mesh 
	# # resolutions at points with the command gmsh.model.occ.mesh.setSize. We start by specifying a distance field 
	# # from the obstacle surface
	indenter_interface = surfaces[0][1]
	distance = gmsh.model.mesh.field.add("Distance")
	gmsh.model.mesh.field.setNumbers(distance, "FacesList", [indenter_interface])

	# # # The next step is to use a threshold function vary the resolution from these surfaces in the following way:
	resolution = r/10
	threshold = gmsh.model.mesh.field.add("Threshold")
	gmsh.model.mesh.field.setNumber(threshold, "IField", distance)
	gmsh.model.mesh.field.setNumber(threshold, "LcMin", resolution)
	gmsh.model.mesh.field.setNumber(threshold, "LcMax", 5*resolution)
	gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.6*r)
	gmsh.model.mesh.field.setNumber(threshold, "DistMax", r)


	# if several fields defined:
	minimum = gmsh.model.mesh.field.add("Min")
	gmsh.model.mesh.field.setNumbers(minimum, "FieldsList", [threshold]) # add other fields in the list if needed
	gmsh.model.mesh.field.setAsBackgroundMesh(minimum)




	gmsh.model.occ.synchronize()
	gmsh.option.setNumber("General.Terminal",1)
	gmsh.option.setNumber("Mesh.Optimize", True)
	gmsh.option.setNumber("Mesh.OptimizeNetgen", True)
	# gmsh.option.setNumber("Mesh.MeshSizeMin", 0.5*lc)
	# gmsh.option.setNumber("Mesh.MeshSizeMax", 5*lc)
	gmsh.model.occ.synchronize()
	# gmsh.option.setNumber("Mesh.MshFileVersion", 2.0)
	gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
	gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
	gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

	gmsh.model.mesh.generate(3)
	gmsh.write("Mouse.msh")
	gmsh.finalize()

exit()