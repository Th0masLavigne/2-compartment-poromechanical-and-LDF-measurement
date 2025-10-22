# Thomas Lavigne: chargement instantanné pose problème
# 02-08-2024
# 
# From https://doi.org/10.1016/j.jmbbm.2023.105902
# Terzaghi Problem
#
# Boundary conditions:
#			   p=0
#			--------
#			|      |
#			|      |
#			|      |
#			|      |
#			|      |
#			|      |
#		ux=0|      |ux=0
#			|      |
#			|      |
#			|      |
#			|      |
#			|      |
#			|      |
#			--------
#			  uy=0
#
#------------------------------------------------------------#
#					Libraries								 #
#------------------------------------------------------------#
# 
try:
	from porous_fenicsx import constitutive_laws, variational_forms_TW, functions
except:
	print("Please run first `python3 -m pip install .` at ./Theoretical_developments_and_service_files/0_Service_files/porous_fenicsx to create the package.")
	exit()
# 
# 
import numpy
import time
import dolfinx
import ufl
import mpi4py
import basix
import petsc4py
# 
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
# 
# 
#------------------------------------------------------------#
#                    User-defined functions                  #
#------------------------------------------------------------#
# 
# 
def terzaghi_p(x):
	"""
	Compute the Exact Terzaghi solution
	Inputs: coordinates
	Outputs: Fluid pressure 
	"""
	kmax = 1e3
	p0,L = pinit, Height
	# Lamé coefficients
	lambda_m, mu = constitutive_laws.Lame_coefficients(E,nu)
	cv   = k_l/mu_l*(lambda_m+2*mu)
	pression = 0
	for k in range(1,int(kmax)):
		pression += p0*4/numpy.pi*(-1)**(k-1)/(2*k-1)*numpy.cos((2*k-1)*0.5*numpy.pi*(x[1]/L))*numpy.exp(-(2*k-1)**2*0.25*numpy.pi**2*cv*t/L**2)
	pl = pression
	return pl
# 
def L2_error_p(mesh,P1,__p):
	"""
	Define the L2_error computation
	Inputs: Mesh, type of element, solution function
	Outputs: L2 error to the analytical solution
	"""
	P1space = dolfinx.fem.functionspace(mesh, P1)
	p_theo  = dolfinx.fem.Function(P1space)
	p_theo.interpolate(terzaghi_p)
	L2_errorp, L2_normp = dolfinx.fem.form(ufl.inner(__p - p_theo, __p - p_theo) * dx), dolfinx.fem.form(ufl.inner(p_theo, p_theo) * dx)
    num_local = dolfinx.fem.assemble_scalar(L2_errorp)
    den_local = dolfinx.fem.assemble_scalar(L2_normp)

    num = mesh.comm.allreduce(num_local, op=mpi4py.MPI.SUM)
    den = mesh.comm.allreduce(den_local, op=mpi4py.MPI.SUM)
    return numpy.sqrt(num / den)
# 
def evaluate_point(mesh, function, contributing_cells, point, output_list, index):
	"""
	Suitable Evaluations functions for Parallel computation
	Inputs: mesh, function to evaluate, contributing cells to the point, point, output list to store the value, index in the list
	Outputs: the evaluated function value is added at the index location in output list
	"""
	from mpi4py            import MPI
	function_eval = None
	if len(contributing_cells) > 0:
		function_eval = function.eval(point, contributing_cells[:1])
	function_eval = mesh.comm.gather(function_eval, root=0)
	# Choose first pressure that is found from the different processors
	if mpi4py.MPI.COMM_WORLD.rank == 0:
		for element in function_eval:
			if element is not None:
				output_list[index]=element[0]
				break
	pass
# Terzaghi analytical solution
def terzaghi(p0,L,cv,y,t,kmax):
	"""
	y as the position, t as the time we are looking to
	p0 the applied pressure
	L the sample's length
	cv the consolidation time
	"""
	pression=0
	for k in range(1,kmax):
		pression += p0*4/numpy.pi*(-1)**(k-1)/(2*k-1)*numpy.cos((2*k-1)*0.5*numpy.pi*(y/L))*numpy.exp(-(2*k-1)**2*0.25*numpy.pi**2*cv*t/L**2)
	pl = pression
	return pl
# 
#------------------------------------------------------------#
#                    Set time counter		                 #
#------------------------------------------------------------#
# 
begin_t = time.time()
# 
#------------------------------------------------------------#
#     		               FE_Mesh 			                 #
#------------------------------------------------------------#
# 
# Create the domain / mesh
Height = 1e-4 #[m]
Width  = 1e-5 #[m]
mesh   = dolfinx.mesh.create_rectangle(mpi4py.MPI.COMM_WORLD, numpy.array([[0,0],[Width, Height]]), [2,100], cell_type=dolfinx.mesh.CellType.triangle)
# 
## Define the boundaries:
# 1 = bottom, 2 = right, 3=top, 4=left
boundaries = [(1, lambda x: numpy.isclose(x[1], 0)),
				(2, lambda x: numpy.isclose(x[0], Width)),
				(3, lambda x: numpy.isclose(x[1], Height)),
				(4, lambda x: numpy.isclose(x[0], 0))]
# 
facet_indices, facet_markers = [], []
fdim = mesh.topology.dim - 1
for (marker, locator) in boundaries:
	facets = dolfinx.mesh.locate_entities_boundary(mesh, fdim, locator)
	facet_indices.append(facets)
	facet_markers.append(numpy.full_like(facets, marker))
# Concatenate and sort the arrays based on facet indices. Left facets marked with 1, right facets with two
facet_indices = numpy.hstack(facet_indices).astype(numpy.int32)
facet_markers = numpy.hstack(facet_markers).astype(numpy.int32)
sorted_facets = numpy.argsort(facet_indices)
facet_tag     = dolfinx.mesh.meshtags(mesh, fdim, facet_indices[sorted_facets], facet_markers[sorted_facets])
# 
# 
with dolfinx.io.XDMFFile(mpi4py.MPI.COMM_WORLD, "_column_incremental_1_mono_tags.xdmf", "w") as xdmf:
	xdmf.write_mesh(mesh)
	xdmf.write_meshtags(facet_tag,mesh.geometry)
	xdmf.close()
# 
#----------------------------------------------------------------------
# Operators
# Create the surfacic element
metadata = {"quadrature_degree": 4}
dx       = ufl.Measure("dx", domain=mesh, metadata=metadata)
ds       = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tag)
# compute the mesh normals to express t^imposed = T.normal
normal   = ufl.FacetNormal(mesh)
# 
#------------------------------------------------------------#
#     		          Evaluation Areas			             #
#------------------------------------------------------------#
# 
# Identify contributing cells to our points of interest for post processing
num_points = 11
# Physical points we want an evaluation in
y_check          = numpy.linspace(0,Height,num_points)
points_for_time  = numpy.array([[Width/2, 0., 0.], [Width/2, Height/2, 0.]])
points_for_space = numpy.zeros((num_points,3))
for ii in range(num_points):
	points_for_space[ii,0] = Width/2
	points_for_space[ii,1] = y_check[ii]
# Create the bounding box tree
tree             = dolfinx.geometry.bb_tree(mesh, mesh.geometry.dim)
points           = numpy.concatenate((points_for_time,points_for_space))
cell_candidates  = dolfinx.geometry.compute_collisions_points(tree, points)
colliding_cells  = dolfinx.geometry.compute_colliding_cells(mesh, cell_candidates, points)
cells_y_0        = colliding_cells.links(0)
cells_y_H_over_2 = colliding_cells.links(1)
# 
#------------------------------------------------------------#
#         Time and Mechanical load parameters                #
#------------------------------------------------------------#
# 
# Time parameters
# Start time [s]
t         = 0             
# End time [s]
Tf        = 6.           
# Number of time steps
num_steps = 1000
# Time step [s]
dt        = dolfinx.default_scalar_type((Tf-t)/num_steps) 
# 
# Mechanical loading 
pinit = 100 #[Pa]
T     = dolfinx.fem.Constant(mesh,dolfinx.default_scalar_type(-pinit))
# 
#------------------------------------------------------------#
#                   Material Definition                      #
#------------------------------------------------------------#
# 
# Solid scaffold
# Young Moduli [Pa]
E            = dolfinx.default_scalar_type(5000)
# Possion's ratios [-]
nu           = dolfinx.default_scalar_type(0.4)
# 
# Porous material
# IF fluid viscosity [Pa.s]
mu_l    = dolfinx.default_scalar_type(1e-2)
# IF Intrinsic permeabilitty [m^2]
k_l     = dolfinx.default_scalar_type(1.8e-15)
# 
#------------------------------------------------------------#
#                   Function Spaces                          #
#------------------------------------------------------------#
# 
# Mixed Space (R,R,R2) -> (pc,plc,u)
P1_v     = basix.ufl.element("P", mesh.topology.cell_name(), degree=1, shape=(mesh.topology.dim,))
P2_v     = basix.ufl.element("P", mesh.topology.cell_name(), degree=2, shape=(mesh.topology.dim,))
P1       = basix.ufl.element("P", mesh.topology.cell_name(), degree=1)
# 
# Function spaces
# Updated mesh space
updated_mesh_space   = dolfinx.fem.functionspace(mesh, mesh.ufl_domain().ufl_coordinate_element())
# 
P1v_space     = dolfinx.fem.functionspace(mesh, P1_v)
P2v_space     = dolfinx.fem.functionspace(mesh, P2_v)
P1_space      = dolfinx.fem.functionspace(mesh, P1)
MS            = dolfinx.fem.functionspace(mesh=mesh, element=basix.ufl.mixed_element([P1,P2_v]))
#  
#------------------------------------------------------------#
#                           Functions                        #
#------------------------------------------------------------#
# 
# Create the function update the nodal positions of the mesh
du_update                = dolfinx.fem.Function(updated_mesh_space)
# Create the function to export a P2 displacement in the xdmf file
displacement_export      = dolfinx.fem.Function(P1v_space)
displacement_export.name = "Displacement"
# 
# Define the solution at the current and previous time steps
solution          = dolfinx.fem.Function(MS)
previous_solution = dolfinx.fem.Function(MS)
# 
# 
# Mapping in the Mixed Space: FunctionSpace, Mapping_in_the_mixed_space = MS.sub(xx).collapse()
Pln_, Pln_to_MS   = MS.sub(0).collapse()
Un_, Un_to_MS     = MS.sub(1).collapse()
#------------------------------------------------------------#
#                   Initial Conditions                       #
#------------------------------------------------------------#
# 
# pl=pinit
previous_solution.x.array[Pln_to_MS] = numpy.full_like(previous_solution.x.array[Pln_to_MS], pinit, dtype=dolfinx.default_scalar_type)
#------------------------------------------------------------#
#                          Expressions                       #
#------------------------------------------------------------#
#
# If updated form
# Computation of the pressure at bottom points
S_bottom = surface_indenter = mesh.comm.allreduce(dolfinx.fem.assemble_scalar(dolfinx.fem.form(1*ds(1))), op=mpi4py.MPI.SUM)
Press_bottom_expr_IF        = dolfinx.fem.form(1/S_bottom*(previous_solution.sub(0)+solution.sub(0))*ds(1))
# Computation of the displacement at top points
S_top = surface_indenter = mesh.comm.allreduce(dolfinx.fem.assemble_scalar(dolfinx.fem.form(1*ds(3)) ), op=mpi4py.MPI.SUM)
Disp_top_expr            = dolfinx.fem.form(1/S_top*(ufl.dot(previous_solution.sub(1)+solution.sub(1),normal))*ds(3))
# 
#------------------------------------------------------------#
#                   Create the export file                   #
#------------------------------------------------------------#
#
# Open file for export
xdmf = dolfinx.io.XDMFFile(mesh.comm, "_column_incremental_1_mono_Result_constant_load.xdmf", "w")
xdmf.write_mesh(mesh)
# 
# Create output lists in time and space for the IF pressure
# time steps to evaluate the pressure in space:
n0, n1, n2 = 200,400,800
# 
displacement_all, pressure_IF_all, time_all = [], [], []
# 
pressure_y_0             = numpy.zeros(num_steps, dtype=dolfinx.default_scalar_type)
pressure_y_Height_over_2 = numpy.zeros(num_steps, dtype=dolfinx.default_scalar_type)
pressure_space0          = numpy.zeros(num_points, dtype=dolfinx.default_scalar_type)
pressure_space1          = numpy.zeros(num_points, dtype=dolfinx.default_scalar_type)
pressure_space2          = numpy.zeros(num_points, dtype=dolfinx.default_scalar_type)
# 
#------------------------------------------------------------#
#                   Initial Conditions                       #
#------------------------------------------------------------#
# None
# 
#------------------------------------------------------------#
#              Dirichlet Boundary Conditions                 #
#------------------------------------------------------------#
#
# 1 = bottom: uy=0, 2 = right: ux=0, 3=top: pl=0 leakage, 4=left: ux=0
bcs    = []
fdim   = mesh.topology.dim - 1
# duy=0
facets = facet_tag.find(1)
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(1).sub(1), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(1).sub(1)))
# dux=0
facets = facet_tag.find(2)
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(1).sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(1).sub(0)))
# dux=0
facets = facet_tag.find(4)
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(1).sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(1).sub(0)))
# leakage dpl=pinit initially then 0
facets = facet_tag.find(3)
p_imp  = dolfinx.fem.Constant(mesh,dolfinx.default_scalar_type(-pinit))
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(p_imp, dofs, MS.sub(0)))
# 
#------------------------------------------------------------#
#                     Variationnal form                      #
#------------------------------------------------------------#
#
# Create the test functions 
ql, w  = ufl.TestFunctions(MS)
# 
F = variational_forms_TW.variational_form_1_mono_updated(previous_solution, solution, ql, w, E, nu, k_l, mu_l, dx, ds, dt, constitutive_laws.Elastic_constitutive_law)
# Add Neuman BCs
F+= - T*ufl.inner(w,normal)*ds(3)
# 
#------------------------------------------------------------#
#                           Solver                           #
#------------------------------------------------------------#
# 
# Non linear problem definition
dsolution = ufl.TrialFunction(MS)
J         = ufl.derivative(F, solution, dsolution)
Problem   = NonlinearProblem(F, solution, bcs = bcs, J = J)
solver    = functions.set_non_linear_solver_parameters(mesh, Problem, 5e-10, 1e-11, "incremental", 15, log_newton=False)
# 
#------------------------------------------------------------#
#                         Computation                        #
#------------------------------------------------------------#
# 
t = 0
L2_p = numpy.zeros(num_steps, dtype=dolfinx.default_scalar_type)
for n in range(num_steps):
	if n == 1:
		p_imp.value = 0
	t+=dt
	# t = round(t,2)
	time_all.append(t)
	try:
		num_its, converged = solver.solve(solution)
	except:
		if mpi4py.MPI.COMM_WORLD.rank == 0:
			print("*************") 
			print("Solver failed")
			print("*************") 
			break
	# 
	solution.x.scatter_forward()
	# 
	# Evaluate the bottom pressure
	PBot = dolfinx.fem.assemble_scalar(Press_bottom_expr_IF)
	pressure_IF_all.append(mesh.comm.allreduce(PBot, op=mpi4py.MPI.SUM))
	# Evaluate the top displacement
	DTop = dolfinx.fem.assemble_scalar(Disp_top_expr)
	displacement_all.append(mesh.comm.allreduce(DTop, op=mpi4py.MPI.SUM))
	# 
	# Update Value for updated variational form
	previous_solution.x.array[:] += solution.x.array[:]
	# 
	previous_solution.x.scatter_forward()
	__p, __u = previous_solution.split()
	__p.name = "Pressure"
	# Displacement to P1 for export
	displacement_export.interpolate(previous_solution.sub(1))
	displacement_export.x.scatter_forward()
	# 
	# Export the results
	xdmf.write_function(displacement_export,t)
	xdmf.write_function(__p,t)
	# 
	# Compute L2 norm for pressure
	error_L2p     = L2_error_p(mesh,P1,__p)
	L2_p[n] = error_L2p
	# 
	# Solve tracking
	if mpi4py.MPI.COMM_WORLD.rank == 0:
		print(f"Time step {n}/{num_steps}, Load {T.value}, L2-error p {L2_p[n]:.2e}") 
		print(f"Pressure IF at bottom points: {pressure_IF_all[-1]:.2e}") 
		print(f"Displacement at top points: {displacement_all[-1]:.2e}") 
	# Evaluate the functions
	# in time
	evaluate_point(mesh, __p, cells_y_0, points[0], pressure_y_0, n)
	evaluate_point(mesh, __p, cells_y_H_over_2, points[1], pressure_y_Height_over_2, n)
	# in space
	if n == n0:
		for ii in range(num_points):
			evaluate_point(mesh, __p, colliding_cells.links(ii+2), points[ii+2], pressure_space0, ii)
		t0 = t
	elif n==n1:
		for ii in range(num_points):
			evaluate_point(mesh, __p, colliding_cells.links(ii+2), points[ii+2], pressure_space1, ii)
		t1 = t
	elif n==n2:
		for ii in range(num_points):
			evaluate_point(mesh, __p, colliding_cells.links(ii+2), points[ii+2], pressure_space2, ii)
		t2 = t
	# Update the mesh coordinates 
	du_update.interpolate(solution.sub(1))
	mesh.geometry.x[:,:mesh.geometry.dim] += du_update.x.array.reshape((-1, mesh.geometry.dim))
xdmf.close()
# 
# 
#------------------------------------------------------------#
#                       Post-Processing                      #
#------------------------------------------------------------#
# 
if mpi4py.MPI.COMM_WORLD.rank == 0:
	print(f"L2 error p, min {numpy.min(L2_p):.2e}, mean {numpy.mean(L2_p):.2e}, max {numpy.max(L2_p):.2e}, std {numpy.std(L2_p):.2e}")
# 
# Evaluate final time
end_t = time.time()
t_hours = int((end_t-begin_t)//3600)
tmin = int(((end_t-begin_t)%3600)//60)
tsec = int(((end_t-begin_t)%3600)%60)
if mpi4py.MPI.COMM_WORLD.rank == 0:
	print(f"FEM operated with {num_steps} iterations, in {t_hours} h {tmin} min {tsec} sec")
	# 
	###################################################
	################ Analytical solutions #############
	###################################################
	# 
	# Lamé coefficients
	lambda_m, mu = constitutive_laws.Lame_coefficients(E,nu)
	cv = k_l/mu_l*(lambda_m+2*mu)
	y=0
	pressure4 = numpy.zeros(num_steps)
	kmax=1e3
	for i in range(num_steps):
		pressure4[i] = terzaghi(pinit,Height,cv,y,time_all[i],int(kmax))
	# 
	y=Height/2
	pressure5 = numpy.zeros(num_steps)
	for i in range(num_steps):
		pressure5[i] = terzaghi(pinit,Height,cv,y,time_all[i],int(kmax))
	# 
	pressure0 = numpy.zeros(num_points)
	for i in range(num_points):
		pressure0[i] = terzaghi(pinit,Height,cv,y_check[i],t0,int(kmax))
	# 
	pressure1 = numpy.zeros(num_points)
	for i in range(num_points):
		pressure1[i] = terzaghi(pinit,Height,cv,y_check[i],t1,int(kmax))
	# 
	pressure2 = numpy.zeros(num_points)
	for i in range(num_points):
		pressure2[i] = terzaghi(pinit,Height,cv,y_check[i],t2,int(kmax))
	# 
	###################################################
	################ Plots ############################
	###################################################
	# 
	import matplotlib.pyplot as plt
	# 
	plt.rcParams.update({'font.size': 15})
	plt.rcParams.update({'legend.loc':'upper right'})
	# 
	fig1, ax1 = plt.subplots()
	ax1.plot(time_all,pressure4,linestyle='-',linewidth=2,label='Analytic y=0',color='powderblue')
	ax1.plot(time_all,pressure5,linestyle='-',linewidth=2,label='Analytic y=h/2',color='bisque')
	ax1.plot(time_all,pressure_IF_all,linestyle='-',linewidth=2,color='navy', label='IF')
	ax1.plot(time_all,pressure_y_0,linestyle=':',linewidth=2,label='FEniCSx y=0',color='cornflowerblue')
	ax1.plot(time_all,pressure_y_Height_over_2,linestyle=':',linewidth=2,label='FEniCSx y=h/2',color='salmon')
	ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
	ax1.set_xlabel('time (s)')
	ax1.set_ylabel('Pressure (Pa)')
	ax1.legend()
	fig1.tight_layout()
	fig1.savefig('_column_incremental_1_mono_Pressure_time.jpg')
	# 
	fig2, ax2 = plt.subplots()
	ax2.plot(pressure0,y_check,linestyle='-',linewidth=2,color='lightgreen',label='Analytic')
	ax2.plot(pressure1,y_check,linestyle='-',linewidth=2,color='lightgreen')
	ax2.plot(pressure2,y_check,marker='',linestyle='-',linewidth=2,color='lightgreen')
	ax2.plot(pressure_space0,y_check,linestyle=':',linewidth=2,color='olivedrab',label='FEniCSx')
	ax2.plot(pressure_space1,y_check,linestyle=':',linewidth=2,color='olivedrab')
	ax2.plot(pressure_space2,y_check,linestyle=':',linewidth=2,color='olivedrab')
	ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
	tt=plt.text(12, 4e-5, f't={numpy.round(t2,1)}s', fontsize = 12 )
	tt=plt.text(35, 4e-5,f't={numpy.round(t1,1)}s', fontsize = 12)
	tt.set_bbox(dict(facecolor='white', alpha=0.7, linewidth=0))
	tt=plt.text(60, 4e-5, f't={numpy.round(t0,1)}s', fontsize = 12)
	tt.set_bbox(dict(facecolor='white', alpha=0.7, linewidth=0))
	ax2.set_xlabel('Pressure (Pa)')
	ax2.set_ylabel('Height (m)')
	ax2.legend()
	fig2.tight_layout()
	fig2.savefig('_column_incremental_1_mono_Pressure_Space.jpg')
	# 
	# 
	# 
	fig1, ax1 = plt.subplots()
	ax1.plot(time_all,displacement_all,linestyle='-',linewidth=2,color='salmon')
	ax1.set_xlabel('time (s)')
	ax1.set_ylabel('u (m)')
	# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
	ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
	fig1.tight_layout()
	fig1.savefig('_column_incremental_1_mono_Displacement_top.jpg')
	# 
	def export_to_csv(data, filename, header=None):
	    import csv
	    try:
	        with open(filename, 'w', newline='') as file:
	            writer = csv.writer(file)
	            if header:
	                writer.writerow(header)
	            writer.writerows(data)
	        print(f"Data exported to {filename} successfully")
	    except Exception as e:
	        print(f"An error occurred while exporting data to {filename}: {e}")
	# 
	export_to_csv([y_check,pressure0,pressure1,L2_p],"_column_incremental_1_mono_Results.csv",["y","pressure0","pressure1","L2P"])
# 
# Evaluate the error
RMSE_p = functions.RMSE(pressure_IF_all,pressure4)
# 
if mpi4py.MPI.COMM_WORLD.rank == 0:
	print(f"RMSE of the pressure at bottom points {RMSE_p} Pa, RMSE/max(p_theo) = {100*RMSE_p/max(pressure4)} %")
# 
# EoF
