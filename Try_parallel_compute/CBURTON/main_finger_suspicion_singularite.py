# Thomas Lavigne
# 18-09-2024
# 
#------------------------------------------------------------#
#                   Libraries 					             #
#------------------------------------------------------------#
# 
import constitutive_laws
import variational_forms
import functions
# 
# FEniCSx Required Libraries
import dolfinx
import ufl
import basix
import mpi4py
from dolfinx.fem.petsc import NonlinearProblem
# 
# Additionnal Libraries
import numpy
import time
# 
# 
#------------------------------------------------------------#
#            Check Version and initialise time counter       #
#------------------------------------------------------------#
# 
# Print in the log the FEniCSx version loaded
if mpi4py.MPI.COMM_WORLD.rank == 0:
	print("Dolfinx version is:",dolfinx.__version__)
# 
# Set time counter
begin_t = time.time()
#
#------------------------------------------------------------#
#                   Load the Geometry from GMSH              #
#------------------------------------------------------------#
# 
# Mesh filename. This mesh is supposed to be a 3D gmsh file
filename = "./Mesh_2_poro.msh"
# 
# Load the mesh and tags in the FEniCSx environment
mesh, cell_tag, facet_tag = dolfinx.io.gmshio.read_from_msh(filename, mpi4py.MPI.COMM_WORLD, 0, gdim=3)
# 
#--------------------------------------------------------------
#______________________________________________________________
#--------------------------------------------------------------
# Identify indices of the cells for each region for material definition
tissue_indices   = [x for x in cell_tag.indices if (cell_tag.values[x] == 1)]
LDF_indices      = [x for x in cell_tag.indices if (cell_tag.values[x] == 2)]
# 
# try :
# 	assert(len(cell_tag.indices) == len(tissue_indices)+len(LDF_indices))
# 	if mpi4py.MPI.COMM_WORLD.rank       == 0:
# 		print("All cell tags have been attributed")
# except:
# 	if mpi4py.MPI.COMM_WORLD.rank       == 0:
# 		print("*************") 
# 		print("Forgotten tags => material badly defined")
# 		print("*************") 
# 		exit()
# 
#--------------------------------------------------------------
#______________________________________________________________
#--------------------------------------------------------------
# Identification of the useful elements from the mesh for the integrals
normal = ufl.FacetNormal(mesh)
# Specify the desired quadrature degree
q_deg  = 4
# attribute the cell_tags to the integral element
dx     = ufl.Measure('dx', metadata={"quadrature_degree":q_deg}, subdomain_data=cell_tag, domain=mesh)
# attribute the facet_tags to the surfacic integral element
ds     = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tag)
# 
#------------------------------------------------------------#
#         Time and Mechanical load parameters                #
#------------------------------------------------------------#
# 
# Initialise mechanical load
load_magnitude = 31.6e3 #[Pa]
# Applied surfacic force for Neumann BCs
T              = dolfinx.fem.Constant(mesh,dolfinx.default_scalar_type(0))
# 
# initial time for imposing the blood pressure gradient gradually
index_ti = 0
ti       = 5
t_ramp   = 20 
t_sust   = 60.
# Time parametrization
t_init   = 0                              # Start time
Tf       = 850         # End time
# Must ensure t%dt==0
dt       = 1.0
t        = t_init
num_steps_pb_imp = int(ti/dt)
num_steps = int((Tf-t_init+ti)/dt)
num_steps = 150
# 
#------------------------------------------------------------#
#                   Material Definition                      #
#------------------------------------------------------------#
# 
# Solid scaffold
# Young Moduli [Pa]
E         = dolfinx.default_scalar_type(1e5) 
# Possion's ratios [-]
nu        = dolfinx.default_scalar_type(0.42)
# Porous material
#
# Extra vascular porosity (initial value)
varepsilon     = dolfinx.default_scalar_type(0.3)
varepsilon_min = dolfinx.default_scalar_type(0.1)
varepsilon_max = dolfinx.default_scalar_type(0.6)
# 
# IF fluid viscosity [Pa.s]
mu_l      = dolfinx.default_scalar_type(1) 
# IF Intrinsic permeabilitty [m^2]
k_l       = dolfinx.default_scalar_type(1e-13)
# 
# Blood fluid viscosity [Pa.s]
mu_b      = dolfinx.default_scalar_type(1)  
# Vascular Intrinsic permeabilitty [m^2]
k_b       = dolfinx.default_scalar_type(1e-9)
# To define a preferential direction for function matrix_vascular_permeability
ratio     = [1, 1, 1]
# exponent to make evolution of kb according to poro
alpha     =  dolfinx.default_scalar_type(5)
# Vessel Compressibility [Pa]
K         = dolfinx.default_scalar_type(numpy.pi/2*1000)
# 
# BLOOD FLOW [m3/sec * 1/m2]
# BFLOW  = dolfinx.default_scalar_type(4.e-4)
# PB_IMP = dolfinx.fem.Constant(mesh,  dolfinx.default_scalar_type(0.2))
PB_IMP = dolfinx.fem.Constant(mesh,  dolfinx.default_scalar_type(4e3))
dPB_IMP = dolfinx.fem.Constant( mesh,  dolfinx.default_scalar_type( PB_IMP.value/num_steps_pb_imp ) )
if mpi4py.MPI.COMM_WORLD.rank == 0:
	print("pression sanguine imposee= ",PB_IMP.value)
# 
# 
# Vascular porosity (initial value)
initial_varepsilonb = 0.08
varepsilonb = dolfinx.default_scalar_type( initial_varepsilonb )
evolution_rule = 'atan'
# 
# 
# 
#------------------------------------------------------------#
#                   Function Spaces                          #
#------------------------------------------------------------#
# 
# Mixed Space (R,R,R2) -> (pl,pb,u)
P1_v     = basix.ufl.element("P", mesh.topology.cell_name(), degree=1, shape=(mesh.topology.dim,))
P2       = basix.ufl.element("P", mesh.topology.cell_name(), degree=2, shape=(mesh.topology.dim,))
P1       = basix.ufl.element("P", mesh.topology.cell_name(), degree=1)
# 
# Function spaces
# Updated mesh space
updated_mesh_space   = dolfinx.fem.functionspace(mesh, mesh.ufl_domain().ufl_coordinate_element())
# 
P1v_space     = dolfinx.fem.functionspace(mesh, P1_v)
P2_space      = dolfinx.fem.functionspace(mesh, P2)
P1_space      = dolfinx.fem.functionspace(mesh, P1)
MS            = dolfinx.fem.functionspace(mesh=mesh, element=basix.ufl.mixed_element([P1,P1,P2]))
#  
#------------------------------------------------------------#
#                           Functions                        #
#------------------------------------------------------------#
# 
# 
# Create the function update the nodal positions of the mesh
du_update                = dolfinx.fem.Function(updated_mesh_space)
# Create the function to export a P2 displacement in the xdmf file
displacement_export      = dolfinx.fem.Function(P1v_space)
displacement_export.name = "Displacement"
# 
# Create the previous solution of the porosity
# porosity_n        = dolfinx.fem.Function(P1_space)
# porosity_n.name   = "porosity"
# Create the previous solution of the vascular porosity
porosity_b_n        = dolfinx.fem.Function(P1_space)
porosity_b_n.name   = "vascular porosity"
# 
# Define the solution at the current and previous time steps
solution          = dolfinx.fem.Function(MS)
previous_solution = dolfinx.fem.Function(MS)
# 
# 
# Mapping in the Mixed Space: FunctionSpace, Mapping_in_the_mixed_space = MS.sub(xx).collapse()
Pln_, Pln_to_MS   = MS.sub(0).collapse()
Pbn_, Pbn_to_MS   = MS.sub(1).collapse()
Un_, Un_to_MS     = MS.sub(2).collapse()
# 
# Create the test functions
ql,qb,w = ufl.TestFunctions(MS)
# 
# post_processing functions
flux_blood      = dolfinx.fem.Function(P1v_space)
flux_blood.name ="blood flux"
#------------------------------------------------------------#
#                          Expressions                       #
#------------------------------------------------------------#
# Useful variables:
# global directional vectors
Nx                  = dolfinx.fem.Constant(mesh, numpy.asarray((1.0,0.0,0.0)))
Ny                  = dolfinx.fem.Constant(mesh, numpy.asarray((0.0,1.0,0.0)))
Nz                  = dolfinx.fem.Constant(mesh, numpy.asarray((0.0,0.0,1.0)))
# 
k_b_current          = constitutive_laws.matrix_vascular_permeability_mono(k_b,ratio,alpha,varepsilonb,previous_solution.sub(0)+solution.sub(0),previous_solution.sub(1)+solution.sub(1),K,evolution_rule)
# 
if evolution_rule == 'linear':
	varepsilon_b_current = constitutive_laws.vascular_porosity_linear_mono(varepsilonb,previous_solution.sub(0)+solution.sub(0),previous_solution.sub(1)+solution.sub(1),K)
elif evolution_rule == 'atan':
	varepsilon_b_current = constitutive_laws.vascular_porosity_atan_mono(varepsilonb,previous_solution.sub(0)+solution.sub(0),previous_solution.sub(1)+solution.sub(1),K)
#
darcy_blood      = -1 * k_b_current /mu_b*ufl.grad(previous_solution.sub(1)+solution.sub(1))
# Evaluation of volume and surfaces to integrate the solutions
Volume_LDF       = mesh.comm.allreduce(dolfinx.fem.assemble_scalar(dolfinx.fem.form(1*dx(2)) ), op=mpi4py.MPI.SUM)
surface_indenter = mesh.comm.allreduce(dolfinx.fem.assemble_scalar(dolfinx.fem.form(1*ds(3)) ), op=mpi4py.MPI.SUM)
# 
# Expressions to be computed before the update of previous solution
# Evaluation of the porosity 
# porosity_n_expr   = dolfinx.fem.Expression(varepsilon_current, P1_space.element.interpolation_points())
# Evaluation of the vascular porosity
porosity_b_n_expr = dolfinx.fem.Expression(varepsilon_b_current, P1_space.element.interpolation_points())
# Evaluation of the blood flux qb = -kb/mu * (vb-vs)
flux_blood_expr   = dolfinx.fem.Expression(darcy_blood , P1v_space.element.interpolation_points())
# vbs = qb/epsb
# speed_blood_expr = dolfinx.fem.Expression(1/varepsilon_b_current *darcy_blood , P1v_space.element.interpolation_points())
# 
# Only in x (in reality q not v just to compare on z direction only)
# LDF_expr_v                    = dolfinx.fem.form(1/Volume_LDF*(ufl.dot(darcy_blood/varepsilon_b_current,Nz))*dx(2))
LDF_expr_v                    = dolfinx.fem.form(1/Volume_LDF*( ufl.sqrt( ufl.dot(darcy_blood,Nz)**2+ufl.dot(darcy_blood,Nz)**2 ) )*dx(2))
LDF_expr_q                    = dolfinx.fem.form(1/Volume_LDF*( ufl.dot(darcy_blood,Nz) )*dx(2))
# Computation of the displacement at top points
Disp_top_expr                 = dolfinx.fem.form(1/surface_indenter*(ufl.dot(previous_solution.sub(2)+solution.sub(2),normal))*ds(3))
# 
#------------------------------------------------------------#
#                   Create the export file                   #
#------------------------------------------------------------#
#
# Open file for export
# xdmf = dolfinx.io.XDMFFile(mesh.comm, "ID_63_.xdmf", "w")
# xdmf.write_mesh(mesh)
# 
#------------------------------------------------------------#
#                   Initial Conditions                       #
#------------------------------------------------------------#
# 
# Porosity
# porosity_n.x.array[:]   = numpy.full_like(porosity_n.x.array[:], varepsilon, dtype=dolfinx.default_scalar_type)
# Vascular Porosity
porosity_b_n.x.array[:] = numpy.full_like(porosity_b_n.x.array[:], varepsilonb, dtype=dolfinx.default_scalar_type)
# 
#------------------------------------------------------------#
#              Dirichlet Boundary Conditions                 #
#------------------------------------------------------------#
# 
bcs    = []
fdim   = mesh.topology.dim - 1
# 
# symmetry_bottom, symmetry, load_ldf, top, left, right, up_cyl, bone = 1,2,3,4,5,6,8,7
# 
# symmetry_bottom
facets = facet_tag.find(1)
# duy=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(2).sub(1), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(2).sub(1)))
# symmetry
facets = facet_tag.find(2)
# dux=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(2).sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(2).sub(0)))
# bone
facets = facet_tag.find(7)
# dux=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(2).sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(2).sub(0)))
# duy=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(2).sub(1), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(2).sub(1)))
# duz=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(2).sub(2), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(2).sub(2)))
# left
facets = facet_tag.find(5)
# duz=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(2).sub(2), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(2).sub(2)))
# dpl=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(0)))
# dpb = dpbimp
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(1), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dPB_IMP, dofs, MS.sub(1)))
# right
facets = facet_tag.find(6)
# duz=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(2).sub(2), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(2).sub(2)))
# dpl=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(0)))
# dpb = dpbimp
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(1), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(1)))
# 
#------------------------------------------------------------#
#                     Variationnal form                      #
#------------------------------------------------------------#
#
# 
F = variational_forms.variational_form_2_mono_updated_anisotropic_k_b(previous_solution, solution, ql, qb, w, E, nu, k_l, mu_l, k_b, ratio, mu_b, varepsilonb, K, alpha, dx, ds, dt, constitutive_laws.Elastic_constitutive_law, evolution_rule)
#
# Add Neuman BCs on indenter surface
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
solver    = functions.set_non_linear_solver_parameters(mesh, Problem, 5e-9, 1e-10, "incremental", 10, log_newton=False)
# 
#------------------------------------------------------------#
#                         Computation                        #
#------------------------------------------------------------#
# Check initial sol
# Split the total solution for export
__pl, __pb, __u = previous_solution.split()
__pl.name              = "IF pressure"
__pb.name              = "Blood pressure"
# Export in the xdmf file 
xdmf.write_function(displacement_export,t)
# xdmf.write_function(porosity_n,t)
xdmf.write_function(porosity_b_n,t)
xdmf.write_function(__pl,t)
xdmf.write_function(__pb,t)
xdmf.write_function(flux_blood,t)
# 
displacement_all, LDF_v_all, LDF_q_all, load_all, time_all = [], [], [], [], []
# 
for n in range(num_steps):
	# Update BCs
	if n == num_steps_pb_imp:
		dPB_IMP.value = 0
	# update time
	t+=dt
	t = round(t,2)
	time_all.append(t)
	if t == ti:
		index_ti = n
  	# update the load value and export it in the list
	T.value = functions.mechanical_load_LDF_2(t,load_magnitude,ti)
	load_all.append(-functions.mechanical_load_LDF_2(t,load_magnitude,ti))
	# 
	# Solve
	if mpi4py.MPI.COMM_WORLD.rank == 0:
		print('Step: ',n+1,'/',num_steps, 'time:', t, ' s',' load: ', T.value, ' Pa')
	try:
		num_its, converged = solver.solve(solution)
	except:
		if mpi4py.MPI.COMM_WORLD.rank == 0:
			print("*************") 
			print("Solver failed")
			print("*************") 
		break
	# Ensure pushing the solution to all processes
	solution.x.scatter_forward()
	# Update the porosity
	# porosity_n.interpolate(porosity_n_expr)
	# porosity_n.x.scatter_forward()
	# Update vascular porosity
	porosity_b_n.interpolate(porosity_b_n_expr)
	porosity_b_n.x.scatter_forward()
	# Update the blood flux
	flux_blood.interpolate(flux_blood_expr)
	flux_blood.x.scatter_forward()
	# 
	# Evaluate the LDF in the cylinder
	ExpLDF_v = dolfinx.fem.assemble_scalar(LDF_expr_v)
	LDF_v_all.append(mesh.comm.allreduce(ExpLDF_v, op=mpi4py.MPI.SUM))
	ExpLDF_q = dolfinx.fem.assemble_scalar(LDF_expr_q)
	LDF_q_all.append(mesh.comm.allreduce(ExpLDF_q, op=mpi4py.MPI.SUM))
	# Evaluate the top displacement
	DTop     = dolfinx.fem.assemble_scalar(Disp_top_expr)
	displacement_all.append(mesh.comm.allreduce(DTop, op=mpi4py.MPI.SUM))
	# 
	# Update previous solution
	previous_solution.x.array[:] += solution.x.array[:]
	previous_solution.x.scatter_forward()
  	# 
  	# Displacement to P1 for export
	displacement_export.interpolate(previous_solution.sub(2))
	displacement_export.x.scatter_forward()
	# 
	__pl, __pb, __u = previous_solution.split()
	# __pl.name              = "IF pressure"
	# __pb.name              = "Blood pressure"
	# 
	# Check error
	# if mpi4py.MPI.COMM_WORLD.rank == 0:
	# 	if t>=ti:
	# 		print(f"LDF_v/Baseline: {100*LDF_v_all[-1]/LDF_v_all[index_ti]}, LDF_q/Baseline: {100*LDF_q_all[-1]/LDF_q_all[index_ti]}")
	# 	print(f"Displacement at top points: {displacement_all[-1]:.2e}") 
	# if (t%0.5)==0:
	# 	# Export in the xdmf file 
	# 	xdmf.write_function(displacement_export,t)
	# 	# xdmf.write_function(porosity_n,t)
	# 	xdmf.write_function(porosity_b_n,t)
	# 	xdmf.write_function(__pl,t)
	# 	xdmf.write_function(__pb,t)
	# 	xdmf.write_function(flux_blood,t)
	# 
	# Update the mesh coordinates 
	du_update.interpolate(solution.sub(2))
	mesh.geometry.x[:,:mesh.geometry.dim] += du_update.x.array.reshape((-1, mesh.geometry.dim))
	# 
# xdmf.close()
# 
# Evaluate final time
# end_t = time.time()
# t_hours = int((end_t-begin_t)//3600)
# tmin = int(((end_t-begin_t)%3600)//60)
# tsec = int(((end_t-begin_t)%3600)%60)
# if mpi4py.MPI.COMM_WORLD.rank == 0:
# 	print(f"FEM operated with {num_steps} iterations, in {t_hours} h {tmin} min {tsec} sec")
# #------------------------------------------------------------#
# #                       Post-Processing                      #
# #------------------------------------------------------------#
# # 
# LDF_baseline_v = [100*x/LDF_v_all[index_ti] for x in LDF_v_all]
# LDF_baseline_q = [100*x/LDF_q_all[index_ti] for x in LDF_q_all]
# # Save raw data
# if mpi4py.MPI.COMM_WORLD.rank == 0:
# 	functions.export_to_csv(numpy.transpose([displacement_all, LDF_v_all, LDF_q_all, LDF_baseline_v, LDF_baseline_q, load_all, time_all]), 'ID_63_export.csv',header=['displacement_all', 'LDF_v_all', 'LDF_q_all', 'LDF_baseline_v', 'LDF_baseline_q', 'load_all', 'time_all'])
# 	# 
# 	# Plot for rapid evaluation
# 	import matplotlib.pyplot as plt
# 	# 
# 	plt.rcParams.update({'font.size': 15})
# 	plt.rcParams.update({'legend.loc':'upper right'})
# 	# 
# 	fig1, ax1 = plt.subplots()
# 	ax1.plot(time_all,LDF_baseline_v,linestyle='-',linewidth=2,color='navy', label='LDF_xz')
# 	ax1.plot(time_all,LDF_baseline_q,linestyle=':',linewidth=2,color='cyan', label='LDF_z')
# 	ax1.set_xlabel('time (s)')
# 	ax1.set_ylabel('LDF (-)')
# 	# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
# 	ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
# 	ax1.legend()
# 	fig1.tight_layout()
# 	fig1.savefig('ID_63_Indentation_LDF_2_comp.jpg')
# 	# 
# 	fig1, ax1 = plt.subplots()
# 	ax1.plot(time_all,load_all,linestyle='-',linewidth=2,color='pink', label='LOAD')
# 	ax1.set_xlabel('time (s)')
# 	ax1.set_ylabel('Load (Pa)')
# 	# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
# 	ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
# 	ax1.legend()
# 	fig1.tight_layout()
# 	fig1.savefig('ID_63_Indentation_Load_2_comp.jpg')
# 	# 
# 	fig1, ax1 = plt.subplots()
# 	ax1.plot(time_all[index_ti:],LDF_baseline_v[index_ti:],linestyle='-',linewidth=2,color='navy', label='LDF_xz')
# 	ax1.plot(time_all[index_ti:],LDF_baseline_q[index_ti:],linestyle=':',linewidth=2,color='cyan', label='LDF_z')
# 	ax1.set_xlabel('time (s)')
# 	ax1.set_ylabel('LDF (-)')
# 	# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
# 	ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
# 	ax1.legend()
# 	fig1.tight_layout()
# 	fig1.savefig('ID_63_Indentation_ZOOM_LDF_comp.jpg')
# 	# 
# 	# 
# 	fig1, ax1 = plt.subplots()
# 	ax1.plot(time_all,displacement_all,linestyle='-',linewidth=2,color='salmon')
# 	ax1.set_xlabel('time (s)')
# 	ax1.set_ylabel('u (m)')
# 	# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
# 	ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
# 	fig1.tight_layout()
# 	fig1.savefig('ID_63_Indentation_Displacement_top_2_comp.jpg')
# 	# 
# 	fig1, ax1 = plt.subplots()
# 	ax1.plot(time_all[index_ti:],displacement_all[index_ti:],linestyle='-',linewidth=2,color='salmon')
# 	ax1.set_xlabel('time (s)')
# 	ax1.set_ylabel('u (m)')
# 	# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
# 	ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
# 	fig1.tight_layout()
# 	fig1.savefig('ID_63_Indentation_ZOOM_Displacement_top_2_comp.jpg')
# 	# 
# 	# EoF