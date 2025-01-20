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
filename = "./Mouse_mesh.msh"
# filename = "./Mouse_local_refine_mesh.msh"
# 
# Load the mesh and tags in the FEniCSx environment
mesh, cell_tag, facet_tag = dolfinx.io.gmshio.read_from_msh(filename, mpi4py.MPI.COMM_WORLD, 0, gdim=3)
# 
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
# load Magnitude
load_magnitude = 2e3 #[Pa]
# Applied surfacic force for Neumann BCs
T              = dolfinx.fem.Constant(mesh,dolfinx.default_scalar_type(0))
# 
# initial time for imposing the blood pressure gradient gradually
index_ti = 0
ti       = 4
t_ramp   = 5 
t_sust   = 25.
# Time parametrization
t_init   = 0                # Start time
Tf       = t_ramp+t_sust    # End time
# Must ensure t%dt==0
dt       = 0.5
t        = t_init
num_steps_pb_imp = int(ti/dt)
num_steps = int((Tf-t_init+ti)/dt)
# 
#------------------------------------------------------------#
#                   Material Definition                      #
#------------------------------------------------------------#
# 
# Solid scaffold
# Young Moduli [Pa]
E         = dolfinx.default_scalar_type(5e3) 
# Possion's ratios [-]
nu        = dolfinx.default_scalar_type(0.2)
# Porous material
#
# Extra vascular porosity (initial value)
varepsilon     = dolfinx.default_scalar_type(0.6)
varepsilon_min = dolfinx.default_scalar_type(0.2)
varepsilon_max = dolfinx.default_scalar_type(0.8)
# varepsilon = dolfinx.default_scalar_type(0.5)
# 
# IF fluid viscosity [Pa.s]
mu_l      = dolfinx.default_scalar_type(1) 
# IF Intrinsic permeabilitty [m^2]
k_l       = dolfinx.default_scalar_type(1e-14)
# 
# Cell fluid viscosity [Pa.s]
# mu_c      = dolfinx.default_scalar_type(0.5)  
# mu_c      = dolfinx.default_scalar_type(20)  
mu_c      = dolfinx.default_scalar_type(5)  
# mu_c      = dolfinx.default_scalar_type(100)  
# Cell Intrinsic permeabilitty [m^2]
k_c       = dolfinx.default_scalar_type(1e-14)
# 
# To get initial% interstitium fluid (other than cell):
varepsilon_l = 0.45
Sc0          = 1-(varepsilon_l/varepsilon)
# 
# Coefficient to compute the saturation related to the tortuosity of the medium
# a = dolfinx.default_scalar_type(600)
a = dolfinx.default_scalar_type(300)
# 
print("Initial capillary pressure: ",constitutive_laws.Initial_plc(Sc0,a)," Pa", " Initial Sc", Sc0, " "," a=",a)
# 
# Blood fluid viscosity [Pa.s]
mu_b      = dolfinx.default_scalar_type(4e-3)  
# Vascular Intrinsic permeabilitty [m^2]
k_b       = dolfinx.default_scalar_type(4e-12)
# To define a preferential direction for function matrix_vascular_permeability
ratio     = [1, 1, 1]
# exponent to make evolution of kb according to poro
alpha     =  dolfinx.default_scalar_type(3)
# Vessel Compressibility [Pa]
K         = dolfinx.default_scalar_type(numpy.pi/2*1000)
# 
# BLOOD FLOW [m3/sec * 1/m2]
BFLOW  = dolfinx.default_scalar_type(4.e-6)
PB_IMP = dolfinx.fem.Constant(mesh,  dolfinx.default_scalar_type(BFLOW * (mu_b/k_b) * 0.006))
dPB_IMP = dolfinx.fem.Constant( mesh,  dolfinx.default_scalar_type( PB_IMP.value/num_steps_pb_imp ) )
# PB_IMP = dolfinx.fem.Constant(mesh,  dolfinx.default_scalar_type(24))
print("pression sanguine imposee= ",PB_IMP.value)
# 
# 
# Vascular porosity (initial value)
initial_varepsilonb = 0.04
varepsilonb = dolfinx.default_scalar_type( initial_varepsilonb )
# 
# 
#------------------------------------------------------------#
#                   Function Spaces                          #
#------------------------------------------------------------#
# 
# Mixed Space (R,R,R2) -> (pc,plc,u)
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
MS            = dolfinx.fem.functionspace(mesh=mesh, element=basix.ufl.mixed_element([P1,P1,P1,P2]))
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
porosity_n        = dolfinx.fem.Function(P1_space)
porosity_n.name   = "porosity"
# Create the previous solution of the vascular porosity
porosity_b_n        = dolfinx.fem.Function(P1_space)
porosity_b_n.name   = "vascular porosity"
# Create the cell pressure function
pc        = dolfinx.fem.Function(P1_space)
pc.name   = "Cell pressure"
# Create the solid pressure function
ps        = dolfinx.fem.Function(P1_space)
ps.name   = "Solid pressure"
# Create the cell saturation function
satur        = dolfinx.fem.Function(P1_space)
satur.name   = "Cell Saturation"
# 
# Define the solution at the current and previous time steps
solution          = dolfinx.fem.Function(MS)
previous_solution = dolfinx.fem.Function(MS)
# 
# 
# Mapping in the Mixed Space: FunctionSpace, Mapping_in_the_mixed_space = MS.sub(xx).collapse()
Pln_, Pln_to_MS   = MS.sub(0).collapse()
Plcn_, Plcn_to_MS = MS.sub(1).collapse()
Pbn_, Pbn_to_MS   = MS.sub(2).collapse()
Un_, Un_to_MS     = MS.sub(3).collapse()
# 
# Create the test functions
qc,ql,qb,w = ufl.TestFunctions(MS)
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
# 
k_b_current          = constitutive_laws.matrix_vascular_permeability_bi(k_b, ratio, alpha, varepsilonb, previous_solution.sub(0)+solution.sub(0), previous_solution.sub(1)+solution.sub(1), previous_solution.sub(2)+solution.sub(2), a, K)
# varepsilon_current   = constitutive_laws.porosity_2_pressure(solution.sub(3),porosity_n,k_b_current,mu_b,previous_solution.sub(2)+solution.sub(2),dt,varepsilon_min,varepsilon_max)
varepsilon_current   = constitutive_laws.porosity_2(solution.sub(3),porosity_n,constitutive_laws.vascular_porosity_atan_bi(varepsilonb,previous_solution.sub(0)+solution.sub(0),previous_solution.sub(1)+solution.sub(1),previous_solution.sub(2)+solution.sub(2),a,K),constitutive_laws.vascular_porosity_atan_bi(varepsilonb,previous_solution.sub(0),previous_solution.sub(1),previous_solution.sub(2),a,K),varepsilon_min,varepsilon_max)
varepsilon_b_current = constitutive_laws.vascular_porosity_atan_bi(varepsilonb,previous_solution.sub(0)+solution.sub(0),previous_solution.sub(1)+solution.sub(1),previous_solution.sub(2)+solution.sub(2),a,K)
#
darcy_blood      = -1 * k_b_current /mu_b*ufl.grad(previous_solution.sub(2)+solution.sub(2))
# Evaluation of volume and surfaces to integrate the solutions
Volume_LDF       = mesh.comm.allreduce(dolfinx.fem.assemble_scalar(dolfinx.fem.form(1*dx(1)) ), op=mpi4py.MPI.SUM)
surface_indenter = mesh.comm.allreduce(dolfinx.fem.assemble_scalar(dolfinx.fem.form(1*ds(1)) ), op=mpi4py.MPI.SUM)
# 
# Expressions to be computed before the update of previous solution
# Evaluation of pc
pc_expr           = dolfinx.fem.Expression( ( (previous_solution.sub(0)+solution.sub(0)) - (previous_solution.sub(1)+solution.sub(1)) ), P1_space.element.interpolation_points())
# Evaluation of ps
ps_expr           = dolfinx.fem.Expression( ( (previous_solution.sub(0)+solution.sub(0)) - constitutive_laws.saturation_cell(solution.sub(1)+previous_solution.sub(1),a)*(previous_solution.sub(1)+solution.sub(1)) ), P1_space.element.interpolation_points())
# Evaluation of the porosity 
porosity_n_expr   = dolfinx.fem.Expression(varepsilon_current, P1_space.element.interpolation_points())
# Evaluation of the vascular porosity
porosity_b_n_expr = dolfinx.fem.Expression(varepsilon_b_current, P1_space.element.interpolation_points())
# Evaluation of the blood flux qb = -kb/mu * (vb-vs)
flux_blood_expr   = dolfinx.fem.Expression(darcy_blood , P1v_space.element.interpolation_points())
# vbs = qb/epsb
# speed_blood_expr = dolfinx.fem.Expression(1/varepsilon_b_current *darcy_blood , P1v_space.element.interpolation_points())
saturation_expr   = dolfinx.fem.Expression(constitutive_laws.saturation_cell(solution.sub(1)+previous_solution.sub(1),a), P1_space.element.interpolation_points())
# 
# Computation of the LDF
# LDF_expr_q                    = dolfinx.fem.form(1/Volume_LDF*(ufl.sqrt(ufl.inner(darcy_blood,darcy_blood))*dx(1)))
# Only in plane
# LDF_expr_v                  = dolfinx.fem.form(1/Volume_LDF*1/varepsilon_b_current*(ufl.sqrt(ufl.dot(darcy_blood,Nx)**2+ufl.dot(darcy_blood,Ny)**2))*dx(1))
# LDF_expr_q                  = dolfinx.fem.form(1/Volume_LDF*(ufl.sqrt(ufl.dot(darcy_blood,Nx)**2+ufl.dot(darcy_blood,Ny)**2))*dx(1))
# Only in x (in reality q not v just to compare xy vs x)
LDF_expr_v                    = dolfinx.fem.form(1/Volume_LDF*(ufl.dot(darcy_blood/varepsilon_b_current,Nx))*dx(1))
LDF_expr_q                    = dolfinx.fem.form(1/Volume_LDF*(ufl.dot(darcy_blood,Nx))*dx(1))
# Computation of the displacement at top points
Disp_top_expr                 = dolfinx.fem.form(1/surface_indenter*(ufl.dot(previous_solution.sub(3)+solution.sub(3),normal))*ds(1))
# 
#------------------------------------------------------------#
#                   Create the export file                   #
#------------------------------------------------------------#
#
# Open file for export
xdmf = dolfinx.io.XDMFFile(mesh.comm, "Indentation_Result_2_comp.xdmf", "w")
xdmf.write_mesh(mesh)
# 
#------------------------------------------------------------#
#                   Initial Conditions                       #
#------------------------------------------------------------#
# 
# Porosity
porosity_n.x.array[:]   = numpy.full_like(porosity_n.x.array[:], varepsilon, dtype=dolfinx.default_scalar_type)
# Vascular Porosity
porosity_b_n.x.array[:] = numpy.full_like(porosity_b_n.x.array[:], varepsilonb, dtype=dolfinx.default_scalar_type)
# Capillary pressure to get an initial value of saturation
previous_solution.x.array[Plcn_to_MS] = numpy.full_like(previous_solution.x.array[Plcn_to_MS], constitutive_laws.Initial_plc(Sc0,a), dtype=dolfinx.default_scalar_type)
# pl to compensate so ps=0
previous_solution.x.array[Pln_to_MS] = numpy.full_like(previous_solution.x.array[Pln_to_MS], Sc0 * constitutive_laws.Initial_plc(Sc0,a), dtype=dolfinx.default_scalar_type)
# saturation in space
satur.x.array[:]=numpy.full_like(satur.x.array[:],Sc0,dtype=dolfinx.default_scalar_type)
# 
#------------------------------------------------------------#
#              Dirichlet Boundary Conditions                 #
#------------------------------------------------------------#
# 
bcs    = []
fdim   = mesh.topology.dim - 1
# 
# bottom_marker, front_marker, back_marker, left_marker, right_marker, top_marker, indenter_marker = 3, 6, 7, 4, 5, 2, 1
# 
# bottom
facets = facet_tag.find(3)
# duz=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(3).sub(2), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(3).sub(2)))
# front
facets = facet_tag.find(6)
# duy=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(3).sub(1), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(3).sub(1)))
# back
facets = facet_tag.find(7)
# duy=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(3).sub(1), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(3).sub(1)))
# dpl=0
dpl_imp_back = dolfinx.fem.Constant(mesh,dolfinx.default_scalar_type(0))
dofs         = dolfinx.fem.locate_dofs_topological(MS.sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dpl_imp_back, dofs, MS.sub(0)))
# dplc=0 
dplc_imp_back = dolfinx.fem.Constant(mesh,dolfinx.default_scalar_type(0))
dofs          = dolfinx.fem.locate_dofs_topological(MS.sub(1), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dplc_imp_back, dofs, MS.sub(1)))
# left
facets = facet_tag.find(4)
# dux=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(3).sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(3).sub(0)))
# dpl=0
dpl_imp_left = dolfinx.fem.Constant(mesh,dolfinx.default_scalar_type(0))
dofs         = dolfinx.fem.locate_dofs_topological(MS.sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dpl_imp_left, dofs, MS.sub(0)))
# dplc=0 
dplc_imp_left = dolfinx.fem.Constant(mesh,dolfinx.default_scalar_type(0))
dofs          = dolfinx.fem.locate_dofs_topological(MS.sub(1), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dplc_imp_left, dofs, MS.sub(1)))
# dpb= 0
dpb_imp_left = dolfinx.fem.Constant(mesh,dolfinx.default_scalar_type(0))
dofs         = dolfinx.fem.locate_dofs_topological(MS.sub(2), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dpb_imp_left, dofs, MS.sub(2)))
# right
facets = facet_tag.find(5)
# dux=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(3).sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(3).sub(0)))
# dpl=0
dpl_imp_right = dolfinx.fem.Constant(mesh,dolfinx.default_scalar_type(0))
dofs          = dolfinx.fem.locate_dofs_topological(MS.sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dpl_imp_right, dofs, MS.sub(0)))
# dplc=0 
dplc_imp_right = dolfinx.fem.Constant(mesh,dolfinx.default_scalar_type(0))
dofs           = dolfinx.fem.locate_dofs_topological(MS.sub(1), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dplc_imp_right, dofs, MS.sub(1)))
# dpb=0 after gradual increase of the pressure to get the good BCs
dpb_imp_right = dPB_IMP
dofs          = dolfinx.fem.locate_dofs_topological(MS.sub(2), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dpb_imp_right, dofs, MS.sub(2)))
# 
#------------------------------------------------------------#
#                     Variationnal form                      #
#------------------------------------------------------------#
#
# 
F = variational_forms.variational_form_2_comp_bi(previous_solution, solution, qc, ql, qb, w, porosity_n, varepsilon_min, varepsilon_max, dt, E, nu, k_l, mu_l, k_c, mu_c, k_b, mu_b, a, K, varepsilonb, alpha, ratio, dx, ds)
#
# Add Neuman BCs on indenter surface
F+= - T*ufl.inner(w,normal)*ds(1)
# 
#------------------------------------------------------------#
#                           Solver                           #
#------------------------------------------------------------#
# 
# Non linear problem definition
dsolution = ufl.TrialFunction(MS)
J         = ufl.derivative(F, solution, dsolution)
Problem   = NonlinearProblem(F, solution, bcs = bcs, J = J)
solver    = functions.set_non_linear_solver_parameters(mesh, Problem, 5e-10, 1e-11, "incremental", 10, log_newton=True)
# 
#------------------------------------------------------------#
#                         Computation                        #
#------------------------------------------------------------#
# Check initial sol
# Split the total solution for export
__pl, __plc, __pb, __u = previous_solution.split()
__pl.name              = "IF pressure"
__plc.name             = "Capillary Pressure"
__pb.name              = "Blood pressure"
# Export in the xdmf file 
xdmf.write_function(displacement_export,t)
xdmf.write_function(porosity_n,t)
xdmf.write_function(porosity_b_n,t)
xdmf.write_function(__pl,t)
xdmf.write_function(__plc,t)
xdmf.write_function(__pb,t)
xdmf.write_function(pc,t)
xdmf.write_function(ps,t)
xdmf.write_function(flux_blood,t)
xdmf.write_function(satur,t)
# 
displacement_all, LDF_v_all, LDF_q_all, load_all, time_all = [], [], [], [], []
# 
for n in range(num_steps):
	# Update BCs
	if n == num_steps_pb_imp:
		dpb_imp_right.value = 0
	# if n==1:
	# 	# dpb_imp_left.value = 0
	# 	dpb_imp_right.value = 0
	# update time
	t+=dt
	t = round(t,2)
	time_all.append(t)
	if t == ti:
		index_ti = n
  	# update the load value and export it in the list
	T.value = functions.mechanical_load(t,ti,t_ramp,load_magnitude)
	load_all.append(-functions.mechanical_load(t,ti,t_ramp,load_magnitude))
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
	porosity_n.interpolate(porosity_n_expr)
	porosity_n.x.scatter_forward()
	# Update vascular porosity
	porosity_b_n.interpolate(porosity_b_n_expr)
	porosity_b_n.x.scatter_forward()
	# Update the cell Pressure
	pc.interpolate(pc_expr)
	pc.x.scatter_forward()
	# Update the solid Pressure
	ps.interpolate(ps_expr)
	ps.x.scatter_forward()
	# Update the blood flux
	flux_blood.interpolate(flux_blood_expr)
	flux_blood.x.scatter_forward
	# Update saturation
	satur.interpolate(saturation_expr)
	satur.x.scatter_forward()
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
	displacement_export.interpolate(previous_solution.sub(3))
	displacement_export.x.scatter_forward()
	# 
	__pl, __plc, __pb, __u = previous_solution.split()
	__pl.name              = "IF pressure"
	__plc.name             = "Capillary Pressure"
	__pb.name              = "Blood pressure"
	# 
	# Check error
	if mpi4py.MPI.COMM_WORLD.rank == 0:
		if t>=ti:
			print(f"LDF_v/Baseline: {100*LDF_v_all[-1]/LDF_v_all[index_ti]}, LDF_q/Baseline: {100*LDF_q_all[-1]/LDF_q_all[index_ti]}")
		print(f"Displacement at top points: {displacement_all[-1]:.2e}") 
	if (t%0.5)==0:
		# Export in the xdmf file 
		xdmf.write_function(displacement_export,t)
		xdmf.write_function(porosity_n,t)
		xdmf.write_function(porosity_b_n,t)
		xdmf.write_function(__pl,t)
		xdmf.write_function(__plc,t)
		xdmf.write_function(__pb,t)
		xdmf.write_function(pc,t)	
		xdmf.write_function(ps,t)	
		xdmf.write_function(flux_blood,t)
		xdmf.write_function(satur,t)
	# 
	# Update the mesh coordinates 
	du_update.interpolate(solution.sub(3))
	mesh.geometry.x[:,:mesh.geometry.dim] += du_update.x.array.reshape((-1, mesh.geometry.dim))
	# 
xdmf.close()
# 
# Evaluate final time
end_t = time.time()
t_hours = int((end_t-begin_t)//3600)
tmin = int(((end_t-begin_t)%3600)//60)
tsec = int(((end_t-begin_t)%3600)%60)
if mpi4py.MPI.COMM_WORLD.rank == 0:
	print(f"FEM operated with {num_steps} iterations, in {t_hours} h {tmin} min {tsec} sec")
#------------------------------------------------------------#
#                       Post-Processing                      #
#------------------------------------------------------------#
# 
LDF_baseline_v = [100*x/LDF_v_all[index_ti] for x in LDF_v_all]
LDF_baseline_q = [100*x/LDF_q_all[index_ti] for x in LDF_q_all]
# Save raw data
functions.export_to_csv(numpy.transpose([displacement_all, LDF_v_all, LDF_q_all, LDF_baseline_v, LDF_baseline_q, load_all, time_all]), 'Indentation_Results_2_comp.csv',header=['displacement_all', 'LDF_v_all', 'LDF_q_all', 'LDF_baseline_v', 'LDF_baseline_q', 'load_all', 'time_all'])
# 
# Plot for rapid evaluation
import matplotlib.pyplot as plt
# 
plt.rcParams.update({'font.size': 15})
plt.rcParams.update({'legend.loc':'upper right'})
# 
fig1, ax1 = plt.subplots()
ax1.plot(time_all,LDF_baseline_v,linestyle='-',linewidth=2,color='navy', label='LDF_v')
ax1.plot(time_all,LDF_baseline_q,linestyle=':',linewidth=2,color='cyan', label='LDF_q')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('LDF (-)')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
ax1.legend()
fig1.tight_layout()
fig1.savefig('Indentation_LDF_2_comp.jpg')
# 
fig1, ax1 = plt.subplots()
ax1.plot(time_all,load_all,linestyle='-',linewidth=2,color='pink', label='LOAD')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('Load (Pa)')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
ax1.legend()
fig1.tight_layout()
fig1.savefig('Indentation_Load_2_comp.jpg')
# 
fig1, ax1 = plt.subplots()
ax1.plot(time_all[index_ti:],LDF_baseline_v[index_ti:],linestyle='-',linewidth=2,color='navy', label='LDF_v')
ax1.plot(time_all[index_ti:],LDF_baseline_q[index_ti:],linestyle=':',linewidth=2,color='cyan', label='LDF_q')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('LDF (-)')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
ax1.legend()
fig1.tight_layout()
fig1.savefig('Indentation_ZOOM_LDF_comp.jpg')
# 
# 
fig1, ax1 = plt.subplots()
ax1.plot(time_all,displacement_all,linestyle='-',linewidth=2,color='salmon')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('u (m)')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
fig1.tight_layout()
fig1.savefig('Indentation_Displacement_top_2_comp.jpg')
# 
fig1, ax1 = plt.subplots()
ax1.plot(time_all[index_ti:],displacement_all[index_ti:],linestyle='-',linewidth=2,color='salmon')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('u (m)')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
fig1.tight_layout()
fig1.savefig('Indentation_ZOOM_Displacement_top_2_comp.jpg')
# 
exit()
# EoF