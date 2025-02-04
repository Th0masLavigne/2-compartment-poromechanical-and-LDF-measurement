# Thomas Lavigne
# 18-09-2024
# 
#------------------------------------------------------------#
#                   Libraries 					             #
#------------------------------------------------------------#
# 
# Functions created alongside the theory
try:
	from porous_fenicsx import constitutive_laws, variational_forms, functions
except:
	print("Please run first `python3 -m pip install .` at ./Theoretical_developments_and_service_files/0_Service_files/porous_fenicsx to create the package.")
	exit()
# 
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
filename = "./Mesh.msh"
# 
# Load the mesh and tags in the FEniCSx environment
mesh, cell_tag, facet_tag = dolfinx.io.gmshio.read_from_msh(filename, mpi4py.MPI.COMM_WORLD, 0, gdim=2)
# 
# Identification of the useful elements from the mesh for the integrals
normal  = ufl.FacetNormal(mesh)
# Specify the desired quadrature degree
q_deg = 4
# attribute the cell_tags to the integral element
dx    = ufl.Measure('dx', metadata={"quadrature_degree":q_deg}, subdomain_data=cell_tag, domain=mesh)
# attribute the facet_tags to the surfacic integral element
ds    = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tag)
# 
#------------------------------------------------------------#
#         Time and Mechanical load parameters                #
#------------------------------------------------------------#
# 
# load Magnitude
load_magnitude = 2e2 #[Pa]
# Applied surfacic force for Neumann BCs
T           = dolfinx.fem.Constant(mesh,dolfinx.default_scalar_type(0))
# 
# 
ti = 1
t_ramp = 0.2   
t_sust = 130.
# Time parametrization
t_init    = 0                # Start time
Tf        = t_ramp+t_sust    # End time
dt = 0.5
t=t_init
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
nu        = dolfinx.default_scalar_type(0.4)
# Porous material
# IF fluid viscosity [Pa.s]
mu_l      = dolfinx.default_scalar_type(1e-2) 
# IF Intrinsic permeabilitty [m^2]
k_l       = dolfinx.default_scalar_type(1e-15)
#
# Cell fluid viscosity [Pa.s]
# mu_c      = dolfinx.default_scalar_type(0.5)  
mu_c      = dolfinx.default_scalar_type(10)  
# mu_c      = dolfinx.default_scalar_type(5e-2)  
k_c       = dolfinx.default_scalar_type(1e-15)
# 
# Extra vascular porosity (initial value)
# varepsilon = dolfinx.default_scalar_type(0.7)
varepsilon     = dolfinx.default_scalar_type(0.6)
varepsilon_min = dolfinx.default_scalar_type(0.2)
varepsilon_max = dolfinx.default_scalar_type(0.8)
# 
# Initial saturation in cells
# Sc0 = 0.6625
# Sc0 = 0.005
# Sc0 = 0.998
Sc0 = 1
# 
# To get 20% IF:
varepsilon_l=0.45
Sc0 = 1-(varepsilon_l/varepsilon)
# 
# Coefficient to compute the saturation related to the tortuosity of the medium
a = 600
# 
print("Initial capillary pressure: ",constitutive_laws.Initial_plc(Sc0,a)," Pa")
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
porosity_n        = dolfinx.fem.Function(P1_space)
porosity_n.name = "porosity"
# Create the cell pressure function
pc        = dolfinx.fem.Function(P1_space)
pc.name = "Cell pressure"
# Create the solid pressure function
ps        = dolfinx.fem.Function(P1_space)
ps.name   = "Solid pressure"
# Create the cell pressure function
satur        = dolfinx.fem.Function(P1_space)
satur.name = "Cell Saturation"
# 
# Define the solution at the current and previous time steps
solution          = dolfinx.fem.Function(MS)
previous_solution = dolfinx.fem.Function(MS)
# 
# 
# Mapping in the Mixed Space: FunctionSpace, Mapping_in_the_mixed_space = MS.sub(xx).collapse()
Pln_, Pln_to_MS   = MS.sub(0).collapse()
Plcn_, Plcn_to_MS = MS.sub(1).collapse()
Un_, Un_to_MS     = MS.sub(2).collapse()
# 
# Create the test functions
qc,ql,w = ufl.TestFunctions(MS)
#------------------------------------------------------------#
#                          Expressions                       #
#------------------------------------------------------------#
#
# Expressions to be computed before the update of previous solution
# Evaluation of pc
pc_expr         = dolfinx.fem.Expression( ( (previous_solution.sub(0)+solution.sub(0)) - (previous_solution.sub(1)+solution.sub(1)) ), P1_space.element.interpolation_points())
# 
ps_expr           = dolfinx.fem.Expression( ( (previous_solution.sub(0)+solution.sub(0)) - constitutive_laws.saturation_cell(solution.sub(1)+previous_solution.sub(1),a)*(previous_solution.sub(1)+solution.sub(1)) ), P1_space.element.interpolation_points())
# Evaluation of the porosity
porosity_n_expr = dolfinx.fem.Expression(constitutive_laws.porosity_1(solution.sub(2),porosity_n,varepsilon_min,varepsilon_max), P1_space.element.interpolation_points())
# Evaluation of the porosity
saturation_expr = dolfinx.fem.Expression(constitutive_laws.saturation_cell(solution.sub(1)+previous_solution.sub(1),a), P1_space.element.interpolation_points())
# 
# Computation of the pressure at bottom points
S_bottom = surface_indenter = mesh.comm.allreduce(dolfinx.fem.assemble_scalar(dolfinx.fem.form(1*ds(1)) ), op=mpi4py.MPI.SUM)
Press_bottom_expr_IF        = dolfinx.fem.form(1/S_bottom*(previous_solution.sub(0)+solution.sub(0))*ds(1))
Press_bottom_expr_s         = dolfinx.fem.form(1/S_bottom*( previous_solution.sub(0)+solution.sub(0)- constitutive_laws.saturation_cell(previous_solution.sub(1)+solution.sub(1),a)*(previous_solution.sub(1)+solution.sub(1)) )*ds(1))
Press_bottom_expr_c         = dolfinx.fem.form(1/S_bottom*( (previous_solution.sub(0)+solution.sub(0)) - (previous_solution.sub(1)+solution.sub(1)) )*ds(1))
# Computation of the displacement at top points
S_top = surface_indenter = mesh.comm.allreduce(dolfinx.fem.assemble_scalar(dolfinx.fem.form(1*ds(4)) ), op=mpi4py.MPI.SUM)
Disp_top_expr            = dolfinx.fem.form(1/S_top*(ufl.dot(previous_solution.sub(2)+solution.sub(2),normal))*ds(4))
# 
#------------------------------------------------------------#
#                   Create the export file                   #
#------------------------------------------------------------#
#
# Open file for export
xdmf = dolfinx.io.XDMFFile(mesh.comm, "_Result.xdmf", "w")
xdmf.write_mesh(mesh)
# 
#------------------------------------------------------------#
#                   Initial Conditions                       #
#------------------------------------------------------------#
# 
# Porosity
porosity_n.x.array[:]=numpy.full_like(porosity_n.x.array[:], varepsilon, dtype=dolfinx.default_scalar_type)
# Capillary pressure to get an initial value of saturation
previous_solution.x.array[Plcn_to_MS]=numpy.full_like(previous_solution.x.array[Plcn_to_MS], constitutive_laws.Initial_plc(Sc0,a), dtype=dolfinx.default_scalar_type)
# pl compensate so ps=0
previous_solution.x.array[Pln_to_MS] = numpy.full_like(previous_solution.x.array[Pln_to_MS], Sc0 * constitutive_laws.Initial_plc(Sc0,a), dtype=dolfinx.default_scalar_type)
# initial cell saturation
satur.x.array[:]=numpy.full_like(satur.x.array[:],Sc0,dtype=dolfinx.default_scalar_type)
#------------------------------------------------------------#
#              Dirichlet Boundary Conditions                 #
#------------------------------------------------------------#
# 
bcs    = []
fdim   = mesh.topology.dim - 1
# 
# bottom_marker, left_marker, right_marker, top_marker = 1,2,3,4
# 
facets = facet_tag.find(1)
# duy=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(2).sub(1), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(2).sub(1)))
# 
facets = facet_tag.find(2)
# dux=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(2).sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(2).sub(0)))
# 
facets = facet_tag.find(3)
# dux=0
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(2).sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), dofs, MS.sub(2).sub(0)))
# 
facets = facet_tag.find(4)
# dpl=0
dpl_imp = dolfinx.fem.Constant(mesh,dolfinx.default_scalar_type(0))
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(0), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dpl_imp, dofs, MS.sub(0)))
# dplc=0 
dplc_imp = dolfinx.fem.Constant(mesh,dolfinx.default_scalar_type(0))
dofs   = dolfinx.fem.locate_dofs_topological(MS.sub(1), fdim, facets)
bcs.append(dolfinx.fem.dirichletbc(dplc_imp, dofs, MS.sub(1)))
# 
#------------------------------------------------------------#
#                     Variationnal form                      #
#------------------------------------------------------------#
#
# 
F = variational_forms.variational_form_1_comp_bi(previous_solution, solution, qc, ql, w, porosity_n, varepsilon_min, varepsilon_max, dt, E, nu, k_l, mu_l, k_c, mu_c, a, dx, ds)
# 
# Add Neuman BCs
F+= - T*ufl.inner(w,normal)*ds(4)
# 
#------------------------------------------------------------#
#                           Solver                           #
#------------------------------------------------------------#
# 
# Non linear problem definition
dsolution = ufl.TrialFunction(MS)
J         = ufl.derivative(F, solution, dsolution)
Problem   = NonlinearProblem(F, solution, bcs = bcs, J = J)
solver    = functions.set_non_linear_solver_parameters(mesh, Problem, 1e-11, 1e-10, "incremental", 10, log_newton=False)
# 
#------------------------------------------------------------#
#                         Computation                        #
#------------------------------------------------------------#
# Check initial sol
# Split the total solution for export
__pl, __plc, __u = previous_solution.split()
__pl.name="IF pressure"
__plc.name = "Capillary Pressure"
# Export in the xdmf dile 
xdmf.write_function(displacement_export,t)
xdmf.write_function(porosity_n,t)
xdmf.write_function(__pl,t)
xdmf.write_function(__plc,t)
xdmf.write_function(pc,t)
xdmf.write_function(ps,t)
xdmf.write_function(satur,t)
# 
displacement_all, pressure_solid_all, pressure_IF_all, pressure_cell_all, load_all, time_all = [], [], [], [], [], []
# 
for n in range(num_steps):
	# update time
	t+=dt
	# t = round(t,2)
	time_all.append(t)
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
	# Update the cell Pressure
	pc.interpolate(pc_expr)
	pc.x.scatter_forward()
	# Update the solid Pressure
	ps.interpolate(ps_expr)
	ps.x.scatter_forward()
	# Update saturation
	satur.interpolate(saturation_expr)
	satur.x.scatter_forward()
	# 
	# Evaluate the bottom pressure
	PBot = dolfinx.fem.assemble_scalar(Press_bottom_expr_IF)
	pressure_IF_all.append(mesh.comm.allreduce(PBot, op=mpi4py.MPI.SUM))
	PBot_s = dolfinx.fem.assemble_scalar(Press_bottom_expr_s)
	pressure_solid_all.append(mesh.comm.allreduce(PBot_s, op=mpi4py.MPI.SUM))
	PBot_c = dolfinx.fem.assemble_scalar(Press_bottom_expr_c)
	pressure_cell_all.append(mesh.comm.allreduce(PBot_c, op=mpi4py.MPI.SUM))
	# Evaluate the top displacement
	DTop = dolfinx.fem.assemble_scalar(Disp_top_expr)
	displacement_all.append(mesh.comm.allreduce(DTop, op=mpi4py.MPI.SUM))
	# Update previous solution
	previous_solution.x.array[:] += solution.x.array[:]
	previous_solution.x.scatter_forward()
  	# 
  	# Displacement to P1 for export
	displacement_export.interpolate(previous_solution.sub(2))
	displacement_export.x.scatter_forward()
	# 
	# Split the total solution for export
	__pl, __plc, __u = previous_solution.split()
	__pl.name="IF pressure"
	__plc.name = "Capillary Pressure"
	# 
	# Check error
	if mpi4py.MPI.COMM_WORLD.rank == 0:
		print(f"Pressure s at bottom points: {pressure_solid_all[-1]:.2e}") 
		print(f"Displacement at top points: {displacement_all[-1]:.2e}") 
	tcheck = round(t,2)
	if (tcheck%1)==0:
		# Export in the xdmf dile 
		xdmf.write_function(displacement_export,t)
		xdmf.write_function(porosity_n,t)
		xdmf.write_function(__pl,t)
		xdmf.write_function(__plc,t)
		xdmf.write_function(pc,t)
		xdmf.write_function(ps,t)
		xdmf.write_function(satur,t)
	# 
	# Update the mesh coordinates 
	du_update.interpolate(solution.sub(2))
	mesh.geometry.x[:,:mesh.geometry.dim] += du_update.x.array.reshape((-1, mesh.geometry.dim))
	# 
	# Compute the new height for the beginning of the Terzaghi comparison
	tcheck = round(t,2)
	if tcheck==ti:
		Height_corrected = 1e-4 + displacement_all[n]
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
lambda_m, mu = constitutive_laws.Lame_coefficients(E,nu)
cvl = k_l/mu_l*(lambda_m+2*mu)
cvc = k_c/mu_c*(lambda_m+2*mu)
y=0
pressurec = numpy.zeros(num_steps)
kmax=1e3
for i in range(num_steps):
	if time_all[i]<ti+t_ramp:
		pressurec[i]=200
	else:
		pressurec[i] = functions.terzaghi(load_magnitude,Height_corrected,cvc,y,time_all[i],ti+t_ramp,int(kmax))

pressurel = numpy.zeros(num_steps)
kmax=1e3
for i in range(num_steps):
	if time_all[i]<ti+t_ramp:
		pressurel[i]=200
		index = i
	else:
		pressurel[i] = functions.terzaghi(load_magnitude,Height_corrected,cvl,y,time_all[i],ti+t_ramp,int(kmax))
# 
# Evaluate the error
# RMSE_p = functions.RMSE(pressure_IF_all[index+1:],pressurel[index+1:])
RMSE_p = functions.RMSE(pressure_solid_all[index+1:],pressurel[index+1:])
# 
if mpi4py.MPI.COMM_WORLD.rank == 0:
	print(f"RMSE of the solid pressure with respect to Terzaghi's fluid solution (i.e comparing our new ps to ps=pl previously) at bottom points {RMSE_p} Pa, RMSE/max(p_theo) = {100*RMSE_p/max(pressurel)} %")
	print("The more a augments, the more the error augments.")
# 
# Save raw data
functions.export_to_csv(numpy.transpose([displacement_all, pressure_solid_all, pressure_IF_all, pressure_cell_all, load_all, time_all ,pressurel]), '_Results.csv',header=['displacement_all', 'pressure_solid_all', 'pressure_IF_all', 'pressure_cell_all', 'load_all', 'time_all', 'pressurel'])
# 
# Plot for rapid evaluation
import matplotlib.pyplot as plt
# 
plt.rcParams.update({'font.size': 15})
plt.rcParams.update({'legend.loc':'upper right'})
# 
fig1, ax1 = plt.subplots()
ax1.plot(time_all,pressure_IF_all,linestyle='-',linewidth=2,color='navy', label='IF')
ax1.plot(time_all,pressure_solid_all,linestyle='-.',linewidth=2,color='cyan', label='S')
ax1.plot(time_all,pressure_cell_all,linestyle='-.',linewidth=2,color='gray', label='C')
ax1.plot(time_all,load_all,linestyle='-',linewidth=2,color='pink', label='LOAD')
# ax1.plot(time_all,pressurec,linestyle=':',linewidth=2,color='red', label='Terz_C')
ax1.plot(time_all,pressurel,linestyle=':',linewidth=2,color='orange', label='TERZ_IF')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('Pressure (Pa)')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
ax1.legend()
fig1.tight_layout()
fig1.savefig('_Pressure_bottom.jpg')
# 
# 
fig1, ax1 = plt.subplots()
ax1.plot(time_all,displacement_all,linestyle='-',linewidth=2,color='salmon')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('u (m)')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
fig1.tight_layout()
fig1.savefig('_Displacement_top.jpg')
# 
exit()
# EoF