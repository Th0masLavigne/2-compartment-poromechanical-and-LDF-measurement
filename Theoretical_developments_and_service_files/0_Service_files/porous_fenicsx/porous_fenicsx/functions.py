# This document holds the functions defined
# by the user for the ease of computation
# such as the solver settings
# 
# Author: Thomas Lavigne
# Date: 18/09/2024
# 
# run first : " python3 -m pip install . " to create the package at ../ location
#------------------------------------------------------------#
#                    Computation tools                       #
#------------------------------------------------------------#
# 
def set_non_linear_solver_parameters(mesh, Problem, atol, rtol, convergence_criterion, max_it, log_newton=True):
	"""
	Configures and returns a non-linear Newton solver for a given mesh and problem.

	Parameters:
	-----------
	mesh : dolfinx.Mesh
	    The computational mesh of the domain. This object provides the communication context for the solver.

	Problem : dolfinx.fem.NonlinearProblem
	    The non-linear problem to be solved. This contains the variational formulation for the system.

	atol : float
	    The absolute tolerance for the solver. Determines when the solution is considered sufficiently converged 
	    regardless of the relative tolerance.

	rtol : float
	    The relative tolerance for the solver. Convergence is determined based on relative changes between 
	    iterations.

	convergence_criterion : str
	    The criterion for convergence, usually set to 'incremental' or 'residual' depending on the desired 
	    stopping criterion.

	max_it : int
	    The maximum number of allowed iterations for the Newton solver before stopping.

	log_newton : bool, optional
	    If set to True (default), Newton convergence logging is enabled to display information during the 
	    solution process. This helps track the solver's progress.

	Returns:
	--------
	solver : dolfinx.nls.petsc.NewtonSolver
	    Configured Newton solver with the specified parameters. The solver is ready to solve the non-linear problem.

	Notes:
	------
	- Uses PETSc Krylov solver for preconditioning, specifically LU factorization with the "mumps" solver type.
	- The Newton solver is initialized with absolute and relative tolerances, as well as a maximum number of iterations.
	- Convergence behavior can be logged using the `log_newton` flag, providing insights into the iteration progress.
	- The Krylov solver (KSP) is configured to use direct factorization (LU) to handle the linear systems at each 
	  Newton step, ensuring robust convergence for difficult problems.

	Exceptions:
	-----------
	If any error occurs during solver configuration, the solver may not initialize correctly, so ensure that the 
	`mesh`, `Problem`, and other arguments are valid.
	"""
	from dolfinx.nls.petsc import NewtonSolver
	import petsc4py
	# 
	if log_newton:
		from dolfinx import log
		log.set_log_level(log.LogLevel.INFO)
	# 
	# set up the non-linear solver
	solver                       = NewtonSolver(mesh.comm, Problem)
	# Absolute tolerance
	solver.atol                  = atol
	# relative tolerance
	solver.rtol                  = rtol
	# Convergence criterion
	solver.convergence_criterion = convergence_criterion
	# Maximum iterations
	solver.max_it                = max_it
	# 
	ksp  = solver.krylov_solver
	opts = petsc4py.PETSc.Options()
	option_prefix = ksp.getOptionsPrefix()
	opts[f"{option_prefix}ksp_type"] = "preonly"
	opts[f"{option_prefix}pc_type"]  = "lu"
	opts[f"{option_prefix}pc_factor_mat_solver_type"] = "mumps"
	ksp.setFromOptions()
	return solver
# 
def export_to_csv(data, filename, header=None):
	"""
	Exports data to a CSV file with an optional header.

	Parameters:
	----------
	data : list of lists
	    The data to be written to the CSV file. Each inner list represents a row.

	filename : str
	    The name of the CSV file to which the data will be written. This file will
	    be created if it doesn't exist, or overwritten if it does.

	header : list, optional
	    A list representing the header row. If provided, this will be written at the
	    top of the CSV file before any data rows.

	Exceptions:
	-----------
	If an error occurs during file writing, an error message will be printed
	with details of the exception.

	Notes:
	------
	- The file is written in write mode ('w'), so any existing file with the same
	  name will be overwritten.
	- The 'newline' parameter is set to an empty string to avoid extra blank lines
	  between rows on Windows systems.
	"""
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
def RMSE(x,xref):
	"""
	Computes the Root Mean Square Error (RMSE) between two numpy arrays.

	Parameters:
	-----------
	x : numpy.ndarray
	    The array to be evaluated. This can represent predicted values, measurements, or any other data.

	xref : numpy.ndarray
	    The reference array of the same size as `x`, representing the true or expected values.

	Returns:
	--------
	rmse_value : float
	    The computed RMSE, which quantifies the average magnitude of the error between `x` and `xref`.

	Notes:
	------
	- RMSE is a commonly used metric in regression and error analysis. It provides a measure of how close predictions 
	  or estimates are to the true values.
	- The formula for RMSE is:
	    RMSE = sqrt((1/N) * Σ (x - xref)²)
	  where N is the number of elements in the arrays, and the sum is taken over all elements.
	- RMSE has the same units as the input arrays and is sensitive to large errors, making it useful for identifying 
	  significant deviations.
	"""
	import numpy as np
	return np.sqrt(np.mean((x-xref)**2))
# 
#------------------------------------------------------------#
#                    User-defined functions                  #
#------------------------------------------------------------#
# 
def mechanical_load(t,ti,t_ramp,magnitude):
	"""
	Computes the temporal evolution of a mechanical load with a smooth ramp-up phase.

	Parameters:
	-----------
	t : float
	    The current time (in seconds) at which to evaluate the load.

	ti : float
	    The time (in seconds) when the load application starts (ramp-up begins).

	t_ramp : float
	    The duration (in seconds) of the ramp-up phase, during which the load smoothly increases from 0 to the full magnitude.

	magnitude : float
	    The maximum magnitude of the load (in Pascals) to be applied once the ramp-up is complete.

	Returns:
	--------
	instantaneous_load : float
	    The instantaneous value of the mechanical load at time `t`. This will be a value between 0 and `magnitude` during the ramp-up phase 
	    and will remain constant at `-magnitude` after the ramp-up phase is complete.

	Notes:
	------
	- The load follows a cosine ramp-up profile, starting from 0 at `t=ti` and smoothly increasing to the specified `magnitude`.
	- Once the ramp-up phase is complete (i.e., after `t = ti + t_ramp`), the load remains constant at the full magnitude.
	- The cosine function ensures a smooth transition, avoiding sudden jumps in the load value.

	Formula:
	--------
	- During the ramp-up phase (ti ≤ t < ti + t_ramp):
	    load = -magnitude * 0.5 * (1 - cos(pi * (t - ti) / t_ramp))

	- After the ramp-up phase (t ≥ ti + t_ramp):
	    load = -magnitude
	"""
	import numpy as np
	if t<ti:
		f1=0
	elif t < (ti+t_ramp):
		tstart=ti
		f1 = 0.5 * (1 - np.cos(np.pi*(t-tstart)/t_ramp))
	else:
		f1 = 1
	return -magnitude*f1
# 
def mechanical_load_LDF(t,t_ramp,magnitude,ti,t_sustained):
	"""
	Computes the temporal evolution of a mechanical load as a function of time.

	Parameters:
	-----------
	t : float
	    The current time (in seconds) at which to evaluate the load.

	t_ramp : float
	    The duration of the ramp-up and ramp-down phases (in seconds), during which the load smoothly increases or decreases.

	magnitude : float
	    The maximum magnitude of the load (in Pascals) applied during the sustained load phases.

	ti : float
	    The start time (in seconds) when the ramp-up begins.

	t_sustained : float
	    The duration (in seconds) of the sustained load phases, during which the load remains constant after each ramp-up.

	Returns:
	--------
	instantaneous_load : float
	    The instantaneous value of the mechanical load at time `t`.

	Notes:
	------
	- The load follows a cyclic pattern with alternating ramp-up, sustained, and ramp-down phases.
	- The function models smooth transitions between no load and the maximum load using a cosine profile during the 
	  ramp-up and ramp-down phases.
	- After several cycles of alternating half-load and full-load, the load transitions to a continuous full-load state 
	  (starting at the 5th cycle), which remains until the end of the process.

	Load Profile:
	-------------
	- **Ramp-up**: The load increases from 0 to half the maximum value or full load.
	- **Sustained phase**: The load is held constant at either half-load or full-load for a specified duration.
	- **Ramp-down**: The load decreases smoothly from its current value back to zero after each sustained phase.
	- This pattern repeats with alternating half and full loads until it eventually reaches full load.
	"""
	import numpy
	if t<ti:
		f1=0
	elif t < (ti+t_ramp):
		tstart=ti
		# the first 0.5 is for half of the magnitude : 20mL or 40 mL
		f1 = 0.5*(0.5 * (1 - numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	elif t < (ti+t_ramp+t_sustained):
		f1 = 0.5*1
	elif t < (ti+2*t_ramp+t_sustained):
		tstart=ti+t_ramp+t_sustained
		f1 = 0.5*(0.5 * (1 + numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	elif t < (ti+2*t_ramp+2*t_sustained):
		f1=0
	elif t < (ti+3*t_ramp+2*t_sustained):
		tstart=ti+2*t_ramp+2*t_sustained
		f1 = 0.5*(0.5 * (1 - numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	elif t < (ti+3*t_ramp+3*t_sustained):
		f1 = 0.5*1
	elif t < (ti+4*t_ramp+3*t_sustained):
		tstart=ti+3*t_ramp+3*t_sustained
		f1 = 0.5*(0.5 * (1 + numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	elif t < (ti+4*t_ramp+4*t_sustained):
		f1=0
	# now it is full load: 40 mL
	elif t < (ti+5*t_ramp+4*t_sustained):
		tstart=ti+4*t_ramp+4*t_sustained
		f1 = (0.5 * (1 - numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	elif t < (ti+5*t_ramp+5*t_sustained):
		f1 = 1
	elif t < (ti+6*t_ramp+5*t_sustained):
		tstart=ti+5*t_ramp+5*t_sustained
		f1 = (0.5 * (1 + numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	elif t < (ti+6*t_ramp+6*t_sustained):
		f1=0
	elif t < (ti+7*t_ramp+6*t_sustained):
		tstart=ti+6*t_ramp+6*t_sustained
		f1 = (0.5 * (1 - numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	elif t < (ti+7*t_ramp+7*t_sustained):
		f1 = 1
	elif t < (ti+8*t_ramp+7*t_sustained):
		tstart=ti+7*t_ramp+7*t_sustained
		f1 = (0.5 * (1 + numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	else:
		f1=0
	return -magnitude*f1
# 
def mechanical_load_LDF_2(t,magnitude,ti):
	"""
	Computes the temporal evolution of a mechanical load as a function of time.

	Parameters:
	-----------
	t : float
	    The current time (in seconds) at which to evaluate the load.

	magnitude : float
	    The maximum magnitude of the load (in Pascals) applied during the sustained load phases.

	ti : float
	    The start time (in seconds) when the ramp-up begins.

	Returns:
	--------
	instantaneous_load : float
	    The instantaneous value of the mechanical load at time `t`.

	"""
	import numpy
	if t<ti:
		f1=0
	elif t < (ti+22):
		tstart=ti
		t_ramp = ti+22-tstart
		# the first 0.5 is for half of the magnitude : 20mL or 40 mL
		f1 = 0.5*(0.5 * (1 - numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	elif t < (ti+1*60+20):
		f1 = 0.5*1
	elif t < (ti+1*60+44):
		tstart=ti+1*60+20
		t_ramp = ti+1*60+44-tstart
		f1 = 0.5*(0.5 * (1 + numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	elif t < (ti+2*60+47):
		f1=0
	elif t < (ti+3*60+7):
		tstart=ti+2*60+47
		t_ramp=ti+3*60+7-tstart
		f1 = 0.5*(0.5 * (1 - numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	elif t < (ti+4*60+9):
		f1 = 0.5*1
	elif t < (ti+4*60+32):
		tstart=ti+4*60+9
		t_ramp = ti+4*60+32-tstart
		f1 = 0.5*(0.5 * (1 + numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	elif t < (ti+5*60+50):
		f1=0
	# now it is full load: 40 mL
	elif t < (ti+6*60+29):
		tstart=ti+5*60+50
		t_ramp=ti+6*60+29-tstart
		f1 = (0.5 * (1 - numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	elif t < (ti+7*60+32):
		f1 = 1
	elif t < (ti+8*60+15):
		tstart=ti+7*60+32
		t_ramp=ti+8*60+15-tstart
		f1 = (0.5 * (1 + numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	elif t < (ti+9*60+17):
		f1=0
	elif t < (ti+9*60+50):
		tstart=ti+9*60+17
		t_ramp=ti+9*60+50-tstart
		f1 = (0.5 * (1 - numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	elif t < (ti+11*60+4):
		f1 = 1
	elif t < (ti+11*60+54):
		tstart=ti+11*60+4
		t_ramp=ti+11*60+54-tstart
		f1 = (0.5 * (1 + numpy.cos(numpy.pi*(t-tstart)/(t_ramp))))
	else:
		f1=0
	return -magnitude*f1
# 
def terzaghi(p0,L,cv,y,t,tstart,kmax):
	"""
	Computes the pore pressure at a given position and time using Terzaghi's one-dimensional consolidation theory.

	Parameters:
	-----------
	p0 : float
	    The initial applied pressure (typically at the surface).
	    
	L : float
	    The length (thickness) of the sample.
	    
	cv : float
	    The coefficient of consolidation, which characterizes the rate of pore pressure dissipation.
	    
	y : float
	    The position within the sample where the pressure is being calculated (0 <= y <= L).
	    
	t : float
	    The time at which the pressure is being evaluated.
	    
	tstart : float
	    The start time of the consolidation process.
	    
	kmax : int
	    The maximum number of terms to include in the summation (affects accuracy).
	    
	Returns:
	--------
	pl : float
	    The pore pressure at position `y` and time `t` based on Terzaghi's consolidation theory.

	Notes:
	------
	- The accuracy of the result depends on the value of `kmax`. A higher `kmax` gives a more accurate solution but requires more computation.
	- This function uses a Fourier series to approximate the solution, which is typical for one-dimensional consolidation problems.
	- The function assumes that the consolidation process begins at `tstart`. For `t < tstart`, the pressure is considered zero.
	"""
	import numpy
	pression=0
	for k in range(1,kmax):
		pression += p0*4/numpy.pi*(-1)**(k-1)/(2*k-1)*numpy.cos((2*k-1)*0.5*numpy.pi*(y/L))*numpy.exp(-(2*k-1)**2*0.25*numpy.pi**2*cv*(t-tstart)/L**2)
	pl = pression
	return pl
# 
#____________________________________________________________#
# 					End of functions
#____________________________________________________________#
# 
# 
if __name__ == "__main__":
    print("Loading of the user-defined functions successfully completed.")
    # EoF