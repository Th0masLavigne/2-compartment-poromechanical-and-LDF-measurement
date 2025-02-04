# This document holds the functions used as part of 
# the poro-mechanical modelling of living tissues 
# according to the theory in ./Theoretical_Developments.pdf
# Theta Wilson method
# Author: Thomas Lavigne
# Date: 18/09/2024
# 
# run first : " python3 -m pip install . " to create the package at ../ location
# 
from porous_fenicsx import constitutive_laws
# 
#------------------------------------------------------------#
#                    Variational forms                       #
#------------------------------------------------------------#
# 
# 
#...........................................................#
# 					 Monophasic
#...........................................................#
# 
#____________________________________________________________#
# 				Single Compartment Monophasic
#____________________________________________________________#
# 
# 
def variational_form_1_mono_total(previous_solution, solution, ql, w, Young, Poisson, k_l, mu_l, dx, ds, dt, material_law):
	"""
	This function provides the total weak form according to Eq. 24 and 25.

	The order of the functions in the solution and previous solution mixed functions must respect (p^l, u^s) in (P1, P2). 
	
	Inputs:
		- The functions at the previous time step:  previous_solution 		(dolfinx functions)
		- The functions increments to be computed:  solution     			(dolfinx functions)
		- The test functions: 					    ql, w 					(dolfinx Testfunctions)
		- The solid parameters: 					Young, Poisson 			(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The IF parameters: 						k_l, mu_l 				(dolfinx scalartype or dolfinx constant or dolfinx function)
		- Elementary integral elements:				dx, ds 					(ufl objects)
		- The time increment: 					    dt 						(dolfinx scalartype or dolfinx constant)
		- The name of the material law to be used: 	material_law 			(function name)
	Output: 
		- The weak form without BCs (Robin or Neumann)
	"""
	import ufl
	# Split the functions
	pl , u     = ufl.split(solution)
	pl_n , u_n = ufl.split(previous_solution)
	theta      = 0.7
	# Mass conservation of the IF, Eq. 24
	MassIF     = theta * k_l/mu_l*ufl.dot(ufl.grad(pl),ufl.grad(ql))*dx + (1-theta)* k_l/mu_l*ufl.dot(ufl.grad(pl_n),ufl.grad(ql))*dx
	MassIF    += 1/dt * ufl.nabla_div(u-u_n)*ql*dx
	# Momentum conservation, Eq. 25
	Momentum   = theta* ufl.inner(ufl.grad(w),material_law(u,Young,Poisson))*dx + (1-theta)*ufl.inner(ufl.grad(w),material_law(u_n,Young,Poisson))*dx 
	Momentum  += - theta * pl * ufl.nabla_div(w)*dx - (1-theta) * pl_n * ufl.nabla_div(w)*dx
	return MassIF + Momentum
# 
def variational_form_1_mono_updated(previous_solution, solution, ql, w, Young, Poisson, k_l, mu_l, dx, ds, dt, material_law):
	"""
	This function provides the incremental weak form according to Eq. 24 and 25.

	The order of the functions in the solution and previous solution mixed functions must respect (p^l, u^s) in (P1, P2). 
	
	Inputs:
		- The functions at the previous time step:  previous_solution 		(dolfinx functions)
		- The functions increments to be computed:  solution     			(dolfinx functions)
		- The test functions: 					    ql, w 					(dolfinx Testfunctions)
		- The solid parameters: 					Young, Poisson 			(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The IF parameters: 						k_l, mu_l 				(dolfinx scalartype or dolfinx constant or dolfinx function)
		- Elementary integral elements:				dx, ds 					(ufl objects)
		- The time increment: 					    dt 						(dolfinx scalartype or dolfinx constant)
		- The name of the material law to be used: 	material_law 			(function name)
	Output: 
		- The weak form without BCs (Robin or Neumann)
	"""
	import ufl
	# Split the functions
	dpl , du   = ufl.split(solution)
	pl_n , u_n = ufl.split(previous_solution)
	theta      = 0.7
	# Mass conservation of the IF, Eq. 24
	MassIF     = theta  * k_l/mu_l*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx + (1-theta) * k_l/mu_l*ufl.dot(ufl.grad(pl_n),ufl.grad(ql))*dx 
	MassIF    += 1/dt * ufl.nabla_div(du)*ql*dx
	# Momentum conservation, Eq. 25
	Momentum   = theta * ufl.inner(ufl.grad(w),material_law(u_n+du,Young,Poisson))*dx + (1-theta) * ufl.inner(ufl.grad(w),material_law(u_n,Young,Poisson))*dx
	Momentum  += - theta * (pl_n+dpl)*ufl.nabla_div(w)*dx - (1-theta) * (pl_n)*ufl.nabla_div(w)*dx 
	return MassIF + Momentum
# 
# 
#____________________________________________________________#
# 				Two Compartment Monophasic
#____________________________________________________________#
# 
# Isotropic K_b
def variational_form_2_mono_updated_isotropic_k_b(previous_solution, solution, ql, qb, w, Young, Poisson, k_l, mu_l, k_b, mu_b, epsb0, K, alpha, dx, ds, dt, material_law,varepsilon_b_law):
	"""
	This function provides the incremental weak form according to Eq. 60, 61 and 62. The vascular permeability is here isotropic.

	The order of the functions in the solution and previous solution mixed functions must respect (p^l, p^b, u^s) in (P1, P1, P2). 
	
	Inputs:
		- The functions at the previous time step:  previous_solution 			(dolfinx functions)
		- The functions increments to be computed:  solution     				(dolfinx functions)
		- The test functions: 					    ql, qb, w					(dolfinx Testfunctions)
		- The solid parameters: 					Young, Poisson 				(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The IF parameters: 						k_l, mu_l 					(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The vascular parameters: 					k_b, mu_b, epsb0, K, alpha	(dolfinx scalartype or dolfinx constant or dolfinx function)
		- Elementary integral elements:				dx, ds 						(ufl objects)
		- The time increment: 					    dt 							(dolfinx scalartype or dolfinx constant)
		- The name of the material law to be used: 	material_law 				(function name)
		- The type of vascular porosity evolution:  varepsilon_b_law       		(string equal to 'linear' or 'atan')
	Output: 
		- The weak form without BCs (Robin or Neumann)
	"""
	import ufl
	# Split the functions
	dpl , dpb , du    = ufl.split(solution)
	pl_n , pb_n , u_n = ufl.split(previous_solution)
	theta             = 0.7
	if varepsilon_b_law == 'linear':
		# Mass conservation of the IF, Eq. 60
		MassIF            = theta      * constitutive_laws.vascular_permeability_scalar_mono(k_b,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)/mu_b*ufl.dot(ufl.grad(pb_n+dpb),ufl.grad(ql))*dx
		MassIF           += (1-theta)  * constitutive_laws.vascular_permeability_scalar_mono(k_b,alpha,epsb0,pl_n,pb_n,K,varepsilon_b_law)/mu_b*ufl.dot(ufl.grad(pb_n),ufl.grad(ql))*dx
		MassIF           += theta      * k_l/mu_l*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
		MassIF           += (1-theta)  * k_l/mu_l*ufl.dot(ufl.grad(pl_n),ufl.grad(ql))*dx
		MassIF           += 1/dt * ufl.nabla_div(du)*ql*dx
		# Mass conservation of the IF, Eq. 61
		MassBlood         =   1/dt    * constitutive_laws.Cep_linear_mono(epsb0,K)*dpl*qb*dx
		MassBlood        += - 1/dt    * constitutive_laws.Cep_linear_mono(epsb0,K)*dpb*qb*dx
		MassBlood        += theta     * constitutive_laws.vascular_permeability_scalar_mono(k_b,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)/mu_b*ufl.dot(ufl.grad(pb_n+dpb),ufl.grad(qb))*dx
		MassBlood        += (1-theta) * constitutive_laws.vascular_permeability_scalar_mono(k_b,alpha,epsb0,pl_n,pb_n,K,varepsilon_b_law)/mu_b*ufl.dot(ufl.grad(pb_n),ufl.grad(qb))*dx 
		MassBlood        +=   1/dt    * constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K) * ufl.nabla_div(du)*qb*dx
		# Momentum conservation, Eq. 62
		Momentum          = theta      * ufl.inner(ufl.grad(w),material_law(u_n+du,Young,Poisson))*dx
		Momentum         += (1-theta)  * ufl.inner(ufl.grad(w),material_law(u_n,Young,Poisson))*dx
		Momentum         += - theta    * (1-constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K) )*(pl_n+dpl)*ufl.nabla_div(w)*dx
		Momentum         += -(1-theta) * (1-constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n,pb_n,K) )*(pl_n)*ufl.nabla_div(w)*dx
		Momentum         += -theta     * constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K)*(pb_n+dpb)*ufl.nabla_div(w)*dx 
		Momentum         += -(1-theta) * constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n,pb_n,K)*(pb_n)*ufl.nabla_div(w)*dx 
	elif varepsilon_b_law == 'atan':
		# Mass conservation of the IF, Eq. 60
		MassIF            = theta     * constitutive_laws.vascular_permeability_scalar_mono(k_b,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)/mu_b*ufl.dot(ufl.grad(pb_n+dpb),ufl.grad(ql))*dx
		MassIF           += (1-theta) * constitutive_laws.vascular_permeability_scalar_mono(k_b,alpha,epsb0,pl_n,pb_n,K,varepsilon_b_law)/mu_b*ufl.dot(ufl.grad(pb_n),ufl.grad(ql))*dx
		MassIF           += theta     * k_l/mu_l*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
		MassIF           += (1-theta) * k_l/mu_l*ufl.dot(ufl.grad(pl_n),ufl.grad(ql))*dx
		MassIF           += 1/dt      * ufl.nabla_div(du)*ql*dx
		# Mass conservation of the IF, Eq. 61
		MassBlood         =  1/dt     * constitutive_laws.Cep_atan_mono(epsb0,K,pl_n+dpl,pb_n+dpb)*dpl*qb*dx
		MassBlood        += -1/dt     * constitutive_laws.Cep_atan_mono(epsb0,K,pl_n+dpl,pb_n+dpb)*dpb*qb*dx
		MassBlood        += theta     * constitutive_laws.vascular_permeability_scalar_mono(k_b,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)/mu_b*ufl.dot(ufl.grad(pb_n+dpb),ufl.grad(qb))*dx
		MassBlood        += (1-theta) * constitutive_laws.vascular_permeability_scalar_mono(k_b,alpha,epsb0,pl_n,pb_n,K,varepsilon_b_law)/mu_b*ufl.dot(ufl.grad(pb_n),ufl.grad(qb))*dx 
		MassBlood        +=  1/dt     * constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K) * ufl.nabla_div(du)*qb*dx
		# Momentum conservation, Eq. 62
		Momentum          = theta     * ufl.inner(ufl.grad(w),material_law(u_n+du,Young,Poisson))*dx
		Momentum         += (1-theta) * ufl.inner(ufl.grad(w),material_law(u_n,Young,Poisson))*dx
		Momentum         += -theta    * (1-constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K) )*(pl_n+dpl)*ufl.nabla_div(w)*dx
		Momentum         += -(1-theta)* (1-constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n,pb_n,K) )*(pl_n)*ufl.nabla_div(w)*dx 
		Momentum         += -theta    * constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K)*(pb_n+dpb)*ufl.nabla_div(w)*dx
		Momentum         += -(1-theta)* constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n,pb_n,K)*(pb_n)*ufl.nabla_div(w)*dx 
	else:
		print("ERROR in the choice of the evolution type. Choose 'linear' or 'atan'.")
	return MassIF + MassBlood + Momentum
# 
# Anisotropic k_b
def variational_form_2_mono_updated_anisotropic_k_b(previous_solution, solution, ql, qb, w, Young, Poisson, k_l, mu_l, k_b, ratio, mu_b, epsb0, K, alpha, dx, ds, dt, material_law,varepsilon_b_law):
	"""
	This function provides the incremental weak form according to Eq. 60, 62 and 62. The vascular permeability is here anisotropic.

	The order of the functions in the solution and previous solution mixed functions must respect (p^l, p^b, u^s) in (P1, P1, P2). 
	
	Inputs:
		- The functions at the previous time step:  previous_solution 				(dolfinx functions)
		- The functions increments to be computed:  solution     					(dolfinx functions)
		- The test functions: 					    ql, qb, w						(dolfinx Testfunctions)
		- The solid parameters: 					Young, Poisson 					(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The IF parameters: 						k_l, mu_l 						(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The vascular parameters: 					k_b, mu_b, epsb0, K, alpha		(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The vascular permability ratio:			ratio 							(list of scalars) 
		- Elementary integral elements:				dx, ds 							(ufl objects)
		- The time increment: 					    dt 								(dolfinx scalartype or dolfinx constant)
		- The name of the material law to be used: 	material_law 					(function name)
		- The type of vascular porosity evolution:  varepsilon_b_law       			(string equal to 'linear' or 'atan')
	Output: 
		- The weak form without BCs (Robin or Neumann)
	"""
	import ufl
	# Split the functions
	dpl , dpb , du    = ufl.split(solution)
	pl_n , pb_n , u_n = ufl.split(previous_solution)
	# 
	theta = 0.7
	if varepsilon_b_law == 'linear':
		# Mass conservation of the IF, Eq. 60
		MassIF            = theta*1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_mono(k_b, ratio,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)*ufl.grad(pb_n+dpb),ufl.grad(ql))*dx
		MassIF           += (1-theta)*1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_mono(k_b, ratio,alpha,epsb0,pl_n,pb_n,K,varepsilon_b_law)*ufl.grad(pb_n),ufl.grad(ql))*dx
		# 
		MassIF           += theta*k_l/mu_l*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
		MassIF           += (1-theta)*k_l/mu_l*ufl.dot(ufl.grad(pl_n),ufl.grad(ql))*dx
		# 
		MassIF           += 1/dt * ufl.nabla_div(du)*ql*dx
		# Mass conservation of the IF, Eq. 61
		MassBlood         = constitutive_laws.Cep_linear_mono(epsb0,K)*dpl/dt*qb*dx
		MassBlood        += -constitutive_laws.Cep_linear_mono(epsb0,K)*dpb/dt*qb*dx
		# 
		MassBlood        += theta*1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_mono(k_b, ratio,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)*ufl.grad(pb_n+dpb),ufl.grad(qb))*dx
		MassBlood        += (1-theta)*1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_mono(k_b, ratio,alpha,epsb0,pl_n,pb_n,K,varepsilon_b_law)*ufl.grad(pb_n),ufl.grad(qb))*dx
		# 
		MassBlood        += 1/dt* constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K) * ufl.nabla_div(du)*qb*dx
		# Momentum conservation, Eq. 62
		Momentum          = theta*ufl.inner(ufl.grad(w),material_law(u_n+du,Young,Poisson))*dx 
		Momentum         += (1-theta)*ufl.inner(ufl.grad(w),material_law(u_n,Young,Poisson))*dx 
		# 
		Momentum         += -theta* (1-constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K) )*(pl_n+dpl)*ufl.nabla_div(w)*dx 
		Momentum         += -(1-theta)* (1-constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n,pb_n,K) )*(pl_n)*ufl.nabla_div(w)*dx 
		# 
		Momentum         += -theta* constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K)*(pb_n+dpb)*ufl.nabla_div(w)*dx 
		Momentum         += -(1-theta)* constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n,pb_n,K)*(pb_n)*ufl.nabla_div(w)*dx 
	elif varepsilon_b_law == 'atan':
		# Mass conservation of the IF, Eq. 60
		MassIF            = theta*1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_mono(k_b, ratio,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)*ufl.grad(pb_n+dpb),ufl.grad(ql))*dx
		MassIF           += (1-theta)*1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_mono(k_b, ratio,alpha,epsb0,pl_n,pb_n,K,varepsilon_b_law)*ufl.grad(pb_n),ufl.grad(ql))*dx
		# 
		MassIF           += theta*k_l/mu_l*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
		MassIF           += (1-theta)*k_l/mu_l*ufl.dot(ufl.grad(pl_n),ufl.grad(ql))*dx
		# 
		MassIF           += 1/dt * ufl.nabla_div(du)*ql*dx
		# Mass conservation of the IF, Eq. 61
		MassBlood         = constitutive_laws.Cep_atan_mono(epsb0,K,pl_n+dpl,pb_n+dpb)*dpl/dt*qb*dx
		MassBlood        += -constitutive_laws.Cep_atan_mono(epsb0,K,pl_n+dpl,pb_n+dpb)*dpb/dt*qb*dx
		# 
		MassBlood        += theta*1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_mono(k_b, ratio,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)*ufl.grad(pb_n+dpb),ufl.grad(qb))*dx
		MassBlood        += (1-theta)*1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_mono(k_b, ratio,alpha,epsb0,pl_n,pb_n,K,varepsilon_b_law)*ufl.grad(pb_n),ufl.grad(qb))*dx
		# 
		MassBlood        += 1/dt* constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K) * ufl.nabla_div(du)*qb*dx
		# Momentum conservation, Eq. 62
		Momentum          = theta*ufl.inner(ufl.grad(w),material_law(u_n+du,Young,Poisson))*dx
		Momentum         += (1-theta)*ufl.inner(ufl.grad(w),material_law(u_n,Young,Poisson))*dx 
		# 
		Momentum         += -theta* (1-constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K) )*(pl_n+dpl)*ufl.nabla_div(w)*dx 
		Momentum         += -(1-theta)* (1-constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n,pb_n,K) )*(pl_n)*ufl.nabla_div(w)*dx 
		# 
		Momentum         += -theta* constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K)*(pb_n+dpb)*ufl.nabla_div(w)*dx 
		Momentum         += -(1-theta)* constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n,pb_n,K)*(pb_n)*ufl.nabla_div(w)*dx 
	else:
		print("ERROR in the choice of the evolution type. Choose 'linear' or 'atan'.")
	return MassIF + MassBlood + Momentum
# 
#...........................................................#
# 					 Biphasic
#...........................................................#
# 
#____________________________________________________________#
# 				Single Compartment Biphasic
#____________________________________________________________#
# 
def variational_form_1_comp_bi(previous_solution, solution, qc, ql, w, poro_n, poro_min, poro_max, dt, Young, Poisson, k_l, mu_l, k_c, mu_c, a, dx, ds):
	"""
	This function provides the incremental weak form according to Eq. 111-113.

	The order of the functions in the solution and previous solution mixed functions must respect (p^l, p^lc, u^s) in (P1, P1, P2). 

	Inputs:
		- The functions at the previous time step:  previous_solution 		(dolfinx functions)
		- The functions increments to be computed:  solution 	 			(dolfinx functions)
		- The test functions: 					    qc, ql, w 				(dolfinx Testfunctions)
		- The porosity at the previous time step:   poro_n 					(dolfinx functions)
		- The time increment: 					    dt 						(dolfinx scalartype or dolfinx constant)
		- The solid parameters: 					Young, Poisson 			(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The IF parameters: 						k_l, mu_l 				(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The cell parameters: 						k_c, mu_c 				(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The parameter for the saturation: 		a 						(dolfinx scalartype or dolfinx constant or dolfinx function)
		- dx, ds resulting from Measure and cell_tag/facet_tag mapping.
	Output: 
		- The weak form without BCs (Robin or Neumann)
	"""
	import ufl
	# Split the functions
	dpl , dplc , du    = ufl.split(solution)
	pl_n , plc_n , u_n = ufl.split(previous_solution)
	theta = 0.7
	# Mass conservation of the cells, Eq. 111
	MassCell  = constitutive_laws.Cmc(poro_n,plc_n+dplc,a) * dplc/dt * qc * dx 
	# 
	MassCell += theta * (k_c/mu_c)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(qc))*dx
	MassCell += (1-theta) * (k_c/mu_c)*ufl.dot(ufl.grad(pl_n),ufl.grad(qc))*dx
	# 
	MassCell += -theta* (k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(qc))*dx 
	MassCell += -(1-theta)* (k_c/mu_c)*ufl.dot(ufl.grad(plc_n),ufl.grad(qc))*dx 
	# 
	MassCell += constitutive_laws.saturation_cell(plc_n+dplc,a)* 1/dt * ufl.nabla_div(du)*qc*dx
	# Mass conservation of the IF, Eq. 112
	MassIF    = theta * (k_c/mu_c + k_l/mu_l)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
	MassIF   += (1-theta) * (k_c/mu_c + k_l/mu_l)*ufl.dot(ufl.grad(pl_n),ufl.grad(ql))*dx
	# 
	MassIF   += - theta* (k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(ql))*dx 
	MassIF   += - (1-theta)* (k_c/mu_c)*ufl.dot(ufl.grad(plc_n),ufl.grad(ql))*dx 
	# 
	MassIF   +=  1/dt * ufl.nabla_div(du)*ql*dx
	# Momentum conservation, Eq. 113
	Momentum  = theta* ufl.inner(ufl.grad(w),constitutive_laws.Elastic_constitutive_law(u_n+du,Young,Poisson))*dx 
	Momentum += (1-theta)* ufl.inner(ufl.grad(w),constitutive_laws.Elastic_constitutive_law(u_n,Young,Poisson))*dx 
	# 
	Momentum += -theta* (pl_n+dpl)*ufl.nabla_div(w)*dx 
	Momentum += -(1-theta)* (pl_n)*ufl.nabla_div(w)*dx 
	# 
	Momentum += theta* constitutive_laws.saturation_cell(plc_n+dplc,a)*(plc_n+dplc)*ufl.nabla_div(w)*dx
	Momentum += (1-theta)* constitutive_laws.saturation_cell(plc_n,a)*(plc_n)*ufl.nabla_div(w)*dx
	return MassCell + MassIF + Momentum
# 
#____________________________________________________________#
# 				Two Compartments Biphasic
#____________________________________________________________#
# 
# Matrix k_b
def variational_form_2_comp_bi(previous_solution, solution, qc, ql, qb, w, poro_n, poro_min, poro_max, dt, Young, Poisson, k_l, mu_l, k_c, mu_c, k_b, mu_b, a, K, epsb0, alpha, ratio, dx, ds):
	"""
	This function provides the incremental weak form according to Eq. 155-158.

	The order of the functions in the solution and previous solution mixed functions must respect (p^l, p^lc, p^b, u^s) in (P1, P1, P1, P2). 

	Inputs:
		- The functions at the previous time step:  previous_solution		(dolfinx functions)
		- The functions increments to be computed:  solution				(dolfinx functions)
		- The test functions: 					    qc, ql, qb, w			(dolfinx Testfunctions)
		- The porosity at the previous time step:   poro_n 					(dolfinx functions)
		- The time increment: 					    dt 						(dolfinx scalartype or dolfinx constant)
		- The solid parameters: 					Young, Poisson 			(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The IF parameters: 						k_l, mu_l 				(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The cell parameters: 						k_c, mu_c 				(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The varscular parameters:					k_b, mu_b 				(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The parameter for the saturation: 		a 						(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The vessel compressibility: 				K 						(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The initial vascular potosity:			K 						(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The exponent from Eq. 63:					alpha					(dolfinx scalartype)
		- Ratio for preferential direction kb:		ratio 					(dolfinx scalartype)
		- Preferential Direction for kb matrix:		direction 				(integer)
		- Dimension of the problem 2D or 3D:		dimension				(integer)
		- dx, ds resulting from Measure and cell_tag/facet_tag mapping
	Output: 
		- The weak form without BCs (Robin or Neumann)
	"""
	import ufl
	# split the solution into: displacement increment, IF pressure increment, Blood Pressure increment
	dpl, dplc, dpb, du    = ufl.split(solution)
	# Previous time steps solutions
	pl_n, plc_n, pb_n, u_n = ufl.split(previous_solution)
	# 
	theta = 0.7
	# 
	k_b_current   = constitutive_laws.matrix_vascular_permeability_bi(k_b, ratio, alpha, epsb0, pl_n+dpl, plc_n+dplc, pb_n+dpb, a, K)
	k_b_current_n = constitutive_laws.matrix_vascular_permeability_bi(k_b, ratio, alpha, epsb0, pl_n, plc_n, pb_n, a, K)
	# 
	S_c           = constitutive_laws.saturation_cell(plc_n+dplc,a)
	S_c_n         = constitutive_laws.saturation_cell(plc_n,a)
	# 
	epsb          = constitutive_laws.vascular_porosity_atan_bi(epsb0,pl_n+dpl,plc_n+dplc,pb_n+dpb,a,K)
	epsb_n        = constitutive_laws.vascular_porosity_atan_bi(epsb0,pl_n,plc_n,pb_n,a,K)
	# 
	coeff_CEP     = constitutive_laws.Cep_atan_bi(epsb0,K,plc_n+dplc,pl_n+dpl,pb_n+dpb,a)
	# 
	# Mass conservation of the cells, Eq. 155
	MassCell   = constitutive_laws.Cmc(poro_n,plc_n+dplc,a) * dplc/dt * qc * dx 
	# 
	MassCell  += theta* S_c*ufl.dot((k_b_current/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(qc))*dx 
	MassCell  += (1-theta)* S_c_n*ufl.dot((k_b_current_n/mu_b)*ufl.grad(pb_n),ufl.grad(qc))*dx 
	# 
	MassCell  += theta* (k_c/mu_c)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(qc))*dx 
	MassCell  += (1-theta)* (k_c/mu_c)*ufl.dot(ufl.grad(pl_n),ufl.grad(qc))*dx 
	# 
	MassCell  += -theta* (k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(qc))*dx 
	MassCell  += -(1-theta)* (k_c/mu_c)*ufl.dot(ufl.grad(plc_n),ufl.grad(qc))*dx 
	# 
	MassCell  += S_c* 1/dt * ufl.nabla_div(du)*qc*dx
	# 
	# Mass conservation of the IF, Eq. 156
	MassIF     = theta* ufl.dot((k_b_current/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(ql))*dx 
	MassIF    += (1-theta)* ufl.dot((k_b_current_n/mu_b)*ufl.grad(pb_n),ufl.grad(ql))*dx 
	# 
	MassIF    += theta* (k_c/mu_c + k_l/mu_l)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx 
	MassIF    += (1-theta)* (k_c/mu_c + k_l/mu_l)*ufl.dot(ufl.grad(pl_n),ufl.grad(ql))*dx 
	# 
	MassIF    += -theta* (k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(ql))*dx 
	MassIF    += -(1-theta)* (k_c/mu_c)*ufl.dot(ufl.grad(plc_n),ufl.grad(ql))*dx 
	# 
	MassIF    += 1/dt * ufl.nabla_div(du)*ql*dx
	# Mass conservation of the blood, Eq. 157
	MassBlood  = coeff_CEP * dpl/dt * qb* dx 
	MassBlood += - coeff_CEP * constitutive_laws.Cstate(plc_n+dplc,a) * dplc/dt * qb * dx 
	MassBlood += - coeff_CEP * dpb/dt * qb * dx  
	# 
	MassBlood +=  theta*ufl.dot((k_b_current/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(qb))*dx 
	MassBlood +=  (1-theta)*ufl.dot((k_b_current_n/mu_b)*ufl.grad(pb_n),ufl.grad(qb))*dx 
	# 
	MassBlood +=  epsb * 1/dt * ufl.nabla_div(du)*qb*dx
	# Momentum conservation, Eq. 158
	Momentum   = theta*ufl.inner(ufl.sym(ufl.grad(w)),constitutive_laws.Elastic_constitutive_law(u_n+du,Young,Poisson))*dx 
	Momentum  += (1-theta)* ufl.inner(ufl.sym(ufl.grad(w)),constitutive_laws.Elastic_constitutive_law(u_n,Young,Poisson))*dx 
	# 
	Momentum  += -theta* (1-epsb)*(pl_n+dpl)*ufl.nabla_div(w)*dx 
	Momentum  += -(1-theta)* (1-epsb_n)*(pl_n)*ufl.nabla_div(w)*dx 
	# 
	Momentum  += theta*(1-epsb)*S_c*(plc_n+dplc)*ufl.nabla_div(w)*dx 
	Momentum  += (1-theta)* (1-epsb_n)*S_c_n*(plc_n)*ufl.nabla_div(w)*dx 
	# 
	Momentum  += - theta * epsb * (pb_n+dpb) * ufl.nabla_div(w)*dx
	Momentum  += - (1-theta) * epsb_n * (pb_n) * ufl.nabla_div(w)*dx
	return MassCell + MassIF + MassBlood + Momentum
# 
# Scalar k_b
def variational_form_2_comp_bi_function_DG0(previous_solution, solution, qc, ql, qb, w, poro_n, poro_min, poro_max, dt, Young, Poisson, k_l, mu_l, k_c, mu_c, k_b, mu_b, a, K, epsb0, alpha, dx, ds):
	"""
	This function provides the incremental weak form according to Eq. 155-158.

	The order of the functions in the solution and previous solution mixed functions must respect (p^l, p^lc, p^b, u^s) in (P1, P1, P1, P2). 

	Inputs:
		- The functions at the previous time step:  previous_solution		(dolfinx functions)
		- The functions increments to be computed:  solution				(dolfinx functions)
		- The test functions: 					    qc, ql, qb, w			(dolfinx Testfunctions)
		- The porosity at the previous time step:   poro_n 					(dolfinx functions)
		- The time increment: 					    dt 						(dolfinx scalartype or dolfinx constant)
		- The solid parameters: 					Young, Poisson 			(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The IF parameters: 						k_l, mu_l 				(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The cell parameters: 						k_c, mu_c 				(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The varscular parameters:					k_b, mu_b 				(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The parameter for the saturation: 		a 						(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The vessel compressibility: 				K 						(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The initial vascular potosity:			K 						(dolfinx scalartype or dolfinx constant or dolfinx function)
		- The exponent from Eq. 63:					alpha					(dolfinx scalartype)
		- Ratio for preferential direction kb:		ratio 					(dolfinx scalartype)
		- Preferential Direction for kb matrix:		direction 				(integer)
		- Dimension of the problem 2D or 3D:		dimension				(integer)
		- dx, ds resulting from Measure and cell_tag/facet_tag mapping
	Output: 
		- The weak form without BCs (Robin or Neumann)
	"""
	import ufl
	# split the solution into: displacement increment, IF pressure increment, Blood Pressure increment
	dpl, dplc, dpb, du    = ufl.split(solution)
	# Previous time steps solutions
	pl_n, plc_n, pb_n, u_n = ufl.split(previous_solution)
	# 
	theta = 0.7
	# 
	k_b_current   = constitutive_laws.vascular_permeability_scalar_bi(k_b, alpha, epsb0, pl_n+dpl, plc_n+dplc, pb_n+dpb, a, K)
	k_b_current_n = constitutive_laws.vascular_permeability_scalar_bi(k_b, alpha, epsb0, pl_n, plc_n, pb_n, a, K)
	# 
	S_c           = constitutive_laws.saturation_cell(plc_n+dplc,a)
	S_c_n         = constitutive_laws.saturation_cell(plc_n,a)
	# 
	epsb          = constitutive_laws.vascular_porosity_atan_bi(epsb0,pl_n+dpl,plc_n+dplc,pb_n+dpb,a,K)
	epsb_n        = constitutive_laws.vascular_porosity_atan_bi(epsb0,pl_n,plc_n,pb_n,a,K)
	# 
	coeff_CEP     = constitutive_laws.Cep_atan_bi(epsb0,K,plc_n+dplc,pl_n+dpl,pb_n+dpb,a)
	# 
	# Mass conservation of the cells, Eq. 155
	MassCell   = constitutive_laws.Cmc(poro_n,plc_n+dplc,a) * dplc/dt * qc * dx 
	# 
	MassCell  += theta* S_c*ufl.dot((k_b_current/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(qc))*dx 
	MassCell  += (1-theta)* S_c_n*ufl.dot((k_b_current_n/mu_b)*ufl.grad(pb_n),ufl.grad(qc))*dx 
	# 
	MassCell  += theta* (k_c/mu_c)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(qc))*dx 
	MassCell  += (1-theta)* (k_c/mu_c)*ufl.dot(ufl.grad(pl_n),ufl.grad(qc))*dx 
	# 
	MassCell  += -theta* (k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(qc))*dx 
	MassCell  += -(1-theta)* (k_c/mu_c)*ufl.dot(ufl.grad(plc_n),ufl.grad(qc))*dx 
	# 
	MassCell  += S_c* 1/dt * ufl.nabla_div(du)*qc*dx
	# 
	# Mass conservation of the IF, Eq. 156
	MassIF     = theta* ufl.dot((k_b_current/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(ql))*dx 
	MassIF    += (1-theta)* ufl.dot((k_b_current_n/mu_b)*ufl.grad(pb_n),ufl.grad(ql))*dx 
	# 
	MassIF    += theta* (k_c/mu_c + k_l/mu_l)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx 
	MassIF    += (1-theta)* (k_c/mu_c + k_l/mu_l)*ufl.dot(ufl.grad(pl_n),ufl.grad(ql))*dx 
	# 
	MassIF    += -theta* (k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(ql))*dx 
	MassIF    += -(1-theta)* (k_c/mu_c)*ufl.dot(ufl.grad(plc_n),ufl.grad(ql))*dx 
	# 
	MassIF    += 1/dt * ufl.nabla_div(du)*ql*dx
	# Mass conservation of the blood, Eq. 157
	MassBlood  = coeff_CEP * dpl/dt * qb* dx 
	MassBlood += - coeff_CEP * constitutive_laws.Cstate(plc_n+dplc,a) * dplc/dt * qb * dx 
	MassBlood += - coeff_CEP * dpb/dt * qb * dx  
	# 
	MassBlood +=  theta*ufl.dot((k_b_current/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(qb))*dx 
	MassBlood +=  (1-theta)*ufl.dot((k_b_current_n/mu_b)*ufl.grad(pb_n),ufl.grad(qb))*dx 
	# 
	MassBlood +=  epsb * 1/dt * ufl.nabla_div(du)*qb*dx
	# Momentum conservation, Eq. 158
	Momentum   = theta*ufl.inner(ufl.sym(ufl.grad(w)),constitutive_laws.Elastic_constitutive_law_DG0(u_n+du,Young,Poisson))*dx 
	Momentum  += (1-theta)* ufl.inner(ufl.sym(ufl.grad(w)),constitutive_laws.Elastic_constitutive_law_DG0(u_n,Young,Poisson))*dx 
	# 
	Momentum  += -theta* (1-epsb)*(pl_n+dpl)*ufl.nabla_div(w)*dx 
	Momentum  += -(1-theta)* (1-epsb_n)*(pl_n)*ufl.nabla_div(w)*dx 
	# 
	Momentum  += theta*(1-epsb)*S_c*(plc_n+dplc)*ufl.nabla_div(w)*dx 
	Momentum  += (1-theta)* (1-epsb_n)*S_c_n*(plc_n)*ufl.nabla_div(w)*dx 
	# 
	Momentum  += - theta * epsb * (pb_n+dpb) * ufl.nabla_div(w)*dx
	Momentum  += - (1-theta) * epsb_n * (pb_n) * ufl.nabla_div(w)*dx
	return MassCell + MassIF + MassBlood + Momentum
#____________________________________________________________#
# 					 End of functions
#____________________________________________________________#
# 
# 
if __name__ == "__main__":
    print("Loading of the constitutive laws successfully completed.")
    # EoF