# This document holds the functions used as part of 
# the poro-mechanical modelling of living tissues 
# according to the theory in ./Theoretical_Developments.pdf
# 
# Author: Thomas Lavigne
# Date: 18/09/2024
# 
import constitutive_laws
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
	# Mass conservation of the IF, Eq. 24
	MassIF     = k_l/mu_l*ufl.dot(ufl.grad(pl),ufl.grad(ql))*dx
	MassIF    += 1/dt * ufl.nabla_div(u-u_n)*ql*dx
	# Momentum conservation, Eq. 25
	Momentum   = ufl.inner(ufl.grad(w),material_law(u,Young,Poisson))*dx 
	Momentum  += - pl*ufl.nabla_div(w)*dx 
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
	# Mass conservation of the IF, Eq. 24
	MassIF     = k_l/mu_l*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
	MassIF    += 1/dt * ufl.nabla_div(du)*ql*dx
	# Momentum conservation, Eq. 25
	Momentum   = ufl.inner(ufl.grad(w),material_law(u_n+du,Young,Poisson))*dx 
	Momentum  += - (pl_n+dpl)*ufl.nabla_div(w)*dx 
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
	if varepsilon_b_law == 'linear':
		# Mass conservation of the IF, Eq. 60
		MassIF            = constitutive_laws.vascular_permeability_scalar_mono(k_b,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)/mu_b*ufl.dot(ufl.grad(pb_n+dpb),ufl.grad(ql))*dx
		MassIF           += k_l/mu_l*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
		MassIF           += 1/dt * ufl.nabla_div(du)*ql*dx
		# Mass conservation of the IF, Eq. 61
		MassBlood         = constitutive_laws.Cep_linear_mono(epsb0,K)*dpl/dt*qb*dx
		MassBlood        += -constitutive_laws.Cep_linear_mono(epsb0,K)*dpb/dt*qb*dx
		MassBlood        += constitutive_laws.vascular_permeability_scalar_mono(k_b,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)/mu_b*ufl.dot(ufl.grad(pb_n+dpb),ufl.grad(qb))*dx
		MassBlood        += 1/dt* constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K) * ufl.nabla_div(du)*qb*dx
		# Momentum conservation, Eq. 62
		Momentum          = ufl.inner(ufl.grad(w),material_law(u_n+du,Young,Poisson))*dx 
		Momentum         += - (1-constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K) )*(pl_n+dpl)*ufl.nabla_div(w)*dx 
		Momentum         += - constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K)*(pb_n+dpb)*ufl.nabla_div(w)*dx 
	elif varepsilon_b_law == 'atan':
		# Mass conservation of the IF, Eq. 60
		MassIF            = constitutive_laws.vascular_permeability_scalar_mono(k_b,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)/mu_b*ufl.dot(ufl.grad(pb_n+dpb),ufl.grad(ql))*dx
		MassIF           += k_l/mu_l*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
		MassIF           += 1/dt * ufl.nabla_div(du)*ql*dx
		# Mass conservation of the IF, Eq. 61
		MassBlood         = constitutive_laws.Cep_atan_mono(epsb0,K,pl_n+dpl,pb_n+dpb)*dpl/dt*qb*dx
		MassBlood        += -constitutive_laws.Cep_atan_mono(epsb0,K,pl_n+dpl,pb_n+dpb)*dpb/dt*qb*dx
		MassBlood        += constitutive_laws.vascular_permeability_scalar_mono(k_b,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)/mu_b*ufl.dot(ufl.grad(pb_n+dpb),ufl.grad(qb))*dx
		MassBlood        += 1/dt* constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K) * ufl.nabla_div(du)*qb*dx
		# Momentum conservation, Eq. 62
		Momentum          = ufl.inner(ufl.grad(w),material_law(u_n+du,Young,Poisson))*dx 
		Momentum         += - (1-constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K) )*(pl_n+dpl)*ufl.nabla_div(w)*dx 
		Momentum         += - constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K)*(pb_n+dpb)*ufl.nabla_div(w)*dx 
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
	if varepsilon_b_law == 'linear':
		# Mass conservation of the IF, Eq. 60
		MassIF            = 1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_mono(k_b, ratio,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)*ufl.grad(pb_n+dpb),ufl.grad(ql))*dx
		MassIF           += k_l/mu_l*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
		MassIF           += 1/dt * ufl.nabla_div(du)*ql*dx
		# Mass conservation of the IF, Eq. 61
		MassBlood         = constitutive_laws.Cep_linear_mono(epsb0,K)*dpl/dt*qb*dx
		MassBlood        += -constitutive_laws.Cep_linear_mono(epsb0,K)*dpb/dt*qb*dx
		MassBlood        += 1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_mono(k_b, ratio,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)*ufl.grad(pb_n+dpb),ufl.grad(qb))*dx
		MassBlood        += 1/dt* constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K) * ufl.nabla_div(du)*qb*dx
		# Momentum conservation, Eq. 62
		Momentum          = ufl.inner(ufl.grad(w),material_law(u_n+du,Young,Poisson))*dx 
		Momentum         += - (1-constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K) )*(pl_n+dpl)*ufl.nabla_div(w)*dx 
		Momentum         += - constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K)*(pb_n+dpb)*ufl.nabla_div(w)*dx 
	elif varepsilon_b_law == 'atan':
		# Mass conservation of the IF, Eq. 60
		MassIF            = 1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_mono(k_b, ratio,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)*ufl.grad(pb_n+dpb),ufl.grad(ql))*dx
		MassIF           += k_l/mu_l*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
		MassIF           += 1/dt * ufl.nabla_div(du)*ql*dx
		# Mass conservation of the IF, Eq. 61
		MassBlood         = constitutive_laws.Cep_atan_mono(epsb0,K,pl_n+dpl,pb_n+dpb)*dpl/dt*qb*dx
		MassBlood        += -constitutive_laws.Cep_atan_mono(epsb0,K,pl_n+dpl,pb_n+dpb)*dpb/dt*qb*dx
		MassBlood        += 1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_mono(k_b, ratio,alpha,epsb0,pl_n+dpl,pb_n+dpb,K,varepsilon_b_law)*ufl.grad(pb_n+dpb),ufl.grad(qb))*dx
		MassBlood        += 1/dt* constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K) * ufl.nabla_div(du)*qb*dx
		# Momentum conservation, Eq. 62
		Momentum          = ufl.inner(ufl.grad(w),material_law(u_n+du,Young,Poisson))*dx 
		Momentum         += - (1-constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K) )*(pl_n+dpl)*ufl.nabla_div(w)*dx 
		Momentum         += - constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K)*(pb_n+dpb)*ufl.nabla_div(w)*dx 
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
	# Mass conservation of the cells, Eq. 111
	MassCell  = constitutive_laws.Cmc(poro_n,plc_n+dplc,a) * dplc/dt * qc * dx 
	MassCell += (k_c/mu_c)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(qc))*dx
	MassCell += - (k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(qc))*dx 
	MassCell += constitutive_laws.saturation_cell(plc_n+dplc,a)* 1/dt * ufl.nabla_div(du)*qc*dx
	# Mass conservation of the IF, Eq. 112
	MassIF    = (k_c/mu_c + k_l/mu_l)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
	MassIF   += - (k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(ql))*dx 
	MassIF   +=  1/dt * ufl.nabla_div(du)*ql*dx
	# Momentum conservation, Eq. 113
	Momentum  = ufl.inner(ufl.grad(w),constitutive_laws.Elastic_constitutive_law(u_n+du,Young,Poisson))*dx 
	Momentum += - (pl_n+dpl)*ufl.nabla_div(w)*dx 
	Momentum += constitutive_laws.saturation_cell(plc_n+dplc,a)*(plc_n+dplc)*ufl.nabla_div(w)*dx
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
	k_b_current = constitutive_laws.matrix_vascular_permeability_bi(k_b, ratio, alpha, epsb0, pl_n+dpl, plc_n+dplc, pb_n+dpb, a, K)
	S_c         = constitutive_laws.saturation_cell(plc_n+dplc,a)
	epsb        = constitutive_laws.vascular_porosity_atan_bi(epsb0,pl_n+dpl,plc_n+dplc,pb_n+dpb,a,K)
	coeff_CEP   = constitutive_laws.Cep_atan_bi(epsb0,K,plc_n+dplc,pl_n+dpl,pb_n+dpb,a)
	# Mass conservation of the cells, Eq. 155
	MassCell   = constitutive_laws.Cmc(poro_n,plc_n+dplc,a) * dplc/dt * qc * dx 
	MassCell  += S_c*ufl.dot((k_b_current/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(qc))*dx 
	MassCell  +=  (k_c/mu_c)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(qc))*dx 
	MassCell  += - (k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(qc))*dx 
	MassCell  += S_c* 1/dt * ufl.nabla_div(du)*qc*dx
	# Mass conservation of the IF, Eq. 156
	MassIF     = ufl.dot((k_b_current/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(ql))*dx 
	MassIF    += (k_c/mu_c + k_l/mu_l)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx 
	MassIF    += - (k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(ql))*dx 
	MassIF    += 1/dt * ufl.nabla_div(du)*ql*dx
	# Mass conservation of the blood, Eq. 157
	MassBlood  = coeff_CEP * dpl/dt * qb* dx 
	MassBlood += - coeff_CEP * constitutive_laws.Cstate(plc_n+dplc,a) * dplc/dt * qb * dx 
	MassBlood += - coeff_CEP * dpb/dt * qb * dx  
	MassBlood +=  ufl.dot((k_b_current/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(qb))*dx 
	MassBlood +=  epsb * 1/dt * ufl.nabla_div(du)*qb*dx
	# Momentum conservation, Eq. 158
	Momentum   = ufl.inner(ufl.sym(ufl.grad(w)),constitutive_laws.Elastic_constitutive_law(u_n+du,Young,Poisson))*dx 
	Momentum  += - (1-epsb)*(pl_n+dpl)*ufl.nabla_div(w)*dx 
	Momentum  += (1-epsb)*S_c*(plc_n+dplc)*ufl.nabla_div(w)*dx 
	Momentum  += - epsb * (pb_n+dpb) * ufl.nabla_div(w)*dx
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
	MassCell  += S_c*ufl.dot((k_b_current/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(qc))*dx 
	# 
	MassCell  += (k_c/mu_c)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(qc))*dx 
	# 
	MassCell  += -(k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(qc))*dx 
	# 
	MassCell  += S_c* 1/dt * ufl.nabla_div(du)*qc*dx
	# 
	# Mass conservation of the IF, Eq. 156
	MassIF     =  ufl.dot((k_b_current/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(ql))*dx 
	# 
	MassIF    += (k_c/mu_c + k_l/mu_l)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx 
	# 
	MassIF    += -(k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(ql))*dx 
	# 
	MassIF    += 1/dt * ufl.nabla_div(du)*ql*dx
	# Mass conservation of the blood, Eq. 157
	MassBlood  = coeff_CEP * dpl/dt * qb* dx 
	MassBlood += - coeff_CEP * constitutive_laws.Cstate(plc_n+dplc,a) * dplc/dt * qb * dx 
	MassBlood += - coeff_CEP * dpb/dt * qb * dx  
	# 
	MassBlood +=  ufl.dot((k_b_current/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(qb))*dx 
	# 
	MassBlood +=  epsb * 1/dt * ufl.nabla_div(du)*qb*dx
	# Momentum conservation, Eq. 158
	Momentum   = ufl.inner(ufl.sym(ufl.grad(w)),constitutive_laws.Elastic_constitutive_law_DG0(u_n+du,Young,Poisson))*dx 
	# 
	Momentum  += - (1-epsb)*(pl_n+dpl)*ufl.nabla_div(w)*dx 
	# 
	Momentum  += (1-epsb)*S_c*(plc_n+dplc)*ufl.nabla_div(w)*dx 
	# 
	Momentum  += - epsb * (pb_n+dpb) * ufl.nabla_div(w)*dx
	return MassCell + MassIF + MassBlood + Momentum
# 
#____________________________________________________________#
# 				Two Compartments Biphasic + O2
#____________________________________________________________#
# # 
# # Correct functions, variables and create main
# # def variational_form_oxygen_constant_cell(plc_n, pl_n, pb_n, u_n, wo2_n, dplc, dpl, dpb, du, dwo2, qc, ql, qb, w, qo2, poro_n, dt, Young, Poisson, k_l, mu_l, rho_l, k_c, mu_c, k_b, mu_b, wo2_b, a, K, alpha, gamma_0, D0_o2, delta, dx, ds):
# #
# def variational_form_oxygen_constant_cell(previous_solution, solution, qc, ql, qb, w, qo2, poro_n, dt, Young, Poisson, k_l, mu_l, rho_l, k_c, mu_c, k_b, mu_b, wo2_b, a, K, alpha, gamma_0, D0_o2, delta, dx, ds):
# 	"""
# 	This function provides the incremental weak form according to Eq. XX-XX.
# 	Inputs:
# 		- The functions at the previous time step:  plc_n, pl_n, pb_n, u_n, wo2_n	(dolfinx functions)
# 		- The functions increments to be computed:  dplc, dpl, dpb, du, dwo2		(dolfinx functions)
# 		- The test functions: 					    qc, ql, qb, w, qo2				(dolfinx Testfunctions)
# 		- The porosity at the previous time step:   poro_n 						(dolfinx functions)
# 		- The time increment: 					    dt 							(dolfinx scalartype or dolfinx constant)
# 		- The solid parameters: 					Young, Poisson 				(dolfinx scalartype or dolfinx constant or dolfinx function)
# 		- The IF parameters: 						k_l, mu_l, rho_l			(dolfinx scalartype or dolfinx constant or dolfinx function)
# 		- The cell parameters: 						k_c, mu_c 					(dolfinx scalartype or dolfinx constant or dolfinx function)
# 		- The varscular parameters:					k_b, mu_b, wo2_b			(dolfinx scalartype or dolfinx constant or dolfinx function)
# 		- The parameter for the saturation: 		a 							(dolfinx scalartype or dolfinx constant or dolfinx function)
# 		- The vessel compressibility: 				K 							(dolfinx scalartype or dolfinx constant or dolfinx function)
# 		- The exponent from Eq. 63:					alpha						(dolfinx scalartype)
# 		- The cell metabolism:						gamma_0						(dolfinx scalartype)
# 		- D0_o2: the diffusion coefficient of oxygen in the bulk interstitial fluid
# 		- delta: coefficient related to the tortuosity of the medium
# 		- dx, ds resulting from Measure and cell_tag/facet_tag mapping
# 	Output: 
# 		- The weak form without BCs (Robin or Neumann)
# 	"""
# 	import ufl
#	# split the solution into: displacement increment, IF pressure increment, Blood Pressure increment
#	dpl, dplc, dpb, du, dwo2      = ufl.split(solution)
#	# Previous time steps solutions
#	pl_n, plc_n, pb_n, u_n, wo2_n = ufl.split(previous_solution)
# 	import ufl 
#	REPARTIR DE CE QUI EST AU DESSUS LAUTRE FORMULATION
# 	KK = vascular_permeability(epsb0,pl_n+dpl,plc_n+dplc,pb_n+dpb,a,K,k_b,ratio,direction,dimension,alpha)
# 	# Mass conservation of the cells, Eq. 109
# 	MassCell   = Cmc(du,poro_n,porosity_2_comp(du,poro_n,KK,mu_b,pb_n+dpb,dt),plc_n+dplc,a) * dplc/dt * qc * dx 
# 	MassCell  += saturation_cell(plc_n+dplc,a)*ufl.dot((KK/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(qc))*dx 
# 	MassCell  +=  (k_c/mu_c)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(qc))*dx 
# 	MassCell  += - (k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(qc))*dx 
# 	MassCell  += saturation_cell(plc_n+dplc,a)* 1/dt * ufl.nabla_div(du)*qc*dx
# 	# Mass conservation of the IF, Eq. 110
# 	MassIF     = ufl.dot((KK/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(ql))*dx 
# 	MassIF    += (k_c/mu_c + k_l/mu_l)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx 
# 	MassIF    += - (k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(ql))*dx 
# 	MassIF    += 1/dt * ufl.nabla_div(du)*ql*dx
# 	# Mass conservation of the blood, Eq. 111
# 	MassBlood  = Cep(epsb0,K,plc_n+dplc,pl_n+dpl,pb_n+dpb,a)* dpl/dt * qb* dx 
# 	MassBlood += - Cep(epsb0,K,plc_n+dplc,pl_n+dpl,pb_n+dpb,a) * Cstate(plc_n+dplc,a) * dplc/dt * qb* dx 
# 	MassBlood += - Cep(epsb0,K,plc_n+dplc,pl_n+dpl,pb_n+dpb,a) * dpb/dt * qb* dx  
# 	MassBlood +=  ufl.dot((KK/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(qb))*dx 
# 	MassBlood +=  vascular_porosity(epsb0,pl_n+dpl,plc_n+dplc,pb_n+dpb,a,K)*1/dt * ufl.nabla_div(du)*qb*dx
# 	# Momentum conservation, Eq. 112
# 	Momentum   = ufl.inner(ufl.grad(w),Elastic_constitutive_law(u_n+du,Young,Poisson))*dx 
# 	Momentum  += - (1-vascular_porosity(epsb0,pl_n+dpl,plc_n+dplc,pb_n+dpb,a,K))*(pl_n+dpl)*ufl.nabla_div(w)*dx 
# 	Momentum  += (1-vascular_porosity(epsb0,pl_n+dpl,plc_n+dplc,pb_n+dpb,a,K))*saturation_cell(plc_n+dplc,a)*(plc_n+dplc)*ufl.nabla_div(w)*dx 
# 	Momentum  += - vascular_porosity(epsb0,pl_n+dpl,plc_n+dplc,pb_n+dpb,a,K) * (pb_n+dpb) * ufl.nabla_div(w)*dx
# 	# Mass conservation of the O2 Eq. 113
# 	Deff       = Deff_O2(D0_o2,delta,pl_n+dpl,plc_n+dpc,pb_n+dpb,a,K,k_b,ratio,direction,dimension,alpha,du,poro_n,mu_b,dt)
# 	MassO2     = (1-saturation_cell(plc_n+dplc,a)) * porosity_2_comp(du,poro_n,KK,mu_b,pb_n+dpb,dt) * dwo2/dt * qo2 * dx
# 	MassO2	  += -k_l/mu_l*ufl.dot(ufl.grad(wo2_n+dwo2),ufl.grad(pl_n+dpl))*qo2*dx
# 	MassO2	  += (1-saturation_cell(plc_n+dplc,a)) * porosity_2_comp(du,poro_n,KK,mu_b,pb_n+dpb,dt) * Deff_O2(D0_o2,delta,pl,plc,pb,a,K,k_b,ratio,direction,dimension,alpha,du,poro_n,mu_b,dt) * ufl.dot(ufl.grad(wo2_n+dwo2),ufl.grad(qo2))*dx
# 	MassO2	  += -1/rho_l*(MO2_b_l(epsb0,pl_n+dpl,plc_n+dpcl,pb_n+dpb,a,K,wo2_b,wo2_n+dwo2) - MO2_l_c_constant(wo2,gamma_0,epsb0,pl,plc,pb,du,a,K,poro_n,k_b,mu_b,ratio,direction,dimension,alpha,dt))*qo2*dx
# 	return MassCell + MassIF + MassBlood + Momentum + MassO2
# # 
# def variational_form_oxygen_conditional_cell(plc_n, pl_n, pb_n, u_n, wo2_n, dplc, dpl, dpb, du, dwo2, qc, ql, qb, w, qo2, poro_n, dt, Young, Poisson, k_l, mu_l, rho_l, k_c, mu_c, k_b, mu_b, wo2_b, a, K, alpha, wcrit, gamma_0, D0_o2, delta, dx, ds):
# 	"""
# 	This function provides the incremental weak form according to Eq. 46-48.
# 	Inputs:
# 		- The functions at the previous time step:  plc_n, pl_n, pb_n, u_n, wo2_n	(dolfinx functions)
# 		- The functions increments to be computed:  dplc, dpl, dpb, du, dwo2		(dolfinx functions)
# 		- The test functions: 					    qc, ql, qb, w, qo2				(dolfinx Testfunctions)
# 		- The porosity at the previous time step:   poro_n 						(dolfinx functions)
# 		- The time increment: 					    dt 							(dolfinx scalartype or dolfinx constant)
# 		- The solid parameters: 					Young, Poisson 				(dolfinx scalartype or dolfinx constant or dolfinx function)
# 		- The IF parameters: 						k_l, mu_l, rho_l			(dolfinx scalartype or dolfinx constant or dolfinx function)
# 		- The cell parameters: 						k_c, mu_c 					(dolfinx scalartype or dolfinx constant or dolfinx function)
# 		- The varscular parameters:					k_b, mu_b, wo2_b			(dolfinx scalartype or dolfinx constant or dolfinx function)
# 		- The parameter for the saturation: 		a 							(dolfinx scalartype or dolfinx constant or dolfinx function)
# 		- The vessel compressibility: 				K 							(dolfinx scalartype or dolfinx constant or dolfinx function)
# 		- The exponent from Eq. 63:					alpha						(dolfinx scalartype)
# 		- The hypoxia threshold:					wcrit						(dolfinx scalartype)
# 		- The cell metabolism:						gamma_0						(dolfinx scalartype)
# 		- D0_o2: the diffusion coefficient of oxygen in the bulk interstitial fluid
# 		- delta: coefficient related to the tortuosity of the medium
# 		- dx, ds resulting from Measure and cell_tag/facet_tag mapping
# 	Output: 
# 		- The weak form without BCs (Robin or Neumann)
# 	"""
# 	import ufl
# 	KK = vascular_permeability(epsb0,pl_n+dpl,plc_n+dplc,pb_n+dpb,a,K,k_b,ratio,direction,dimension,alpha)
# 	# Mass conservation of the cells, Eq. 109
# 	MassCell   = Cmc(du,poro_n,porosity_2_comp(du,poro_n,KK,mu_b,pb_n+dpb,dt),plc_n+dplc,a) * dplc/dt * qc * dx 
# 	MassCell  += saturation_cell(plc_n+dplc,a)*ufl.dot((KK/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(qc))*dx 
# 	MassCell  +=  (k_c/mu_c)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(qc))*dx 
# 	MassCell  += - (k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(qc))*dx 
# 	MassCell  += saturation_cell(plc_n+dplc,a)* 1/dt * ufl.nabla_div(du)*qc*dx
# 	# Mass conservation of the IF, Eq. 110
# 	MassIF     = ufl.dot((KK/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(ql))*dx 
# 	MassIF    += (k_c/mu_c + k_l/mu_l)*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx 
# 	MassIF    += - (k_c/mu_c)*ufl.dot(ufl.grad(plc_n+dplc),ufl.grad(ql))*dx 
# 	MassIF    += 1/dt * ufl.nabla_div(du)*ql*dx
# 	# Mass conservation of the blood, Eq. 111
# 	MassBlood  = Cep(epsb0,K,plc_n+dplc,pl_n+dpl,pb_n+dpb,a)* dpl/dt * qb* dx 
# 	MassBlood += - Cep(epsb0,K,plc_n+dplc,pl_n+dpl,pb_n+dpb,a) * Cstate(plc_n+dplc,a) * dplc/dt * qb* dx 
# 	MassBlood += - Cep(epsb0,K,plc_n+dplc,pl_n+dpl,pb_n+dpb,a) * dpb/dt * qb* dx  
# 	MassBlood +=  ufl.dot((KK/mu_b)*ufl.grad(pb_n+dpb),ufl.grad(qb))*dx 
# 	MassBlood +=  vascular_porosity(epsb0,pl_n+dpl,plc_n+dplc,pb_n+dpb,a,K)*1/dt * ufl.nabla_div(du)*qb*dx
# 	# Momentum conservation, Eq. 112
# 	Momentum   = ufl.inner(ufl.grad(w),Elastic_constitutive_law(u_n+du,Young,Poisson))*dx 
# 	Momentum  += - (1-vascular_porosity(epsb0,pl_n+dpl,plc_n+dplc,pb_n+dpb,a,K))*(pl_n+dpl)*ufl.nabla_div(w)*dx 
# 	Momentum  += (1-vascular_porosity(epsb0,pl_n+dpl,plc_n+dplc,pb_n+dpb,a,K))*saturation_cell(plc_n+dplc,a)*(plc_n+dplc)*ufl.nabla_div(w)*dx 
# 	Momentum  += - vascular_porosity(epsb0,pl_n+dpl,plc_n+dplc,pb_n+dpb,a,K) * (pb_n+dpb) * ufl.nabla_div(w)*dx
# 	# Mass conservation of the O2 Eq. 113
# 	Deff       = Deff_O2(D0_o2,delta,pl_n+dpl,plc_n+dpc,pb_n+dpb,a,K,k_b,ratio,direction,dimension,alpha,du,poro_n,mu_b,dt)
# 	MassO2     = (1-saturation_cell(plc_n+dplc,a)) * porosity_2_comp(du,poro_n,KK,mu_b,pb_n+dpb,dt) * dwo2/dt * qo2 * dx
# 	MassO2	  += -k_l/mu_l*ufl.dot(ufl.grad(wo2_n+dwo2),ufl.grad(pl_n+dpl))*qo2*dx
# 	MassO2	  += (1-saturation_cell(plc_n+dplc,a)) * porosity_2_comp(du,poro_n,KK,mu_b,pb_n+dpb,dt) * Deff_O2(D0_o2,delta,pl,plc,pb,a,K,k_b,ratio,direction,dimension,alpha,du,poro_n,mu_b,dt) * ufl.dot(ufl.grad(wo2_n+dwo2),ufl.grad(qo2))*dx
# 	MassO2	  += -1/rho_l*(MO2_b_l(epsb0,pl_n+dpl,plc_n+dpcl,pb_n+dpb,a,K,wo2_b,wo2_n+dwo2) - MO2_l_c(wcrit,wo2,gamma_0,epsb0,pl,plc,pb,du,a,K,poro_n,k_b,mu_b,ratio,direction,dimension,alpha,dt))*qo2*dx
# 	return MassCell + MassIF + MassBlood + Momentum + MassO2
# # 
#------------------------------------------------------------#
#              Debug Variational forms                       #
#------------------------------------------------------------#
# 
#____________________________________________________________#
# 				Two Compartment Monophasic debug forms
#____________________________________________________________#
# 
def variational_form_2_mono_updated_anisotropic_k_b_constant(previous_solution, solution, ql, qb, w, Young, Poisson, k_l, mu_l, k_b, ratio, mu_b, epsb0, K, alpha, dx, ds, dt, material_law,varepsilon_b_law):
	"""
	This function provides the incremental weak form according to Eq. 60, 62 and 62. The vascular permeability is here anisotropic and constant to assess if it is the implementation with (epsb/ebsb0)**alpha that creates weird result.

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
	if varepsilon_b_law == 'linear':
		# Mass conservation of the IF, Eq. 60
		MassIF            = 1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_constant(k_b, ratio)*ufl.grad(pb_n+dpb),ufl.grad(ql))*dx
		MassIF           += k_l/mu_l*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
		MassIF           += 1/dt * ufl.nabla_div(du)*ql*dx
		# Mass conservation of the IF, Eq. 61
		MassBlood         = constitutive_laws.Cep_linear_mono(epsb0,K)*dpl/dt*qb*dx
		MassBlood        += -constitutive_laws.Cep_linear_mono(epsb0,K)*dpb/dt*qb*dx
		MassBlood        += 1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_constant(k_b, ratio)*ufl.grad(pb_n+dpb),ufl.grad(qb))*dx
		MassBlood        += 1/dt* constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K) * ufl.nabla_div(du)*qb*dx
		# Momentum conservation, Eq. 62
		Momentum          = ufl.inner(ufl.grad(w),material_law(u_n+du,Young,Poisson))*dx 
		Momentum         += - (1-constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K) )*(pl_n+dpl)*ufl.nabla_div(w)*dx 
		Momentum         += - constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K)*(pb_n+dpb)*ufl.nabla_div(w)*dx 
	elif varepsilon_b_law == 'atan':
		# Mass conservation of the IF, Eq. 60
		MassIF            = 1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_constant(k_b, ratio)*ufl.grad(pb_n+dpb),ufl.grad(ql))*dx
		MassIF           += k_l/mu_l*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
		MassIF           += 1/dt * ufl.nabla_div(du)*ql*dx
		# Mass conservation of the IF, Eq. 61
		MassBlood         = constitutive_laws.Cep_atan_mono(epsb0,K,pl_n+dpl,pb_n+dpb)*dpl/dt*qb*dx
		MassBlood        += -constitutive_laws.Cep_atan_mono(epsb0,K,pl_n+dpl,pb_n+dpb)*dpb/dt*qb*dx
		MassBlood        += 1/mu_b*ufl.dot(constitutive_laws.matrix_vascular_permeability_constant(k_b, ratio)*ufl.grad(pb_n+dpb),ufl.grad(qb))*dx
		MassBlood        += 1/dt* constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K) * ufl.nabla_div(du)*qb*dx
		# Momentum conservation, Eq. 62
		Momentum          = ufl.inner(ufl.grad(w),material_law(u_n+du,Young,Poisson))*dx 
		Momentum         += - (1-constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K) )*(pl_n+dpl)*ufl.nabla_div(w)*dx 
		Momentum         += - constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K)*(pb_n+dpb)*ufl.nabla_div(w)*dx 
	else:
		print("ERROR in the choice of the evolution type. Choose 'linear' or 'atan'.")
	return MassIF + MassBlood + Momentum
# 
def variational_form_2_mono_updated_anisotropic_k_b_product_outside_dot(previous_solution, solution, ql, qb, w, Young, Poisson, k_l, mu_l, k_b, ratio, mu_b, epsb0, K, alpha, dx, ds, dt, material_law,varepsilon_b_law):
	"""
	This function provides the incremental weak form according to Eq. 60, 62 and 62. The vascular permeability is here anisotropic and constant to assess if it is the implementation with (epsb/ebsb0)**alpha that creates weird result.

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
	if varepsilon_b_law == 'linear':
		# Mass conservation of the IF, Eq. 60
		MassIF            = 1/mu_b*((constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K)/epsb0)**alpha) * ufl.dot(constitutive_laws.matrix_vascular_permeability_constant(k_b, ratio)*ufl.grad(pb_n+dpb),ufl.grad(ql))*dx
		MassIF           += k_l/mu_l*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
		MassIF           += 1/dt * ufl.nabla_div(du)*ql*dx
		# Mass conservation of the IF, Eq. 61
		MassBlood         = constitutive_laws.Cep_linear_mono(epsb0,K)*dpl/dt*qb*dx
		MassBlood        += -constitutive_laws.Cep_linear_mono(epsb0,K)*dpb/dt*qb*dx
		MassBlood        += 1/mu_b*((constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K)/epsb0)**alpha) *ufl.dot(constitutive_laws.matrix_vascular_permeability_constant(k_b, ratio)*ufl.grad(pb_n+dpb),ufl.grad(qb))*dx
		MassBlood        += 1/dt* constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K) * ufl.nabla_div(du)*qb*dx
		# Momentum conservation, Eq. 62
		Momentum          = ufl.inner(ufl.grad(w),material_law(u_n+du,Young,Poisson))*dx 
		Momentum         += - (1-constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K) )*(pl_n+dpl)*ufl.nabla_div(w)*dx 
		Momentum         += - constitutive_laws.vascular_porosity_linear_mono(epsb0,pl_n+dpl,pb_n+dpb,K)*(pb_n+dpb)*ufl.nabla_div(w)*dx 
	elif varepsilon_b_law == 'atan':
		# Mass conservation of the IF, Eq. 60
		MassIF            = 1/mu_b*((constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K)/epsb0)**alpha) *ufl.dot(constitutive_laws.matrix_vascular_permeability_constant(k_b, ratio)*ufl.grad(pb_n+dpb),ufl.grad(ql))*dx
		MassIF           += k_l/mu_l*ufl.dot(ufl.grad(pl_n+dpl),ufl.grad(ql))*dx
		MassIF           += 1/dt * ufl.nabla_div(du)*ql*dx
		# Mass conservation of the IF, Eq. 61
		MassBlood         = constitutive_laws.Cep_atan_mono(epsb0,K,pl_n+dpl,pb_n+dpb)*dpl/dt*qb*dx
		MassBlood        += -constitutive_laws.Cep_atan_mono(epsb0,K,pl_n+dpl,pb_n+dpb)*dpb/dt*qb*dx
		MassBlood        += 1/mu_b*((constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K)/epsb0)**alpha) *ufl.dot(constitutive_laws.matrix_vascular_permeability_constant(k_b, ratio)*ufl.grad(pb_n+dpb),ufl.grad(qb))*dx
		MassBlood        += 1/dt* constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K) * ufl.nabla_div(du)*qb*dx
		# Momentum conservation, Eq. 62
		Momentum          = ufl.inner(ufl.grad(w),material_law(u_n+du,Young,Poisson))*dx 
		Momentum         += - (1-constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K) )*(pl_n+dpl)*ufl.nabla_div(w)*dx 
		Momentum         += - constitutive_laws.vascular_porosity_atan_mono(epsb0,pl_n+dpl,pb_n+dpb,K)*(pb_n+dpb)*ufl.nabla_div(w)*dx 
	else:
		print("ERROR in the choice of the evolution type. Choose 'linear' or 'atan'.")
	return MassIF + MassBlood + Momentum
# 
#____________________________________________________________#
# 					 End of functions
#____________________________________________________________#
# 
# 
if __name__ == "__main__":
    print("Loading of the constitutive laws successfully completed.")
    # EoF