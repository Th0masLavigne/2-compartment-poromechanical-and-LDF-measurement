# This document holds the functions used as part of 
# the poro-mechanical modelling of living tissues 
# according to the theory in ./Theoretical_Developments.pdf
# 
# Author: Thomas Lavigne
# Date: 18/09/2024
# run first : " python3 -m pip install . " to create the package at ../ location
# 
#------------------------------------------------------------#
#                    Constitutive Laws                       #
#------------------------------------------------------------#
# 
#____________________________________________________________#
# 					 Solid Scaffold
#____________________________________________________________#
# 
def Lame_coefficients(Young,Poisson):
	"""
	Computes the Lamé coefficients from the Young's modulus and Poisson's ratio.

	Parameters:
	-----------
	Young : dolfinx.default_scalar_type or dolfinx.fem.Constant 
	    The Young's modulus, which represents the material's stiffness.

	Poisson : dolfinx.default_scalar_type or dolfinx.fem.Constant
	    The Poisson's ratio, which characterizes the material's compressibility.

	Returns:
	--------
	lambda : dolfinx.default_scalar_type
	    The first Lamé parameter, λ, which relates to compressive strength.

	modulus : dolfinx.default_scalar_type
	    The shear modulus, μ, which characterizes the material's resistance to shear deformations.

	Notes:
	------
	- The Lamé coefficients are commonly used in linear elasticity to describe the relationship 
	  between stress and strain in isotropic materials.
	- λ and μ are derived using the relations:
	    λ = Young * Poisson / ((1 + Poisson) * (1 - 2 * Poisson))
	    μ = Young / (2 * (1 + Poisson))
	"""
	import dolfinx
	return dolfinx.default_scalar_type(Young*Poisson/((1+Poisson)*(1-2*Poisson))), dolfinx.default_scalar_type(Young/(2*(1+Poisson)))
# 
# 
def Lame_coefficients_functions(Young,Poisson):
	"""
	Computes the Lamé coefficients from the Young's modulus and Poisson's ratio.

	Parameters:
	-----------
	Young : dolfinx.fem.Function
	    The Young's modulus, which represents the material's stiffness.

	Poisson : dolfinx.fem.Function
	    The Poisson's ratio, which characterizes the material's compressibility.

	Returns:
	--------
	lambda : dolfinx.default_scalar_type
	    The first Lamé parameter, λ, which relates to compressive strength.

	modulus : dolfinx.default_scalar_type
	    The shear modulus, μ, which characterizes the material's resistance to shear deformations.

	Notes:
	------
	- The Lamé coefficients are commonly used in linear elasticity to describe the relationship 
	  between stress and strain in isotropic materials.
	- λ and μ are derived using the relations:
	    λ = Young * Poisson / ((1 + Poisson) * (1 - 2 * Poisson))
	    μ = Young / (2 * (1 + Poisson))
	"""
	return Young*Poisson/((1+Poisson)*(1-2*Poisson)), Young/(2*(1+Poisson))
# 
def Local_deformation(u):
	"""
	Computes the local deformation (strain) tensor for a linear elastic model.

	Parameters:
	-----------
	u : dolfinx.fem.Function
	    The displacement field, representing the displacements of material points in the domain.

	Returns:
	--------
	deformation_tensor : ufl object
	    The strain tensor,  sym(grad(u)), which measures the local deformation of the material.

	Notes:
	------
	- The deformation tensor is symmetric and describes how the material stretches or compresses locally.
	- In a linear elastic model, the deformation tensor is directly related to the stress tensor.
	"""
	import ufl
	return ufl.sym(ufl.grad(u))
# 
def Elastic_constitutive_law(u,Young,Poisson):
	"""
	Computes the stress tensor for a linear elastic material using the constitutive law (Hooke's Law).
	This law corresponds to t^eff Eq. 12.

	Parameters:
	-----------
	u : dolfinx.fem.Function
	    The displacement field, representing the displacements of material points in the domain.

	Young : dolfinx.default_scalar_type or dolfinx.fem.Constant or dolfinx.fem.Function
	    The Young's modulus, characterizing the material's stiffness.

	Poisson : dolfinx.default_scalar_type or dolfinx.fem.Constant or dolfinx.fem.Function
	    The Poisson's ratio, characterizing the material's compressibility.

	Returns:
	--------
	stress_tensor : ufl object
	    The stress tensor, σ(u), representing internal forces per unit area in the material.

	Notes:
	------
	- The stress tensor is computed using Hooke's law for isotropic materials: 
	    σ = 2 * μ * sym(∇u) + λ * (∇·u) * I
	  where sym(∇u) is the strain tensor, λ and μ are the Lamé coefficients, and I is the identity matrix.
	- This function combines the displacement field with material properties to compute the internal stresses.
	"""
	import ufl
	lmbda, mu = Lame_coefficients(Young,Poisson)
	return 2 * mu * Local_deformation(u) + lmbda * ufl.nabla_div(u) * ufl.Identity(len(u))
# 
def Elastic_constitutive_law_DG0(u,Young,Poisson):
	"""
	Computes the stress tensor for a linear elastic material using the constitutive law (Hooke's Law).
	This law corresponds to t^eff Eq. 12.

	Parameters:
	-----------
	u : dolfinx.fem.Function
	    The displacement field, representing the displacements of material points in the domain.

	Young : dolfinx.fem.Function
	    The Young's modulus, characterizing the material's stiffness.

	Poisson : dolfinx.fem.Function
	    The Poisson's ratio, characterizing the material's compressibility.

	Returns:
	--------
	stress_tensor : ufl object
	    The stress tensor, σ(u), representing internal forces per unit area in the material.

	Notes:
	------
	- The stress tensor is computed using Hooke's law for isotropic materials: 
	    σ = 2 * μ * sym(∇u) + λ * (∇·u) * I
	  where sym(∇u) is the strain tensor, λ and μ are the Lamé coefficients, and I is the identity matrix.
	- This function combines the displacement field with material properties to compute the internal stresses.
	"""
	import ufl
	lmbda, mu = Lame_coefficients_functions(Young,Poisson)
	return 2 * mu * Local_deformation(u) + lmbda * ufl.nabla_div(u) * ufl.Identity(len(u))
# 
def Neo_Hooke(u,Young,Poisson):
	"""
	Computes the Cauchy stress tensor from a compressible neo-Hookean formulation.

	The strain energy function is given by:
	    W = (μ / 2) * (Ic - tr(I)) - μ * ln(J) + (λ_m / 2) * (ln(J))^2
	where Ic is the first invariant of the right Cauchy-Green tensor, J is the determinant of the deformation gradient, and λ_m, μ are the Lamé coefficients.

	Parameters:
	-----------
	Young : dolfinx.fem.Function
	    The Young's modulus, characterizing the material's stiffness.

	Poisson : dolfinx.fem.Function
	    The Poisson's ratio, characterizing the material's compressibility.

	u : dolfinx.fem.Function
	    The displacement field, representing the deformation of material points.

	Returns:
	--------
	cauchy_stress : ufl object
	    The Cauchy stress tensor, computed using Nanson’s formula and the derivative of W with respect to F.
	    
	Notes:
	------
	- The First Piola-Kirchhoff stress tensor is obtained from the derivative of W with respect to the deformation gradient F.
	- The Cauchy stress tensor follows from the Piola transformation using Nanson's formula.
	"""
	import ufl 
	lmbda, mu = Lame_coefficients_functions(Young,Poisson)
	## Deformation gradient
	F  = ufl.variable(ufl.Identity(3) + ufl.grad(u))
	J  = ufl.variable(ufl.det(F))
	## Right Cauchy-Green tensor
	C  = ufl.variable(F.T * F)
	##Invariants of deformation tensors
	Ic = ufl.variable(ufl.tr(C))
	## Strain energy density function
	W  = (mu / 2) * (Ic - ufl.tr(ufl.Identity(3))) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J))**2
	return (1/J)*ufl.diff(W, F)*F.T
# 
#____________________________________________________________#
# 					 Extra-vascular fluids
#____________________________________________________________#
# 
# 
#...........................................................#
# 					 Single comportment
#...........................................................#
# 
def porosity_1(du,poro_n,poro_min,poro_max):
	"""
	Computes the updated porosity based on the displacement increment field and the previous porosity value, according to Eq. 16.

	Parameters:
	-----------
	du : dolfinx.fem.Function
	    The displacement increment field, representing the change in displacement over the current time step.

	poro_n : dolfinx.fem.Function
	    The porosity at the previous time step, representing the fraction of fluid space in the material.

	poro_min/poro_max : dolfinx.default_scalartype
	    The maximum and minimum of physically relevant porosity to ensure avoiding non-physical values (for instance negative). 
	    This might suggest a wrong choice in spatial and/or temporal discretizing if the bounds are reached.

	Returns:
	--------
	porosity : ufl object
	    The updated porosity, ε(t), representing the current void fraction in the material after the deformation.

	Notes:
	------
	- The updated porosity is computed using the formula:
	    ε = (poro_n + ∇·du) / (1 + ∇·du)
	  where `poro_n` is the previous porosity, and `∇·du` is the volumetric strain.
	- This function is useful in modeling flow through porous media where the porosity evolves due to deformation.
	"""
	import ufl
	import dolfinx
	return ufl.min_value(ufl.max_value( (poro_n + ufl.div(du))/(1+ufl.div(du)), dolfinx.default_scalar_type(poro_min) ), dolfinx.default_scalar_type(poro_max) )
# 
#...........................................................#
# 					 Two-compartment
#...........................................................#
def porosity_2(du,poro_n,poro_b,poro_b_n,poro_min,poro_max):
	"""
	Computes the updated porosity based on displacement, vascular porosity, and vascular parameters as described in Eq. 147.

	Parameters:
	-----------
	du : dolfinx.function or ufl expression
	    The displacement increment field representing the changes in displacement over time.

	poro_n : dolfinx.function or ufl expression
	    The porosity at the previous time step.

	poro_b : dolfinx.function or ufl expression
	    The current vascular porosity.

	poro_b_n : dolfinx.function or ufl expression
	    The previous vascular porosity.

	poro_min/poro_max : dolfinx.default_scalartype
	    The maximum and minimum of physically relevant porosity to ensure avoiding non-physical values (for instance negative). 
	    This might suggest a wrong choice in spatial and/or temporal discretizing if the bounds are reached.

	Returns:
	--------
	updated_porosity : ufl expression
	    The expression representing the updated porosity at the current time step, accounting for mechanical deformation, fluid flow, 
	    and the time evolution of the system.
	"""
	import ufl
	import dolfinx
	return ufl.min_value( ufl.max_value( ( poro_n + poro_b_n - poro_b * ( 1 + ufl.div(du) ) + ufl.div(du) )/( 1+ufl.div(du) ), dolfinx.default_scalar_type(poro_min) ) , dolfinx.default_scalar_type(poro_max) )
# 
# 
def porosity_2_pressure(du,poro_n,k_b,mu_b,pb,dt,poro_min,poro_max):
	"""
	Computes the updated porosity based on displacement, pressure gradients, and vascular parameters as described in Eq. 54.

	Parameters:
	-----------
	du : dolfinx.function or ufl expression
	    The displacement increment field representing the changes in displacement over time.

	poro_n : dolfinx.function or ufl expression
	    The porosity at the previous time step.

	k_b : dolfinx.function or ufl expression
	    The current vascular permeability, which influences the fluid flow through the tissue.

	mu_b : dolfinx.function or ufl expression
	    The viscosity of the blood, affecting the resistance to fluid flow.

	pb : dolfinx.function or ufl expression
	    The blood pressure, driving the fluid flow within the vascular system.

	dt : float
	    The time step (in seconds) for updating the porosity.

	poro_min/poro_max : dolfinx.default_scalartype
	    The maximum and minimum of physically relevant porosity to ensure avoiding non-physical values (for instance negative). 
	    This might suggest a wrong choice in spatial and/or temporal discretizing if the bounds are reached.

	Returns:
	--------
	updated_porosity : ufl expression
	    The expression representing the updated porosity at the current time step, accounting for mechanical deformation, fluid flow, 
	    and the time evolution of the system.

	Notes:
	------
	- The function updates the porosity by accounting for both mechanical deformation (through the displacement field `du`) 
	  and fluid flow (through the gradient of blood pressure `pb`).
	- The term `(k_b / mu_b) * grad(pb)` represents the Darcy velocity, which accounts for the flow of blood through the porous 
	  medium due to pressure gradients.
	- The final expression divides the sum of the previous porosity, mechanical deformation, and fluid effects by `(1 + div(du))` 
	  to ensure the updated porosity respects volume conservation.

	Formula:
	--------
	- The updated porosity is computed as:
	  updated_porosity = (poro_n + div(du) - div((k_b / mu_b) * grad(pb)) * dt) / (1 + div(du))
	"""
	import ufl
	import dolfinx
	return ufl.min_value( ufl.max_value( ( poro_n + ufl.div(du) - ufl.div( (k_b/mu_b)*ufl.grad(pb) ) * dt )/(1+ufl.div(du)), dolfinx.default_scalar_type(poro_min) ) , dolfinx.default_scalar_type(poro_max) )
# 
#...........................................................#
# 					 Bi-phasic specific laws
#...........................................................#
# 
def saturation_cell(plc,a):
	"""
	Computes the saturation in cells of the extra-vascular space based on capillary pressure (Eq. 70).

	Parameters:
	-----------
	plc : ufl expression or dolfinx.function
	    The capillary pressure within the extra-vascular space.

	a : float
	    A constant parameter that reflects the characteristics of the connective tissue microstructure. It modulates the relationship between capillary pressure and saturation.

	Returns:
	--------
	Sc : ufl expression
	    The saturation in cells of the extra-vascular space, ranging from 0 (fully unsaturated) to 1 (fully saturated).

	Notes:
	------
	- The function implements a smooth approximation of the saturation based on the arctangent function. As `plc` increases, saturation `Sc` approaches 1, indicating a higher degree of fluid saturation in the tissue.
	- The parameter `a` controls how quickly saturation changes with pressure. A smaller `a` results in a steeper transition from unsaturated to saturated states.
	"""
	import numpy
	import ufl
	return 1 - ( 2/numpy.pi * ufl.atan( plc / a ) )
# 
def derivative_Sc_plc(plc,a):
	"""
	Computes the derivative of the saturation in cells of the extra-vascular space with respect to capillary pressure (Eq. 105).

	Parameters:
	-----------
	plc : float or ufl expression
	    The capillary pressure within the extra-vascular space.

	a : float
	    A constant parameter that characterizes the connective tissue microstructure, influencing how saturation changes with pressure.

	Returns:
	--------
	dScdplc : float or ufl expression
	    The derivative of the saturation with respect to capillary pressure, indicating how sensitive the saturation is to changes in capillary pressure.

	Notes:
	------
	- The derivative is calculated based on the formula:
	  \\[
	  \\frac{dSc}{dplc} = - \\frac{2}{a \\cdot \\pi} \\cdot \\frac{1}{1 + \\left(\\frac{plc}{a}\\right)^2}
	  \\]
	- This derivative provides insight into the rate of change of saturation in response to variations in capillary pressure. A steep gradient indicates high sensitivity, while a flatter gradient suggests lower sensitivity.
	"""
	import numpy
	Coeff = - 2 / ( a * numpy.pi )
	return Coeff * 1 / ( 1 + ( plc / a ) **2 )
# 
def Initial_plc(Sc0,a):
	"""
	Computes the initial capillary pressure needed to achieve a specified initial saturation in the extra-vascular space (Eq. derived from saturation formula).

	Parameters:
	-----------
	Sc0 : float
	    The expected initial saturation in the cells of the extra-vascular space. This value should be between 0 and 1.

	a : float
	    A constant parameter that characterizes the connective tissue microstructure and influences how capillary pressure affects saturation.

	Returns:
	--------
	plc_initial : float
	    The initial capillary pressure required to achieve the specified initial saturation `Sc0`.

	Notes:
	------
	- The function is derived from inverting the saturation equation:
	  \\[
	  plc\\_initial = a \\cdot \\tan\\left(\\frac{\\pi}{2} \\cdot (1 - Sc0)\\right)
	  \\]
	- This formula is based on the assumption that the saturation function follows a particular relationship with capillary pressure, parameterized by `a`.
	"""
	import numpy
	return a * numpy.tan( numpy.pi / 2 * ( 1 - Sc0 ) )
# 
# 
def Cmc(poro_n,plc,a):
	"""
	Computes the coefficient Cmc as introduced in Eq. 106, which relates the porosity and capillary pressure in the extra-vascular space.

	Parameters:
	-----------
	poro_n : ufl function
	    The porosity at the previous time step, which serves as the reference state for the current calculation.

	plc : ufl function or dolfinx constant
	    The capillary pressure in the tissue, which affects the saturation and, in turn, the porosity.

	a : float
	    A constant parameter that characterizes the connective tissue microstructure, affecting how the capillary pressure influences saturation.

	Returns:
	--------
	Cmc_value : ufl expression
	    The coefficient Cmc, as defined in Eq. 106. This coefficient is used in further calculations related to tissue mechanics and fluid transport.

	Notes:
	------
	- The function calculates the product of the current porosity and the derivative of the saturation with respect to capillary pressure.
	- The derivative of the saturation with respect to capillary pressure, computed via `derivative_Sc_plc(plc, a)`, reflects the sensitivity of the saturation to changes in capillary pressure.
	"""
	return poro_n*derivative_Sc_plc(plc,a)
# 
def Cstate(plc,a):
	"""
	Computes the Cstate coefficient as described in Eq. 122, which is a combination of the saturation in the extra-vascular space and its derivative with respect to capillary pressure.

	Parameters:
	-----------
	plc : ufl function or dolfinx constant
	    The capillary pressure in the tissue. This variable significantly influences the saturation in the extra-vascular space.

	a : float
	    A constant parameter characterizing the connective tissue microstructure. It affects how the capillary pressure influences the saturation.

	Returns:
	--------
	Cstate_value : ufl expression
	    The Cstate coefficient, calculated using the saturation and its derivative with respect to capillary pressure, as defined in Eq. 122.

	Notes:
	------
	- The `Cstate` coefficient integrates both the saturation in the extra-vascular space (`saturation_cell(plc, a)`) and its sensitivity to capillary pressure (`plc * derivative_Sc_plc(plc, a)`).
	- This coefficient is crucial for understanding the mechanical behavior of the tissue in response to changes in capillary pressure.
	"""
	return saturation_cell(plc,a)+plc*derivative_Sc_plc(plc,a)
# 
#____________________________________________________________#
# 					 Vascular fluid
#____________________________________________________________#
# 
# 
#...........................................................#
# 					 Monophasic
#...........................................................#
# 
def vascular_porosity_atan_mono(epsb0,pl,pb,K):
	"""
	Computes the vascular porosity based on pressure differences and tissue properties, following the equation described in Eq. 118.

	Parameters:
	-----------
	epsb0 : float or ufl expression
	    The initial vascular porosity, representing the porosity at the start of the process.

	pl : dolfinx.function or ufl expression
	    The pressure of the non-wetting phase (e.g., interstitial fluid pressure).

	pb : dolfinx.function or ufl expression
	    The blood pressure, driving the flow within the vascular system.

	K : float
	    A parameter related to vessel compressibility, which governs how responsive the vascular porosity is to pressure differences.

	Returns:
	--------
	vascular_porosity : ufl expression
	    The updated vascular porosity, which accounts for the interplay between fluid phase pressure,  
	    and blood pressure, as described by Eq. 33 but with the atan expression.

	Notes:
	------
	- This function computes the vascular porosity `varepsilon_b` based on the initial porosity `epsb0` and pressure conditions in the system.
	- The pressure difference `pl  - pb` governs how much the fluid phase and blood pressure influence the porosity.
	- The `atan` function is used to smoothly transition between pressure states, with the constant `K` controlling the vessel's compressibility and response to pressure changes.
	
	Formula:
	--------
	- The vascular porosity is computed as:
	  vascular_porosity = epsb0 * (1 - (2 / pi) * atan((pl - pb) / K))
	"""
	import numpy
	import ufl
	return epsb0*( 1-( 2/numpy.pi*ufl.atan( ( pl - pb ) / K ) ) )
# 
def Cep_atan_mono(epsb0,K,pl,pb):
	"""
	Computes the Cep coefficient, which represents the sensitivity of vascular porosity to pressure changes, following Eq. 123.

	Parameters:
	-----------
	epsb0 : float
	    The initial vascular porosity, representing the baseline porosity before any pressure changes.

	K : float
	    A parameter related to the compressibility of the vessel, determining how sensitive the porosity is to pressure variations.

	pl : dolfinx.function or ufl expression
	    The pressure of the non-wetting phase (e.g., interstitial fluid pressure).

	pb : dolfinx.function or ufl expression
	    The blood pressure, which drives fluid flow in the vascular system.

	Returns:
	--------
	Cep_coefficient : ufl expression
	    The coefficient `Cep`, which represents the sensitivity of vascular porosity to pressure differences. This value is derived from the difference between 
	    non-wetting phase pressure, blood pressure, and the compressibility of the vessel system.

	Notes:
	------
	- This function computes the coefficient `Cep`, which adjusts the porosity based on pressure gradients and vessel compressibility.
	- The denominator `(1 + ((pl - pb) / K)^2)` describes how the pressure difference is smoothed by the compressibility parameter `K`, 
	  ensuring a stable transition between pressure states.

	Formula:
	--------
	- The Cep coefficient is computed as:
	  Cep_coefficient = -(2 * epsb0) / (K * pi) * 1 / (1 + ((pl - pb) / K)^2)
	"""
	import numpy
	Coeff = -(2*epsb0)/(K*numpy.pi)
	return Coeff * 1/(1+((pl-pb)/K)**2)
# 
def vascular_porosity_linear_mono(epsb0,pl,pb,K):
	"""
	Computes the vascular porosity based on pressure differences and tissue properties, following the equation described in Eq. 33.

	Parameters:
	-----------
	epsb0 : float or ufl expression
	    The initial vascular porosity, representing the porosity at the start of the process.

	pl : dolfinx.function or ufl expression
	    The pressure of the non-wetting phase (e.g., interstitial fluid pressure).

	pb : dolfinx.function or ufl expression
	    The blood pressure, driving the flow within the vascular system.

	K : float
	    A parameter related to vessel compressibility, which governs how responsive the vascular porosity is to pressure differences.

	Returns:
	--------
	vascular_porosity : ufl expression
	    The updated vascular porosity, which accounts for the interplay between fluid phase pressure, and blood pressure.

	Notes:
	------
	- This function computes the vascular porosity `varepsilon_b` based on the initial porosity `epsb0` and pressure conditions in the system.
	
	Formula:
	--------
	- The vascular porosity is computed as:
	  vascular_porosity = epsb0 * (1 - (  pl - pb ) / K)
	"""
	return epsb0*( 1- ( pl - pb ) / K )
# 
def Cep_linear_mono(epsb0,K):
	"""
	Computes the Cep coefficient, which represents the sensitivity of vascular porosity to pressure changes, following Eq. 36.

	Parameters:
	-----------
	epsb0 : float
	    The initial vascular porosity, representing the baseline porosity before any pressure changes.

	K : float
	    A parameter related to the compressibility of the vessel, determining how sensitive the porosity is to pressure variations.

	pl : dolfinx.function or ufl expression
	    The pressure of the non-wetting phase (e.g., interstitial fluid pressure).

	pb : dolfinx.function or ufl expression
	    The blood pressure, which drives fluid flow in the vascular system.

	Returns:
	--------
	Cep_coefficient : ufl expression
	    The coefficient `Cep`, which represents the sensitivity of vascular porosity to pressure differences. This value is derived from the difference between 
	    non-wetting phase pressure, blood pressure, and the compressibility of the vessel system.

	Notes:
	------
	- This function computes the coefficient `Cep`, which adjusts the porosity based on pressure gradients and vessel compressibility.
	
	Formula:
	--------
	- The Cep coefficient is computed as:
	  Cep_coefficient = - epsb0 / K 
	"""
	return - epsb0 / K 
# 
def vascular_permeability_scalar_mono(k_b,alpha,epsb0,pl,pb,K,evolution_type):
	"""
	Computes the current vascular permeability scalar based on pressure conditions and a specified evolution model, as described in Eq. 39.

	Parameters:
	-----------
	k_b : float
	    The initial vascular permeability in the preferential direction.

	alpha : float
	    The exponent used to control the permeability evolution (alpha >= 2) as per Eq. 39.

	epsb0 : float
	    The initial vascular porosity, representing the baseline porosity before any pressure changes.

	pl : dolfinx.function or ufl expression
	    The pressure of the non-wetting phase (e.g., interstitial or tissue pressure).

	pb : dolfinx.function or ufl expression
	    The blood pressure, which drives fluid flow in the vascular system.

	K : float
	    A parameter related to the compressibility of the vessel, determining how sensitive the porosity and permeability are to pressure variations.

	evolution_type : str
	    Specifies the type of evolution model to compute vascular porosity. Choose between:
	    - 'linear': Uses a linear evolution model.
	    - 'atan': Uses an arctangent evolution model for porosity.

	Returns:
	--------
	current_k_b : float
	    The updated vascular permeability scalar based on the specified evolution model and input pressures, computed as per Eq. 39.

	Notes:
	------
	- The permeability `current_k_b` evolves according to the vascular porosity, which is dependent on the pressure difference between `pl` (non-wetting phase pressure) 
	  and `pb` (blood pressure), and on the compressibility factor `K`.
	- The evolution can follow two models:
	  1. **Linear model**: Uses `vascular_porosity_linear_mono`.
	  2. **Arctangent model**: Uses `vascular_porosity_atan_mono`.
	- If an invalid `evolution_type` is provided, the function will print an error message and return `None`.

	Formula:
	--------
	The permeability evolution follows the equation:
	current_k_b = k_b * (vascular_porosity / epsb0) ** alpha
	"""
	if evolution_type == 'linear':
		current_k_b = k_b*(vascular_porosity_linear_mono(epsb0,pl,pb,K)/epsb0)**alpha
	elif evolution_type == 'atan':
		current_k_b = k_b*(vascular_porosity_atan_mono(epsb0,pl,pb,K)/epsb0)**alpha
	else:
		print("ERROR in the choice of the evolution type. Choose 'linear' or 'atan'.")
	return current_k_b
# 
#
def matrix_vascular_permeability_constant_mono(k_b, ratio):
	"""
	Computes the initial vascular permeability matrix, accounting for different permeability ratios along each spatial direction, as described in Eq. 40.

	Parameters:
	-----------
	k_b : float
	    The vascular permeability in the preferential direction, which serves as the base permeability value.

	ratio : list of floats
	    A list representing the scaling ratios for each spatial direction (e.g., [ratio_x, ratio_y, (ratio_z)]), where each ratio defines how permeability 
	    in a given direction compares to the preferential direction. The ratio in the preferential direction should be 1.

	Returns:
	--------
	permeability_matrix : ufl tensor
	    The permeability matrix in tensor form. For example, in a 3D case with directions x, y, and z, it will return a matrix like:
	    [[k_b * ratio_x, 0, 0],
	     [0, k_b * ratio_y, 0],
	     [0, 0, k_b * ratio_z]]
	    In lower dimensions (e.g., 2D), it returns a corresponding 2x2 matrix.

	Notes:
	------
	- This function computes a diagonal permeability matrix based on the input base permeability `k_b` and direction-specific scaling factors `ratio`.
	- The permeability in the preferential direction (e.g., x) is `k_b`, while in other directions, it is scaled by the corresponding ratio values in the `ratio` list.
	- The matrix is symmetric and diagonal, reflecting anisotropic permeability in different directions.

	Example:
	--------
	For a 2D problem with base permeability `k_b = 1.0`, and ratios `[1.0, 0.5]` for x and y directions, respectively:

	>>> permeability_matrix = matrix_vascular_permeability(k_b=1.0, ratio=[1.0, 0.5])
	>>> print(permeability_matrix)

	This will output:
	[[1.0, 0.0],
	 [0.0, 0.5]]

	For a 3D case with `k_b = 1.0` and ratios `[1.0, 0.5, 0.3]`:

	>>> permeability_matrix = matrix_vascular_permeability(k_b=1.0, ratio=[1.0, 0.5, 0.3])
	>>> print(permeability_matrix)

	This will output:
	[[1.0, 0.0, 0.0],
	 [0.0, 0.5, 0.0],
	 [0.0, 0.0, 0.3]]
	"""
	import ufl
	import numpy
	import dolfinx
	# Initialize a zero matrix of the same dimension as the ratio list
	dimension = len(ratio)
	Matrix = numpy.zeros((dimension, dimension),dtype=dolfinx.default_scalar_type)
	#
	# Fill the matrix based on the ratio values
	for i in range(dimension):
	    Matrix[i][i] = dolfinx.default_scalar_type(k_b * ratio[i])
	#    
	# Display the permeability matrix for debugging
	print("The permeability matrix looks like:", Matrix)
	#
	return ufl.as_tensor(Matrix)
# 
def matrix_vascular_permeability_mono(k_b, ratio, alpha, epsb0, pl, pb, K, evolution_type):
	"""
	Computes the current vascular permeability matrix, adjusting for different permeability ratios in each spatial direction, and accounts for permeability evolution over time based on pressure conditions (Eq. 40).

	Parameters:
	-----------
	k_b : float
	    The base vascular permeability in the preferential direction (e.g., along the x-axis).

	ratio : list of floats
	    A list of scaling ratios for each spatial direction, where the permeability in the preferential direction is scaled by 1.0, and other directions are scaled by their respective ratios. 
	    For example, in 2D: [ratio_x, ratio_y], and in 3D: [ratio_x, ratio_y, ratio_z].

	alpha : float
	    The exponent that governs the evolution of permeability based on vascular porosity (as per Eq. 39, typically alpha >= 2).

	epsb0 : float
	    The initial vascular porosity, serving as a reference for porosity evolution due to changes in pressure.

	pl : dolfinx.function or ufl expression
	    The non-wetting phase pressure (e.g., interstitial or tissue pressure).

	pb : dolfinx.function or ufl expression
	    The blood pressure driving vascular fluid flow.

	K : float
	    A constant controlling how vessel compressibility affects porosity and permeability changes.

	evolution_type : str
	    Specifies the model to use for porosity evolution. Choose from:
	    - 'linear': Uses a linear porosity evolution model.
	    - 'atan': Uses an arctangent porosity evolution model.

	Returns:
	--------
	current_Matrix : ufl tensor
	    A tensor representing the current permeability matrix based on the given base permeability, directional scaling ratios, and porosity evolution.
	    For a 2D case, this is a 2x2 tensor, and for a 3D case, a 3x3 tensor:
	    - In 2D: [[k_b * ratio_x, 0], [0, k_b * ratio_y]]
	    - In 3D: [[k_b * ratio_x, 0, 0], [0, k_b * ratio_y, 0], [0, 0, k_b * ratio_z]]

	Notes:
	------
	- The permeability matrix is diagonal, with permeability in each direction scaled by the respective ratio. In the preferential direction (e.g., x), the ratio is 1.0.
	- This function computes the matrix and then adjusts it based on porosity evolution, which depends on the pressure difference between `pl` and `pb`.
	- Two evolution models are supported:
	  1. **Linear model**: Uses `vascular_porosity_linear_mono`.
	  2. **Arctangent model**: Uses `vascular_porosity_atan_mono`.
	"""
	import ufl
	import numpy
	import dolfinx
	# Initialize a zero matrix of the same dimension as the ratio list
	dimension = len(ratio)
	Matrix = numpy.zeros((dimension, dimension),dtype=dolfinx.default_scalar_type)
	#
	# Fill the matrix based on the ratio values
	for i in range(dimension):
	    Matrix[i][i] = dolfinx.default_scalar_type(k_b * ratio[i])
	#    
	# Display the permeability matrix for debugging
	print("The permeability matrix looks like:", Matrix)
	#
	# Return the permeability matrix as a UFL tensor
	if evolution_type == 'linear':
		current_Matrix = ((vascular_porosity_linear_mono(epsb0,pl,pb,K)/epsb0)**alpha) * ufl.as_tensor(Matrix)
	elif evolution_type == 'atan':
		current_Matrix = ((vascular_porosity_atan_mono(epsb0,pl,pb,K)/epsb0)**alpha) * ufl.as_tensor(Matrix)
	else:
		print("ERROR in the choice of the evolution type. Choose 'linear' or 'atan'.")
	return current_Matrix
# 
#...........................................................#
# 					 Bi-phasic specific laws
#...........................................................#
# 
def vascular_porosity_atan_bi(epsb0,pl,plc,pb,a,K):
	"""
	Computes the vascular porosity based on pressure differences and tissue properties, following the equation described in Eq. 118.

	Parameters:
	-----------
	epsb0 : float or ufl expression
	    The initial vascular porosity, representing the porosity at the start of the process.

	pl : dolfinx.function or ufl expression
	    The pressure of the non-wetting phase (e.g., interstitial or tissue pressure).

	plc : dolfinx.function or ufl expression
	    The capillary pressure, which controls the fluid saturation in the connective tissue.

	pb : dolfinx.function or ufl expression
	    The blood pressure, driving the flow within the vascular system.

	a : float
	    A constant parameter that depends on the microstructure of the connective tissue, affecting how the capillary pressure influences fluid exchange.

	K : float
	    A parameter related to vessel compressibility, which governs how responsive the vascular porosity is to pressure differences.

	Returns:
	--------
	vascular_porosity : ufl expression
	    The updated vascular porosity, which accounts for the interplay between non-wetting phase pressure, capillary pressure, 
	    and blood pressure, as described by Eq. 55.

	Notes:
	------
	- This function computes the vascular porosity `varepsilon_b` based on the initial porosity `epsb0` and pressure conditions in the system.
	- The pressure difference `(pl - saturation_cell(plc, a) * plc) - pb` governs how much the non-wetting phase and blood pressure influence the porosity.
	- The `atan` function is used to smoothly transition between pressure states, with the constant `K` controlling the vessel's compressibility and response to pressure changes.
	- `saturation_cell(plc, a)` is assumed to be a function that calculates the saturation level based on the capillary pressure `plc` and tissue parameter `a`.

	Formula:
	--------
	- The vascular porosity is computed as:
	  vascular_porosity = epsb0 * (1 - (2 / pi) * atan(((pl - saturation_cell(plc, a) * plc) - pb) / K))
	"""
	import numpy
	import ufl
	return epsb0*( 1-( 2/numpy.pi*ufl.atan( ( ( pl - saturation_cell(plc,a) * plc ) - pb )/K ) ) )
# 
def Cep_atan_bi(epsb0,K,plc,pl,pb,a):
	"""
	Computes the Cep coefficient, which represents the sensitivity of vascular porosity to pressure changes, following Eq. 123.

	Parameters:
	-----------
	epsb0 : float
	    The initial vascular porosity, representing the baseline porosity before any pressure changes.

	K : float
	    A parameter related to the compressibility of the vessel, determining how sensitive the porosity is to pressure variations.

	plc : dolfinx.function or ufl expression
	    The capillary pressure, which affects fluid exchange and porosity in the connective tissue.

	pl : dolfinx.function or ufl expression
	    The pressure of the non-wetting phase (e.g., interstitial or tissue pressure).

	pb : dolfinx.function or ufl expression
	    The blood pressure, which drives fluid flow in the vascular system.

	a : float
	    A constant parameter that depends on the microstructure of the connective tissue, influencing the effect of capillary pressure on porosity.

	Returns:
	--------
	Cep_coefficient : ufl expression
	    The coefficient `Cep`, which represents the sensitivity of vascular porosity to pressure differences. This value is derived from the difference between 
	    non-wetting phase pressure, blood pressure, and the compressibility of the vessel system, as described in Eq. 60.

	Notes:
	------
	- This function computes the coefficient `Cep`, which adjusts the porosity based on pressure gradients and vessel compressibility.
	- The term `saturation_cell(plc, a)` is assumed to compute the fluid saturation in the tissue as a function of capillary pressure `plc` and microstructural parameter `a`.
	- The denominator `(1 + ((pl - saturation_cell(plc, a) * plc - pb) / K)^2)` describes how the pressure difference is smoothed by the compressibility parameter `K`, 
	  ensuring a stable transition between pressure states.

	Formula:
	--------
	- The Cep coefficient is computed as:
	  Cep_coefficient = -(2 * epsb0) / (K * pi) * 1 / (1 + ((pl - saturation_cell(plc, a) * plc - pb) / K)^2)
	"""
	import numpy
	Coeff = -(2*epsb0)/(K*numpy.pi)
	return Coeff * 1/(1+((pl-saturation_cell(plc,a)*plc-pb)/K)**2)
# 
# 
def matrix_vascular_permeability_bi(k_b, ratio, alpha, epsb0, pl, plc, pb, a, K):
	"""
	Computes the current vascular permeability matrix, adjusting for different permeability ratios in each spatial direction, and accounts for permeability evolution over time based on pressure conditions (Eq. 40).

	Parameters:
	-----------
	k_b : float
	    The base vascular permeability in the preferential direction (e.g., along the x-axis).

	ratio : list of floats
	    A list of scaling ratios for each spatial direction, where the permeability in the preferential direction is scaled by 1.0, and other directions are scaled by their respective ratios. 
	    For example, in 2D: [ratio_x, ratio_y], and in 3D: [ratio_x, ratio_y, ratio_z].

	alpha : float
	    The exponent that governs the evolution of permeability based on vascular porosity (as per Eq. 39, typically alpha >= 2).

	epsb0 : float
	    The initial vascular porosity, serving as a reference for porosity evolution due to changes in pressure.

	pl : dolfinx.function or ufl expression
	    The non-wetting phase pressure (e.g., interstitial or tissue pressure).

	plc : dolfinx.function or ufl expression
	    The capillary pressure, which affects fluid exchange and porosity in the connective tissue.

	pb : dolfinx.function or ufl expression
	    The blood pressure driving vascular fluid flow.

	a : float
	    A constant parameter that depends on the microstructure of the connective tissue, influencing the effect of capillary pressure on porosity.

	K : float
	    A constant controlling how vessel compressibility affects porosity and permeability changes.

	Returns:
	--------
	current_Matrix : ufl tensor
	    A tensor representing the current permeability matrix based on the given base permeability, directional scaling ratios, and porosity evolution.
	    For a 2D case, this is a 2x2 tensor, and for a 3D case, a 3x3 tensor:
	    - In 2D: [[k_b * ratio_x, 0], [0, k_b * ratio_y]]
	    - In 3D: [[k_b * ratio_x, 0, 0], [0, k_b * ratio_y, 0], [0, 0, k_b * ratio_z]]

	Notes:
	------
	- The permeability matrix is diagonal, with permeability in each direction scaled by the respective ratio. In the preferential direction (e.g., x), the ratio is 1.0.
	- This function computes the matrix and then adjusts it based on porosity evolution, which depends on the pressure difference between `ps` and `pb`.
	- A single evolution models are supported:
	  1. **Arctangent model**: Uses `vascular_porosity_atan_bi`.
	"""
	import ufl
	import numpy
	import dolfinx
	# Initialize a zero matrix of the same dimension as the ratio list
	dimension = len(ratio)
	Matrix = numpy.zeros((dimension, dimension),dtype=dolfinx.default_scalar_type)
	#
	# Fill the matrix based on the ratio values
	for i in range(dimension):
	    Matrix[i][i] = dolfinx.default_scalar_type(k_b * ratio[i])
	#    
	# Display the permeability matrix for debugging
	print("The permeability matrix looks like:", Matrix)
	#
	current_Matrix = ((vascular_porosity_atan_bi(epsb0,pl,plc,pb,a,K)/epsb0)**alpha) * ufl.as_tensor(Matrix)
	# 
	return current_Matrix
# 
def vascular_permeability_scalar_bi(k_b, alpha, epsb0, pl, plc, pb, a, K):
	"""
	Computes the current vascular permeability scalar based on pressure conditions and a specified evolution model, as described in Eq. 39.

	Parameters:
	-----------
	k_b : float
	    The initial vascular permeability in the preferential direction.

	alpha : float
	    The exponent used to control the permeability evolution (alpha >= 2) as per Eq. 39.

	epsb0 : float
	    The initial vascular porosity, representing the baseline porosity before any pressure changes.

	pl : dolfinx.function or ufl expression
	    The pressure of the non-wetting phase (e.g., interstitial or tissue pressure).

	pb : dolfinx.function or ufl expression
	    The blood pressure, which drives fluid flow in the vascular system.

	K : float
	    A parameter related to the compressibility of the vessel, determining how sensitive the porosity and permeability are to pressure variations.

	evolution_type : str
	    Specifies the type of evolution model to compute vascular porosity. Choose between:
	    - 'linear': Uses a linear evolution model.
	    - 'atan': Uses an arctangent evolution model for porosity.

	Returns:
	--------
	current_k_b : float
	    The updated vascular permeability scalar based on the specified evolution model and input pressures, computed as per Eq. 39.

	Notes:
	------
	- The permeability `current_k_b` evolves according to the vascular porosity, which is dependent on the pressure difference between `pl` (non-wetting phase pressure) 
	  and `pb` (blood pressure), and on the compressibility factor `K`.
	- The evolution can follow two models:
	  1. **Linear model**: Uses `vascular_porosity_linear_mono`.
	  2. **Arctangent model**: Uses `vascular_porosity_atan_mono`.
	- If an invalid `evolution_type` is provided, the function will print an error message and return `None`.

	Formula:
	--------
	The permeability evolution follows the equation:
	current_k_b = k_b * (vascular_porosity / epsb0) ** alpha
	"""
	current_k_b = k_b*(vascular_porosity_atan_bi(epsb0,pl,plc,pb,a,K)/epsb0)**alpha
	return current_k_b
# 
#____________________________________________________________#
# 					 Biological species (O2)
#____________________________________________________________#
# 
# Make it only for bi because of hypotheses and keep varepsilon >0  + check variable unuseful anymore and update equation number and do help nicely
# 
# def MO2_b_l(epsb0,pl,plc,pb,a,K,wo2_b,wo2):
# 	"""
# 	Mass exchange between the blood and the IF Eq. XX.
# 	Input:
# 		- epsb0: The initial vascular porosity,
# 		- pl: non-wetting phase pressure,
# 		- plc: capillary pressure,
# 		- pb: blood pressure,
# 		- a: constant parameter depending on the connective tissue microstructure,
# 		- K: parameter related to the vessel compressibility,
# 		- hv: coefficient representative of the vessel wall permeability,
# 		- wo2b: O2 fraction in the blood,
# 		- wo2: O2 mass fraction in the fluid.
# 	Output:
# 		- the mass exchange term Eq. XX.
# 	"""
# 	return hv*vascular_porosity_atan_bi(epsb0,pl,plc,pb,a,K)*(wo2_b-wo2)
# # 
# def MO2_l_c(wcrit,wo2,gamma_0,epsb0,pl,plc,pb,du,a,K,poro_n,k_b,mu_b,ratio,direction,dimension,alpha,dt,poro_min,poro_max):
# 	"""
# 	O2 consuption from the cells Eq. XX-XX.
# 	Input:
# 		- wcrit: hypoxia threshold,
# 		- wo2: O2 mass fraction in the fluid,
# 		- gamma_0: cell metabolism,
# 		- epsb0: The initial vascular porosity,
# 		- pl: non-wetting phase pressure,
# 		- plc: capillary pressure,
# 		- pb: blood pressure,
# 		- du: displacement increment,
# 		- a: constant parameter depending on the connective tissue microstructure,
# 		- K: parameter related to the vessel compressibility,
# 		- poro_n: the porosity at the previous time step,
# 		- k_b: Vascular permeability of the preferential direction,
# 		- mu_b:  blood viscosity,
# 		- ratio: k_b/ratio for other directions,
# 		- direction: principal direction (0,1 (,2 if 3D)),
# 		- dimension: dimension of the problem (2 or 3),
# 		- alpha: exponent >= 2 according to Eq. 63,
# 		- dt: time increment.
# 	Output:
# 		- the mass exchange term Eq. XX.
# 	"""
# 	import ufl
# 	import numpy
# 	Sc         = saturation_cell(plc,a)
# 	KK         = matrix_vascular_permeability_bi(k_b, ratio, alpha, epsb0, pl, plc, pb, a, K)
# 	varepsilon = porosity_2(du,poro_n,KK,mu_b,pb,dt,poro_min,poro_max)
# 	truevalue  = Sc*varepsilon*gamma_0*0.5*(1 - ufl.cos(numpy.pi * wo2/wcrit ) ) 
# 	falsevalue = Sc*varepsilon*gamma_0
# 	return ufl.conditional(wo2<=wcrit,truevalue,falsevalue)
# # 
# def MO2_l_c_constant(wo2,gamma_0,epsb0,pl,plc,pb,du,a,K,poro_n,k_b,mu_b,ratio,direction,dimension,alpha,dt,poro_min,poro_max):
# 	"""
# 	O2 consuption from the cells Eq. XX.
# 	Input:
# 		- wcrit: hypoxia threshold,
# 		- wo2: O2 mass fraction in the fluid,
# 		- gamma_0: cell metabolism,
# 		- epsb0: The initial vascular porosity,
# 		- pl: non-wetting phase pressure,
# 		- plc: capillary pressure,
# 		- pb: blood pressure,
# 		- du: displacement increment,
# 		- a: constant parameter depending on the connective tissue microstructure,
# 		- K: parameter related to the vessel compressibility,
# 		- poro_n: the porosity at the previous time step,
# 		- k_b: Vascular permeability of the preferential direction,
# 		- mu_b:  blood viscosity,
# 		- ratio: k_b/ratio for other directions,
# 		- direction: principal direction (0,1 (,2 if 3D)),
# 		- dimension: dimension of the problem (2 or 3),
# 		- alpha: exponent >= 2 according to Eq. XX,
# 		- dt: time increment.
# 	Output:
# 		- the mass exchange term Eq. XX.
# 	"""
# 	import ufl
# 	import numpy
# 	Sc         = saturation_cell(plc,a)
# 	KK         = constitutive_laws.matrix_vascular_permeability_bi(k_b, ratio, alpha, epsb0, pl, plc, pb, a, K)
# 	varepsilon = porosity_2(du,poro_n,KK,mu_b,pb,dt,poro_min,poro_max)
# 	return Sc*varepsilon*gamma_0
# # 
# def Deff_O2(D0_o2,delta,pl,plc,pb,a,K,k_b,ratio,direction,dimension,alpha,du,poro_n,mu_b,dt,poro_min,poro_max):
# 	"""
# 	Compute the effective diffusion coefficient of the oxygen Eq. XX.
# 	Input:
# 		- D0_o2: the diffusion coefficient of oxygen in the bulk interstitial fluid
# 		- delta: coefficient related to the tortuosity of the medium
# 		- pl: non-wetting phase pressure,
# 		- plc: capillary pressure,
# 		- pb: blood pressure,
# 		- a: constant parameter depending on the connective tissue microstructure,
# 		- K: parameter related to the vessel compressibility,
# 		- k_b: Vascular permeability of the preferential direction,
# 		- ratio: k_b/ratio for other directions,
# 		- direction: principal direction (0,1 (,2 if 3D)),
# 		- dimension: dimension of the problem (2 or 3),
# 		- alpha: exponent >= 2 according to Eq. XX,
# 		- du: displacement increment,
# 		- poro_n: the porosity at the previous time step,
# 		- mu_b:  blood viscosity,
# 		- dt: time increment.
# 	Output:
# 		- Effective diffusion according to Eq. XX. 
# 	"""
# 	KK = matrix_vascular_permeability_bi(k_b, ratio, alpha, epsb0, pl, plc, pb, a, K)
# 	return D0_o2*( (1-saturation_cell(plc,a)) * porosity_2(du,poro_n,KK,mu_b,pb,dt,poro_min,poro_max) )**delta
# 
#____________________________________________________________#
# 					 End of functions
#____________________________________________________________#
# 
if __name__ == "__main__":
    print("Loading of the constitutive laws successfully completed.")
    # EoF