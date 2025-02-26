# Introduction of Multi-compartment poromechanics for living tissues

This repository contains the functions created in FEniCSx to model the biomechanical response of skin as part of Thomas Lavigne PhD. The here-after functions are consistent with a version 0.9.0 of FEniCSx.


Please run first `python3 -m pip install .`  at ./0_Service_files/porous_fenicsx to create the package.


## Description of the repository

Here-after is a graphical representation of the contentes of this repository:

```
+ ./
|
+ Readme.md 
...........*This is the readme file you are currently reading*
|
+ ----+ ./0_Service_files/
|     |
|     + Lavigne_et_al_2023.pdf
|     ...........*Tutorial article in version 0.5.2 on the implementation of poromechanics in the FEniCSx framework*
|     |
|     + Theoretical_Developments.pdf
|     ...........*All the theoretical definitions of the constitutive laws and variational forms required for the different models*
|     |
|     + pyproject.toml -+
|     ...........*Package porous_fenicsx setup file*
|     |
|     + ./porous_fenicsx/constitutive_laws.py -+
|     ...........*All the constitutive laws according to the theoretical developments*
|     |                       + Lame_coefficients()
|     |                       + Local_deformation()
|     |                       + Elastic_constitutive_law()
|	  |                       + Neo_Hooke()
|     |                       + porosity_1()
|     |                       + porosity_2()
|     |                       + porosity_2_pressure()
|     |                       + saturation_cell()
|     |                       + derivative_Sc_plc()
|     |                       + Initial_plc()
|     |                       + Cmc()
|     |                       + Cstate()
|     |                       + vascular_porosity_atan_mono()
|     |                       + Cep_atan_mono()
|     |                       + vascular_porosity_linear_mono()
|     |                       + Cep_linear_mono()
|     |                       + vascular_permeability_scalar_mono()
|     |                       + matrix_vascular_permeability_constant_mono()
|     |                       + matrix_vascular_permeability_mono()
|     |                       + vascular_porosity_atan_bi()
|     |                       + Cep_atan_bi()
|     |                       + matrix_vascular_permeability_bi()
|     |                       + MO2_b_l()
|     |                       + MO2_l_c()
|     |                       + MO2_l_c_constant()
|     |                       + Deff_O2()
|     |
|     |
|     + ./porous_fenicsx/variational_forms.py -+
|     ...........*All the variational forms according to the theoretical developments using backward Euler*
|     |                       + variational_form_1_mono_total()
|     |                       + variational_form_1_mono_updated()
|     |                       + variational_form_2_mono_updated_isotropic_k_b()
|     |                       + variational_form_2_mono_updated_anisotropic_k_b()
|     |                       + variational_form_1_comp_bi()
|     |                       + variational_form_2_comp_bi()
|     |                       + variational_form_oxygen_constant_cell()
|     |                       + variational_form_oxygen_conditional_cell()
|     |                       + variational_form_2_mono_updated_anisotropic_k_b_constant()
|     |                       + variational_form_2_mono_updated_anisotropic_k_b_product_outside_dot()
|     |
|     + ./porous_fenicsx/variational_forms_Neo_Hooke.py -+
|     ...........*All the variational forms according to the theoretical developments using backward Euler and Neo-Hooke hyperelastic solid*
|     |                       + variational_form_1_mono_total()
|     |                       + variational_form_1_mono_updated()
|     |                       + variational_form_2_mono_updated_isotropic_k_b()
|     |                       + variational_form_2_mono_updated_anisotropic_k_b()
|     |                       + variational_form_1_comp_bi()
|     |                       + variational_form_2_comp_bi()
|     |                       + variational_form_oxygen_constant_cell()
|     |                       + variational_form_oxygen_conditional_cell()
|     |                       + variational_form_2_mono_updated_anisotropic_k_b_constant()
|     |                       + variational_form_2_mono_updated_anisotropic_k_b_product_outside_dot()
|     |
|	  |
|     |
|     + ./porous_fenicsx/variational_forms_TW.py -+
|     ...........*All the variational forms according to the theoretical developments using Theta Wilson method*
|     |                       + variational_form_1_mono_total()
|     |                       + variational_form_1_mono_updated()
|     |                       + variational_form_2_mono_updated_isotropic_k_b()
|     |                       + variational_form_2_mono_updated_anisotropic_k_b()
|     |                       + variational_form_1_comp_bi()
|     |                       + variational_form_2_comp_bi()
|     |                       + variational_form_oxygen_constant_cell()
|     |                       + variational_form_oxygen_conditional_cell()
|     |                       + variational_form_2_mono_updated_anisotropic_k_b_constant()
|     |                       + variational_form_2_mono_updated_anisotropic_k_b_product_outside_dot()
|     |
|     + ./porous_fenicsx/functions.py -+
|     ...........*Additional required functions*
|     |                       + set_non_linear_solver_parameters()
|     |                       + export_to_csv()
|     |                       + RMSE()
|     |                       + mechanical_load()
|     |                       + mechanical_load_LDF()
|     |                       + terzaghi()
|	  + __init__.py
|	  ...........*Package's modules declaration*
|
|
+ ----+ ./1_Mesh_functions/
|     ...........*Contains the .msh files and python functions to generate them*
|     
|
+ ----+ ./2_main_files/
|     ...........*Contains the main_XXX.py files to run the codes. 1 and 2 stand for single and two compartments; mono and bi stands for monophasic or biphasic interstitium*
|
|
+ ----+ ./3_Results/
|     ...........*Contains the results of the main_XXX.py codes.*
|
|
```

## Acknowledgments

This research was funded in whole or in part by the Luxembourg National Research Fund (FNR), grant reference No. 17013182. For the purpose of open access, the author has applied a Creative Commons Attribution 4.0 International (CC BY 4.0) licence to any Author-Accredited manuscript version arising from this submission.

The authors also acknowledge the use of the Cassiopee Arts et Métiers Institute of Technology HPC Centre made available for conducting the research reported in this article.


