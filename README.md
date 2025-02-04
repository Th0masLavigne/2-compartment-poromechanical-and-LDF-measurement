# A 2-compartment poromechanical approach to account for the vascular network in human skin: Theoretical developments and comparison to a *in vivo* Laser Doppler Flowmetry measurement.

This repository contains the codes used to generate the results presented in *Lavigne et al.*[^1]. This paper proposes a proof of concept for a poroelastic model to account to simulate the mechanical behaviour and micro-circulation of human skin. A 2-compartment model has been considered and a sensitivity analysis has been performed to study the dominant parameters of the model to reproduce the micro-circulation of *in vivo* human skin, using a porous media approach. The simulation outputs were compared to an experimental campaign (ethical agreement PREVAILXXX). The here-after functions are consistent with a version 0.9.0 of FEniCSx.


Please run first `python3 setup.py build` then `sudo python3 setup.py install` at ./Theoretical_developments_and_service_files/0_Service_files/porous_fenicsx to create the package. This allows not to copy paste the service files (modules) but use them directly as a package.

## Acknowledgments

This research was funded in whole or in part by the Luxembourg National Research Fund (FNR), grant reference No. 17013182. For the purpose of open access, the author has applied a Creative Commons Attribution 4.0 International (CC BY 4.0) licence to any Author-Accredited manuscript version arising from this submission.

The authors also acknowledge the use of the Cassiopee Arts et Métiers Institute of Technology HPC Centre made available for conducting the research reported in this article.


## References
[^1]: *Lavigne et al*, 202X, A 2-compartment poromechanical approach to account for the vascular network in human skin: Theoretical developments and comparison to a *in vivo* Laser Doppler Flowmetry measurement.


## Archive organization

Here-after is a graphical representation of the contentes of this repository:

```
+ ./
|
+-------+ ./Experimental_Data/
|		...........*Contains the excel file of the LDF signal for all the subjects*
|		|
|		+ Subject_ID.xlsx
|		:
|		:
|		+ Subject_ID.xlsx
|		|
|		+ LDF_post_treatment.py
|		...........*Computes the metrics and creates the exeperimental plots.*
|		|
|
|
|-------+ ./1D_Column/
|		...........*Contains the main code and service files (same as the ones of ./Theoretical_developments_and_service_files) to obtain the 1D consolidation test case*
|		|
|		+ constitutive_laws.py
|		...........*All the constitutive laws according to the theoretical developments*
|		|
|		+ functions.py
|		...........*Additional required functions*
|		|
|		+ main_column_2_bi.py
|		...........*Main file to compute the finite element solution for the 1D consolidation test on a column*
|		|
|		+ variational_forms.py
|		...........*All the variational forms according to the theoretical developments using backward Euler*
|		|
|		+ Mesh.msh
|		...........*The Mesh File*
|		|
|		+-----+ ./post-treatment/
|		|     ...........*Contains the code to plot the results*
|		|
|		
|
|-------+ ./3D_Phalanx/
|		...........*Contains the main code and service files (same as the ones of ./Theoretical_developments_and_service_files) to obtain the 3D phalanx indented case*
|		|
|		+-------+ ./Mesh/
|		|		+ New_geom_finger_gmsh.py
|		|		...........*GMSH API code to generate the mesh. Handles local refinement*
|		|		|
|		|		+ Mesh_2_poro.msh
|		|		...........*The Mesh File*
|		|		|
|		|		+ Mesh_2_poro_local_refine.msh
|		|		...........*The locally refined Mesh File*
|		|		|
|		|		|
|		|
|		+-------+ ./Finite_Element_Files/
|		|		+ constitutive_laws.py
|		|		...........*All the constitutive laws according to the theoretical developments*
|		|		|
|		|		+ functions.py
|		|		...........*Additional required functions*
|		|		|
|		|		+ main_finger.py
|		|		...........*Main file to compute the finite element solution for the 3D indentation test on a mid phalanx*
|		|		|
|		|		+ variational_forms.py
|		|		...........*All the variational forms according to the theoretical developments using backward Euler*
|		|		|
|		|		+ New_geom_finger_gmsh.py
|		|		...........*GMSH API code to generate the mesh. Handles local refinement*
|		|		|
|		|		+ Mesh_2_poro.msh
|		|		...........*The Mesh File*
|		|		|
|		|		+ Mesh_2_poro_local_refine.msh
|		|		...........*The locally refined Mesh File*
|		|		|
|		|		|
|		|
|		+-------+ ./Sensitivity_analysis/
|		|		+ sensitivity_analysis.py
|		|		...........*Compute the Sobol indices of the model parameters: first and second order*
|		|		|
|		|		+ README.md
|		|		...........*Provides the parameters used for each computation*
|		|		|
|		|		+ ID_1_export.csv
|		|		...........*Export quantities of the model for the parameters ID_1*
|		|		:
|		|		+ ID_XX_export.csv
|		|		...........*Export quantities of the model for the parameters ID_XX*
|		|		:
|		|		|
|		|		+ appendix_plot_sensitivity.py
|		|		...........*Contains the code to plot the results for the appendix*
|		|		|
|		|		|
|		|
|
|
+-------+ ./Theoretical_developments_and_service_files/
|		|
|		+ Readme.md 
|		...........*This is the readme file you are currently reading*
|		| 
|		+ .gitignore
|		...........*Avoid adding the simulation results to the repository*
|		|
|		+-----+ ./0_Service_files/
|		|     |
|		|     + Lavigne_et_al_2023.pdf
|		|     ...........*Tutorial article in version 0.5.2 on the implementation of poromechanics in the FEniCSx framework*
|		|     |
|		|     + Theoretical_Developments.pdf
|		|     ...........*All the theoretical definitions of the constitutive laws and variational forms required for the different models*
|		|     |
|		|     + constitutive_laws.py -+
|		|     ...........*All the constitutive laws according to the theoretical developments*
|		|     |                       + Lame_coefficients()
|		|     |                       + Local_deformation()
|		|     |                       + Elastic_constitutive_law()
|		|     |                       + porosity_1()
|		|     |                       + porosity_2()
|		|     |                       + porosity_2_pressure()
|		|     |                       + saturation_cell()
|		|     |                       + derivative_Sc_plc()
|		|     |                       + Initial_plc()
|		|     |                       + Cmc()
|		|     |                       + Cstate()
|		|     |                       + vascular_porosity_atan_mono()
|		|     |                       + Cep_atan_mono()
|		|     |                       + vascular_porosity_linear_mono()
|		|     |                       + Cep_linear_mono()
|		|     |                       + vascular_permeability_scalar_mono()
|		|     |                       + matrix_vascular_permeability_constant_mono()
|		|     |                       + matrix_vascular_permeability_mono()
|		|     |                       + vascular_porosity_atan_bi()
|		|     |                       + Cep_atan_bi()
|		|     |                       + matrix_vascular_permeability_bi()
|		|     |                       + MO2_b_l()
|		|     |                       + MO2_l_c()
|		|     |                       + MO2_l_c_constant()
|		|     |                       + Deff_O2()
|		|     |
|		|     |
|		|     + variational_forms.py -+
|		|     ...........*All the variational forms according to the theoretical developments using backward Euler*
|		|     |                       + variational_form_1_mono_total()
|		|     |                       + variational_form_1_mono_updated()
|		|     |                       + variational_form_2_mono_updated_isotropic_k_b()
|		|     |                       + variational_form_2_mono_updated_anisotropic_k_b()
|		|     |                       + variational_form_1_comp_bi()
|		|     |                       + variational_form_2_comp_bi()
|		|     |                       + variational_form_oxygen_constant_cell()
|		|     |                       + variational_form_oxygen_conditional_cell()
|		|     |                       + variational_form_2_mono_updated_anisotropic_k_b_constant()
|		|     |                       + variational_form_2_mono_updated_anisotropic_k_b_product_outside_dot()
|		|     |
|		|     |
|		|     + variational_forms_TW.py -+
|		|     ...........*All the variational forms according to the theoretical developments using Theta Wilson*
|		|     |                       + variational_form_1_mono_total()
|		|     |                       + variational_form_1_mono_updated()
|		|     |                       + variational_form_2_mono_updated_isotropic_k_b()
|		|     |                       + variational_form_2_mono_updated_anisotropic_k_b()
|		|     |                       + variational_form_1_comp_bi()
|		|     |                       + variational_form_2_comp_bi()
|		|     |                       + variational_form_oxygen_constant_cell()
|		|     |                       + variational_form_oxygen_conditional_cell()
|		|     |                       + variational_form_2_mono_updated_anisotropic_k_b_constant()
|		|     |                       + variational_form_2_mono_updated_anisotropic_k_b_product_outside_dot()
|		|     |
|		|     + functions.py -+
|		|     ...........*Additional required functions*
|		|     |                       + set_non_linear_solver_parameters()
|		|     |                       + export_to_csv()
|		|     |                       + RMSE()
|		|     |                       + mechanical_load()
|		|     |                       + mechanical_load_LDF()
|		|     |                       + terzaghi()
|		|
|		|
|		+-----+ ./1_Mesh_functions/
|		|     ...........*Contains the .msh files and python functions to generate them*
|		|     
|		|
|		+-----+ ./2_main_files/
|		|     ...........*Contains the main_XXX.py files to run the codes. 1 and 2 stand for single and two compartments; mono and bi stands for monophasic or biphasic interstitium*
|		|
|		|
|		+-----+ ./3_Results/
|		|     ...........*Contains the results of the main_XXX.py codes.*
|		|
|		|
|
|
```

