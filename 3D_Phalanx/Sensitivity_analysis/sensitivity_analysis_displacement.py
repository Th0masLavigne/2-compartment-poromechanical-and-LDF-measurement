import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

def difference(x,xref):
	"""
	help
	"""
	import numpy as np
	return 1/len(x)*np.sum([x[i]/xref[i] for i in range(len(xref))])

def compute_alpha(p,pref):
	"""
	help: courbe passe par 0 -> 1 
	"""
	return (p-pref)/pref


# Load the Excel file
file_path = "Sensitivity.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")

# # Display data from each sheet
# for sheet_name, data in all_sheets.items():
#     print(f"Sheet Name: {sheet_name}")
#     print(data.head())  # Print the first 5 rows of each sheet
#     print("\n" + "="*50 + "\n")  # Separator for readability

# Optional: Extract specific columns for all sheets
column_names = ["displacement_all",	"LDF_v_all",	"LDF_q_all",	"LDF_baseline_v",	"LDF_baseline_q",	"load_all",	"time_all"]

# Extract data from each sheet into a structured dictionary
model = {}
sheet_names = []
for sheet_name, data in all_sheets.items():
    model[sheet_name] = {col: data[col] for col in column_names if col in data.columns}
    sheet_names.append(sheet_name)

# Compute the relative difference sum
Size_list = len(model["REF_export"]["displacement_all"])
first_ID = 4
all_difference_sheets = []
for ii in range(len(sheet_names[first_ID:])):
	all_difference_sheets.append(difference(model[sheet_names[first_ID+ii]]["displacement_all"],model["REF_export"]["displacement_all"]))

# Alphas for sobol
alphas = []
pointeur_E = []
pointeur_kl = []
pointeur_kb = []
pointeur_alpha = []
pointeur_K = []
pointeur_Ekl = []
pointeur_Ekb = []
pointeur_Ealpha = []
pointeur_EK = []
pointeur_klkb = []
pointeur_klalpha = []
pointeur_klK = []
pointeur_kbalpha = []
pointeur_kbK = []
pointeur_alphaK = []


# ID_1: kl
alphas.append(compute_alpha(1e-15,1e-14))
pointeur_kl.append(0)

# ID_2 : kl
alphas.append(compute_alpha(1e-13,1e-14))
pointeur_kl.append(1)

# Broken code

# ID_4 : kb
alphas.append(compute_alpha(1e-8,1e-9))
pointeur_kb.append(2)

#ID_5 : 
alphas.append(compute_alpha(2,3))
pointeur_alpha.append(3)

#ID_6 : 
alphas.append(compute_alpha(4,3))
pointeur_alpha.append(4)

#ID_7 : 
alphas.append(compute_alpha(5000,1000))
pointeur_K.append(5)

#ID_8 : 
alphas.append(compute_alpha(500,1000))
pointeur_K.append(6)

#ID_9 : 
alphas.append(compute_alpha(5e4,2e5))
pointeur_E.append(7)

#ID_10 : 
alphas.append([compute_alpha(5e4,2e5),compute_alpha(1e-15,1e-14)])
pointeur_Ekl.append(8)

#ID_11 : 
alphas.append([compute_alpha(5e4,2e5),compute_alpha(1e-13,1e-14)])
pointeur_Ekl.append(9)

# Broken code

#ID_13 : 
alphas.append([compute_alpha(5e4,2e5),compute_alpha(1e-8,1e-9)])
pointeur_Ekb.append(10)

#ID_14 : 
alphas.append([compute_alpha(5e4,2e5),compute_alpha(2,3)])
pointeur_Ealpha.append(11)

#ID_15 : 
alphas.append([compute_alpha(5e4,2e5),compute_alpha(4,3)])
pointeur_Ealpha.append(12)

#ID_16 : 
alphas.append([compute_alpha(5e4,2e5),compute_alpha(5000,1000)])
pointeur_EK.append(13)

#ID_17 : 
alphas.append([compute_alpha(5e4,2e5),compute_alpha(500,1000)])
pointeur_EK.append(14)

#ID_18: 
alphas.append(compute_alpha(5e5,2e5))
pointeur_E.append(15)

#ID_19 : 
alphas.append([compute_alpha(5e5,2e5),compute_alpha(1e-15,1e-14)])
pointeur_Ekl.append(16)

#ID_20 : 
alphas.append([compute_alpha(5e5,2e5),compute_alpha(1e-13,1e-14)])
pointeur_Ekl.append(17)

# Broken code

#ID_22 : 
alphas.append([compute_alpha(5e5,2e5),compute_alpha(1e-8,1e-9)])
pointeur_Ekb.append(18)

#ID_23 : 
alphas.append([compute_alpha(5e5,2e5),compute_alpha(2,3)])
pointeur_Ealpha.append(19)

#ID_24 : 
alphas.append([compute_alpha(5e5,2e5),compute_alpha(4,3)])
pointeur_Ealpha.append(20)

#ID_25 : 
alphas.append([compute_alpha(5e5,2e5),compute_alpha(5000,1000)])
pointeur_EK.append(21)

#ID_26 : 
alphas.append([compute_alpha(5e5,2e5),compute_alpha(500,1000)])
pointeur_EK.append(22)

#ID_27 : 
alphas.append([compute_alpha(1e-15,1e-14),compute_alpha(1e-10,1e-9)])
pointeur_klkb.append(23)

#ID_28 : 
alphas.append([compute_alpha(1e-15,1e-14),compute_alpha(1e-8,1e-9)])
pointeur_klkb.append(24)

#ID_29 : 
alphas.append([compute_alpha(1e-15,1e-14),compute_alpha(2,3)])
pointeur_klalpha.append(25)

#ID_30 : 
alphas.append([compute_alpha(1e-15,1e-14),compute_alpha(3,4)])
pointeur_klalpha.append(26)

#ID_31 : 
alphas.append([compute_alpha(1e-15,1e-14),compute_alpha(5000,1000)])
pointeur_klK.append(27)

#ID_32 : 
alphas.append([compute_alpha(1e-15,1e-14),compute_alpha(500,1000)])
pointeur_klK.append(28)

# Broken code

#ID_34 : 
alphas.append([compute_alpha(1e-13,1e-14),compute_alpha(1e-8,1e-9)])
pointeur_klkb.append(29)

#ID_35 : 
alphas.append([compute_alpha(1e-13,1e-14),compute_alpha(2,3)])
pointeur_klalpha.append(30)

#ID_36 : 
alphas.append([compute_alpha(1e-13,1e-14),compute_alpha(3,4)])
pointeur_klalpha.append(31)

#ID_37 : 
alphas.append([compute_alpha(1e-13,1e-14),compute_alpha(5000,1000)])
pointeur_klK.append(32)

#ID_38 : 
alphas.append([compute_alpha(1e-13,1e-14),compute_alpha(500,1000)])
pointeur_klK.append(33)

# Broken code

# Broken code

# Broken code

# Broken code

#ID_43 : 
alphas.append([compute_alpha(1e-8,1e-9),compute_alpha(2,3)])
pointeur_kbalpha.append(34)

#ID_44 : 
alphas.append([compute_alpha(1e-8,1e-9),compute_alpha(4,3)])
pointeur_kbalpha.append(35)

#ID_45 : 
alphas.append([compute_alpha(1e-8,1e-9),compute_alpha(5000,10000)])
pointeur_kbK.append(36)

#ID_46 : 
alphas.append([compute_alpha(1e-8,1e-9),compute_alpha(500,10000)])
pointeur_kbK.append(37)

#ID_47 : 
alphas.append([compute_alpha(2,3),compute_alpha(5000,10000)])
pointeur_alphaK.append(38)

#ID_48 : 
alphas.append([compute_alpha(2,3),compute_alpha(500,10000)])
pointeur_alphaK.append(39)

#ID_49 : 
alphas.append([compute_alpha(4,3),compute_alpha(5000,10000)])
pointeur_alphaK.append(40)

#ID_50 : 
alphas.append([compute_alpha(4,3),compute_alpha(500,10000)])
pointeur_alphaK.append(41)

#ID_51: 
alphas.append(compute_alpha(3e5,2e5))
pointeur_E.append(42)

#ID_52: 
alphas.append(compute_alpha(1e6,2e5))
pointeur_E.append(43)

# ID_53: kl
alphas.append(compute_alpha(5e-14,1e-14))
pointeur_kl.append(44)

# ID_54 : kl
alphas.append(compute_alpha(5e-15,1e-14))
pointeur_kl.append(45)

# ID_55 : kb
alphas.append(compute_alpha(5e-9,1e-9))
pointeur_kb.append(46)

# ID_56 : kb
alphas.append(compute_alpha(5e-9,5e-10))
pointeur_kb.append(47)

#ID_57 : 
alphas.append(compute_alpha(6,3))
pointeur_alpha.append(48)

#ID_58 : 
alphas.append(compute_alpha(6,3))
pointeur_alpha.append(49)

#ID_59 : 
alphas.append(compute_alpha(2000,1000))
pointeur_K.append(50)

#ID_60 : 
alphas.append(compute_alpha(800,1000))
pointeur_K.append(51)

#ID_61 : 
alphas.append([compute_alpha(1e-13,1e-14),compute_alpha(5,3)])
pointeur_klalpha.append(52)

#ID_62 : 
alphas.append([compute_alpha(1e-13,1e-14),compute_alpha(6,3)])
pointeur_klalpha.append(53)

# skip code

# skip code

# skip code

# skip code

#ID_67 : 
# alphas.append(compute_alpha(3e5,2e5)*compute_alpha(6,3))
# pointeur_Ealpha.append(54) # 58 dans difference list

def function_minimize_E(thetai):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_E])
	expression = np.asarray([1 + thetai * alphas[i] for i in pointeur_E])
	return np.sum(objectif - expression)

theta_cherche_E = scipy.optimize.least_squares(function_minimize_E,1)
theta_E = theta_cherche_E.x[0]
# Linalg .

def function_minimize_kl(thetai):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_kl])
	expression = np.asarray([1 + thetai * alphas[i] for i in pointeur_kl])
	return np.sum(objectif - expression)

theta_cherche_kl = scipy.optimize.least_squares(function_minimize_kl,1)
theta_kl = theta_cherche_kl.x[0]

def function_minimize_kb(thetai):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_kb])
	expression = np.asarray([1 + thetai * alphas[i] for i in pointeur_kb])
	return np.sum(objectif - expression)

theta_cherche_kb = scipy.optimize.least_squares(function_minimize_kb,1)
theta_kb = theta_cherche_kb.x[0]

def function_minimize_alpha(thetai):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_alpha])
	expression = np.asarray([1 + thetai * alphas[i] for i in pointeur_alpha])
	return np.sum(objectif - expression)

theta_cherche_alpha = scipy.optimize.least_squares(function_minimize_alpha,1)
theta_alpha = theta_cherche_alpha.x[0]


def function_minimize_K(thetai):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_K])
	expression = np.asarray([1 + thetai * alphas[i] for i in pointeur_K])
	return np.sum(objectif - expression)

theta_cherche_K = scipy.optimize.least_squares(function_minimize_K,1)
theta_K = theta_cherche_K.x[0]

def function_minimize_Ekl(thetaij):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_Ekl])
	expression = np.asarray([1 + theta_E * alphas[i][0] + thetaij * alphas[i][0]* alphas[i][1] for i in pointeur_Ekl])
	return np.sum(objectif - expression)

theta_cherche_Ekl = scipy.optimize.least_squares(function_minimize_Ekl,1)
theta_Ekl = theta_cherche_Ekl.x[0]

def function_minimize_Ekb(thetaij):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_Ekb])
	expression = np.asarray([1 + theta_E * alphas[i][0] + thetaij * alphas[i][0]* alphas[i][1] for i in pointeur_Ekb])
	return np.sum(objectif - expression)

theta_cherche_Ekb = scipy.optimize.least_squares(function_minimize_Ekb,1)
theta_Ekb = theta_cherche_Ekb.x[0]


def function_minimize_Ealpha(thetaij):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_Ealpha[1:]])
	expression = np.asarray([1 + theta_E * alphas[i][0] + thetaij * alphas[i][0]* alphas[i][1] for i in pointeur_Ealpha[1:]])
	return np.sum(objectif - expression)

theta_cherche_Ealpha = scipy.optimize.least_squares(function_minimize_Ealpha,1)
# Removed first one because code was blocked
theta_Ealpha = theta_cherche_Ealpha.x[0]

def function_minimize_EK(thetaij):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_EK])
	expression = np.asarray([1 + theta_E * alphas[i][0] + thetaij * alphas[i][0]* alphas[i][1] for i in pointeur_EK])
	return np.sum(objectif - expression)

theta_cherche_EK = scipy.optimize.least_squares(function_minimize_EK,1)
theta_EK = theta_cherche_EK.x[0]

def function_minimize_klkb(thetaij):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_klkb])
	expression = np.asarray([1 + theta_kl * alphas[i][0] + thetaij * alphas[i][0]* alphas[i][1] for i in pointeur_klkb])
	return np.sum(objectif - expression)

theta_cherche_klkb = scipy.optimize.least_squares(function_minimize_klkb,1)
theta_klkb = theta_cherche_klkb.x[0]

def function_minimize_klalpha(thetaij):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_klalpha])
	expression = np.asarray([1 + theta_kl * alphas[i][0] + thetaij * alphas[i][0]* alphas[i][1] for i in pointeur_klalpha])
	return np.sum(objectif - expression)

theta_cherche_klalpha = scipy.optimize.least_squares(function_minimize_klalpha,1)
theta_klalpha = theta_cherche_klalpha.x[0]

def function_minimize_klK(thetaij):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_klK])
	expression = np.asarray([1 + theta_kl * alphas[i][0] + thetaij * alphas[i][0]* alphas[i][1] for i in pointeur_klK])
	return np.sum(objectif - expression)

theta_cherche_klK = scipy.optimize.least_squares(function_minimize_klK,1)
theta_klK = theta_cherche_klK.x[0]

def function_minimize_kbalpha(thetaij):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_kbalpha[1:]])
	expression = np.asarray([1 + theta_kb * alphas[i][0] + thetaij * alphas[i][0]* alphas[i][1] for i in pointeur_kbalpha[1:]])
	return np.sum(objectif - expression)

theta_cherche_kbalpha = scipy.optimize.least_squares(function_minimize_kbalpha,1)
theta_kbalpha = theta_cherche_kbalpha.x[0]

def function_minimize_kbK(thetaij):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_kbK])
	expression = np.asarray([1 + theta_kb * alphas[i][0] + thetaij * alphas[i][0]* alphas[i][1] for i in pointeur_kbK])
	return np.sum(objectif - expression)

theta_cherche_kbK = scipy.optimize.least_squares(function_minimize_kbK,1)
theta_kbK = theta_cherche_kbK.x[0]

def function_minimize_alphaK(thetaij):
	objectif = np.asarray([all_difference_sheets[i] for i in pointeur_alphaK[1:]])
	expression = np.asarray([1 + theta_alpha * alphas[i][0] + thetaij * alphas[i][0]* alphas[i][1] for i in pointeur_alphaK[1:]])
	return np.sum(objectif - expression)

theta_cherche_alphaK = scipy.optimize.least_squares(function_minimize_alphaK,1)
theta_alphaK = theta_cherche_alphaK.x[0]


theta_i = [theta_E, theta_kl, theta_kb, theta_alpha, theta_K]
print("[theta_E, theta_kl, theta_kb, theta_alpha, theta_K]")
print(theta_i)

theta_ij = [theta_Ekl, theta_Ekb, theta_Ealpha, theta_EK, theta_klkb, theta_klalpha, theta_klK, theta_kbalpha, theta_kbK, theta_alphaK]
print("[theta_Ekl, theta_Ekb, theta_Ealpha, theta_EK, theta_klkb, theta_klalpha, theta_klK, theta_kbalpha, theta_kbK, theta_alphaK]")
print(theta_ij)

def flatten(nested_list):
    flattened = []
    for item in nested_list:
        if isinstance(item, list):
            flattened.extend(flatten(item))
        else:
            flattened.append(item)
    return flattened

print("min%",100*min(flatten(alphas)))
print("max%",100*max(flatten(alphas)))
print("len alphas", len(alphas))

sumsquare = np.sum([x**2 for x in theta_i]) + np.sum([x**2 for x in theta_ij])

Si = [x**2/sumsquare for x in theta_i]
Sij = [x**2/sumsquare for x in theta_ij]

print("[S_E, S_kl, S_kb, S_alpha, S_K]")
print(Si)
print("[S_Ekl, S_Ekb, S_Ealpha, S_EK, S_klkb, S_klalpha, S_klK, S_kbalpha, S_kbK, S_alphaK]")
print(Sij)

print("100*[S_E, S_kl, S_kb, S_alpha, S_K]")
print([100*x for x in Si])
print("100*[S_Ekl, S_Ekb, S_Ealpha, S_EK, S_klkb, S_klalpha, S_klK, S_kbalpha, S_kbK, S_alphaK]")
print([100 * x for x in Sij])

print("First-order indices")
sumsquare = np.sum([x**2 for x in theta_i])
Si = [100*x**2/sumsquare for x in theta_i]
print("100*[S_E, S_kl, S_kb, S_alpha, S_K]")
print([x for x in Si])