import pandas as pd
import numpy as np
from scipy.signal import butter, filtfilt
import matplotlib.pyplot as plt


# Load the Excel file
file_path = "Sensitivity.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["displacement_all",	"LDF_v_all",	"LDF_q_all",	"LDF_baseline_v",	"LDF_baseline_q",	"load_all",	"time_all"]


# Extract data from each sheet into a structured dictionary
model = {}
for sheet_name, data in all_sheets.items():
    model[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

plt.rcParams.update({'font.size': 25})


plt.figure(figsize=(16, 9))
plt.plot(model["ID_9_export"]["displacement_all"][:794], model["ID_9_export"]["LDF_baseline_q"][6:800], color='g', label="5e4 Pa", linewidth=2)
plt.plot(model["REF_export"]["displacement_all"][:794], model["REF_export"]["LDF_baseline_q"][6:800], color='r', label="2e5 Pa", linewidth=2)
plt.plot(model["ID_18_export"]["displacement_all"][:794], model["ID_18_export"]["LDF_baseline_q"][6:800], color='b', label="5e5 Pa", linewidth=2)
plt.plot(model["ID_52_export"]["displacement_all"][:794], model["ID_52_export"]["LDF_baseline_q"][6:800], color='k', label="1e6 Pa", linewidth=2)
plt.legend()
plt.legend(loc='upper left')
plt.xlabel("Displacement [m]")
plt.ylabel("LDF [%]")
plt.grid()
plt.title("$E$")
plt.ylim([0,220])
plt.tight_layout()
plt.savefig('BFU_E_sensitivity.jpg')


plt.figure(figsize=(16, 9))
plt.plot(model["ID_9_export"]["time_all"][:794], model["ID_9_export"]["LDF_baseline_q"][6:800], color='g', label="5e4 Pa", linewidth=2)
plt.plot(model["REF_export"]["time_all"][:794], model["REF_export"]["LDF_baseline_q"][6:800], color='r', label="2e5 Pa", linewidth=2)
plt.plot(model["ID_18_export"]["time_all"][:794], model["ID_18_export"]["LDF_baseline_q"][6:800], color='b', label="5e5 Pa", linewidth=2)
plt.plot(model["ID_52_export"]["time_all"][:794], model["ID_52_export"]["LDF_baseline_q"][6:800], color='k', label="1e6 Pa", linewidth=2)
plt.legend()
plt.legend(loc='upper left')
plt.xlabel("Time [s]")
plt.ylabel("LDF [%]")
plt.grid()
plt.title("$E$")
plt.xlim([0,800])
plt.ylim([0,220])
plt.tight_layout()
plt.savefig('BFT_E_sensitivity.jpg')




plt.figure(figsize=(16, 9))
plt.plot(model["ID_1_export"]["displacement_all"][:794], model["ID_1_export"]["LDF_baseline_q"][6:800], color='g', label="1e-15 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.plot(model["REF_export"]["displacement_all"][:794], model["REF_export"]["LDF_baseline_q"][6:800], color='r', label="1e-14 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.plot(model["ID_53_export"]["displacement_all"][:794], model["ID_53_export"]["LDF_baseline_q"][6:800], color='b', label="5e-14 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.plot(model["ID_2_export"]["displacement_all"][:794], model["ID_2_export"]["LDF_baseline_q"][6:800], color='k', label="1e-13 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.legend()
plt.legend(loc='upper left')
plt.xlabel("Displacement [m]")
plt.ylabel("LDF [%]")
plt.grid()
plt.title("$L^l$")
plt.xlim([-0.0003,0.00005])
plt.ylim([0,220])
plt.tight_layout()
plt.savefig('BFU_Ll_sensitivity.jpg')


plt.figure(figsize=(16, 9))
plt.plot(model["ID_1_export"]["time_all"][:794], model["ID_1_export"]["LDF_baseline_q"][6:800], color='g', label="1e-15 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.plot(model["REF_export"]["time_all"][:794], model["REF_export"]["LDF_baseline_q"][6:800], color='r', label="1e-14 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.plot(model["ID_53_export"]["time_all"][:794], model["ID_53_export"]["LDF_baseline_q"][6:800], color='b', label="5e-14 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.plot(model["ID_2_export"]["time_all"][:794], model["ID_2_export"]["LDF_baseline_q"][6:800], color='k', label="1e6 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.legend()
plt.legend(loc='upper left')
plt.xlabel("Time [s]")
plt.ylabel("LDF [%]")
plt.grid()
plt.title("$L^l$")
plt.xlim([0,800])
plt.ylim([0,220])
plt.tight_layout()
plt.savefig('BFT_Ll_sensitivity.jpg')


plt.figure(figsize=(16, 9))
plt.plot(model["ID_56_export"]["displacement_all"][:794], model["ID_56_export"]["LDF_baseline_q"][6:800], color='g', label="5e-10 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.plot(model["REF_export"]["displacement_all"][:794], model["REF_export"]["LDF_baseline_q"][6:800], color='r', label="1e-9 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.plot(model["ID_55_export"]["displacement_all"][:794], model["ID_55_export"]["LDF_baseline_q"][6:800], color='b', label="5e-9 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.plot(model["ID_4_export"]["displacement_all"][:794], model["ID_4_export"]["LDF_baseline_q"][6:800], color='k', label="1e-8 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.legend()
plt.legend(loc='upper left')
plt.xlabel("Displacement [m]")
plt.ylabel("LDF [%]")
plt.grid()
plt.title("$L^b$")
plt.xlim([-0.0003,0.00005])
plt.ylim([0,220])
plt.tight_layout()
plt.savefig('BFU_Lb_sensitivity.jpg')


plt.figure(figsize=(16, 9))
plt.plot(model["ID_56_export"]["time_all"][:794], model["ID_56_export"]["LDF_baseline_q"][6:800], color='g', label="1e-15 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.plot(model["REF_export"]["time_all"][:794], model["REF_export"]["LDF_baseline_q"][6:800], color='r', label="2e5 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.plot(model["ID_55_export"]["time_all"][:794], model["ID_55_export"]["LDF_baseline_q"][6:800], color='b', label="5e5 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.plot(model["ID_4_export"]["time_all"][:794], model["ID_4_export"]["LDF_baseline_q"][6:800], color='k', label="1e6 $m^2.Pa^{-1}.s^{-1}$", linewidth=2)
plt.legend()
plt.legend(loc='upper left')
plt.xlabel("Time [s]")
plt.ylabel("LDF [%]")
plt.grid()
plt.title("$L^b$")
plt.xlim([0,800])
plt.ylim([0,220])
plt.tight_layout()
plt.savefig('BFT_Lb_sensitivity.jpg')


plt.figure(figsize=(16, 9))
plt.plot(model["ID_5_export"]["displacement_all"][:794], model["ID_5_export"]["LDF_baseline_q"][6:800], color='g', label="2", linewidth=2)
plt.plot(model["REF_export"]["displacement_all"][:794], model["REF_export"]["LDF_baseline_q"][6:800], color='r', label="3", linewidth=2)
plt.plot(model["ID_55_export"]["displacement_all"][:794], model["ID_55_export"]["LDF_baseline_q"][6:800], color='b', label="4", linewidth=2)
plt.plot(model["ID_58_export"]["displacement_all"][:794], model["ID_58_export"]["LDF_baseline_q"][6:800], color='k', label="6", linewidth=2)
plt.legend()
plt.legend(loc='upper left')
plt.xlabel("Displacement [m]")
plt.ylabel("LDF [%]")
plt.grid()
plt.title("$\\alpha$")
plt.xlim([-0.0003,0.00005])
plt.ylim([0,220])
plt.tight_layout()
plt.savefig('BFU_alpha_sensitivity.jpg')


plt.figure(figsize=(16, 9))
plt.plot(model["ID_5_export"]["time_all"][:794], model["ID_5_export"]["LDF_baseline_q"][6:800], color='g', label="2", linewidth=2)
plt.plot(model["REF_export"]["time_all"][:794], model["REF_export"]["LDF_baseline_q"][6:800], color='r', label="3", linewidth=2)
plt.plot(model["ID_6_export"]["time_all"][:794], model["ID_6_export"]["LDF_baseline_q"][6:800], color='b', label="4", linewidth=2)
plt.plot(model["ID_58_export"]["time_all"][:794], model["ID_58_export"]["LDF_baseline_q"][6:800], color='k', label="6", linewidth=2)
plt.legend()
plt.legend(loc='upper left')
plt.xlabel("Time [s]")
plt.ylabel("LDF [%]")
plt.grid()
plt.title("$\\alpha$")
plt.xlim([0,800])
plt.ylim([0,220])
plt.tight_layout()
plt.savefig('BFT_alpha_sensitivity.jpg')


plt.figure(figsize=(16, 9))
plt.plot(model["ID_8_export"]["displacement_all"][:794], model["ID_8_export"]["LDF_baseline_q"][6:800], color='g', label="500 Pa", linewidth=2)
plt.plot(model["REF_export"]["displacement_all"][:794], model["REF_export"]["LDF_baseline_q"][6:800], color='r', label="1000 Pa", linewidth=2)
plt.plot(model["ID_59_export"]["displacement_all"][:794], model["ID_59_export"]["LDF_baseline_q"][6:800], color='b', label="2000 Pa", linewidth=2)
plt.plot(model["ID_7_export"]["displacement_all"][:794], model["ID_7_export"]["LDF_baseline_q"][6:800], color='k', label="5000 Pa", linewidth=2)
plt.legend()
plt.legend(loc='upper left')
plt.xlabel("Displacement [m]")
plt.ylabel("LDF [%]")
plt.grid()
plt.title("$K$")
plt.xlim([-0.0003,0.00005])
plt.ylim([0,220])
plt.tight_layout()
plt.savefig('BFU_K_sensitivity.jpg')


plt.figure(figsize=(16, 9))
plt.plot(model["ID_8_export"]["time_all"][:794], model["ID_8_export"]["LDF_baseline_q"][6:800], color='g', label="500 Pa", linewidth=2)
plt.plot(model["REF_export"]["time_all"][:794], model["REF_export"]["LDF_baseline_q"][6:800], color='r', label="1000 Pa", linewidth=2)
plt.plot(model["ID_59_export"]["time_all"][:794], model["ID_59_export"]["LDF_baseline_q"][6:800], color='b', label="2000 Pa", linewidth=2)
plt.plot(model["ID_7_export"]["time_all"][:794], model["ID_7_export"]["LDF_baseline_q"][6:800], color='k', label="5000 Pa", linewidth=2)
plt.legend()
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.legend(loc='upper left')
plt.xlabel("Time [s]")
plt.ylabel("LDF [%]")
plt.grid()
plt.title("$K$")
plt.xlim([0,800])
plt.ylim([0,220])
plt.tight_layout()
plt.savefig('BFT_K.jpg')
