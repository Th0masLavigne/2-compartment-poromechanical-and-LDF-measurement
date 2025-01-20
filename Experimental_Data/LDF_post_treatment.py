import pandas as pd
import numpy as np
from scipy.signal import butter, filtfilt
import matplotlib.pyplot as plt

# Butterworth low-pass filter
def butter_lowpass_filter(data, cutoff, fs, order=5):
    """
    Apply a Butterworth low-pass filter to the signal.
    :param data: Input signal (numpy array)
    :param cutoff: Cutoff frequency (Hz)
    :param fs: Sampling frequency (Hz)
    :param order: Filter order (higher means sharper cutoff)
    :return: Filtered signal
    """
    nyquist = 0.5 * fs  # Nyquist frequency
    normal_cutoff = cutoff / nyquist  # Normalized cutoff frequency
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)  # Apply filter with zero phase shift
    return y

# Load the Excel file
file_path = "JOL_PREVAIL_LBM_008_for_plot.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")

# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
JOL_PREVAIL_LBM_008 = {}
for sheet_name, data in all_sheets.items():
    JOL_PREVAIL_LBM_008[sheet_name] = {col: data[col] for col in column_names if col in data.columns}


index_11 = 0
for ii in range(len(JOL_PREVAIL_LBM_008["j1_1"]["temps (s)"])):
	if JOL_PREVAIL_LBM_008["j1_1"]["temps (s)"][ii] == 0:
		index_11 = ii 
		break
index_12 = 0
for ii in range(len(JOL_PREVAIL_LBM_008["j1_2"]["temps (s)"])):
	if JOL_PREVAIL_LBM_008["j1_2"]["temps (s)"][ii] == 0:
		index_12 = ii 
		break
index_21 = 0
for ii in range(len(JOL_PREVAIL_LBM_008["j2_1"]["temps (s)"])):
	if JOL_PREVAIL_LBM_008["j2_1"]["temps (s)"][ii] == 0:
		index_21 = ii 
		break
index_22 = 0
for ii in range(len(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"])):
	if JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][ii] == 0:
		index_22 = ii 
		break		

length_list = min(len(JOL_PREVAIL_LBM_008["j1_1"]["temps (s)"][index_11:]) , len(JOL_PREVAIL_LBM_008["j1_2"]["temps (s)"][index_12:]) , len(JOL_PREVAIL_LBM_008["j2_1"]["temps (s)"][index_21:]) , len(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22:]))

JOL_Temperaturemean_signal=[]
JOL_Temperaturemps = []
JOL_Temperaturemms = []
JOL_Temperaturestd_signal = []
JOL_mean_signal=[]
JOL_mps = []
JOL_mms = []
JOL_std_signal = []

for i in range(length_list):
	JOL_mean_signal.append(np.mean([JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i]]))
	JOL_std_signal.append(np.std([JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i]]))
	JOL_mms.append(np.mean([JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i]])-np.std([JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i]]))
	JOL_mps.append(np.mean([JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i]])+np.std([JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i]]))
	JOL_Temperaturemean_signal.append(np.mean([JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i]]))
	JOL_Temperaturestd_signal.append(np.std([JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i]]))
	JOL_Temperaturemms.append(np.mean([JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i]])-np.std([JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i]]))
	JOL_Temperaturemps.append(np.mean([JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i]])+np.std([JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i]]))

# Filter parameters
fs = 1000  # Sampling frequency in Hz
cutoff = 25  # Cutoff frequency in Hz
JOL_filtered_signal = butter_lowpass_filter(JOL_mean_signal, cutoff, fs)

# Load the Excel file
file_path = "ASE_PREVAIL_LBM_009_for_plot.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
ASE_PREVAIL_LBM_009 = {}
for sheet_name, data in all_sheets.items():
    ASE_PREVAIL_LBM_009[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11_ASE = 0
for ii in range(len(ASE_PREVAIL_LBM_009["j1_1"]["temps (s)"])):
	if ASE_PREVAIL_LBM_009["j1_1"]["temps (s)"][ii] == 0:
		index_11_ASE = ii 
		break
index_12_ASE = 0
for ii in range(len(ASE_PREVAIL_LBM_009["j1_2"]["temps (s)"])):
	if ASE_PREVAIL_LBM_009["j1_2"]["temps (s)"][ii] == 0:
		index_12_ASE = ii 
		break
index_21_ASE = 0
for ii in range(len(ASE_PREVAIL_LBM_009["j2_1"]["temps (s)"])):
	if ASE_PREVAIL_LBM_009["j2_1"]["temps (s)"][ii] == 0:
		index_21_ASE = ii 
		break
index_22_ASE = 0
for ii in range(len(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"])):
	if ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][ii] == 0:
		index_22_ASE = ii 
		break		

length_list_ASE = min(len(ASE_PREVAIL_LBM_009["j1_1"]["temps (s)"][index_11_ASE:]) , len(ASE_PREVAIL_LBM_009["j1_2"]["temps (s)"][index_12_ASE:]) , len(ASE_PREVAIL_LBM_009["j2_1"]["temps (s)"][index_21_ASE:]) , len(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE:]))

ASE_Temperaturemean_signal=[]
ASE_Temperaturemps = []
ASE_Temperaturemms = []
ASE_Temperaturestd_signal = []
ASE_mean_signal=[]
ASE_mps = []
ASE_mms = []
ASE_std_signal = []

for i in range(length_list_ASE):
	ASE_mean_signal.append(np.mean([ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i]]))
	ASE_std_signal.append(np.std([ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i]]))
	ASE_mms.append(np.mean([ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i]])-np.std([ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i]]))
	ASE_mps.append(np.mean([ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i]])+np.std([ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i]]))
	ASE_Temperaturemean_signal.append(np.mean([ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i]]))
	ASE_Temperaturestd_signal.append(np.std([ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i]]))
	ASE_Temperaturemms.append(np.mean([ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i]])-np.std([ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i]]))
	ASE_Temperaturemps.append(np.mean([ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i]])+np.std([ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i]]))

# Filter parameters
fs = 1000  # Sampling frequency in Hz
cutoff = 25  # Cutoff frequency in Hz
ASE_filtered_signal = butter_lowpass_filter(ASE_mean_signal, cutoff, fs)

JOL_Temperaturemean_signal_raw=[]
JOL_Temperaturemps_raw = []
JOL_Temperaturemms_raw = []
JOL_Temperaturestd_signal_raw = []
JOL_mean_signal_raw=[]
JOL_mps_raw = []
JOL_mms_raw = []
JOL_std_signal_raw = []

for i in range(length_list+130):
	JOL_mean_signal_raw.append(np.mean([JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130]]))
	JOL_std_signal_raw.append(np.std([JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130]]))
	JOL_mms_raw.append(np.mean([JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130]])-np.std([JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130]]))
	JOL_mps_raw.append(np.mean([JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130]])+np.std([JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130]]))
	JOL_Temperaturemean_signal_raw.append(np.mean([JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-130],JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-130],JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-130],JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-130]]))
	JOL_Temperaturestd_signal_raw.append(np.std([JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-130],JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-130],JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-130],JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-130]]))
	JOL_Temperaturemms_raw.append(np.mean([JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-130],JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-130],JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-130],JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-130]])-np.std([JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-130],JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-130],JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-130],JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-130]]))
	JOL_Temperaturemps_raw.append(np.mean([JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-130],JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-130],JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-130],JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-130]])+np.std([JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-130],JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-130],JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-130],JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-130]]))

ASE_Temperaturemean_signal_raw=[]
ASE_Temperaturemps_raw = []
ASE_Temperaturemms_raw = []
ASE_Temperaturestd_signal_raw = []
ASE_mean_signal_raw=[]
ASE_mps_raw = []
ASE_mms_raw = []
ASE_std_signal_raw = []

for i in range(length_list_ASE+130):
	ASE_mean_signal_raw.append(np.mean([ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130]]))
	ASE_std_signal_raw.append(np.std([ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130]]))
	ASE_mms_raw.append(np.mean([ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130]])-np.std([ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130]]))
	ASE_mps_raw.append(np.mean([ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130]])+np.std([ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130]]))
	ASE_Temperaturemean_signal_raw.append(np.mean([ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-130],ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-130],ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-130],ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-130]]))
	ASE_Temperaturestd_signal_raw.append(np.std([ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-130],ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-130],ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-130],ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-130]]))
	ASE_Temperaturemms_raw.append(np.mean([ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-130],ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-130],ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-130],ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-130]])-np.std([ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-130],ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-130],ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-130],ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-130]]))
	ASE_Temperaturemps_raw.append(np.mean([ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-130],ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-130],ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-130],ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-130]])+np.std([ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-130],ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-130],ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-130],ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-130]]))


# Raw
all_means=[]
all_std=[]
for i in range(min(length_list_ASE+130,length_list+130)):
	all_means.append(np.mean([ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130],JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130]]))
	all_std.append(np.std([ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130],JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130]]))

all_mms = [all_means[i]-all_std[i] for i in range(len(all_means))]
all_mps = [all_means[i]+all_std[i] for i in range(len(all_means))]

all_means_pc=[]
all_std_pc=[]
# Percentage
for i in range(min(length_list_ASE,length_list)):
	all_means_pc.append(np.mean([ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i],JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i]]))
	all_std_pc.append(np.std([ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i],JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i]]))

all_mms_pc = [all_means_pc[i]-all_std_pc[i] for i in range(len(all_means_pc))]
all_mps_pc = [all_means_pc[i]+all_std_pc[i] for i in range(len(all_means_pc))]



first_occlusions = np.mean( [ all_means[170:270],all_means[520:620] ] )
first_occlusionsold = np.mean( [ np.mean(ASE_mean_signal_raw[175:270]),np.mean(JOL_mean_signal_raw[170:270]),np.mean(ASE_mean_signal_raw[500:620]) ,np.mean(JOL_mean_signal_raw[510:620]) ] )
print(first_occlusions,first_occlusionsold)

std_fo = np.std([np.mean(all_means[175:270]),np.mean(all_means[500:620]) ] )
std_foold = np.std([np.mean(ASE_mean_signal_raw[175:270]),np.mean(JOL_mean_signal_raw[170:270]),np.mean(ASE_mean_signal_raw[500:620]),np.mean(JOL_mean_signal_raw[510:620]) ] )
print(first_occlusions,std_fo,std_foold)

last_occlusions = np.mean( [ np.mean(ASE_mean_signal_raw[875:1000]),np.mean(JOL_mean_signal_raw[875:1000]),np.mean(ASE_mean_signal_raw[500:620]) ,np.mean(JOL_mean_signal_raw[1300:1460]) ] )
print(last_occlusions)

std_lo = np.std([np.mean(ASE_mean_signal_raw[875:1000]),np.mean(JOL_mean_signal_raw[875:1000]),np.mean(ASE_mean_signal_raw[500:620]),np.mean(JOL_mean_signal_raw[1300:1460]) ] )
print(last_occlusions,std_lo)

idex0 = 1300
idex = 1460
t = np.linspace(idex0,idex,idex-idex0)

plt.figure(figsize=(10, 6))
plt.plot(t, ASE_mean_signal_raw[idex0:idex], label="Original Signal", alpha=0.5)
plt.savefig('check_time_range_metrics.png')