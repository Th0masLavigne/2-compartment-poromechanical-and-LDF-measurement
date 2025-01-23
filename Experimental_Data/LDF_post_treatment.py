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
file_path = "JOL_PREVAIL_LBM_008.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")

# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
JOL_PREVAIL_LBM_008 = {}
for sheet_name, data in all_sheets.items():
    JOL_PREVAIL_LBM_008[sheet_name] = {col: data[col] for col in column_names if col in data.columns}


index_11 = 0
for ii in range(len(JOL_PREVAIL_LBM_008["j1_1"]["temps (s)"] ) ):
	if JOL_PREVAIL_LBM_008["j1_1"]["temps (s)"][ii] == 0:
		index_11 = ii 
		break
index_12 = 0
for ii in range(len(JOL_PREVAIL_LBM_008["j1_2"]["temps (s)"] ) ):
	if JOL_PREVAIL_LBM_008["j1_2"]["temps (s)"][ii] == 0:
		index_12 = ii 
		break
index_21 = 0
for ii in range(len(JOL_PREVAIL_LBM_008["j2_1"]["temps (s)"] ) ):
	if JOL_PREVAIL_LBM_008["j2_1"]["temps (s)"][ii] == 0:
		index_21 = ii 
		break
index_22 = 0
for ii in range(len(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"] ) ):
	if JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][ii] == 0:
		index_22 = ii 
		break		

length_list = min(len(JOL_PREVAIL_LBM_008["j1_1"]["temps (s)"][index_11:]) , len(JOL_PREVAIL_LBM_008["j1_2"]["temps (s)"][index_12:]) , len(JOL_PREVAIL_LBM_008["j2_1"]["temps (s)"][index_21:]) , len(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22:] ) )

JOL_Temperaturemean_signal=[]
JOL_Temperaturemps = []
JOL_Temperaturemms = []
JOL_Temperaturestd_signal = []
JOL_mean_signal=[]
JOL_mps = []
JOL_mms = []
JOL_std_signal = []

for i in range(length_list):
	JOL_mean_signal.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
									   JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
									   JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
									   JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i] ] ) )
	# 
	JOL_std_signal.append( np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
									 JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
									 JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
									 JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i] ] ) )
	# 
	JOL_mms.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
							   JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
							   JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
							   JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i] ] ) 
					- 1.96 * np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
								JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
								JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
								JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i] ] ) )
	# 
	JOL_mps.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
								JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
								JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
								JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i] ] )
					+ 1.96* np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
										JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
										JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
										JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i] ] ) )
	# 
	JOL_Temperaturemean_signal.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],
												  JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],
												  JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],
												  JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i] ] ) )
	# 
	JOL_Temperaturestd_signal.append( np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],
												JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],
												JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],
												JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i] ] ) )
	# 
	JOL_Temperaturemms.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],
										  JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],
										  JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],
										  JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i] ] )
								- 1.96* np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],
									JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],
									JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],
									JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i] ] ) )
	# 
	JOL_Temperaturemps.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],
										  JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],
										  JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],
										  JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i] ] ) 
								+ 1.96 * np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],
													JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],
													JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],
													JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i] ] ) )

# Filter parameters
fs = 1000  # Sampling frequency in Hz
cutoff = 25  # Cutoff frequency in Hz
JOL_filtered_signal = butter_lowpass_filter(JOL_mean_signal, cutoff, fs)


JOL_Temperaturemean_signal_raw=[]
JOL_Temperaturemps_raw = []
JOL_Temperaturemms_raw = []
JOL_Temperaturestd_signal_raw = []
JOL_mean_signal_raw=[]
JOL_mps_raw = []
JOL_mms_raw = []
JOL_std_signal_raw = []

for i in range(length_list+130):
	JOL_mean_signal_raw.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],
											JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],
											JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],
											JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130] ] ) )
	# 
	JOL_std_signal_raw.append( np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],
										JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],
										JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],
										JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130] ] ) )
	# 
	JOL_mms_raw.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],
									JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],
									JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],
									JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130] ] ) 
						- np.std( [JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],
										JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],
										JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],
										JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130] ] ) )
	# 
	JOL_mps_raw.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],
									JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],
									JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],
									JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130] ] )
						+ np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],
											JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],
											JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],
											JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130] ] ) )
	# 
	JOL_Temperaturemean_signal_raw.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-130],
													JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-130],
													JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-130],
													JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-130] ] ) )
	# 
	JOL_Temperaturestd_signal_raw.append( np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-130],
													JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-130],
													JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-130],
													JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-130] ] ) )
	# 
	JOL_Temperaturemms_raw.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-130],
											  JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-130],
											  JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-130],
											  JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-130] ] )
									-np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-130],
										JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-130],
										JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-130],
										JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-130] ] ) )
	# 
	JOL_Temperaturemps_raw.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-130],
												JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-130],
												JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-130],
												JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-130] ] )
									+np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-130],
														JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-130],
														JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-130],
														JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-130] ] ) )


plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_mean_signal_raw, color='b', label="P$_1$", alpha=0.5)
plt.fill_between(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_mms_raw, JOL_mps_raw, color='b', alpha=0.2)
# 
# plt.plot(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_mean_signal_raw, color='r', label="P$_2$", alpha=0.5)
# plt.fill_between(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_mms_raw, ASE_mps_raw, color='r', alpha=0.2)
# # 
# plt.plot(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_mean_signal_raw, color='g', label="P$_3$", alpha=0.5)
# plt.fill_between(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_mms_raw, MDI_mps_raw, color='g', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 250])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('superp_mean_std_raw_P1.jpg', bbox_inches='tight')


plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_Temperaturemean_signal_raw, color='b', label="P$_1$", alpha=0.5)
plt.fill_between(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_Temperaturemms_raw, JOL_Temperaturemps_raw, color='b', alpha=0.2)
# plt.plot(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_Temperaturemean_signal_raw, color='r', label="P$_2$", alpha=0.5)
# plt.fill_between(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_Temperaturemms_raw, ASE_Temperaturemps_raw, color='r', alpha=0.2)
# # 
# plt.plot(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_Temperaturemean_signal_raw, color='g', label="P$_3$", alpha=0.5)
# plt.fill_between(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_Temperaturemms_raw, MDI_Temperaturemps_raw, color='g', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('superp_mean_std_temp_raw_P1.jpg', bbox_inches='tight')
















# Load the Excel file
file_path = "ASE_PREVAIL_LBM_009.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
ASE_PREVAIL_LBM_009 = {}
for sheet_name, data in all_sheets.items():
    ASE_PREVAIL_LBM_009[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11_ASE = 0
for ii in range(len(ASE_PREVAIL_LBM_009["j1_1"]["temps (s)"] ) ):
	if ASE_PREVAIL_LBM_009["j1_1"]["temps (s)"][ii] == 0:
		index_11_ASE = ii 
		break
index_12_ASE = 0
for ii in range(len(ASE_PREVAIL_LBM_009["j1_2"]["temps (s)"] ) ):
	if ASE_PREVAIL_LBM_009["j1_2"]["temps (s)"][ii] == 0:
		index_12_ASE = ii 
		break
index_21_ASE = 0
for ii in range(len(ASE_PREVAIL_LBM_009["j2_1"]["temps (s)"] ) ):
	if ASE_PREVAIL_LBM_009["j2_1"]["temps (s)"][ii] == 0:
		index_21_ASE = ii 
		break
index_22_ASE = 0
for ii in range(len(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"] ) ):
	if ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][ii] == 0:
		index_22_ASE = ii 
		break		

length_list_ASE = min(len(ASE_PREVAIL_LBM_009["j1_1"]["temps (s)"][index_11_ASE:]) , len(ASE_PREVAIL_LBM_009["j1_2"]["temps (s)"][index_12_ASE:]) , len(ASE_PREVAIL_LBM_009["j2_1"]["temps (s)"][index_21_ASE:]) , len(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE:] ) )

ASE_Temperaturemean_signal=[]
ASE_Temperaturemps = []
ASE_Temperaturemms = []
ASE_Temperaturestd_signal = []
ASE_mean_signal=[]
ASE_mps = []
ASE_mms = []
ASE_std_signal = []

for i in range(length_list_ASE):
	ASE_mean_signal.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],
										ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],
										ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],
										ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i] ] ) )
	# 
	ASE_std_signal.append( np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],
										ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],
										ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],
										ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i] ] ) )
	# 
	ASE_mms.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],
								ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],
								ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],
								ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i] ] )
					- 1.96* np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],
								ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],
								ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],
								ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i] ] ) )
	# 
	ASE_mps.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],
								ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],
								ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],
								ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i] ] )
					+1.96*np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],
								ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],
								ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],
								ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i] ] ) )
	# 
	ASE_Temperaturemean_signal.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i],
													ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i],
													ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i],
													ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i] ] ) )
	# 
	ASE_Temperaturestd_signal.append( np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i],
												ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i],
												ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i],
												ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i] ] ) )
	# 
	ASE_Temperaturemms.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i],
											ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i],
											ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i],
											ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i] ] )
								-1.96*np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i],
											ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i],
											ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i],
											ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i] ] ) )
							# 
	ASE_Temperaturemps.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i],
											ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i],
											ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i],
											ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i] ] )
							+1.96*np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i],
											ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i],
											ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i],
											ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i] ] ) )

# Filter parameters
fs = 1000  # Sampling frequency in Hz
cutoff = 25  # Cutoff frequency in Hz
ASE_filtered_signal = butter_lowpass_filter(ASE_mean_signal, cutoff, fs)



ASE_Temperaturemean_signal_raw=[]
ASE_Temperaturemps_raw = []
ASE_Temperaturemms_raw = []
ASE_Temperaturestd_signal_raw = []
ASE_mean_signal_raw=[]
ASE_mps_raw = []
ASE_mms_raw = []
ASE_std_signal_raw = []

for i in range(length_list_ASE+130):
	ASE_mean_signal_raw.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130] ] ) )
	# 
	ASE_std_signal_raw.append( np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130] ] ) )
	# 
	ASE_mms_raw.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130] ] )
	-np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130] ] ) )
	# 
	ASE_mps_raw.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130] ] )
	+np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130] ] ) )
	# 
	ASE_Temperaturemean_signal_raw.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-130],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-130] ] ) )
	# 
	ASE_Temperaturestd_signal_raw.append( np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-130],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-130] ] ) )
	# 
	ASE_Temperaturemms_raw.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-130],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-130] ] )
	-np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-130],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-130] ] ) )
	# 
	ASE_Temperaturemps_raw.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-130],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-130] ] )
	+np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-130],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-130] ] ) )





plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_mean_signal_raw, color='b', label="P$_1$", alpha=0.5)
# plt.fill_between(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_mms_raw, JOL_mps_raw, color='b', alpha=0.2)
# 
plt.plot(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_mean_signal_raw, color='r', label="P$_2$", alpha=0.5)
plt.fill_between(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_mms_raw, ASE_mps_raw, color='r', alpha=0.2)
# # 
# plt.plot(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_mean_signal_raw, color='g', label="P$_3$", alpha=0.5)
# plt.fill_between(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_mms_raw, MDI_mps_raw, color='g', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 250])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('superp_mean_std_raw_P2.jpg', bbox_inches='tight')


plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_Temperaturemean_signal_raw, color='b', label="P$_1$", alpha=0.5)
# plt.fill_between(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_Temperaturemms_raw, JOL_Temperaturemps_raw, color='b', alpha=0.2)
plt.plot(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_Temperaturemean_signal_raw, color='r', label="P$_2$", alpha=0.5)
plt.fill_between(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_Temperaturemms_raw, ASE_Temperaturemps_raw, color='r', alpha=0.2)
# # 
# plt.plot(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_Temperaturemean_signal_raw, color='g', label="P$_3$", alpha=0.5)
# plt.fill_between(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_Temperaturemms_raw, MDI_Temperaturemps_raw, color='g', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('superp_mean_std_temp_raw_P2.jpg', bbox_inches='tight')










# Load the Excel file
file_path = "MDI_PREVAIL_LBM_026.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
MDI_PREVAIL_LBM_026 = {}
for sheet_name, data in all_sheets.items():
    MDI_PREVAIL_LBM_026[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11_MDI = 0
for ii in range(len(MDI_PREVAIL_LBM_026["j1_1"]["temps (s)"] ) ):
	if MDI_PREVAIL_LBM_026["j1_1"]["temps (s)"][ii] == 0:
		index_11_MDI = ii 
		break
index_12_MDI = 0
for ii in range(len(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"] ) ):
	if MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][ii] == 0:
		index_12_MDI = ii 
		break
index_21_MDI = 0
for ii in range(len(MDI_PREVAIL_LBM_026["j2_1"]["temps (s)"] ) ):
	if MDI_PREVAIL_LBM_026["j2_1"]["temps (s)"][ii] == 0:
		index_21_MDI = ii 
		break
index_22_MDI = 0
for ii in range(len(MDI_PREVAIL_LBM_026["j2_2"]["temps (s)"] ) ):
	if MDI_PREVAIL_LBM_026["j2_2"]["temps (s)"][ii] == 0:
		index_22_MDI = ii 
		break		

length_list_MDI = min(len(MDI_PREVAIL_LBM_026["j1_1"]["temps (s)"][index_11_MDI:]) , len(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI:]) )# , len(MDI_PREVAIL_LBM_026["j2_1"]["temps (s)"][index_21_MDI:]) , len(MDI_PREVAIL_LBM_026["j2_2"]["temps (s)"][index_22_MDI:] ) )

MDI_Temperaturemean_signal=[]
MDI_Temperaturemps = []
MDI_Temperaturemms = []
MDI_Temperaturestd_signal = []
MDI_mean_signal=[]
MDI_mps = []
MDI_mms = []
MDI_std_signal = []

for i in range(length_list_MDI):
	MDI_mean_signal.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22_MDI+i] ] ) )
	# 
	MDI_std_signal.append( np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22_MDI+i] ] ) )
	# 
	MDI_mms.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22_MDI+i] ] )
	- 1.96* np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22_MDI+i] ] ) )
	# 
	MDI_mps.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22_MDI+i] ] )
	+1.96*np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22_MDI+i] ] ) )
	# 
	MDI_Temperaturemean_signal.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i] ] ) )
	# 
	MDI_Temperaturestd_signal.append( np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i] ] ) )
	# 
	MDI_Temperaturemms.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i] ] )
	-1.96*np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i] ] ) )
	# 
	MDI_Temperaturemps.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i] ] )
	+1.96*np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i] ] ) )
	
MDI_Temperaturemean_signal_raw=[]
MDI_Temperaturemps_raw = []
MDI_Temperaturemms_raw = []
MDI_Temperaturestd_signal_raw = []
MDI_mean_signal_raw=[]
MDI_mps_raw = []
MDI_mms_raw = []
MDI_std_signal_raw = []

for i in range(length_list_MDI+130):
	MDI_mean_signal_raw.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-130],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-130] ] ) )
	# 
	MDI_std_signal_raw.append( np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-130],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-130] ] ) )
	# 
	MDI_mms_raw.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-130],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-130] ] )
	-np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-130],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-130] ] ) )
	# 
	MDI_mps_raw.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-130],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-130] ,
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-130] ] )
	+np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-130],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-130] ] ) )
	# 
	MDI_Temperaturemean_signal_raw.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-130],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-130] ] ) )
	# 
	MDI_Temperaturestd_signal_raw.append( np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-130],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-130] ] ) )
	# 
	MDI_Temperaturemms_raw.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-130],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-130] ] )
	-np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-130],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-130] ] ) )
	# 
	MDI_Temperaturemps_raw.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-130],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-130] ] )
	+np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-130],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-130] ] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_mean_signal_raw, color='b', label="P$_1$", alpha=0.5)
# plt.fill_between(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_mms_raw, JOL_mps_raw, color='b', alpha=0.2)
# # 
# plt.plot(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_mean_signal_raw, color='r', label="P$_2$", alpha=0.5)
# plt.fill_between(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_mms_raw, ASE_mps_raw, color='r', alpha=0.2)
# # 
plt.plot(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_mean_signal_raw, color='g', label="P$_3$", alpha=0.5)
plt.fill_between(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_mms_raw, MDI_mps_raw, color='g', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 250])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('superp_mean_std_raw_P3.jpg', bbox_inches='tight')




plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_Temperaturemean_signal_raw, color='b', label="P$_1$", alpha=0.5)
# plt.fill_between(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_Temperaturemms_raw, JOL_Temperaturemps_raw, color='b', alpha=0.2)
# plt.plot(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_Temperaturemean_signal_raw, color='r', label="P$_2$", alpha=0.5)
# plt.fill_between(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_Temperaturemms_raw, ASE_Temperaturemps_raw, color='r', alpha=0.2)
# # 
plt.plot(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_Temperaturemean_signal_raw, color='g', label="P$_3$", alpha=0.5)
plt.fill_between(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_Temperaturemms_raw, MDI_Temperaturemps_raw, color='g', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('superp_mean_std_temp_raw_P3.jpg', bbox_inches='tight')






















# Raw

# A compléter

all_means=[]
all_std=[]
for i in range(min(length_list_ASE+130,length_list+130,length_list_MDI+130)):
	all_means.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130],
		JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],
		JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],
		JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],
		JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130],
		MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-130],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-130] ] ) )
	# 
	all_std.append( np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-130],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-130],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-130],
		JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-130],
		JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-130],
		JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-130],
		JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-130],
		MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-130],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-130],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-130] ] ) )

all_mms = [all_means[i]-all_std[i] for i in range(len(all_means))]
all_mps = [all_means[i]+all_std[i] for i in range(len(all_means))]

all_means_pc=[]
all_std_pc=[]
# Percentage
for i in range(min(length_list_ASE,length_list)):
	all_means_pc.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i],JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i],MDI_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11_MDI+i],MDI_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12_MDI+i],MDI_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21_MDI+i],MDI_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22_MDI+i]                                  ] ) )
	all_std_pc.append( np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i],JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i],MDI_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11_MDI+i],MDI_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12_MDI+i],MDI_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21_MDI+i],MDI_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22_MDI+i]                                     ] ) )

# 68% confidence
all_mms_pc = [all_means_pc[i] - 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]
all_mps_pc = [all_means_pc[i] + 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]

# # 95% Confidence
# all_mms_pc = [all_means_pc[i] - 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]
# all_mps_pc = [all_means_pc[i] + 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]




initial_baseline = np.mean( all_means[:130] )
std_initial = np.mean( all_std[:1130] )
print(initial_baseline, '$\\pm$', std_initial)



first_occlusions = np.mean( [ all_means[60+130:160+130],all_means[350+130:450+130] ] )
std_fo = np.mean( [ np.mean(all_std[60+130:160+130]),np.mean(all_std[60+130:150+130]) ] )
print(first_occlusions, '$\\pm$', std_fo)

last_occlusions = np.mean( [ all_means[720+130:820+130],all_means[1125+130:1225+130] ] )
std_lo = np.mean( [ np.mean(all_std[720+130:820+130]),np.mean(all_std[1125+130:1225+130]) ] )
print(last_occlusions, '$\\pm$', std_lo)



idex0 = 1300
idex = 1460
t = np.linspace(idex0,idex,idex-idex0)

plt.figure(figsize=(10, 6))
plt.plot(t, ASE_mean_signal_raw[idex0:idex], label="Original Signal", alpha=0.5)
plt.savefig('check_time_range_metrics.png')






plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_mean_signal_raw, color='b', label="P$_1$", alpha=0.5, linewidth=3)
# plt.fill_between(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_mms_raw, JOL_mps_raw, color='b', alpha=0.2)
# 
plt.plot(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_mean_signal_raw, color='r', label="P$_2$", alpha=0.5, linewidth=3)
# plt.fill_between(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_mms_raw, ASE_mps_raw, color='r', alpha=0.2)
# 
plt.plot(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_mean_signal_raw, color='g', label="P$_3$", alpha=0.5, linewidth=3)
# plt.fill_between(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_mms_raw, MDI_mps_raw, color='g', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 250])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('superp_mean_std_raw.jpg', bbox_inches='tight')

plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_Temperaturemean_signal_raw, color='b', label="P$_1$", alpha=0.5, linewidth=3)
# plt.fill_between(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-130:index_22+length_list], JOL_Temperaturemms_raw, JOL_Temperaturemps_raw, color='b', alpha=0.2)
plt.plot(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_Temperaturemean_signal_raw, color='r', label="P$_2$", alpha=0.5, linewidth=3)
# plt.fill_between(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-130:index_22_ASE+length_list_ASE], ASE_Temperaturemms_raw, ASE_Temperaturemps_raw, color='r', alpha=0.2)
# 
plt.plot(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_Temperaturemean_signal_raw, color='g', label="P$_3$", alpha=0.5, linewidth=3)
# plt.fill_between(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-130:index_12_MDI+length_list_MDI], MDI_Temperaturemms_raw, MDI_Temperaturemps_raw, color='g', alpha=0.2)
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlim([-60, 780])
plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('superp_mean_std_temp_raw.jpg', bbox_inches='tight')




# # 

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
# 
plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE:index_22_ASE+length_list_ASE], ASE_mean_signal, color='r', label="P$_2$", alpha=0.6)
plt.fill_between(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE:index_22_ASE+length_list_ASE], ASE_mms, ASE_mps, color='r', alpha=0.2)
# 
plt.plot(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI:index_12_MDI+length_list_MDI], MDI_mean_signal , linestyle='-', color='g', label="P$_3$", alpha=0.6)
plt.fill_between(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI:index_12_MDI+length_list_MDI], MDI_mms, MDI_mps , linestyle=':', color='g', alpha=0.2)
# plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22:index_22+length_list], all_means_pc, color='k', label="Exp", alpha=0.5, linewidth=3)
# plt.fill_between(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22:index_22+length_list], all_mms_pc, all_mps_pc, color='k', alpha=0.2)
# plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22:index_22+length_list], all_means_pc, color='k', label="Exp", alpha=0.5, linewidth=3)
# plt.fill_between(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22:index_22+length_list], MDI_mms[:length_list], ASE_mps[:length_list], color='k', alpha=0.2)

# 
plt.plot(model["REF_export"]["time_all"][:774], model["REF_export"]["LDF_baseline_q"][6:780], color='k', linestyle='--', label="Set$_1$", linewidth=3)
plt.plot(model["ID_63_export"]["time_all"][:774], model["ID_63_export"]["LDF_baseline_q"][6:780], color='k', linestyle='-', label="Set$_2$", linewidth=3)
plt.plot(model["ID_66_export"]["time_all"][:774], model["ID_66_export"]["LDF_baseline_q"][6:780], color='k', linestyle=':', label="Set$_3$", linewidth=3)
# plt.legend()
plt.xlim([0, 750])
plt.ylim([0, 850])
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Time [s]")
plt.ylabel("LDF [%]")
plt.grid()
plt.savefig('all_mean_std.jpg', bbox_inches='tight')



