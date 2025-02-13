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

for i in range(length_list+110):
	JOL_mean_signal_raw.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
											JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
											JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
											JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110] ] ) )
	# 
	JOL_std_signal_raw.append( np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
										JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
										JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
										JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110] ] ) )
	# 
	JOL_mms_raw.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
									JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
									JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
									JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110] ] ) 
						- np.std( [JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
										JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
										JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
										JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110] ] ) )
	# 
	JOL_mps_raw.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
									JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
									JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
									JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110] ] )
						+ np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
											JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
											JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
											JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110] ] ) )
	# 
	JOL_Temperaturemean_signal_raw.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
													JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
													JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
													JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110] ] ) )
	# 
	JOL_Temperaturestd_signal_raw.append( np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
													JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
													JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
													JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110] ] ) )
	# 
	JOL_Temperaturemms_raw.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
											  JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
											  JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
											  JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110] ] )
									-np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
										JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
										JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
										JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110] ] ) )
	# 
	JOL_Temperaturemps_raw.append( np.mean( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
												JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
												JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
												JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110] ] )
									+np.std( [ JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
														JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
														JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
														JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110] ] ) )


plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-110:index_22+length_list], JOL_mean_signal_raw, color='black', label="P$_1$", alpha=0.5)
plt.fill_between(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-110:index_22+length_list], JOL_mms_raw, JOL_mps_raw, color='black', alpha=0.2)
# 
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 250])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_raw_P1.jpg', bbox_inches='tight')


plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-110:index_22+length_list], JOL_Temperaturemean_signal_raw, color='black', label="P$_1$", alpha=0.5)
plt.fill_between(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-110:index_22+length_list], JOL_Temperaturemms_raw, JOL_Temperaturemps_raw, color='black', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P1.jpg', bbox_inches='tight')



plt.close()












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

for i in range(length_list_ASE+110):
	ASE_mean_signal_raw.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-110] ] ) )
	# 
	ASE_std_signal_raw.append( np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-110] ] ) )
	# 
	ASE_mms_raw.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-110] ] )
	-np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-110] ] ) )
	# 
	ASE_mps_raw.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-110] ] )
	+np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-110] ] ) )
	# 
	ASE_Temperaturemean_signal_raw.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-110] ] ) )
	# 
	ASE_Temperaturestd_signal_raw.append( np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-110] ] ) )
	# 
	ASE_Temperaturemms_raw.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-110] ] )
	-np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-110] ] ) )
	# 
	ASE_Temperaturemps_raw.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-110] ] )
	+np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-110] ] ) )





plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# 
plt.plot(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-110:index_22_ASE+length_list_ASE], ASE_mean_signal_raw, color='b', label="P$_2$", alpha=0.5)
plt.fill_between(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-110:index_22_ASE+length_list_ASE], ASE_mms_raw, ASE_mps_raw, color='b', alpha=0.2)
# # 
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 250])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_raw_P2.jpg', bbox_inches='tight')


plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-110:index_22_ASE+length_list_ASE], ASE_Temperaturemean_signal_raw, color='b', label="P$_2$", alpha=0.5)
plt.fill_between(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-110:index_22_ASE+length_list_ASE], ASE_Temperaturemms_raw, ASE_Temperaturemps_raw, color='b', alpha=0.2)
# # 
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P2.jpg', bbox_inches='tight')



plt.close()








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

length_list_MDI = min(len(MDI_PREVAIL_LBM_026["j1_1"]["temps (s)"][index_11_MDI:]) , len(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI:]) , len(MDI_PREVAIL_LBM_026["j2_1"]["temps (s)"][index_21_MDI:]) , len(MDI_PREVAIL_LBM_026["j2_2"]["temps (s)"][index_22_MDI:] ) )

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

for i in range(length_list_MDI+110):
	MDI_mean_signal_raw.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-110] ] ) )
	# 
	MDI_std_signal_raw.append( np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-110] ] ) )
	# 
	MDI_mms_raw.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-110] ] )
	-np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-110] ] ) )
	# 
	MDI_mps_raw.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-110] ,
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-110] ] )
	+np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-110] ] ) )
	# 
	MDI_Temperaturemean_signal_raw.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-110] ] ) )
	# 
	MDI_Temperaturestd_signal_raw.append( np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-110] ] ) )
	# 
	MDI_Temperaturemms_raw.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-110] ] )
	-np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-110] ] ) )
	# 
	MDI_Temperaturemps_raw.append( np.mean( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-110] ] )
	+np.std( [ MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-110] ] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-110:index_12_MDI+length_list_MDI], MDI_mean_signal_raw, color='r', label="P$_3$", alpha=0.5)
plt.fill_between(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-110:index_12_MDI+length_list_MDI], MDI_mms_raw, MDI_mps_raw, color='r', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 100])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_raw_P3.jpg', bbox_inches='tight')




plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-110:index_12_MDI+length_list_MDI], MDI_Temperaturemean_signal_raw, color='r', label="P$_3$", alpha=0.5)
plt.fill_between(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-110:index_12_MDI+length_list_MDI], MDI_Temperaturemms_raw, MDI_Temperaturemps_raw, color='r', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P3.jpg', bbox_inches='tight')





plt.close()









































# Load the Excel file
file_path = "MMA_PREVAIL_LBM_028.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
MMA_PREVAIL_LBM_028 = {}
for sheet_name, data in all_sheets.items():
    MMA_PREVAIL_LBM_028[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11_MMA = 0
for ii in range(len(MMA_PREVAIL_LBM_028["j1_1"]["temps (s)"] ) ):
	if MMA_PREVAIL_LBM_028["j1_1"]["temps (s)"][ii] == 0:
		index_11_MMA = ii 
		break
index_12_MMA = 0
for ii in range(len(MMA_PREVAIL_LBM_028["j1_2"]["temps (s)"] ) ):
	if MMA_PREVAIL_LBM_028["j1_2"]["temps (s)"][ii] == 0:
		index_12_MMA = ii 
		break
index_21_MMA = 0
for ii in range(len(MMA_PREVAIL_LBM_028["j2_1"]["temps (s)"] ) ):
	if MMA_PREVAIL_LBM_028["j2_1"]["temps (s)"][ii] == 0:
		index_21_MMA = ii 
		break
index_22_MMA = 0
for ii in range(len(MMA_PREVAIL_LBM_028["j2_2"]["temps (s)"] ) ):
	if MMA_PREVAIL_LBM_028["j2_2"]["temps (s)"][ii] == 0:
		index_22_MMA = ii 
		break		

length_list_MMA = min(len(MMA_PREVAIL_LBM_028["j1_1"]["temps (s)"][index_11_MMA:]) , len(MMA_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12_MMA:]) , len(MMA_PREVAIL_LBM_028["j2_1"]["temps (s)"][index_21_MMA:]) , len(MMA_PREVAIL_LBM_028["j2_2"]["temps (s)"][index_22_MMA:] ) )

MMA_Temperaturemean_signal=[]
MMA_Temperaturemps = []
MMA_Temperaturemms = []
MMA_Temperaturestd_signal = []
MMA_mean_signal=[]
MMA_mps = []
MMA_mms = []
MMA_std_signal = []

for i in range(length_list_MMA):
	MMA_mean_signal.append( np.mean( [ MMA_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22_MMA+i] 
		] ) )
	# 
	MMA_std_signal.append( np.std( [ MMA_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22_MMA+i] 
		] ) )
	# 
	MMA_mms.append( np.mean( [ MMA_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22_MMA+i] 
		] )
	- 1.96* np.std( [ MMA_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22_MMA+i] 
		] ) )
	# 
	MMA_mps.append( np.mean( [ MMA_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22_MMA+i] 
		] )
	+1.96*np.std( [ MMA_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22_MMA+i] 
		] ) )
	# 
	MMA_Temperaturemean_signal.append( np.mean( [ MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i] 
		] ) )
	# 
	MMA_Temperaturestd_signal.append( np.std( [ MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i] 
		] ) )
	# 
	MMA_Temperaturemms.append( np.mean( [ MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i] 
		] )
	-1.96*np.std( [ MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i] 
		] ) )
	# 
	MMA_Temperaturemps.append( np.mean( [ MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i] 
		] )
	+1.96*np.std( [ MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i] 
		] ) )
	
MMA_Temperaturemean_signal_raw=[]
MMA_Temperaturemps_raw = []
MMA_Temperaturemms_raw = []
MMA_Temperaturestd_signal_raw = []
MMA_mean_signal_raw=[]
MMA_mps_raw = []
MMA_mms_raw = []
MMA_std_signal_raw = []

for i in range(length_list_MMA+110):
	MMA_mean_signal_raw.append( np.mean( [ MMA_PREVAIL_LBM_028["j1_1"]["PU"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["PU"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["PU"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["PU"][index_22_MMA+i-110] 
		] ) )
	# 
	MMA_std_signal_raw.append( np.std( [ MMA_PREVAIL_LBM_028["j1_1"]["PU"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["PU"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["PU"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["PU"][index_22_MMA+i-110] 
		] ) )
	# 
	MMA_mms_raw.append( np.mean( [ MMA_PREVAIL_LBM_028["j1_1"]["PU"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["PU"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["PU"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["PU"][index_22_MMA+i-110] 
		] )
	-np.std( [ MMA_PREVAIL_LBM_028["j1_1"]["PU"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["PU"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["PU"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["PU"][index_22_MMA+i-110] 
		] ) )
	# 
	MMA_mps_raw.append( np.mean( [ MMA_PREVAIL_LBM_028["j1_1"]["PU"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["PU"][index_12_MMA+i-110] ,
		MMA_PREVAIL_LBM_028["j2_1"]["PU"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["PU"][index_22_MMA+i-110] 
		] )
	+np.std( [ MMA_PREVAIL_LBM_028["j1_1"]["PU"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["PU"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["PU"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["PU"][index_22_MMA+i-110] 
		] ) )
	# 
	MMA_Temperaturemean_signal_raw.append( np.mean( [ MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i-110] 
		] ) )
	# 
	MMA_Temperaturestd_signal_raw.append( np.std( [ MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i-110] 
		] ) )
	# 
	MMA_Temperaturemms_raw.append( np.mean( [ MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i-110] 
		] )
	-np.std( [ MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i-110] 
		] ) )
	# 
	MMA_Temperaturemps_raw.append( np.mean( [ MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i-110] 
		] )
	+np.std( [ MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(MMA_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12_MMA-110:index_12_MMA+length_list_MMA], MMA_mean_signal_raw,  linestyle='-',color='black', label="P$_4$", alpha=0.5)
plt.fill_between(MMA_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12_MMA-110:index_12_MMA+length_list_MMA], MMA_mms_raw, MMA_mps_raw, color='black', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 40])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_raw_P4.jpg', bbox_inches='tight')




plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(MMA_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12_MMA-110:index_12_MMA+length_list_MMA], MMA_Temperaturemean_signal_raw, linestyle='-', color='black', label="P$_4$", alpha=0.5)
plt.fill_between(MMA_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12_MMA-110:index_12_MMA+length_list_MMA], MMA_Temperaturemms_raw, MMA_Temperaturemps_raw, color='black', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P4.jpg', bbox_inches='tight')









plt.close()


























# Load the Excel file
file_path = "CBO_PREVAIL_LBM_029.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
CBO_PREVAIL_LBM_029 = {}
for sheet_name, data in all_sheets.items():
    CBO_PREVAIL_LBM_029[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11_CBO = 0
for ii in range(len(CBO_PREVAIL_LBM_029["j1_1"]["temps (s)"] ) ):
	if CBO_PREVAIL_LBM_029["j1_1"]["temps (s)"][ii] == 0:
		index_11_CBO = ii 
		break
index_12_CBO = 0
for ii in range(len(CBO_PREVAIL_LBM_029["j1_2"]["temps (s)"] ) ):
	if CBO_PREVAIL_LBM_029["j1_2"]["temps (s)"][ii] == 0:
		index_12_CBO = ii 
		break
index_21_CBO = 0
for ii in range(len(CBO_PREVAIL_LBM_029["j2_1"]["temps (s)"] ) ):
	if CBO_PREVAIL_LBM_029["j2_1"]["temps (s)"][ii] == 0:
		index_21_CBO = ii 
		break
index_22_CBO = 0
for ii in range(len(CBO_PREVAIL_LBM_029["j2_2"]["temps (s)"] ) ):
	if CBO_PREVAIL_LBM_029["j2_2"]["temps (s)"][ii] == 0:
		index_22_CBO = ii 
		break		

length_list_CBO = min(len(CBO_PREVAIL_LBM_029["j1_1"]["temps (s)"][index_11_CBO:]) , len(CBO_PREVAIL_LBM_029["j1_2"]["temps (s)"][index_12_CBO:]) , len(CBO_PREVAIL_LBM_029["j2_1"]["temps (s)"][index_21_CBO:]) , len(CBO_PREVAIL_LBM_029["j2_2"]["temps (s)"][index_22_CBO:] ) )

CBO_Temperaturemean_signal=[]
CBO_Temperaturemps = []
CBO_Temperaturemms = []
CBO_Temperaturestd_signal = []
CBO_mean_signal=[]
CBO_mps = []
CBO_mms = []
CBO_std_signal = []

for i in range(length_list_CBO):
	CBO_mean_signal.append( np.mean( [ CBO_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22_CBO+i] 
		] ) )
	# 
	CBO_std_signal.append( np.std( [ CBO_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22_CBO+i] 
		] ) )
	# 
	CBO_mms.append( np.mean( [ CBO_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22_CBO+i] 
		] )
	- 1.96* np.std( [ CBO_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22_CBO+i] 
		] ) )
	# 
	CBO_mps.append( np.mean( [ CBO_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22_CBO+i] 
		] )
	+1.96*np.std( [ CBO_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22_CBO+i] 
		] ) )
	# 
	CBO_Temperaturemean_signal.append( np.mean( [ CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i] 
		] ) )
	# 
	CBO_Temperaturestd_signal.append( np.std( [ CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i] 
		] ) )
	# 
	CBO_Temperaturemms.append( np.mean( [ CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i] 
		] )
	-1.96*np.std( [ CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i] 
		] ) )
	# 
	CBO_Temperaturemps.append( np.mean( [ CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i] 
		] )
	+1.96*np.std( [ CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i] 
		] ) )
	
CBO_Temperaturemean_signal_raw=[]
CBO_Temperaturemps_raw = []
CBO_Temperaturemms_raw = []
CBO_Temperaturestd_signal_raw = []
CBO_mean_signal_raw=[]
CBO_mps_raw = []
CBO_mms_raw = []
CBO_std_signal_raw = []

for i in range(length_list_CBO+110):
	CBO_mean_signal_raw.append( np.mean( [ CBO_PREVAIL_LBM_029["j1_1"]["PU"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["PU"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["PU"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["PU"][index_22_CBO+i-110] 
		] ) )
	# 
	CBO_std_signal_raw.append( np.std( [ CBO_PREVAIL_LBM_029["j1_1"]["PU"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["PU"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["PU"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["PU"][index_22_CBO+i-110] 
		] ) )
	# 
	CBO_mms_raw.append( np.mean( [ CBO_PREVAIL_LBM_029["j1_1"]["PU"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["PU"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["PU"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["PU"][index_22_CBO+i-110] 
		] )
	-np.std( [ CBO_PREVAIL_LBM_029["j1_1"]["PU"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["PU"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["PU"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["PU"][index_22_CBO+i-110] 
		] ) )
	# 
	CBO_mps_raw.append( np.mean( [ CBO_PREVAIL_LBM_029["j1_1"]["PU"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["PU"][index_12_CBO+i-110] ,
		CBO_PREVAIL_LBM_029["j2_1"]["PU"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["PU"][index_22_CBO+i-110] 
		] )
	+np.std( [ CBO_PREVAIL_LBM_029["j1_1"]["PU"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["PU"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["PU"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["PU"][index_22_CBO+i-110] 
		] ) )
	# 
	CBO_Temperaturemean_signal_raw.append( np.mean( [ CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i-110] 
		] ) )
	# 
	CBO_Temperaturestd_signal_raw.append( np.std( [ CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i-110] 
		] ) )
	# 
	CBO_Temperaturemms_raw.append( np.mean( [ CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i-110] 
		] )
	-np.std( [ CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i-110] 
		] ) )
	# 
	CBO_Temperaturemps_raw.append( np.mean( [ CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i-110] 
		] )
	+np.std( [ CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(CBO_PREVAIL_LBM_029["j1_2"]["temps (s)"][index_12_CBO-110:index_12_CBO+length_list_CBO], CBO_mean_signal_raw, color='darkgreen', label="P$_5$", alpha=0.5)
plt.fill_between(CBO_PREVAIL_LBM_029["j1_2"]["temps (s)"][index_12_CBO-110:index_12_CBO+length_list_CBO], CBO_mms_raw, CBO_mps_raw, color='darkgreen', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 60])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_raw_P5.jpg', bbox_inches='tight')




plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(CBO_PREVAIL_LBM_029["j1_2"]["temps (s)"][index_12_CBO-110:index_12_CBO+length_list_CBO], CBO_Temperaturemean_signal_raw, color='darkgreen', label="P$_5$", alpha=0.5)
plt.fill_between(CBO_PREVAIL_LBM_029["j1_2"]["temps (s)"][index_12_CBO-110:index_12_CBO+length_list_CBO], CBO_Temperaturemms_raw, CBO_Temperaturemps_raw, color='darkgreen', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P5.jpg', bbox_inches='tight')








plt.close()






































# Load the Excel file
file_path = "PLA_PREVAIL_LBM_030.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
PLA_PREVAIL_LBM_030 = {}
for sheet_name, data in all_sheets.items():
    PLA_PREVAIL_LBM_030[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11_PLA = 0
for ii in range(len(PLA_PREVAIL_LBM_030["j1_1"]["temps (s)"] ) ):
	if PLA_PREVAIL_LBM_030["j1_1"]["temps (s)"][ii] == 0:
		index_11_PLA = ii 
		break
index_12_PLA = 0
for ii in range(len(PLA_PREVAIL_LBM_030["j1_2"]["temps (s)"] ) ):
	if PLA_PREVAIL_LBM_030["j1_2"]["temps (s)"][ii] == 0:
		index_12_PLA = ii 
		break
index_21_PLA = 0
for ii in range(len(PLA_PREVAIL_LBM_030["j2_1"]["temps (s)"] ) ):
	if PLA_PREVAIL_LBM_030["j2_1"]["temps (s)"][ii] == 0:
		index_21_PLA = ii 
		break
index_22_PLA = 0
for ii in range(len(PLA_PREVAIL_LBM_030["j2_2"]["temps (s)"] ) ):
	if PLA_PREVAIL_LBM_030["j2_2"]["temps (s)"][ii] == 0:
		index_22_PLA = ii 
		break		

length_list_PLA = min(len(PLA_PREVAIL_LBM_030["j1_1"]["temps (s)"][index_11_PLA:]) , len(PLA_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12_PLA:]) , len(PLA_PREVAIL_LBM_030["j2_1"]["temps (s)"][index_21_PLA:]) , len(PLA_PREVAIL_LBM_030["j2_2"]["temps (s)"][index_22_PLA:] ) )

PLA_Temperaturemean_signal=[]
PLA_Temperaturemps = []
PLA_Temperaturemms = []
PLA_Temperaturestd_signal = []
PLA_mean_signal=[]
PLA_mps = []
PLA_mms = []
PLA_std_signal = []

for i in range(length_list_PLA):
	PLA_mean_signal.append( np.mean( [ PLA_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22_PLA+i] 
		] ) )
	# 
	PLA_std_signal.append( np.std( [ PLA_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22_PLA+i] 
		] ) )
	# 
	PLA_mms.append( np.mean( [ PLA_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22_PLA+i] 
		] )
	- 1.96* np.std( [ PLA_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22_PLA+i] 
		] ) )
	# 
	PLA_mps.append( np.mean( [ PLA_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22_PLA+i] 
		] )
	+1.96*np.std( [ PLA_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22_PLA+i] 
		] ) )
	# 
	PLA_Temperaturemean_signal.append( np.mean( [ PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i] 
		] ) )
	# 
	PLA_Temperaturestd_signal.append( np.std( [ PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i] 
		] ) )
	# 
	PLA_Temperaturemms.append( np.mean( [ PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i] 
		] )
	-1.96*np.std( [ PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i] 
		] ) )
	# 
	PLA_Temperaturemps.append( np.mean( [ PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i] 
		] )
	+1.96*np.std( [ PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i] 
		] ) )
	
PLA_Temperaturemean_signal_raw=[]
PLA_Temperaturemps_raw = []
PLA_Temperaturemms_raw = []
PLA_Temperaturestd_signal_raw = []
PLA_mean_signal_raw=[]
PLA_mps_raw = []
PLA_mms_raw = []
PLA_std_signal_raw = []

for i in range(length_list_PLA+110):
	PLA_mean_signal_raw.append( np.mean( [ PLA_PREVAIL_LBM_030["j1_1"]["PU"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["PU"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["PU"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["PU"][index_22_PLA+i-110] 
		] ) )
	# 
	PLA_std_signal_raw.append( np.std( [ PLA_PREVAIL_LBM_030["j1_1"]["PU"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["PU"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["PU"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["PU"][index_22_PLA+i-110] 
		] ) )
	# 
	PLA_mms_raw.append( np.mean( [ PLA_PREVAIL_LBM_030["j1_1"]["PU"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["PU"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["PU"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["PU"][index_22_PLA+i-110] 
		] )
	-np.std( [ PLA_PREVAIL_LBM_030["j1_1"]["PU"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["PU"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["PU"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["PU"][index_22_PLA+i-110] 
		] ) )
	# 
	PLA_mps_raw.append( np.mean( [ PLA_PREVAIL_LBM_030["j1_1"]["PU"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["PU"][index_12_PLA+i-110] ,
		PLA_PREVAIL_LBM_030["j2_1"]["PU"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["PU"][index_22_PLA+i-110] 
		] )
	+np.std( [ PLA_PREVAIL_LBM_030["j1_1"]["PU"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["PU"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["PU"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["PU"][index_22_PLA+i-110] 
		] ) )
	# 
	PLA_Temperaturemean_signal_raw.append( np.mean( [ PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i-110] 
		] ) )
	# 
	PLA_Temperaturestd_signal_raw.append( np.std( [ PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i-110] 
		] ) )
	# 
	PLA_Temperaturemms_raw.append( np.mean( [ PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i-110] 
		] )
	-np.std( [ PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i-110] 
		] ) )
	# 
	PLA_Temperaturemps_raw.append( np.mean( [ PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i-110] 
		] )
	+np.std( [ PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(PLA_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12_PLA-110:index_12_PLA+length_list_PLA], PLA_mean_signal_raw, color='gold', label="P$_6$", alpha=0.5)
plt.fill_between(PLA_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12_PLA-110:index_12_PLA+length_list_PLA], PLA_mms_raw, PLA_mps_raw, color='gold', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 100])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_raw_P6.jpg', bbox_inches='tight')




plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(PLA_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12_PLA-110:index_12_PLA+length_list_PLA], PLA_Temperaturemean_signal_raw, color='gold', label="P$_6$", alpha=0.5)
plt.fill_between(PLA_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12_PLA-110:index_12_PLA+length_list_PLA], PLA_Temperaturemms_raw, PLA_Temperaturemps_raw, color='gold', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([29, 32])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P6.jpg', bbox_inches='tight')








plt.close()























# Load the Excel file
file_path = "JPE_PREVAIL_LBM_031.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
JPE_PREVAIL_LBM_031 = {}
for sheet_name, data in all_sheets.items():
    JPE_PREVAIL_LBM_031[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11_JPE = 0
for ii in range(len(JPE_PREVAIL_LBM_031["j1_1"]["temps (s)"] ) ):
	if JPE_PREVAIL_LBM_031["j1_1"]["temps (s)"][ii] == 0:
		index_11_JPE = ii 
		break
index_12_JPE = 0
for ii in range(len(JPE_PREVAIL_LBM_031["j1_2"]["temps (s)"] ) ):
	if JPE_PREVAIL_LBM_031["j1_2"]["temps (s)"][ii] == 0:
		index_12_JPE = ii 
		break
index_21_JPE = 0
for ii in range(len(JPE_PREVAIL_LBM_031["j2_1"]["temps (s)"] ) ):
	if JPE_PREVAIL_LBM_031["j2_1"]["temps (s)"][ii] == 0:
		index_21_JPE = ii 
		break
index_22_JPE = 0
for ii in range(len(JPE_PREVAIL_LBM_031["j2_2"]["temps (s)"] ) ):
	if JPE_PREVAIL_LBM_031["j2_2"]["temps (s)"][ii] == 0:
		index_22_JPE = ii 
		break		

length_list_JPE = min(len(JPE_PREVAIL_LBM_031["j1_1"]["temps (s)"][index_11_JPE:]) , len(JPE_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12_JPE:]) , len(JPE_PREVAIL_LBM_031["j2_1"]["temps (s)"][index_21_JPE:]) , len(JPE_PREVAIL_LBM_031["j2_2"]["temps (s)"][index_22_JPE:] ) )

JPE_Temperaturemean_signal=[]
JPE_Temperaturemps = []
JPE_Temperaturemms = []
JPE_Temperaturestd_signal = []
JPE_mean_signal=[]
JPE_mps = []
JPE_mms = []
JPE_std_signal = []

for i in range(length_list_JPE):
	JPE_mean_signal.append( np.mean( [ JPE_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22_JPE+i] 
		] ) )
	# 
	JPE_std_signal.append( np.std( [ JPE_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22_JPE+i] 
		] ) )
	# 
	JPE_mms.append( np.mean( [ JPE_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22_JPE+i] 
		] )
	- 1.96* np.std( [ JPE_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22_JPE+i] 
		] ) )
	# 
	JPE_mps.append( np.mean( [ JPE_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22_JPE+i] 
		] )
	+1.96*np.std( [ JPE_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22_JPE+i] 
		] ) )
	# 
	JPE_Temperaturemean_signal.append( np.mean( [ JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i] 
		] ) )
	# 
	JPE_Temperaturestd_signal.append( np.std( [ JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i] 
		] ) )
	# 
	JPE_Temperaturemms.append( np.mean( [ JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i] 
		] )
	-1.96*np.std( [ JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i] 
		] ) )
	# 
	JPE_Temperaturemps.append( np.mean( [ JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i] 
		] )
	+1.96*np.std( [ JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i] 
		] ) )
	
JPE_Temperaturemean_signal_raw=[]
JPE_Temperaturemps_raw = []
JPE_Temperaturemms_raw = []
JPE_Temperaturestd_signal_raw = []
JPE_mean_signal_raw=[]
JPE_mps_raw = []
JPE_mms_raw = []
JPE_std_signal_raw = []

for i in range(length_list_JPE+110):
	JPE_mean_signal_raw.append( np.mean( [ JPE_PREVAIL_LBM_031["j1_1"]["PU"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["PU"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["PU"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["PU"][index_22_JPE+i-110] 
		] ) )
	# 
	JPE_std_signal_raw.append( np.std( [ JPE_PREVAIL_LBM_031["j1_1"]["PU"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["PU"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["PU"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["PU"][index_22_JPE+i-110] 
		] ) )
	# 
	JPE_mms_raw.append( np.mean( [ JPE_PREVAIL_LBM_031["j1_1"]["PU"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["PU"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["PU"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["PU"][index_22_JPE+i-110] 
		] )
	-np.std( [ JPE_PREVAIL_LBM_031["j1_1"]["PU"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["PU"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["PU"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["PU"][index_22_JPE+i-110] 
		] ) )
	# 
	JPE_mps_raw.append( np.mean( [ JPE_PREVAIL_LBM_031["j1_1"]["PU"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["PU"][index_12_JPE+i-110] ,
		JPE_PREVAIL_LBM_031["j2_1"]["PU"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["PU"][index_22_JPE+i-110] 
		] )
	+np.std( [ JPE_PREVAIL_LBM_031["j1_1"]["PU"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["PU"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["PU"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["PU"][index_22_JPE+i-110] 
		] ) )
	# 
	JPE_Temperaturemean_signal_raw.append( np.mean( [ JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i-110] 
		] ) )
	# 
	JPE_Temperaturestd_signal_raw.append( np.std( [ JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i-110] 
		] ) )
	# 
	JPE_Temperaturemms_raw.append( np.mean( [ JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i-110] 
		] )
	-np.std( [ JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i-110] 
		] ) )
	# 
	JPE_Temperaturemps_raw.append( np.mean( [ JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i-110] 
		] )
	+np.std( [ JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(JPE_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12_JPE-110:index_12_JPE+length_list_JPE], JPE_mean_signal_raw, linestyle='-', color='blue', label="P$_7$", alpha=0.5)
plt.fill_between(JPE_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12_JPE-110:index_12_JPE+length_list_JPE], JPE_mms_raw, JPE_mps_raw, color='blue', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 150])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_raw_P7.jpg', bbox_inches='tight')




plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(JPE_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12_JPE-110:index_12_JPE+length_list_JPE], JPE_Temperaturemean_signal_raw, linestyle='-', color='blue', label="P$_7$", alpha=0.5)
plt.fill_between(JPE_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12_JPE-110:index_12_JPE+length_list_JPE], JPE_Temperaturemms_raw, JPE_Temperaturemps_raw, color='blue', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([29.8, 31.8])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P7.jpg', bbox_inches='tight')








plt.close()
















# Load the Excel file
file_path = "MLM_PREVAIL_LBM_032.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
MLM_PREVAIL_LBM_032 = {}
for sheet_name, data in all_sheets.items():
    MLM_PREVAIL_LBM_032[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11_MLM = 0
for ii in range(len(MLM_PREVAIL_LBM_032["j1_1"]["temps (s)"] ) ):
	if MLM_PREVAIL_LBM_032["j1_1"]["temps (s)"][ii] == 0:
		index_11_MLM = ii 
		break
index_12_MLM = 0
for ii in range(len(MLM_PREVAIL_LBM_032["j1_2"]["temps (s)"] ) ):
	if MLM_PREVAIL_LBM_032["j1_2"]["temps (s)"][ii] == 0:
		index_12_MLM = ii 
		break
# index_21_MLM = 0
# for ii in range(len(MLM_PREVAIL_LBM_032["j2_1"]["temps (s)"] ) ):
# 	if MLM_PREVAIL_LBM_032["j2_1"]["temps (s)"][ii] == 0:
# 		index_21_MLM = ii 
# 		break
# index_22_MLM = 0
# for ii in range(len(MLM_PREVAIL_LBM_032["j2_2"]["temps (s)"] ) ):
# 	if MLM_PREVAIL_LBM_032["j2_2"]["temps (s)"][ii] == 0:
# 		index_22_MLM = ii 
# 		break		

length_list_MLM = min(len(MLM_PREVAIL_LBM_032["j1_1"]["temps (s)"][index_11_MLM:]) , len(MLM_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12_MLM:]) )# , len(MLM_PREVAIL_LBM_032["j2_1"]["temps (s)"][index_21_MLM:]) , len(MLM_PREVAIL_LBM_032["j2_2"]["temps (s)"][index_22_MLM:] ) )

MLM_Temperaturemean_signal=[]
MLM_Temperaturemps = []
MLM_Temperaturemms = []
MLM_Temperaturestd_signal = []
MLM_mean_signal=[]
MLM_mps = []
MLM_mms = []
MLM_std_signal = []

for i in range(length_list_MLM):
	MLM_mean_signal.append( np.mean( [ MLM_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22_MLM+i] 
		] ) )
	# 
	MLM_std_signal.append( np.std( [ MLM_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22_MLM+i] 
		] ) )
	# 
	MLM_mms.append( np.mean( [ MLM_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22_MLM+i] 
		] )
	- 1.96* np.std( [ MLM_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22_MLM+i] 
		] ) )
	# 
	MLM_mps.append( np.mean( [ MLM_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22_MLM+i] 
		] )
	+1.96*np.std( [ MLM_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22_MLM+i] 
		] ) )
	# 
	MLM_Temperaturemean_signal.append( np.mean( [ MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_MLM+i] 
		] ) )
	# 
	MLM_Temperaturestd_signal.append( np.std( [ MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_MLM+i] 
		] ) )
	# 
	MLM_Temperaturemms.append( np.mean( [ MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_MLM+i] 
		] )
	-1.96*np.std( [ MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_MLM+i] 
		] ) )
	# 
	MLM_Temperaturemps.append( np.mean( [ MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_MLM+i] 
		] )
	+1.96*np.std( [ MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_MLM+i]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_MLM+i] 
		] ) )
	
MLM_Temperaturemean_signal_raw=[]
MLM_Temperaturemps_raw = []
MLM_Temperaturemms_raw = []
MLM_Temperaturestd_signal_raw = []
MLM_mean_signal_raw=[]
MLM_mps_raw = []
MLM_mms_raw = []
MLM_std_signal_raw = []

for i in range(length_list_MLM+110):
	MLM_mean_signal_raw.append( np.mean( [ MLM_PREVAIL_LBM_032["j1_1"]["PU"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["PU"][index_12_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["PU"][index_21_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["PU"][index_22_MLM+i-110] 
		] ) )
	# 
	MLM_std_signal_raw.append( np.std( [ MLM_PREVAIL_LBM_032["j1_1"]["PU"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["PU"][index_12_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["PU"][index_21_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["PU"][index_22_MLM+i-110] 
		] ) )
	# 
	MLM_mms_raw.append( np.mean( [ MLM_PREVAIL_LBM_032["j1_1"]["PU"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["PU"][index_12_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["PU"][index_21_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["PU"][index_22_MLM+i-110] 
		] )
	-np.std( [ MLM_PREVAIL_LBM_032["j1_1"]["PU"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["PU"][index_12_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["PU"][index_21_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["PU"][index_22_MLM+i-110] 
		] ) )
	# 
	MLM_mps_raw.append( np.mean( [ MLM_PREVAIL_LBM_032["j1_1"]["PU"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["PU"][index_12_MLM+i-110] #,
		# MLM_PREVAIL_LBM_032["j2_1"]["PU"][index_21_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["PU"][index_22_MLM+i-110] 
		] )
	+np.std( [ MLM_PREVAIL_LBM_032["j1_1"]["PU"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["PU"][index_12_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["PU"][index_21_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["PU"][index_22_MLM+i-110] 
		] ) )
	# 
	MLM_Temperaturemean_signal_raw.append( np.mean( [ MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_MLM+i-110] 
		] ) )
	# 
	MLM_Temperaturestd_signal_raw.append( np.std( [ MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_MLM+i-110] 
		] ) )
	# 
	MLM_Temperaturemms_raw.append( np.mean( [ MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_MLM+i-110] 
		] )
	-np.std( [ MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_MLM+i-110] 
		] ) )
	# 
	MLM_Temperaturemps_raw.append( np.mean( [ MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_MLM+i-110] 
		] )
	+np.std( [ MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_MLM+i-110]#,
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_MLM+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(MLM_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12_MLM-110:index_12_MLM+length_list_MLM], MLM_mean_signal_raw, linestyle='-', color='red', label="P$_8$", alpha=0.5)
plt.fill_between(MLM_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12_MLM-110:index_12_MLM+length_list_MLM], MLM_mms_raw, MLM_mps_raw, color='red', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 150])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_raw_P8.jpg', bbox_inches='tight')




plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(MLM_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12_MLM-110:index_12_MLM+length_list_MLM], MLM_Temperaturemean_signal_raw, linestyle='-', color='red', label="P$_8$", alpha=0.5)
plt.fill_between(MLM_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12_MLM-110:index_12_MLM+length_list_MLM], MLM_Temperaturemms_raw, MLM_Temperaturemps_raw, color='red', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P8.jpg', bbox_inches='tight')









plt.close()


































# Load the Excel file
file_path = "ALA_PREVAIL_LBM_033.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
ALA_PREVAIL_LBM_033 = {}
for sheet_name, data in all_sheets.items():
    ALA_PREVAIL_LBM_033[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11_ALA = 0
for ii in range(len(ALA_PREVAIL_LBM_033["j1_1"]["temps (s)"] ) ):
	if ALA_PREVAIL_LBM_033["j1_1"]["temps (s)"][ii] == 0:
		index_11_ALA = ii 
		break
index_12_ALA = 0
for ii in range(len(ALA_PREVAIL_LBM_033["j1_2"]["temps (s)"] ) ):
	if ALA_PREVAIL_LBM_033["j1_2"]["temps (s)"][ii] == 0:
		index_12_ALA = ii 
		break
index_21_ALA = 0
for ii in range(len(ALA_PREVAIL_LBM_033["j2_1"]["temps (s)"] ) ):
	if ALA_PREVAIL_LBM_033["j2_1"]["temps (s)"][ii] == 0:
		index_21_ALA = ii 
		break
index_22_ALA = 0
for ii in range(len(ALA_PREVAIL_LBM_033["j2_2"]["temps (s)"] ) ):
	if ALA_PREVAIL_LBM_033["j2_2"]["temps (s)"][ii] == 0:
		index_22_ALA = ii 
		break		

length_list_ALA = min(len(ALA_PREVAIL_LBM_033["j1_1"]["temps (s)"][index_11_ALA:]) , len(ALA_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12_ALA:]), len(ALA_PREVAIL_LBM_033["j2_1"]["temps (s)"][index_21_ALA:]) , len(ALA_PREVAIL_LBM_033["j2_2"]["temps (s)"][index_22_ALA:] ) )

ALA_Temperaturemean_signal=[]
ALA_Temperaturemps = []
ALA_Temperaturemms = []
ALA_Temperaturestd_signal = []
ALA_mean_signal=[]
ALA_mps = []
ALA_mms = []
ALA_std_signal = []

for i in range(length_list_ALA):
	ALA_mean_signal.append( np.mean( [ ALA_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21_ALA+i],
		ALA_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22_ALA+i] 
		] ) )
	# 
	ALA_std_signal.append( np.std( [ ALA_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21_ALA+i],
		ALA_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22_ALA+i] 
		] ) )
	# 
	ALA_mms.append( np.mean( [ ALA_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21_ALA+i],
		ALA_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22_ALA+i] 
		] )
	- 1.96* np.std( [ ALA_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21_ALA+i],
		ALA_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22_ALA+i] 
		] ) )
	# 
	ALA_mps.append( np.mean( [ ALA_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21_ALA+i],
		ALA_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22_ALA+i] 
		] )
	+1.96*np.std( [ ALA_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21_ALA+i],
		ALA_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22_ALA+i] 
		] ) )
	# 
	ALA_Temperaturemean_signal.append( np.mean( [ ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i] 
		] ) )
	# 
	ALA_Temperaturestd_signal.append( np.std( [ ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i] 
		] ) )
	# 
	ALA_Temperaturemms.append( np.mean( [ ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i] 
		] )
	-1.96*np.std( [ ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i] 
		] ) )
	# 
	ALA_Temperaturemps.append( np.mean( [ ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i] 
		] )
	+1.96*np.std( [ ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i] 
		] ) )
	
ALA_Temperaturemean_signal_raw=[]
ALA_Temperaturemps_raw = []
ALA_Temperaturemms_raw = []
ALA_Temperaturestd_signal_raw = []
ALA_mean_signal_raw=[]
ALA_mps_raw = []
ALA_mms_raw = []
ALA_std_signal_raw = []

for i in range(length_list_ALA+110):
	ALA_mean_signal_raw.append( np.mean( [ ALA_PREVAIL_LBM_033["j1_1"]["PU"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["PU"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["PU"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["PU"][index_22_ALA+i-110] 
		] ) )
	# 
	ALA_std_signal_raw.append( np.std( [ ALA_PREVAIL_LBM_033["j1_1"]["PU"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["PU"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["PU"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["PU"][index_22_ALA+i-110] 
		] ) )
	# 
	ALA_mms_raw.append( np.mean( [ ALA_PREVAIL_LBM_033["j1_1"]["PU"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["PU"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["PU"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["PU"][index_22_ALA+i-110] 
		] )
	-np.std( [ ALA_PREVAIL_LBM_033["j1_1"]["PU"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["PU"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["PU"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["PU"][index_22_ALA+i-110] 
		] ) )
	# 
	ALA_mps_raw.append( np.mean( [ ALA_PREVAIL_LBM_033["j1_1"]["PU"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["PU"][index_12_ALA+i-110] ,
		ALA_PREVAIL_LBM_033["j2_1"]["PU"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["PU"][index_22_ALA+i-110] 
		] )
	+np.std( [ ALA_PREVAIL_LBM_033["j1_1"]["PU"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["PU"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["PU"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["PU"][index_22_ALA+i-110] 
		] ) )
	# 
	ALA_Temperaturemean_signal_raw.append( np.mean( [ ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i-110] 
		] ) )
	# 
	ALA_Temperaturestd_signal_raw.append( np.std( [ ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i-110] 
		] ) )
	# 
	ALA_Temperaturemms_raw.append( np.mean( [ ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i-110] 
		] )
	-np.std( [ ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i-110] 
		] ) )
	# 
	ALA_Temperaturemps_raw.append( np.mean( [ ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i-110] 
		] )
	+np.std( [ ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(ALA_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12_ALA-110:index_12_ALA+length_list_ALA], ALA_mean_signal_raw, linestyle='-', color='fuchsia', label="P$_9$", alpha=0.5)
plt.fill_between(ALA_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12_ALA-110:index_12_ALA+length_list_ALA], ALA_mms_raw, ALA_mps_raw, color='fuchsia', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 120])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_raw_P9.jpg', bbox_inches='tight')



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(ALA_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12_ALA-110:index_12_ALA+length_list_ALA], ALA_Temperaturemean_signal_raw, linestyle='-', color='fuchsia', label="P$_9$", alpha=0.5)
plt.fill_between(ALA_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12_ALA-110:index_12_ALA+length_list_ALA], ALA_Temperaturemms_raw, ALA_Temperaturemps_raw, color='fuchsia', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P9.jpg', bbox_inches='tight')




plt.close()



































# # Load the Excel file
# file_path = "AES_PREVAIL_LBM_034.xlsx"  # Change this to the path of your Excel file

# # Read all sheets into a dictionary of DataFrames
# all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# # Optional: Extract specific columns for all sheets
# column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# # Extract data from each sheet into a structured dictionary
# AES_PREVAIL_LBM_034 = {}
# for sheet_name, data in all_sheets.items():
#     AES_PREVAIL_LBM_034[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

# index_11_AES = 0
# for ii in range(len(AES_PREVAIL_LBM_034["j1_1"]["temps (s)"] ) ):
# 	if AES_PREVAIL_LBM_034["j1_1"]["temps (s)"][ii] == 0:
# 		index_11_AES = ii 
# 		break
# index_12_AES = 0
# for ii in range(len(AES_PREVAIL_LBM_034["j1_2"]["temps (s)"] ) ):
# 	if AES_PREVAIL_LBM_034["j1_2"]["temps (s)"][ii] == 0:
# 		index_12_AES = ii 
# 		break
# # index_21_AES = 0
# # for ii in range(len(AES_PREVAIL_LBM_034["j2_1"]["temps (s)"] ) ):
# # 	if AES_PREVAIL_LBM_034["j2_1"]["temps (s)"][ii] == 0:
# # 		index_21_AES = ii 
# # 		break
# # index_22_AES = 0
# # for ii in range(len(AES_PREVAIL_LBM_034["j2_2"]["temps (s)"] ) ):
# # 	if AES_PREVAIL_LBM_034["j2_2"]["temps (s)"][ii] == 0:
# # 		index_22_AES = ii 
# # 		break		

# length_list_AES = min(len(AES_PREVAIL_LBM_034["j1_1"]["temps (s)"][index_11_AES:]) , len(AES_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12_AES:]) )# , len(AES_PREVAIL_LBM_034["j2_1"]["temps (s)"][index_21_AES:]) , len(AES_PREVAIL_LBM_034["j2_2"]["temps (s)"][index_22_AES:] ) )

# AES_Temperaturemean_signal=[]
# AES_Temperaturemps = []
# AES_Temperaturemms = []
# AES_Temperaturestd_signal = []
# AES_mean_signal=[]
# AES_mps = []
# AES_mms = []
# AES_std_signal = []

# for i in range(length_list_AES):
# 	AES_mean_signal.append( np.mean( [ AES_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11_AES+i],
# 		AES_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22_AES+i] 
# 		] ) )
# 	# 
# 	AES_std_signal.append( np.std( [ AES_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11_AES+i],
# 		AES_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22_AES+i] 
# 		] ) )
# 	# 
# 	AES_mms.append( np.mean( [ AES_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11_AES+i],
# 		AES_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22_AES+i] 
# 		] )
# 	- 1.96* np.std( [ AES_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11_AES+i],
# 		AES_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22_AES+i] 
# 		] ) )
# 	# 
# 	AES_mps.append( np.mean( [ AES_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11_AES+i],
# 		AES_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22_AES+i] 
# 		] )
# 	+1.96*np.std( [ AES_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11_AES+i],
# 		AES_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22_AES+i] 
# 		] ) )
# 	# 
# 	AES_Temperaturemean_signal.append( np.mean( [ AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i],
# 		AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i] 
# 		] ) )
# 	# 
# 	AES_Temperaturestd_signal.append( np.std( [ AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i],
# 		AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i] 
# 		] ) )
# 	# 
# 	AES_Temperaturemms.append( np.mean( [ AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i],
# 		AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i] 
# 		] )
# 	-1.96*np.std( [ AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i],
# 		AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i] 
# 		] ) )
# 	# 
# 	AES_Temperaturemps.append( np.mean( [ AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i],
# 		AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i] 
# 		] )
# 	+1.96*np.std( [ AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i],
# 		AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i] 
# 		] ) )
	
# AES_Temperaturemean_signal_raw=[]
# AES_Temperaturemps_raw = []
# AES_Temperaturemms_raw = []
# AES_Temperaturestd_signal_raw = []
# AES_mean_signal_raw=[]
# AES_mps_raw = []
# AES_mms_raw = []
# AES_std_signal_raw = []

# for i in range(length_list_AES+110):
# 	AES_mean_signal_raw.append( np.mean( [ AES_PREVAIL_LBM_034["j1_1"]["PU"][index_11_AES+i-110],
# 		AES_PREVAIL_LBM_034["j1_2"]["PU"][index_12_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["PU"][index_21_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["PU"][index_22_AES+i-110] 
# 		] ) )
# 	# 
# 	AES_std_signal_raw.append( np.std( [ AES_PREVAIL_LBM_034["j1_1"]["PU"][index_11_AES+i-110],
# 		AES_PREVAIL_LBM_034["j1_2"]["PU"][index_12_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["PU"][index_21_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["PU"][index_22_AES+i-110] 
# 		] ) )
# 	# 
# 	AES_mms_raw.append( np.mean( [ AES_PREVAIL_LBM_034["j1_1"]["PU"][index_11_AES+i-110],
# 		AES_PREVAIL_LBM_034["j1_2"]["PU"][index_12_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["PU"][index_21_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["PU"][index_22_AES+i-110] 
# 		] )
# 	-np.std( [ AES_PREVAIL_LBM_034["j1_1"]["PU"][index_11_AES+i-110],
# 		AES_PREVAIL_LBM_034["j1_2"]["PU"][index_12_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["PU"][index_21_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["PU"][index_22_AES+i-110] 
# 		] ) )
# 	# 
# 	AES_mps_raw.append( np.mean( [ AES_PREVAIL_LBM_034["j1_1"]["PU"][index_11_AES+i-110],
# 		AES_PREVAIL_LBM_034["j1_2"]["PU"][index_12_AES+i-110] #,
# 		# AES_PREVAIL_LBM_034["j2_1"]["PU"][index_21_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["PU"][index_22_AES+i-110] 
# 		] )
# 	+np.std( [ AES_PREVAIL_LBM_034["j1_1"]["PU"][index_11_AES+i-110],
# 		AES_PREVAIL_LBM_034["j1_2"]["PU"][index_12_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["PU"][index_21_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["PU"][index_22_AES+i-110] 
# 		] ) )
# 	# 
# 	AES_Temperaturemean_signal_raw.append( np.mean( [ AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i-110],
# 		AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i-110] 
# 		] ) )
# 	# 
# 	AES_Temperaturestd_signal_raw.append( np.std( [ AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i-110],
# 		AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i-110] 
# 		] ) )
# 	# 
# 	AES_Temperaturemms_raw.append( np.mean( [ AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i-110],
# 		AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i-110] 
# 		] )
# 	-np.std( [ AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i-110],
# 		AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i-110] 
# 		] ) )
# 	# 
# 	AES_Temperaturemps_raw.append( np.mean( [ AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i-110],
# 		AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i-110] 
# 		] )
# 	+np.std( [ AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i-110],
# 		AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i-110]#,
# 		# AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i-110] 
# 		] ) )



# plt.rcParams.update({'font.size': 25})
# plt.figure(figsize=(20, 12))
# # # 
# plt.plot(AES_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12_AES-110:index_12_AES+length_list_AES], AES_mean_signal_raw, linestyle='-', color='darkgreen', label="P$_10$", alpha=0.5)
# plt.fill_between(AES_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12_AES-110:index_12_AES+length_list_AES], AES_mms_raw, AES_mps_raw, color='darkgreen', alpha=0.2)
# plt.legend()
# plt.xlim([-60, 780])
# plt.ylim([0, 250])
# plt.xlabel("Time [s]")
# plt.ylabel("LDF [AU]")
# plt.grid()
# plt.savefig('./Figures/superp_mean_std_raw_P10.jpg', bbox_inches='tight')

# plt.rcParams.update({'font.size': 25})
# plt.figure(figsize=(20, 12))
# # # 
# plt.plot(AES_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12_AES-110:index_12_AES+length_list_AES], AES_Temperaturemean_signal_raw, linestyle='-', color='darkgreen', label="P$_10$", alpha=0.5)
# plt.fill_between(AES_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12_AES-110:index_12_AES+length_list_AES], AES_Temperaturemms_raw, AES_Temperaturemps_raw, color='darkgreen', alpha=0.2)
# plt.legend()
# plt.xlim([-60, 780])
# # plt.ylim([27, 33])
# plt.xlabel("Time [s]")
# plt.ylabel("Temperature [°C]")
# plt.grid()
# plt.savefig('./Figures/superp_mean_std_temp_raw_P10.jpg', bbox_inches='tight')






























# Load the Excel file
file_path = "CDE_PREVAIL_LBM_035.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
CDE_PREVAIL_LBM_035 = {}
for sheet_name, data in all_sheets.items():
    CDE_PREVAIL_LBM_035[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11_CDE = 0
for ii in range(len(CDE_PREVAIL_LBM_035["j1_1"]["temps (s)"] ) ):
	if CDE_PREVAIL_LBM_035["j1_1"]["temps (s)"][ii] == 0:
		index_11_CDE = ii 
		break
index_12_CDE = 0
for ii in range(len(CDE_PREVAIL_LBM_035["j1_2"]["temps (s)"] ) ):
	if CDE_PREVAIL_LBM_035["j1_2"]["temps (s)"][ii] == 0:
		index_12_CDE = ii 
		break
index_21_CDE = 0
for ii in range(len(CDE_PREVAIL_LBM_035["j2_1"]["temps (s)"] ) ):
	if CDE_PREVAIL_LBM_035["j2_1"]["temps (s)"][ii] == 0:
		index_21_CDE = ii 
		break
index_22_CDE = 0
for ii in range(len(CDE_PREVAIL_LBM_035["j2_2"]["temps (s)"] ) ):
	if CDE_PREVAIL_LBM_035["j2_2"]["temps (s)"][ii] == 0:
		index_22_CDE = ii 
		break		

length_list_CDE = min(len(CDE_PREVAIL_LBM_035["j1_1"]["temps (s)"][index_11_CDE:]) , len(CDE_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12_CDE:]) , len(CDE_PREVAIL_LBM_035["j2_1"]["temps (s)"][index_21_CDE:]) , len(CDE_PREVAIL_LBM_035["j2_2"]["temps (s)"][index_22_CDE:] ) )

CDE_Temperaturemean_signal=[]
CDE_Temperaturemps = []
CDE_Temperaturemms = []
CDE_Temperaturestd_signal = []
CDE_mean_signal=[]
CDE_mps = []
CDE_mms = []
CDE_std_signal = []

for i in range(length_list_CDE):
	CDE_mean_signal.append( np.mean( [ CDE_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22_CDE+i] 
		] ) )
	# 
	CDE_std_signal.append( np.std( [ CDE_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22_CDE+i] 
		] ) )
	# 
	CDE_mms.append( np.mean( [ CDE_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22_CDE+i] 
		] )
	- 1.96* np.std( [ CDE_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22_CDE+i] 
		] ) )
	# 
	CDE_mps.append( np.mean( [ CDE_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22_CDE+i] 
		] )
	+1.96*np.std( [ CDE_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22_CDE+i] 
		] ) )
	# 
	CDE_Temperaturemean_signal.append( np.mean( [ CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i] 
		] ) )
	# 
	CDE_Temperaturestd_signal.append( np.std( [ CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i] 
		] ) )
	# 
	CDE_Temperaturemms.append( np.mean( [ CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i] 
		] )
	-1.96*np.std( [ CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i] 
		] ) )
	# 
	CDE_Temperaturemps.append( np.mean( [ CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i] 
		] )
	+1.96*np.std( [ CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i] 
		] ) )
	
CDE_Temperaturemean_signal_raw=[]
CDE_Temperaturemps_raw = []
CDE_Temperaturemms_raw = []
CDE_Temperaturestd_signal_raw = []
CDE_mean_signal_raw=[]
CDE_mps_raw = []
CDE_mms_raw = []
CDE_std_signal_raw = []

for i in range(length_list_CDE+110):
	CDE_mean_signal_raw.append( np.mean( [ CDE_PREVAIL_LBM_035["j1_1"]["PU"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["PU"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["PU"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["PU"][index_22_CDE+i-110] 
		] ) )
	# 
	CDE_std_signal_raw.append( np.std( [ CDE_PREVAIL_LBM_035["j1_1"]["PU"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["PU"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["PU"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["PU"][index_22_CDE+i-110] 
		] ) )
	# 
	CDE_mms_raw.append( np.mean( [ CDE_PREVAIL_LBM_035["j1_1"]["PU"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["PU"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["PU"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["PU"][index_22_CDE+i-110] 
		] )
	-np.std( [ CDE_PREVAIL_LBM_035["j1_1"]["PU"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["PU"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["PU"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["PU"][index_22_CDE+i-110] 
		] ) )
	# 
	CDE_mps_raw.append( np.mean( [ CDE_PREVAIL_LBM_035["j1_1"]["PU"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["PU"][index_12_CDE+i-110] ,
		CDE_PREVAIL_LBM_035["j2_1"]["PU"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["PU"][index_22_CDE+i-110] 
		] )
	+np.std( [ CDE_PREVAIL_LBM_035["j1_1"]["PU"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["PU"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["PU"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["PU"][index_22_CDE+i-110] 
		] ) )
	# 
	CDE_Temperaturemean_signal_raw.append( np.mean( [ CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i-110] 
		] ) )
	# 
	CDE_Temperaturestd_signal_raw.append( np.std( [ CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i-110] 
		] ) )
	# 
	CDE_Temperaturemms_raw.append( np.mean( [ CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i-110] 
		] )
	-np.std( [ CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i-110] 
		] ) )
	# 
	CDE_Temperaturemps_raw.append( np.mean( [ CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i-110] 
		] )
	+np.std( [ CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(CDE_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12_CDE-110:index_12_CDE+length_list_CDE], CDE_mean_signal_raw, linestyle='-', color='gold', label="P$_{11}$", alpha=0.5)
plt.fill_between(CDE_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12_CDE-110:index_12_CDE+length_list_CDE], CDE_mms_raw, CDE_mps_raw, color='gold', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 60])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_raw_P11.jpg', bbox_inches='tight')




plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(CDE_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12_CDE-110:index_12_CDE+length_list_CDE], CDE_Temperaturemean_signal_raw, linestyle='-', color='gold', label="P$_{11}$", alpha=0.5)
plt.fill_between(CDE_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12_CDE-110:index_12_CDE+length_list_CDE], CDE_Temperaturemms_raw, CDE_Temperaturemps_raw, color='gold', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P11.jpg', bbox_inches='tight')






plt.close()








































































# Raw

# A compléter

all_means=[]
all_means_temp = []
all_std=[]
all_std_temp =[]
# for i in range(min(length_list_ASE+110,length_list+110,length_list_MDI+110,length_list_MMA+110,length_list_CBO+110,length_list_PLA+110,length_list_JPE+110,length_list_MLM+110,length_list_ALA+110,length_list_AES+110,length_list_CDE+110)):
for i in range(min(length_list_ASE+110,length_list+110,length_list_MDI+110,length_list_MMA+110,length_list_CBO+110,length_list_PLA+110,length_list_JPE+110,length_list_MLM+110,length_list_ALA+110,length_list_CDE+110)):
	all_means.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-110],
		# 
		JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
		JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
		JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
		JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110],
		# 
		MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-110],
		# 
		MMA_PREVAIL_LBM_028["j1_1"]["PU"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["PU"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["PU"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["PU"][index_22_MMA+i-110],  
		# 
		CBO_PREVAIL_LBM_029["j1_1"]["PU"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["PU"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["PU"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["PU"][index_22_CBO+i-110], 
		# 
		PLA_PREVAIL_LBM_030["j1_1"]["PU"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["PU"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["PU"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["PU"][index_22_PLA+i-110] ,
		# 
		JPE_PREVAIL_LBM_031["j1_1"]["PU"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["PU"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["PU"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["PU"][index_22_JPE+i-110] , 
		# 
		MLM_PREVAIL_LBM_032["j1_1"]["PU"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["PU"][index_12_MLM+i-110],
		# MLM_PREVAIL_LBM_032["j2_1"]["PU"][index_21_PLA+i-110],
		# MLM_PREVAIL_LBM_032["j2_2"]["PU"][index_22_PLA+i-110] ,
		# #
		ALA_PREVAIL_LBM_033["j1_1"]["PU"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["PU"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["PU"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["PU"][index_22_ALA+i-110] ,
		# # 
		# AES_PREVAIL_LBM_034["j1_1"]["PU"][index_11_AES+i-110],
		# AES_PREVAIL_LBM_034["j1_2"]["PU"][index_12_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_1"]["PU"][index_21_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_2"]["PU"][index_22_AES+i-110] ,
		# # 
		CDE_PREVAIL_LBM_035["j1_1"]["PU"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["PU"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["PU"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["PU"][index_22_CDE+i-110] 
		# 
		] ) )
	# 
	all_std.append( np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-110],
		# 
		JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
		JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
		JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
		JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110],
		# 
		MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-110] ,
		# 
		MMA_PREVAIL_LBM_028["j1_1"]["PU"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["PU"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["PU"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["PU"][index_22_MMA+i-110],  
		# 
		CBO_PREVAIL_LBM_029["j1_1"]["PU"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["PU"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["PU"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["PU"][index_22_CBO+i-110], 
		# 
		PLA_PREVAIL_LBM_030["j1_1"]["PU"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["PU"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["PU"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["PU"][index_22_PLA+i-110] ,
		# 
		JPE_PREVAIL_LBM_031["j1_1"]["PU"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["PU"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["PU"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["PU"][index_22_JPE+i-110] ,
		# 
		# 
		MLM_PREVAIL_LBM_032["j1_1"]["PU"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["PU"][index_12_MLM+i-110],
		# MLM_PREVAIL_LBM_032["j2_1"]["PU"][index_21_PLA+i-110],
		# MLM_PREVAIL_LBM_032["j2_2"]["PU"][index_22_PLA+i-110] ,
		# #
		ALA_PREVAIL_LBM_033["j1_1"]["PU"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["PU"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["PU"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["PU"][index_22_ALA+i-110] ,
		# # 
		# AES_PREVAIL_LBM_034["j1_1"]["PU"][index_11_AES+i-110],
		# AES_PREVAIL_LBM_034["j1_2"]["PU"][index_12_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_1"]["PU"][index_21_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_2"]["PU"][index_22_AES+i-110] ,
		# # 
		CDE_PREVAIL_LBM_035["j1_1"]["PU"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["PU"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["PU"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["PU"][index_22_CDE+i-110] 
		] ) )
	# 
	all_means_temp.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-110],
		# 
		JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
		JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
		JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
		JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110],
		# 
		MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-110],
		# 
		MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i-110],  
		# 
		CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i-110], 
		# 
		PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i-110] ,
		# 
		JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i-110] ,
		# 
		# 
		MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i-110],
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_PLA+i-110],
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_PLA+i-110] ,
		# #
		ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i-110] ,
		# # 
		# AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i-110],
		# AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i-110] ,
		# # 
		CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i-110] 
		] ) )
	# 
	all_std_temp.append( np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-110],
		# 
		JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
		JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
		JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
		JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110],
		# 
		MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-110] ,
		# 
		MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i-110],  
		# 
		CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i-110], 
		# 
		PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i-110] ,
		# 
		JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i-110] ,
		# 
		MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i-110],
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_PLA+i-110],
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_PLA+i-110] ,
		# #
		ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i-110] ,
		# # 
		# AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i-110],
		# AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i-110] ,
		# # 
		CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i-110] 
		# 
		] ) )

all_mms = [all_means[i]-all_std[i] for i in range(len(all_means))]
all_mps = [all_means[i]+all_std[i] for i in range(len(all_means))]

all_means_pc=[]
all_std_pc=[]
# Percentage
# for i in range(min(length_list_ASE,length_list,length_list_MDI,length_list_MMA,length_list_CBO,length_list_PLA,length_list_JPE,length_list_MLM,length_list_ALA,length_list_AES,length_list_CDE)):
for i in range(min(length_list_ASE,length_list,length_list_MDI,length_list_MMA,length_list_CBO,length_list_PLA,length_list_JPE,length_list_MLM+110,length_list_ALA+110,length_list_CDE+110)):
	all_means_pc.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],
		ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],
		ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],
		ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i],
		# 
		JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
		JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
		JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
		JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i],
		# 
		MDI_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22_MDI+i]  ,
		# 
		MMA_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22_MMA+i],  
		# 
		CBO_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22_CBO+i], 
		# 
		PLA_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22_PLA+i] ,
		# 
		JPE_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22_JPE+i] ,
		# 
		MLM_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12_MLM+i],
		# MLM_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21_PLA+i],
		# MLM_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22_PLA+i] ,
		# #
		ALA_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21_ALA+i],
		ALA_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22_ALA+i] ,
		# # 
		# AES_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11_AES+i],
		# AES_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12_AES+i],
		# # AES_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21_AES+i],
		# # AES_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22_AES+i] ,
		# # 
		CDE_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22_CDE+i] 
		#                               
		] ) )
	all_std_pc.append( np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],
		ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],
		ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],
		ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i],
		# 
		JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
		JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
		JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
		JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i],
		# 
		MDI_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22_MDI+i] ,
		# 
		MMA_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22_MMA+i],  
		# 
		CBO_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22_CBO+i], 
		# 
		PLA_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22_PLA+i] ,
		# 
		JPE_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22_JPE+i] ,
		# 
		MLM_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12_MLM+i],
		# MLM_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21_PLA+i],
		# MLM_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22_PLA+i] ,
		# #
		ALA_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21_ALA+i],
		ALA_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22_ALA+i] ,
		# # 
		# AES_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11_AES+i],
		# AES_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12_AES+i],
		# # AES_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21_AES+i],
		# # AES_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22_AES+i] ,
		# # 
		CDE_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22_CDE+i] 
		#                                    
		 ] ) )

# 68% confidence
all_mms_pc = [all_means_pc[i] -  all_std_pc[i] for i in range(len(all_means_pc))]
all_mps_pc = [all_means_pc[i] +  all_std_pc[i] for i in range(len(all_means_pc))]

# # 95% Confidence
# all_mms_pc = [all_means_pc[i] - 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]
# all_mps_pc = [all_means_pc[i] + 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]




initial_baseline = np.mean( all_means[:110] )
std_initial = np.mean( all_std[:110] )
print('Baseline LDF AU:',initial_baseline, '$\\pm$', std_initial)

initial_temp = np.mean( all_means_temp[:110] )
std_initial_temp = np.mean( all_std_temp[:110] )
print('Initial Temperature °C:',initial_temp, '$\\pm$', std_initial_temp)


first_occlusions = np.mean( [ all_means[60+110:160+110],all_means[350+110:450+110] ] )
std_fo = np.mean( [ np.mean(all_std[60+110:160+110]),np.mean(all_std[60+110:150+110]) ] )
print('First occlusion LDF AU:',first_occlusions, '$\\pm$', std_fo)

last_occlusions = np.mean( [ all_means[720+110:820+110],all_means[1125+110:1225+110] ] )
std_lo = np.mean( [ np.mean(all_std[720+110:820+110]),np.mean(all_std[1125+110:1225+110]) ] )
print('Last occlusion LDF AU:',last_occlusions, '$\\pm$', std_lo)


first_occlusions_pc = np.mean( [ all_means_pc[60:160],all_means_pc[350:450] ] )
std_fo_pc = np.mean( [ np.mean(all_std_pc[60:160]),np.mean(all_std_pc[60:150]) ] )
print('First occlusion LDF %:',first_occlusions_pc, '$\\pm$', std_fo_pc)

last_occlusions_pc = np.mean( [ all_means_pc[720:820],all_means_pc[1125:1225] ] )
std_lo_pc = np.mean( [ np.mean(all_std_pc[720:820]),np.mean(all_std_pc[1125:1225]) ] )
print('Last occlusion LDF %:',last_occlusions_pc, '$\\pm$', std_lo_pc)

idex0 = 1100
idex = 1460
t = np.linspace(idex0,idex,idex-idex0)

plt.figure(figsize=(10, 6))
plt.plot(t, ASE_mean_signal_raw[idex0:idex], label="Original Signal", alpha=0.5)
plt.savefig('./Figures/check_time_range_metrics.png')






plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-110:index_22+length_list], JOL_mean_signal_raw, color='black', label="P$_1$", alpha=1, linewidth=3)
# 
plt.plot(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-110:index_22_ASE+length_list_ASE], ASE_mean_signal_raw, color='blue', label="P$_2$", alpha=1, linewidth=3)
# 
plt.plot(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-110:index_12_MDI+length_list_MDI], MDI_mean_signal_raw, color='red', label="P$_3$", alpha=1, linewidth=3)
plt.plot(MMA_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12_CBO-110:index_12_CBO+length_list_CBO], CBO_mean_signal_raw, color='darkgreen', label="P$_5$", alpha=1, linewidth=3)
plt.plot(PLA_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12_PLA-110:index_12_PLA+length_list_PLA], PLA_mean_signal_raw, color='gold', label="P$_6$", alpha=1, linewidth=3)
plt.plot(ALA_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12_ALA-110:index_12_ALA+length_list_ALA], ALA_mean_signal_raw, color='fuchsia', label="P$_9$", alpha=0.5, linewidth=3)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 250])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_raw_h.jpg', bbox_inches='tight')
plt.close()

plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(MMA_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12_MMA-110:index_12_MMA+length_list_MMA], MMA_mean_signal_raw, linestyle='-', color='black', label="P$_4$", alpha=1, linewidth=3)
plt.plot(JPE_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12_JPE-110:index_12_JPE+length_list_JPE], JPE_mean_signal_raw, linestyle='-', color='blue', label="P$_7$", alpha=1, linewidth=3)
plt.plot(MLM_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12_MLM-110:index_12_MLM+length_list_MLM], MLM_mean_signal_raw, linestyle='-', color='red', label="P$_8$", alpha=1, linewidth=3)
# plt.plot(AES_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12_AES-110:index_12_AES+length_list_AES], AES_mean_signal_raw, linestyle='-', color='darkgreen', label="P$_10$", alpha=1, linewidth=3)
plt.plot(CDE_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12_CDE-110:index_12_CDE+length_list_CDE], CDE_mean_signal_raw, linestyle='-', color='gold', label="P$_{11}$", alpha=1, linewidth=3)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 250])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_raw_f.jpg', bbox_inches='tight')
plt.close()


plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22:index_22+length_list], JOL_mean_signal, color='black', label="P$_1$", alpha=1, linewidth=3)
# 
plt.plot(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE:index_22_ASE+length_list_ASE], ASE_mean_signal, color='blue', label="P$_2$", alpha=1, linewidth=3)
# 
plt.plot(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI:index_12_MDI+length_list_MDI], MDI_mean_signal, color='red', label="P$_3$", alpha=1, linewidth=3)
plt.plot(MMA_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12_CBO:index_12_CBO+length_list_CBO], CBO_mean_signal, color='darkgreen', label="P$_5$", alpha=1, linewidth=3)
plt.plot(PLA_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12_PLA:index_12_PLA+length_list_PLA], PLA_mean_signal, color='gold', label="P$_6$", alpha=1, linewidth=3)
plt.plot(ALA_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12_ALA:index_12_ALA+length_list_ALA], ALA_mean_signal, color='fuchsia', label="P$_9$", alpha=0.5, linewidth=3)
plt.legend()
plt.xlim([0, 780])
plt.ylim([0, 650])
plt.xlabel("Time [s]")
plt.ylabel("LDF [%]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_pc_h.jpg', bbox_inches='tight')
plt.close()


plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(MMA_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12_MMA:index_12_MMA+length_list_MMA], MMA_mean_signal, linestyle='-', color='black', label="P$_4$", alpha=1, linewidth=3)
plt.plot(JPE_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12_JPE:index_12_JPE+length_list_JPE], JPE_mean_signal, linestyle='-', color='blue', label="P$_7$", alpha=1, linewidth=3)
plt.plot(MLM_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12_MLM:index_12_MLM+length_list_MLM], MLM_mean_signal, linestyle='-', color='red', label="P$_8$", alpha=1, linewidth=3)
# plt.plot(AES_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12_AES:index_12_JPE+length_list_AES], AES_mean_signal, linestyle='-', color='darkgreen', label="P$_10$", alpha=1, linewidth=3)
plt.plot(CDE_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12_CDE:index_12_CDE+length_list_CDE], CDE_mean_signal, linestyle='-', color='gold', label="P$_{11}$", alpha=1, linewidth=3)
plt.legend()
plt.xlim([0, 780])
plt.ylim([0, 650])
plt.xlabel("Time [s]")
plt.ylabel("LDF [%]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_pc_f.jpg', bbox_inches='tight')
plt.close()





plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(JOL_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-110:index_22+length_list], JOL_Temperaturemean_signal_raw, color='black', label="P$_1$", alpha=1, linewidth=3)
# 
plt.plot(ASE_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22_ASE-110:index_22_ASE+length_list_ASE], ASE_Temperaturemean_signal_raw, color='blue', label="P$_2$", alpha=1, linewidth=3)
# 
plt.plot(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI-110:index_12_MDI+length_list_MDI], MDI_Temperaturemean_signal_raw, color='red', label="P$_3$", alpha=1, linewidth=3)
plt.plot(MMA_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12_CBO-110:index_12_CBO+length_list_CBO], CBO_Temperaturemean_signal_raw, color='darkgreen', label="P$_5$", alpha=1, linewidth=3)
plt.plot(PLA_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12_PLA-110:index_12_PLA+length_list_PLA], PLA_Temperaturemean_signal_raw, color='gold', label="P$_6$", alpha=1, linewidth=3)
plt.plot(ALA_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12_ALA-110:index_12_ALA+length_list_ALA], ALA_Temperaturemean_signal_raw, color='fuchsia', label="P$_9$", alpha=0.5, linewidth=3)
# plt.legend()
plt.xlim([-60, 780])
plt.ylim([24, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_h.jpg', bbox_inches='tight')
plt.close()


plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(MMA_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12_MMA-110:index_12_MMA+length_list_MMA], MMA_Temperaturemean_signal_raw, linestyle='-', color='black', label="P$_4$", alpha=1, linewidth=3)
plt.plot(JPE_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12_JPE-110:index_12_JPE+length_list_JPE], JPE_Temperaturemean_signal_raw, linestyle='-', color='blue', label="P$_7$", alpha=1, linewidth=3)
plt.plot(MLM_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12_MLM-110:index_12_MLM+length_list_MLM], MLM_Temperaturemean_signal_raw, linestyle='-', color='red', label="P$_8$", alpha=1, linewidth=3)
# plt.plot(AES_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12_AES-110:index_12_AES+length_list_AES], AES_Temperaturemean_signal_raw, linestyle='-', color='darkgreen', label="P$_10$", alpha=1, linewidth=3)
plt.plot(CDE_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12_CDE-110:index_12_CDE+length_list_CDE], CDE_Temperaturemean_signal_raw, linestyle='-', color='gold', label="P$_{11}$", alpha=1, linewidth=3)
# plt.legend()
plt.xlim([-60, 780])
plt.ylim([24, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_f.jpg', bbox_inches='tight')
plt.close()


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
plt.fill_between(MDI_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12_MDI:index_12_MDI+length_list_MDI], MDI_mms, MDI_mps , linestyle='-', color='g', alpha=0.2)
# 
plt.plot(CBO_PREVAIL_LBM_029["j1_2"]["temps (s)"][index_12_CBO:index_12_CBO+length_list_CBO], CBO_mean_signal , linestyle='-', color='turquoise', label="P$_5$", alpha=0.6)
plt.fill_between(CBO_PREVAIL_LBM_029["j1_2"]["temps (s)"][index_12_CBO:index_12_CBO+length_list_CBO], CBO_mms, CBO_mps , linestyle='-', color='turquoise', alpha=0.2)
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
plt.savefig('./Figures/all_mean_std.jpg', bbox_inches='tight')
plt.close()














# Raw

# A compléter

all_means=[]
all_means_temp = []
all_std=[]
all_std_temp =[]
# for i in range(min(length_list_ASE+110,length_list+110,length_list_MDI+110,length_list_MMA+110,length_list_CBO+110,length_list_PLA+110,length_list_JPE+110,length_list_MLM+110,length_list_ALA+110,length_list_AES+110,length_list_CDE+110)):
for i in range(min(length_list_ASE+110,length_list+110,length_list_MDI+110,length_list_CBO+110,length_list_PLA+110,length_list_ALA+110)):
	all_means.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-110],
		# 
		JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
		JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
		JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
		JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110],
		# 
		MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-110],
		# 
		# 
		CBO_PREVAIL_LBM_029["j1_1"]["PU"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["PU"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["PU"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["PU"][index_22_CBO+i-110], 
		# 
		PLA_PREVAIL_LBM_030["j1_1"]["PU"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["PU"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["PU"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["PU"][index_22_PLA+i-110] ,
		# 
		# 
		# #
		ALA_PREVAIL_LBM_033["j1_1"]["PU"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["PU"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["PU"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["PU"][index_22_ALA+i-110] 
		# 
		] ) )
	# 
	all_std.append( np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["PU"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["PU"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["PU"][index_22_ASE+i-110],
		# 
		JOL_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
		JOL_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
		JOL_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
		JOL_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110],
		# 
		MDI_PREVAIL_LBM_026["j1_1"]["PU"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["PU"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["PU"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["PU"][index_22_MDI+i-110] ,
		# 
		CBO_PREVAIL_LBM_029["j1_1"]["PU"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["PU"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["PU"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["PU"][index_22_CBO+i-110], 
		# 
		PLA_PREVAIL_LBM_030["j1_1"]["PU"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["PU"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["PU"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["PU"][index_22_PLA+i-110] ,
		# 
		JPE_PREVAIL_LBM_031["j1_1"]["PU"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["PU"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["PU"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["PU"][index_22_JPE+i-110] ,
		# 
		# #
		ALA_PREVAIL_LBM_033["j1_1"]["PU"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["PU"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["PU"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["PU"][index_22_ALA+i-110] 
		] ) )
	# 
	all_means_temp.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-110],
		# 
		JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
		JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
		JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
		JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110],
		# 
		MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-110],
		# 
		CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i-110], 
		# 
		PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i-110] ,
		# 
		# 
		# #
		ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i-110] 
		] ) )
	# 
	all_std_temp.append( np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11_ASE+i-110],
		ASE_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21_ASE+i-110],
		ASE_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22_ASE+i-110],
		# 
		JOL_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
		JOL_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
		JOL_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
		JOL_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110],
		# 
		MDI_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11_MDI+i-110],
		MDI_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21_MDI+i-110],
		MDI_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22_MDI+i-110] ,
		# 
		# 
		CBO_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11_CBO+i-110],
		CBO_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21_CBO+i-110],
		CBO_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22_CBO+i-110], 
		# 
		PLA_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11_PLA+i-110],
		PLA_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21_PLA+i-110],
		PLA_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22_PLA+i-110] ,
		# #
		ALA_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11_ALA+i-110],
		ALA_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21_ALA+i-110],
		ALA_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22_ALA+i-110] 
		# 
		] ) )

all_mms = [all_means[i]-all_std[i] for i in range(len(all_means))]
all_mps = [all_means[i]+all_std[i] for i in range(len(all_means))]

all_means_pc=[]
all_std_pc=[]
# Percentage
# for i in range(min(length_list_ASE,length_list,length_list_MDI,length_list_MMA,length_list_CBO,length_list_PLA,length_list_JPE,length_list_MLM,length_list_ALA,length_list_AES,length_list_CDE)):
for i in range(min(length_list_ASE,length_list,length_list_MDI,length_list_CBO,length_list_PLA,length_list_ALA+110)):
	all_means_pc.append( np.mean( [ ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],
		ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],
		ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],
		ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i],
		# 
		JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
		JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
		JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
		JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i],
		# 
		MDI_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22_MDI+i]  ,
		# 
		CBO_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22_CBO+i], 
		# 
		PLA_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22_PLA+i] ,
		# 
		# #
		ALA_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21_ALA+i],
		# ALA_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22_ALA+i] 
		#                               
		] ) )
	all_std_pc.append( np.std( [ ASE_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11_ASE+i],
		ASE_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12_ASE+i],
		ASE_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21_ASE+i],
		ASE_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22_ASE+i],
		# 
		JOL_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
		JOL_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
		JOL_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
		JOL_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i],
		# 
		MDI_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11_MDI+i],
		MDI_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12_MDI+i],
		MDI_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21_MDI+i],
		MDI_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22_MDI+i] ,
		# 
		# 
		CBO_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11_CBO+i],
		CBO_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12_CBO+i],
		CBO_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21_CBO+i],
		CBO_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22_CBO+i], 
		# 
		PLA_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11_PLA+i],
		PLA_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12_PLA+i],
		PLA_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21_PLA+i],
		PLA_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22_PLA+i] ,
		# 
		# 
		# #
		ALA_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11_ALA+i],
		ALA_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12_ALA+i],
		ALA_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21_ALA+i],
		# ALA_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22_ALA+i] 
		# # 
		#                                    
		 ] ) )

# 68% confidence
all_mms_pc = [all_means_pc[i] -  all_std_pc[i] for i in range(len(all_means_pc))]
all_mps_pc = [all_means_pc[i] +  all_std_pc[i] for i in range(len(all_means_pc))]

# # 95% Confidence
# all_mms_pc = [all_means_pc[i] - 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]
# all_mps_pc = [all_means_pc[i] + 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]




initial_baseline = np.mean( all_means[:110] )
std_initial = np.mean( all_std[:110] )
print('Baseline LDF AU: M',initial_baseline, '$\\pm$', std_initial)

initial_temp = np.mean( all_means_temp[:110] )
std_initial_temp = np.mean( all_std_temp[:110] )
print('Initial Temperature °C: M',initial_temp, '$\\pm$', std_initial_temp)


first_occlusions = np.mean( [ all_means[60+110:160+110],all_means[350+110:450+110] ] )
std_fo = np.mean( [ np.mean(all_std[60+110:160+110]),np.mean(all_std[60+110:150+110]) ] )
print('First occlusion LDF AU: M',first_occlusions, '$\\pm$', std_fo)

last_occlusions = np.mean( [ all_means[720+110:820+110],all_means[1125+110:1225+110] ] )
std_lo = np.mean( [ np.mean(all_std[720+110:820+110]),np.mean(all_std[1125+110:1225+110]) ] )
print('Last occlusion LDF AU: M',last_occlusions, '$\\pm$', std_lo)


first_occlusions_pc = np.mean( [ all_means_pc[60:160],all_means_pc[350:450] ] )
std_fo_pc = np.mean( [ np.mean(all_std_pc[60:160]),np.mean(all_std_pc[60:150]) ] )
print('First occlusion LDF M %:',first_occlusions_pc, '$\\pm$', std_fo_pc)

last_occlusions_pc = np.mean( [ all_means_pc[720:820],all_means_pc[1125:1225] ] )
std_lo_pc = np.mean( [ np.mean(all_std_pc[720:820]),np.mean(all_std_pc[1125:1225]) ] )
print('Last occlusion LDF M %:',last_occlusions_pc, '$\\pm$', std_lo_pc)











# Raw

# A compléter

all_means=[]
all_means_temp = []
all_std=[]
all_std_temp =[]
# for i in range(min(length_list_ASE+110,length_list+110,length_list_MDI+110,length_list_MMA+110,length_list_CBO+110,length_list_PLA+110,length_list_JPE+110,length_list_MLM+110,length_list_ALA+110,length_list_AES+110,length_list_CDE+110)):
for i in range(min(length_list_MMA+110,length_list_JPE+110,length_list_MLM+110,length_list_CDE+110)):
	all_means.append( np.mean( [ 
		MMA_PREVAIL_LBM_028["j1_1"]["PU"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["PU"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["PU"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["PU"][index_22_MMA+i-110],  
		# 
		JPE_PREVAIL_LBM_031["j1_1"]["PU"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["PU"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["PU"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["PU"][index_22_JPE+i-110] , 
		# 
		MLM_PREVAIL_LBM_032["j1_1"]["PU"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["PU"][index_12_MLM+i-110],
		# MLM_PREVAIL_LBM_032["j2_1"]["PU"][index_21_PLA+i-110],
		# MLM_PREVAIL_LBM_032["j2_2"]["PU"][index_22_PLA+i-110] ,
		# #
		# # 
		# AES_PREVAIL_LBM_034["j1_1"]["PU"][index_11_AES+i-110],
		# AES_PREVAIL_LBM_034["j1_2"]["PU"][index_12_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_1"]["PU"][index_21_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_2"]["PU"][index_22_AES+i-110] ,
		# # 
		CDE_PREVAIL_LBM_035["j1_1"]["PU"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["PU"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["PU"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["PU"][index_22_CDE+i-110] 
		# 
		] ) )
	# 
	all_std.append( np.std( [ 
		# 
		MMA_PREVAIL_LBM_028["j1_1"]["PU"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["PU"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["PU"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["PU"][index_22_MMA+i-110],  
		# 
		JPE_PREVAIL_LBM_031["j1_1"]["PU"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["PU"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["PU"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["PU"][index_22_JPE+i-110] ,
		# 
		# 
		MLM_PREVAIL_LBM_032["j1_1"]["PU"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["PU"][index_12_MLM+i-110],
		# MLM_PREVAIL_LBM_032["j2_1"]["PU"][index_21_PLA+i-110],
		# MLM_PREVAIL_LBM_032["j2_2"]["PU"][index_22_PLA+i-110] ,
		# #
		# # 
		# AES_PREVAIL_LBM_034["j1_1"]["PU"][index_11_AES+i-110],
		# AES_PREVAIL_LBM_034["j1_2"]["PU"][index_12_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_1"]["PU"][index_21_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_2"]["PU"][index_22_AES+i-110] ,
		# # 
		CDE_PREVAIL_LBM_035["j1_1"]["PU"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["PU"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["PU"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["PU"][index_22_CDE+i-110] 
		] ) )
	# 
	all_means_temp.append( np.mean( [ 
		# 
		MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i-110],  
		# 
		# 
		MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i-110],
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_PLA+i-110],
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_PLA+i-110] ,
		# #
		# # 
		# AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i-110],
		# AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i-110] ,
		# # 
		CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i-110] 
		] ) )
	# 
	all_std_temp.append( np.std( [ 
		# 
		MMA_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11_MMA+i-110],
		MMA_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21_MMA+i-110],
		MMA_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22_MMA+i-110],  
		# 
		# 
		JPE_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11_JPE+i-110],
		JPE_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21_JPE+i-110],
		JPE_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22_JPE+i-110] ,
		# 
		MLM_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11_MLM+i-110],
		MLM_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12_MLM+i-110],
		# MLM_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21_PLA+i-110],
		# MLM_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22_PLA+i-110] ,
		# #
		# # 
		# AES_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11_AES+i-110],
		# AES_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21_AES+i-110],
		# # AES_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22_AES+i-110] ,
		# # 
		CDE_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11_CDE+i-110],
		CDE_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21_CDE+i-110],
		CDE_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22_CDE+i-110] 
		# 
		] ) )

all_mms = [all_means[i]-all_std[i] for i in range(len(all_means))]
all_mps = [all_means[i]+all_std[i] for i in range(len(all_means))]

all_means_pc=[]
all_std_pc=[]
# Percentage
# for i in range(min(length_list_ASE,length_list,length_list_MDI,length_list_MMA,length_list_CBO,length_list_PLA,length_list_JPE,length_list_MLM,length_list_ALA,length_list_AES,length_list_CDE)):
for i in range(min(length_list_MMA,length_list_JPE,length_list_MLM,length_list_CDE)):
	all_means_pc.append( np.mean( [ 
		# 
		MMA_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22_MMA+i],  
		# 
		JPE_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22_JPE+i] ,
		# 
		MLM_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12_MLM+i],
		# MLM_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21_PLA+i],
		# MLM_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22_PLA+i] ,
		# #
		# # 
		# AES_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11_AES+i],
		# AES_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12_AES+i],
		# # AES_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21_AES+i],
		# # AES_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22_AES+i] ,
		# # 
		CDE_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22_CDE+i] 
		#                               
		] ) )
	all_std_pc.append( np.std( [ 
		# 
		# 
		# 
		MMA_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11_MMA+i],
		MMA_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12_MMA+i],
		MMA_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21_MMA+i],
		MMA_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22_MMA+i],  
		# 
		# 
		# 
		JPE_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11_JPE+i],
		JPE_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12_JPE+i],
		JPE_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21_JPE+i],
		JPE_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22_JPE+i] ,
		# 
		MLM_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11_MLM+i],
		MLM_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12_MLM+i],
		# MLM_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21_PLA+i],
		# MLM_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22_PLA+i] ,
		# #
		# # 
		# AES_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11_AES+i],
		# AES_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12_AES+i],
		# # AES_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21_AES+i],
		# # AES_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22_AES+i] ,
		# # 
		CDE_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11_CDE+i],
		CDE_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12_CDE+i],
		CDE_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21_CDE+i],
		CDE_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22_CDE+i] 
		#                                    
		 ] ) )

# 68% confidence
all_mms_pc = [all_means_pc[i] -  all_std_pc[i] for i in range(len(all_means_pc))]
all_mps_pc = [all_means_pc[i] +  all_std_pc[i] for i in range(len(all_means_pc))]

# # 95% Confidence
# all_mms_pc = [all_means_pc[i] - 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]
# all_mps_pc = [all_means_pc[i] + 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]




initial_baseline = np.mean( all_means[:110] )
std_initial = np.mean( all_std[:110] )
print('Baseline LDF AU F:',initial_baseline, '$\\pm$', std_initial)

initial_temp = np.mean( all_means_temp[:110] )
std_initial_temp = np.mean( all_std_temp[:110] )
print('Initial Temperature °C F:',initial_temp, '$\\pm$', std_initial_temp)


first_occlusions = np.mean( [ all_means[60+110:160+110],all_means[350+110:450+110] ] )
std_fo = np.mean( [ np.mean(all_std[60+110:160+110]),np.mean(all_std[60+110:150+110]) ] )
print('First occlusion LDF AU F:',first_occlusions, '$\\pm$', std_fo)

last_occlusions = np.mean( [ all_means[720+110:820+110],all_means[1125+110:1225+110] ] )
std_lo = np.mean( [ np.mean(all_std[720+110:820+110]),np.mean(all_std[1125+110:1225+110]) ] )
print('Last occlusion LDF AU F:',last_occlusions, '$\\pm$', std_lo)


first_occlusions_pc = np.mean( [ all_means_pc[60:160],all_means_pc[350:450] ] )
std_fo_pc = np.mean( [ np.mean(all_std_pc[60:160]),np.mean(all_std_pc[60:150]) ] )
print('First occlusion LDF F % :',first_occlusions_pc, '$\\pm$', std_fo_pc)

last_occlusions_pc = np.mean( [ all_means_pc[720:820],all_means_pc[1125:1225] ] )
std_lo_pc = np.mean( [ np.mean(all_std_pc[720:820]),np.mean(all_std_pc[1125:1225]) ] )
print('Last occlusion LDF F %:',last_occlusions_pc, '$\\pm$', std_lo_pc)