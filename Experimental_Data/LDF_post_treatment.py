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
    y = filtfilt(b, a, data)  # Apply filter with zero ph_009 shift
    return y

# Load the Excel file
file_path = "008_PREVAIL_LBM_008.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")

# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
_008_PREVAIL_LBM_008 = {}
for sheet_name, data in all_sheets.items():
    _008_PREVAIL_LBM_008[sheet_name] = {col: data[col] for col in column_names if col in data.columns}


index_11 = 0
for ii in range(len(_008_PREVAIL_LBM_008["j1_1"]["temps (s)"] ) ):
	if _008_PREVAIL_LBM_008["j1_1"]["temps (s)"][ii] == 0:
		index_11 = ii 
		break
index_12 = 0
for ii in range(len(_008_PREVAIL_LBM_008["j1_2"]["temps (s)"] ) ):
	if _008_PREVAIL_LBM_008["j1_2"]["temps (s)"][ii] == 0:
		index_12 = ii 
		break
index_21 = 0
for ii in range(len(_008_PREVAIL_LBM_008["j2_1"]["temps (s)"] ) ):
	if _008_PREVAIL_LBM_008["j2_1"]["temps (s)"][ii] == 0:
		index_21 = ii 
		break
index_22 = 0
for ii in range(len(_008_PREVAIL_LBM_008["j2_2"]["temps (s)"] ) ):
	if _008_PREVAIL_LBM_008["j2_2"]["temps (s)"][ii] == 0:
		index_22 = ii 
		break		

length_list = min(len(_008_PREVAIL_LBM_008["j1_1"]["temps (s)"][index_11:]) , len(_008_PREVAIL_LBM_008["j1_2"]["temps (s)"][index_12:]) , len(_008_PREVAIL_LBM_008["j2_1"]["temps (s)"][index_21:]) , len(_008_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22:] ) )

_008_Temperaturemean_signal=[]
_008_Temperaturemps = []
_008_Temperaturemms = []
_008_Temperaturestd_signal = []
_008_mean_signal=[]
_008_mps = []
_008_mms = []
_008_std_signal = []

for i in range(length_list):
	_008_mean_signal.append( np.mean( [ _008_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
									   _008_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
									   _008_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
									   _008_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i] ] ) )
	# 
	_008_std_signal.append( np.std( [ _008_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
									 _008_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
									 _008_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
									 _008_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i] ] ) )
	# 
	_008_mms.append( np.mean( [ _008_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
							   _008_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
							   _008_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
							   _008_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i] ] ) 
					- 1.96 * np.std( [ _008_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
								_008_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
								_008_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
								_008_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i] ] ) )
	# 
	_008_mps.append( np.mean( [ _008_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
								_008_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
								_008_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
								_008_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i] ] )
					+ 1.96* np.std( [ _008_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
										_008_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
										_008_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
										_008_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i] ] ) )
	# 
	_008_Temperaturemean_signal.append( np.mean( [ _008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],
												  _008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],
												  _008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],
												  _008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i] ] ) )
	# 
	_008_Temperaturestd_signal.append( np.std( [ _008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],
												_008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],
												_008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],
												_008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i] ] ) )
	# 
	_008_Temperaturemms.append( np.mean( [ _008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],
										  _008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],
										  _008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],
										  _008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i] ] )
								- 1.96* np.std( [ _008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],
									_008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],
									_008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],
									_008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i] ] ) )
	# 
	_008_Temperaturemps.append( np.mean( [ _008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],
										  _008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],
										  _008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],
										  _008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i] ] ) 
								+ 1.96 * np.std( [ _008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i],
													_008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i],
													_008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i],
													_008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i] ] ) )

# Filter parameters
fs = 1000  # Sampling frequency in Hz
cutoff = 25  # Cutoff frequency in Hz
_008_filtered_signal = butter_lowpass_filter(_008_mean_signal, cutoff, fs)


_008_Temperaturemean_signal_raw=[]
_008_Temperaturemps_raw = []
_008_Temperaturemms_raw = []
_008_Temperaturestd_signal_raw = []
_008_mean_signal_raw=[]
_008_mps_raw = []
_008_mms_raw = []
_008_std_signal_raw = []

for i in range(length_list+110):
	_008_mean_signal_raw.append( np.mean( [ _008_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
											_008_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
											_008_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
											_008_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110] ] ) )
	# 
	_008_std_signal_raw.append( np.std( [ _008_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
										_008_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
										_008_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
										_008_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110] ] ) )
	# 
	_008_mms_raw.append( np.mean( [ _008_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
									_008_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
									_008_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
									_008_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110] ] ) 
						- np.std( [_008_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
										_008_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
										_008_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
										_008_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110] ] ) )
	# 
	_008_mps_raw.append( np.mean( [ _008_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
									_008_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
									_008_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
									_008_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110] ] )
						+ np.std( [ _008_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
											_008_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
											_008_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
											_008_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110] ] ) )
	# 
	_008_Temperaturemean_signal_raw.append( np.mean( [ _008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
													_008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
													_008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
													_008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110] ] ) )
	# 
	_008_Temperaturestd_signal_raw.append( np.std( [ _008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
													_008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
													_008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
													_008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110] ] ) )
	# 
	_008_Temperaturemms_raw.append( np.mean( [ _008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
											  _008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
											  _008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
											  _008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110] ] )
									-np.std( [ _008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
										_008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
										_008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
										_008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110] ] ) )
	# 
	_008_Temperaturemps_raw.append( np.mean( [ _008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
												_008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
												_008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
												_008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110] ] )
									+np.std( [ _008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
														_008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
														_008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
														_008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110] ] ) )


plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(_008_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-110:index_22+length_list], _008_mean_signal_raw, color='black', label="P$_1$", alpha=0.5)
plt.fill_between(_008_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-110:index_22+length_list], _008_mms_raw, _008_mps_raw, color='black', alpha=0.2)
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
plt.plot(_008_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-110:index_22+length_list], _008_Temperaturemean_signal_raw, color='black', label="P$_1$", alpha=0.5)
plt.fill_between(_008_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-110:index_22+length_list], _008_Temperaturemms_raw, _008_Temperaturemps_raw, color='black', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P1.jpg', bbox_inches='tight')



plt.close()












# Load the Excel file
file_path = "009_PREVAIL_LBM_009.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
_009_PREVAIL_LBM_009 = {}
for sheet_name, data in all_sheets.items():
    _009_PREVAIL_LBM_009[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11__009 = 0
for ii in range(len(_009_PREVAIL_LBM_009["j1_1"]["temps (s)"] ) ):
	if _009_PREVAIL_LBM_009["j1_1"]["temps (s)"][ii] == 0:
		index_11__009 = ii 
		break
index_12__009 = 0
for ii in range(len(_009_PREVAIL_LBM_009["j1_2"]["temps (s)"] ) ):
	if _009_PREVAIL_LBM_009["j1_2"]["temps (s)"][ii] == 0:
		index_12__009 = ii 
		break
index_21__009 = 0
for ii in range(len(_009_PREVAIL_LBM_009["j2_1"]["temps (s)"] ) ):
	if _009_PREVAIL_LBM_009["j2_1"]["temps (s)"][ii] == 0:
		index_21__009 = ii 
		break
index_22__009 = 0
for ii in range(len(_009_PREVAIL_LBM_009["j2_2"]["temps (s)"] ) ):
	if _009_PREVAIL_LBM_009["j2_2"]["temps (s)"][ii] == 0:
		index_22__009 = ii 
		break		

length_list__009 = min(len(_009_PREVAIL_LBM_009["j1_1"]["temps (s)"][index_11__009:]) , len(_009_PREVAIL_LBM_009["j1_2"]["temps (s)"][index_12__009:]) , len(_009_PREVAIL_LBM_009["j2_1"]["temps (s)"][index_21__009:]) , len(_009_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22__009:] ) )

_009_Temperaturemean_signal=[]
_009_Temperaturemps = []
_009_Temperaturemms = []
_009_Temperaturestd_signal = []
_009_mean_signal=[]
_009_mps = []
_009_mms = []
_009_std_signal = []

for i in range(length_list__009):
	_009_mean_signal.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11__009+i],
										_009_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12__009+i],
										_009_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21__009+i],
										_009_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22__009+i] ] ) )
	# 
	_009_std_signal.append( np.std( [ _009_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11__009+i],
										_009_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12__009+i],
										_009_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21__009+i],
										_009_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22__009+i] ] ) )
	# 
	_009_mms.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11__009+i],
								_009_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12__009+i],
								_009_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21__009+i],
								_009_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22__009+i] ] )
					- 1.96* np.std( [ _009_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11__009+i],
								_009_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12__009+i],
								_009_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21__009+i],
								_009_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22__009+i] ] ) )
	# 
	_009_mps.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11__009+i],
								_009_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12__009+i],
								_009_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21__009+i],
								_009_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22__009+i] ] )
					+1.96*np.std( [ _009_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11__009+i],
								_009_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12__009+i],
								_009_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21__009+i],
								_009_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22__009+i] ] ) )
	# 
	_009_Temperaturemean_signal.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i],
													_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i],
													_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i],
													_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i] ] ) )
	# 
	_009_Temperaturestd_signal.append( np.std( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i],
												_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i],
												_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i],
												_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i] ] ) )
	# 
	_009_Temperaturemms.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i],
											_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i],
											_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i],
											_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i] ] )
								-1.96*np.std( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i],
											_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i],
											_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i],
											_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i] ] ) )
							# 
	_009_Temperaturemps.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i],
											_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i],
											_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i],
											_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i] ] )
							+1.96*np.std( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i],
											_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i],
											_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i],
											_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i] ] ) )

# Filter parameters
fs = 1000  # Sampling frequency in Hz
cutoff = 25  # Cutoff frequency in Hz
_009_filtered_signal = butter_lowpass_filter(_009_mean_signal, cutoff, fs)



_009_Temperaturemean_signal_raw=[]
_009_Temperaturemps_raw = []
_009_Temperaturemms_raw = []
_009_Temperaturestd_signal_raw = []
_009_mean_signal_raw=[]
_009_mps_raw = []
_009_mms_raw = []
_009_std_signal_raw = []

for i in range(length_list__009+110):
	_009_mean_signal_raw.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["PU"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["PU"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["PU"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["PU"][index_22__009+i-110] ] ) )
	# 
	_009_std_signal_raw.append( np.std( [ _009_PREVAIL_LBM_009["j1_1"]["PU"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["PU"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["PU"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["PU"][index_22__009+i-110] ] ) )
	# 
	_009_mms_raw.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["PU"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["PU"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["PU"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["PU"][index_22__009+i-110] ] )
	-np.std( [ _009_PREVAIL_LBM_009["j1_1"]["PU"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["PU"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["PU"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["PU"][index_22__009+i-110] ] ) )
	# 
	_009_mps_raw.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["PU"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["PU"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["PU"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["PU"][index_22__009+i-110] ] )
	+np.std( [ _009_PREVAIL_LBM_009["j1_1"]["PU"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["PU"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["PU"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["PU"][index_22__009+i-110] ] ) )
	# 
	_009_Temperaturemean_signal_raw.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i-110] ] ) )
	# 
	_009_Temperaturestd_signal_raw.append( np.std( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i-110] ] ) )
	# 
	_009_Temperaturemms_raw.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i-110] ] )
	-np.std( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i-110] ] ) )
	# 
	_009_Temperaturemps_raw.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i-110] ] )
	+np.std( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i-110] ] ) )





plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# 
plt.plot(_009_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22__009-110:index_22__009+length_list__009], _009_mean_signal_raw, color='b', label="P$_2$", alpha=0.5)
plt.fill_between(_009_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22__009-110:index_22__009+length_list__009], _009_mms_raw, _009_mps_raw, color='b', alpha=0.2)
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
plt.plot(_009_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22__009-110:index_22__009+length_list__009], _009_Temperaturemean_signal_raw, color='b', label="P$_2$", alpha=0.5)
plt.fill_between(_009_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22__009-110:index_22__009+length_list__009], _009_Temperaturemms_raw, _009_Temperaturemps_raw, color='b', alpha=0.2)
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
file_path = "026_PREVAIL_LBM_026.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
_026_PREVAIL_LBM_026 = {}
for sheet_name, data in all_sheets.items():
    _026_PREVAIL_LBM_026[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11__026 = 0
for ii in range(len(_026_PREVAIL_LBM_026["j1_1"]["temps (s)"] ) ):
	if _026_PREVAIL_LBM_026["j1_1"]["temps (s)"][ii] == 0:
		index_11__026 = ii 
		break
index_12__026 = 0
for ii in range(len(_026_PREVAIL_LBM_026["j1_2"]["temps (s)"] ) ):
	if _026_PREVAIL_LBM_026["j1_2"]["temps (s)"][ii] == 0:
		index_12__026 = ii 
		break
index_21__026 = 0
for ii in range(len(_026_PREVAIL_LBM_026["j2_1"]["temps (s)"] ) ):
	if _026_PREVAIL_LBM_026["j2_1"]["temps (s)"][ii] == 0:
		index_21__026 = ii 
		break
index_22__026 = 0
for ii in range(len(_026_PREVAIL_LBM_026["j2_2"]["temps (s)"] ) ):
	if _026_PREVAIL_LBM_026["j2_2"]["temps (s)"][ii] == 0:
		index_22__026 = ii 
		break		

length_list__026 = min(len(_026_PREVAIL_LBM_026["j1_1"]["temps (s)"][index_11__026:]) , len(_026_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12__026:]) , len(_026_PREVAIL_LBM_026["j2_1"]["temps (s)"][index_21__026:]) , len(_026_PREVAIL_LBM_026["j2_2"]["temps (s)"][index_22__026:] ) )

_026_Temperaturemean_signal=[]
_026_Temperaturemps = []
_026_Temperaturemms = []
_026_Temperaturestd_signal = []
_026_mean_signal=[]
_026_mps = []
_026_mms = []
_026_std_signal = []

for i in range(length_list__026):
	_026_mean_signal.append( np.mean( [ _026_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22__026+i] ] ) )
	# 
	_026_std_signal.append( np.std( [ _026_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22__026+i] ] ) )
	# 
	_026_mms.append( np.mean( [ _026_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22__026+i] ] )
	- 1.96* np.std( [ _026_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22__026+i] ] ) )
	# 
	_026_mps.append( np.mean( [ _026_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22__026+i] ] )
	+1.96*np.std( [ _026_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22__026+i] ] ) )
	# 
	_026_Temperaturemean_signal.append( np.mean( [ _026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i] ] ) )
	# 
	_026_Temperaturestd_signal.append( np.std( [ _026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i] ] ) )
	# 
	_026_Temperaturemms.append( np.mean( [ _026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i] ] )
	-1.96*np.std( [ _026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i] ] ) )
	# 
	_026_Temperaturemps.append( np.mean( [ _026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i] ] )
	+1.96*np.std( [ _026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i] ] ) )
	
_026_Temperaturemean_signal_raw=[]
_026_Temperaturemps_raw = []
_026_Temperaturemms_raw = []
_026_Temperaturestd_signal_raw = []
_026_mean_signal_raw=[]
_026_mps_raw = []
_026_mms_raw = []
_026_std_signal_raw = []

for i in range(length_list__026+110):
	_026_mean_signal_raw.append( np.mean( [ _026_PREVAIL_LBM_026["j1_1"]["PU"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["PU"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["PU"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["PU"][index_22__026+i-110] ] ) )
	# 
	_026_std_signal_raw.append( np.std( [ _026_PREVAIL_LBM_026["j1_1"]["PU"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["PU"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["PU"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["PU"][index_22__026+i-110] ] ) )
	# 
	_026_mms_raw.append( np.mean( [ _026_PREVAIL_LBM_026["j1_1"]["PU"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["PU"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["PU"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["PU"][index_22__026+i-110] ] )
	-np.std( [ _026_PREVAIL_LBM_026["j1_1"]["PU"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["PU"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["PU"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["PU"][index_22__026+i-110] ] ) )
	# 
	_026_mps_raw.append( np.mean( [ _026_PREVAIL_LBM_026["j1_1"]["PU"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["PU"][index_12__026+i-110] ,
		_026_PREVAIL_LBM_026["j2_1"]["PU"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["PU"][index_22__026+i-110] ] )
	+np.std( [ _026_PREVAIL_LBM_026["j1_1"]["PU"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["PU"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["PU"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["PU"][index_22__026+i-110] ] ) )
	# 
	_026_Temperaturemean_signal_raw.append( np.mean( [ _026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i-110] ] ) )
	# 
	_026_Temperaturestd_signal_raw.append( np.std( [ _026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i-110] ] ) )
	# 
	_026_Temperaturemms_raw.append( np.mean( [ _026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i-110] ] )
	-np.std( [ _026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i-110] ] ) )
	# 
	_026_Temperaturemps_raw.append( np.mean( [ _026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i-110] ] )
	+np.std( [ _026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i-110] ] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(_026_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12__026-110:index_12__026+length_list__026], _026_mean_signal_raw, color='r', label="P$_3$", alpha=0.5)
plt.fill_between(_026_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12__026-110:index_12__026+length_list__026], _026_mms_raw, _026_mps_raw, color='r', alpha=0.2)
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
plt.plot(_026_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12__026-110:index_12__026+length_list__026], _026_Temperaturemean_signal_raw, color='r', label="P$_3$", alpha=0.5)
plt.fill_between(_026_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12__026-110:index_12__026+length_list__026], _026_Temperaturemms_raw, _026_Temperaturemps_raw, color='r', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P3.jpg', bbox_inches='tight')





plt.close()









































# Load the Excel file
file_path = "028_PREVAIL_LBM_028.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
_028_PREVAIL_LBM_028 = {}
for sheet_name, data in all_sheets.items():
    _028_PREVAIL_LBM_028[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11__028 = 0
for ii in range(len(_028_PREVAIL_LBM_028["j1_1"]["temps (s)"] ) ):
	if _028_PREVAIL_LBM_028["j1_1"]["temps (s)"][ii] == 0:
		index_11__028 = ii 
		break
index_12__028 = 0
for ii in range(len(_028_PREVAIL_LBM_028["j1_2"]["temps (s)"] ) ):
	if _028_PREVAIL_LBM_028["j1_2"]["temps (s)"][ii] == 0:
		index_12__028 = ii 
		break
index_21__028 = 0
for ii in range(len(_028_PREVAIL_LBM_028["j2_1"]["temps (s)"] ) ):
	if _028_PREVAIL_LBM_028["j2_1"]["temps (s)"][ii] == 0:
		index_21__028 = ii 
		break
index_22__028 = 0
for ii in range(len(_028_PREVAIL_LBM_028["j2_2"]["temps (s)"] ) ):
	if _028_PREVAIL_LBM_028["j2_2"]["temps (s)"][ii] == 0:
		index_22__028 = ii 
		break		

length_list__028 = min(len(_028_PREVAIL_LBM_028["j1_1"]["temps (s)"][index_11__028:]) , len(_028_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12__028:]) , len(_028_PREVAIL_LBM_028["j2_1"]["temps (s)"][index_21__028:]) , len(_028_PREVAIL_LBM_028["j2_2"]["temps (s)"][index_22__028:] ) )

_028_Temperaturemean_signal=[]
_028_Temperaturemps = []
_028_Temperaturemms = []
_028_Temperaturestd_signal = []
_028_mean_signal=[]
_028_mps = []
_028_mms = []
_028_std_signal = []

for i in range(length_list__028):
	_028_mean_signal.append( np.mean( [ _028_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22__028+i] 
		] ) )
	# 
	_028_std_signal.append( np.std( [ _028_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22__028+i] 
		] ) )
	# 
	_028_mms.append( np.mean( [ _028_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22__028+i] 
		] )
	- 1.96* np.std( [ _028_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22__028+i] 
		] ) )
	# 
	_028_mps.append( np.mean( [ _028_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22__028+i] 
		] )
	+1.96*np.std( [ _028_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22__028+i] 
		] ) )
	# 
	_028_Temperaturemean_signal.append( np.mean( [ _028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i] 
		] ) )
	# 
	_028_Temperaturestd_signal.append( np.std( [ _028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i] 
		] ) )
	# 
	_028_Temperaturemms.append( np.mean( [ _028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i] 
		] )
	-1.96*np.std( [ _028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i] 
		] ) )
	# 
	_028_Temperaturemps.append( np.mean( [ _028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i] 
		] )
	+1.96*np.std( [ _028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i] 
		] ) )
	
_028_Temperaturemean_signal_raw=[]
_028_Temperaturemps_raw = []
_028_Temperaturemms_raw = []
_028_Temperaturestd_signal_raw = []
_028_mean_signal_raw=[]
_028_mps_raw = []
_028_mms_raw = []
_028_std_signal_raw = []

for i in range(length_list__028+110):
	_028_mean_signal_raw.append( np.mean( [ _028_PREVAIL_LBM_028["j1_1"]["PU"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["PU"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["PU"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["PU"][index_22__028+i-110] 
		] ) )
	# 
	_028_std_signal_raw.append( np.std( [ _028_PREVAIL_LBM_028["j1_1"]["PU"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["PU"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["PU"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["PU"][index_22__028+i-110] 
		] ) )
	# 
	_028_mms_raw.append( np.mean( [ _028_PREVAIL_LBM_028["j1_1"]["PU"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["PU"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["PU"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["PU"][index_22__028+i-110] 
		] )
	-np.std( [ _028_PREVAIL_LBM_028["j1_1"]["PU"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["PU"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["PU"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["PU"][index_22__028+i-110] 
		] ) )
	# 
	_028_mps_raw.append( np.mean( [ _028_PREVAIL_LBM_028["j1_1"]["PU"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["PU"][index_12__028+i-110] ,
		_028_PREVAIL_LBM_028["j2_1"]["PU"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["PU"][index_22__028+i-110] 
		] )
	+np.std( [ _028_PREVAIL_LBM_028["j1_1"]["PU"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["PU"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["PU"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["PU"][index_22__028+i-110] 
		] ) )
	# 
	_028_Temperaturemean_signal_raw.append( np.mean( [ _028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i-110] 
		] ) )
	# 
	_028_Temperaturestd_signal_raw.append( np.std( [ _028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i-110] 
		] ) )
	# 
	_028_Temperaturemms_raw.append( np.mean( [ _028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i-110] 
		] )
	-np.std( [ _028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i-110] 
		] ) )
	# 
	_028_Temperaturemps_raw.append( np.mean( [ _028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i-110] 
		] )
	+np.std( [ _028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(_028_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12__028-110:index_12__028+length_list__028], _028_mean_signal_raw,  linestyle='-',color='black', label="P$_4$", alpha=0.5)
plt.fill_between(_028_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12__028-110:index_12__028+length_list__028], _028_mms_raw, _028_mps_raw, color='black', alpha=0.2)
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
plt.plot(_028_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12__028-110:index_12__028+length_list__028], _028_Temperaturemean_signal_raw, linestyle='-', color='black', label="P$_4$", alpha=0.5)
plt.fill_between(_028_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12__028-110:index_12__028+length_list__028], _028_Temperaturemms_raw, _028_Temperaturemps_raw, color='black', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P4.jpg', bbox_inches='tight')









plt.close()


























# Load the Excel file
file_path = "029_PREVAIL_LBM_029.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
_029_PREVAIL_LBM_029 = {}
for sheet_name, data in all_sheets.items():
    _029_PREVAIL_LBM_029[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11__029 = 0
for ii in range(len(_029_PREVAIL_LBM_029["j1_1"]["temps (s)"] ) ):
	if _029_PREVAIL_LBM_029["j1_1"]["temps (s)"][ii] == 0:
		index_11__029 = ii 
		break
index_12__029 = 0
for ii in range(len(_029_PREVAIL_LBM_029["j1_2"]["temps (s)"] ) ):
	if _029_PREVAIL_LBM_029["j1_2"]["temps (s)"][ii] == 0:
		index_12__029 = ii 
		break
index_21__029 = 0
for ii in range(len(_029_PREVAIL_LBM_029["j2_1"]["temps (s)"] ) ):
	if _029_PREVAIL_LBM_029["j2_1"]["temps (s)"][ii] == 0:
		index_21__029 = ii 
		break
index_22__029 = 0
for ii in range(len(_029_PREVAIL_LBM_029["j2_2"]["temps (s)"] ) ):
	if _029_PREVAIL_LBM_029["j2_2"]["temps (s)"][ii] == 0:
		index_22__029 = ii 
		break		

length_list__029 = min(len(_029_PREVAIL_LBM_029["j1_1"]["temps (s)"][index_11__029:]) , len(_029_PREVAIL_LBM_029["j1_2"]["temps (s)"][index_12__029:]) , len(_029_PREVAIL_LBM_029["j2_1"]["temps (s)"][index_21__029:]) , len(_029_PREVAIL_LBM_029["j2_2"]["temps (s)"][index_22__029:] ) )

_029_Temperaturemean_signal=[]
_029_Temperaturemps = []
_029_Temperaturemms = []
_029_Temperaturestd_signal = []
_029_mean_signal=[]
_029_mps = []
_029_mms = []
_029_std_signal = []

for i in range(length_list__029):
	_029_mean_signal.append( np.mean( [ _029_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22__029+i] 
		] ) )
	# 
	_029_std_signal.append( np.std( [ _029_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22__029+i] 
		] ) )
	# 
	_029_mms.append( np.mean( [ _029_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22__029+i] 
		] )
	- 1.96* np.std( [ _029_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22__029+i] 
		] ) )
	# 
	_029_mps.append( np.mean( [ _029_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22__029+i] 
		] )
	+1.96*np.std( [ _029_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22__029+i] 
		] ) )
	# 
	_029_Temperaturemean_signal.append( np.mean( [ _029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i] 
		] ) )
	# 
	_029_Temperaturestd_signal.append( np.std( [ _029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i] 
		] ) )
	# 
	_029_Temperaturemms.append( np.mean( [ _029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i] 
		] )
	-1.96*np.std( [ _029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i] 
		] ) )
	# 
	_029_Temperaturemps.append( np.mean( [ _029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i] 
		] )
	+1.96*np.std( [ _029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i] 
		] ) )
	
_029_Temperaturemean_signal_raw=[]
_029_Temperaturemps_raw = []
_029_Temperaturemms_raw = []
_029_Temperaturestd_signal_raw = []
_029_mean_signal_raw=[]
_029_mps_raw = []
_029_mms_raw = []
_029_std_signal_raw = []

for i in range(length_list__029+110):
	_029_mean_signal_raw.append( np.mean( [ _029_PREVAIL_LBM_029["j1_1"]["PU"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["PU"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["PU"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["PU"][index_22__029+i-110] 
		] ) )
	# 
	_029_std_signal_raw.append( np.std( [ _029_PREVAIL_LBM_029["j1_1"]["PU"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["PU"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["PU"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["PU"][index_22__029+i-110] 
		] ) )
	# 
	_029_mms_raw.append( np.mean( [ _029_PREVAIL_LBM_029["j1_1"]["PU"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["PU"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["PU"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["PU"][index_22__029+i-110] 
		] )
	-np.std( [ _029_PREVAIL_LBM_029["j1_1"]["PU"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["PU"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["PU"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["PU"][index_22__029+i-110] 
		] ) )
	# 
	_029_mps_raw.append( np.mean( [ _029_PREVAIL_LBM_029["j1_1"]["PU"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["PU"][index_12__029+i-110] ,
		_029_PREVAIL_LBM_029["j2_1"]["PU"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["PU"][index_22__029+i-110] 
		] )
	+np.std( [ _029_PREVAIL_LBM_029["j1_1"]["PU"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["PU"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["PU"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["PU"][index_22__029+i-110] 
		] ) )
	# 
	_029_Temperaturemean_signal_raw.append( np.mean( [ _029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i-110] 
		] ) )
	# 
	_029_Temperaturestd_signal_raw.append( np.std( [ _029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i-110] 
		] ) )
	# 
	_029_Temperaturemms_raw.append( np.mean( [ _029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i-110] 
		] )
	-np.std( [ _029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i-110] 
		] ) )
	# 
	_029_Temperaturemps_raw.append( np.mean( [ _029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i-110] 
		] )
	+np.std( [ _029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(_029_PREVAIL_LBM_029["j1_2"]["temps (s)"][index_12__029-110:index_12__029+length_list__029], _029_mean_signal_raw, color='darkgreen', label="P$_5$", alpha=0.5)
plt.fill_between(_029_PREVAIL_LBM_029["j1_2"]["temps (s)"][index_12__029-110:index_12__029+length_list__029], _029_mms_raw, _029_mps_raw, color='darkgreen', alpha=0.2)
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
plt.plot(_029_PREVAIL_LBM_029["j1_2"]["temps (s)"][index_12__029-110:index_12__029+length_list__029], _029_Temperaturemean_signal_raw, color='darkgreen', label="P$_5$", alpha=0.5)
plt.fill_between(_029_PREVAIL_LBM_029["j1_2"]["temps (s)"][index_12__029-110:index_12__029+length_list__029], _029_Temperaturemms_raw, _029_Temperaturemps_raw, color='darkgreen', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P5.jpg', bbox_inches='tight')








plt.close()






































# Load the Excel file
file_path = "030_PREVAIL_LBM_030.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
_030_PREVAIL_LBM_030 = {}
for sheet_name, data in all_sheets.items():
    _030_PREVAIL_LBM_030[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11__030 = 0
for ii in range(len(_030_PREVAIL_LBM_030["j1_1"]["temps (s)"] ) ):
	if _030_PREVAIL_LBM_030["j1_1"]["temps (s)"][ii] == 0:
		index_11__030 = ii 
		break
index_12__030 = 0
for ii in range(len(_030_PREVAIL_LBM_030["j1_2"]["temps (s)"] ) ):
	if _030_PREVAIL_LBM_030["j1_2"]["temps (s)"][ii] == 0:
		index_12__030 = ii 
		break
index_21__030 = 0
for ii in range(len(_030_PREVAIL_LBM_030["j2_1"]["temps (s)"] ) ):
	if _030_PREVAIL_LBM_030["j2_1"]["temps (s)"][ii] == 0:
		index_21__030 = ii 
		break
index_22__030 = 0
for ii in range(len(_030_PREVAIL_LBM_030["j2_2"]["temps (s)"] ) ):
	if _030_PREVAIL_LBM_030["j2_2"]["temps (s)"][ii] == 0:
		index_22__030 = ii 
		break		

length_list__030 = min(len(_030_PREVAIL_LBM_030["j1_1"]["temps (s)"][index_11__030:]) , len(_030_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12__030:]) , len(_030_PREVAIL_LBM_030["j2_1"]["temps (s)"][index_21__030:]) , len(_030_PREVAIL_LBM_030["j2_2"]["temps (s)"][index_22__030:] ) )

_030_Temperaturemean_signal=[]
_030_Temperaturemps = []
_030_Temperaturemms = []
_030_Temperaturestd_signal = []
_030_mean_signal=[]
_030_mps = []
_030_mms = []
_030_std_signal = []

for i in range(length_list__030):
	_030_mean_signal.append( np.mean( [ _030_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22__030+i] 
		] ) )
	# 
	_030_std_signal.append( np.std( [ _030_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22__030+i] 
		] ) )
	# 
	_030_mms.append( np.mean( [ _030_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22__030+i] 
		] )
	- 1.96* np.std( [ _030_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22__030+i] 
		] ) )
	# 
	_030_mps.append( np.mean( [ _030_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22__030+i] 
		] )
	+1.96*np.std( [ _030_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22__030+i] 
		] ) )
	# 
	_030_Temperaturemean_signal.append( np.mean( [ _030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i] 
		] ) )
	# 
	_030_Temperaturestd_signal.append( np.std( [ _030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i] 
		] ) )
	# 
	_030_Temperaturemms.append( np.mean( [ _030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i] 
		] )
	-1.96*np.std( [ _030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i] 
		] ) )
	# 
	_030_Temperaturemps.append( np.mean( [ _030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i] 
		] )
	+1.96*np.std( [ _030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i] 
		] ) )
	
_030_Temperaturemean_signal_raw=[]
_030_Temperaturemps_raw = []
_030_Temperaturemms_raw = []
_030_Temperaturestd_signal_raw = []
_030_mean_signal_raw=[]
_030_mps_raw = []
_030_mms_raw = []
_030_std_signal_raw = []

for i in range(length_list__030+110):
	_030_mean_signal_raw.append( np.mean( [ _030_PREVAIL_LBM_030["j1_1"]["PU"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["PU"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["PU"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["PU"][index_22__030+i-110] 
		] ) )
	# 
	_030_std_signal_raw.append( np.std( [ _030_PREVAIL_LBM_030["j1_1"]["PU"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["PU"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["PU"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["PU"][index_22__030+i-110] 
		] ) )
	# 
	_030_mms_raw.append( np.mean( [ _030_PREVAIL_LBM_030["j1_1"]["PU"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["PU"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["PU"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["PU"][index_22__030+i-110] 
		] )
	-np.std( [ _030_PREVAIL_LBM_030["j1_1"]["PU"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["PU"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["PU"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["PU"][index_22__030+i-110] 
		] ) )
	# 
	_030_mps_raw.append( np.mean( [ _030_PREVAIL_LBM_030["j1_1"]["PU"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["PU"][index_12__030+i-110] ,
		_030_PREVAIL_LBM_030["j2_1"]["PU"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["PU"][index_22__030+i-110] 
		] )
	+np.std( [ _030_PREVAIL_LBM_030["j1_1"]["PU"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["PU"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["PU"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["PU"][index_22__030+i-110] 
		] ) )
	# 
	_030_Temperaturemean_signal_raw.append( np.mean( [ _030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i-110] 
		] ) )
	# 
	_030_Temperaturestd_signal_raw.append( np.std( [ _030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i-110] 
		] ) )
	# 
	_030_Temperaturemms_raw.append( np.mean( [ _030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i-110] 
		] )
	-np.std( [ _030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i-110] 
		] ) )
	# 
	_030_Temperaturemps_raw.append( np.mean( [ _030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i-110] 
		] )
	+np.std( [ _030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(_030_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12__030-110:index_12__030+length_list__030], _030_mean_signal_raw, color='gold', label="P$_6$", alpha=0.5)
plt.fill_between(_030_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12__030-110:index_12__030+length_list__030], _030_mms_raw, _030_mps_raw, color='gold', alpha=0.2)
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
plt.plot(_030_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12__030-110:index_12__030+length_list__030], _030_Temperaturemean_signal_raw, color='gold', label="P$_6$", alpha=0.5)
plt.fill_between(_030_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12__030-110:index_12__030+length_list__030], _030_Temperaturemms_raw, _030_Temperaturemps_raw, color='gold', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([29, 32])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P6.jpg', bbox_inches='tight')








plt.close()























# Load the Excel file
file_path = "031_PREVAIL_LBM_031.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
_031_PREVAIL_LBM_031 = {}
for sheet_name, data in all_sheets.items():
    _031_PREVAIL_LBM_031[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11__031 = 0
for ii in range(len(_031_PREVAIL_LBM_031["j1_1"]["temps (s)"] ) ):
	if _031_PREVAIL_LBM_031["j1_1"]["temps (s)"][ii] == 0:
		index_11__031 = ii 
		break
index_12__031 = 0
for ii in range(len(_031_PREVAIL_LBM_031["j1_2"]["temps (s)"] ) ):
	if _031_PREVAIL_LBM_031["j1_2"]["temps (s)"][ii] == 0:
		index_12__031 = ii 
		break
index_21__031 = 0
for ii in range(len(_031_PREVAIL_LBM_031["j2_1"]["temps (s)"] ) ):
	if _031_PREVAIL_LBM_031["j2_1"]["temps (s)"][ii] == 0:
		index_21__031 = ii 
		break
index_22__031 = 0
for ii in range(len(_031_PREVAIL_LBM_031["j2_2"]["temps (s)"] ) ):
	if _031_PREVAIL_LBM_031["j2_2"]["temps (s)"][ii] == 0:
		index_22__031 = ii 
		break		

length_list__031 = min(len(_031_PREVAIL_LBM_031["j1_1"]["temps (s)"][index_11__031:]) , len(_031_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12__031:]) , len(_031_PREVAIL_LBM_031["j2_1"]["temps (s)"][index_21__031:]) , len(_031_PREVAIL_LBM_031["j2_2"]["temps (s)"][index_22__031:] ) )

_031_Temperaturemean_signal=[]
_031_Temperaturemps = []
_031_Temperaturemms = []
_031_Temperaturestd_signal = []
_031_mean_signal=[]
_031_mps = []
_031_mms = []
_031_std_signal = []

for i in range(length_list__031):
	_031_mean_signal.append( np.mean( [ _031_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22__031+i] 
		] ) )
	# 
	_031_std_signal.append( np.std( [ _031_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22__031+i] 
		] ) )
	# 
	_031_mms.append( np.mean( [ _031_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22__031+i] 
		] )
	- 1.96* np.std( [ _031_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22__031+i] 
		] ) )
	# 
	_031_mps.append( np.mean( [ _031_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22__031+i] 
		] )
	+1.96*np.std( [ _031_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22__031+i] 
		] ) )
	# 
	_031_Temperaturemean_signal.append( np.mean( [ _031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i] 
		] ) )
	# 
	_031_Temperaturestd_signal.append( np.std( [ _031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i] 
		] ) )
	# 
	_031_Temperaturemms.append( np.mean( [ _031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i] 
		] )
	-1.96*np.std( [ _031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i] 
		] ) )
	# 
	_031_Temperaturemps.append( np.mean( [ _031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i] 
		] )
	+1.96*np.std( [ _031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i] 
		] ) )
	
_031_Temperaturemean_signal_raw=[]
_031_Temperaturemps_raw = []
_031_Temperaturemms_raw = []
_031_Temperaturestd_signal_raw = []
_031_mean_signal_raw=[]
_031_mps_raw = []
_031_mms_raw = []
_031_std_signal_raw = []

for i in range(length_list__031+110):
	_031_mean_signal_raw.append( np.mean( [ _031_PREVAIL_LBM_031["j1_1"]["PU"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["PU"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["PU"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["PU"][index_22__031+i-110] 
		] ) )
	# 
	_031_std_signal_raw.append( np.std( [ _031_PREVAIL_LBM_031["j1_1"]["PU"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["PU"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["PU"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["PU"][index_22__031+i-110] 
		] ) )
	# 
	_031_mms_raw.append( np.mean( [ _031_PREVAIL_LBM_031["j1_1"]["PU"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["PU"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["PU"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["PU"][index_22__031+i-110] 
		] )
	-np.std( [ _031_PREVAIL_LBM_031["j1_1"]["PU"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["PU"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["PU"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["PU"][index_22__031+i-110] 
		] ) )
	# 
	_031_mps_raw.append( np.mean( [ _031_PREVAIL_LBM_031["j1_1"]["PU"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["PU"][index_12__031+i-110] ,
		_031_PREVAIL_LBM_031["j2_1"]["PU"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["PU"][index_22__031+i-110] 
		] )
	+np.std( [ _031_PREVAIL_LBM_031["j1_1"]["PU"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["PU"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["PU"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["PU"][index_22__031+i-110] 
		] ) )
	# 
	_031_Temperaturemean_signal_raw.append( np.mean( [ _031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i-110] 
		] ) )
	# 
	_031_Temperaturestd_signal_raw.append( np.std( [ _031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i-110] 
		] ) )
	# 
	_031_Temperaturemms_raw.append( np.mean( [ _031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i-110] 
		] )
	-np.std( [ _031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i-110] 
		] ) )
	# 
	_031_Temperaturemps_raw.append( np.mean( [ _031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i-110] 
		] )
	+np.std( [ _031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(_031_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12__031-110:index_12__031+length_list__031], _031_mean_signal_raw, linestyle='-', color='blue', label="P$_7$", alpha=0.5)
plt.fill_between(_031_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12__031-110:index_12__031+length_list__031], _031_mms_raw, _031_mps_raw, color='blue', alpha=0.2)
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
plt.plot(_031_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12__031-110:index_12__031+length_list__031], _031_Temperaturemean_signal_raw, linestyle='-', color='blue', label="P$_7$", alpha=0.5)
plt.fill_between(_031_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12__031-110:index_12__031+length_list__031], _031_Temperaturemms_raw, _031_Temperaturemps_raw, color='blue', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([29.8, 31.8])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P7.jpg', bbox_inches='tight')








plt.close()
















# Load the Excel file
file_path = "032_PREVAIL_LBM_032.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
_032_PREVAIL_LBM_032 = {}
for sheet_name, data in all_sheets.items():
    _032_PREVAIL_LBM_032[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11__032 = 0
for ii in range(len(_032_PREVAIL_LBM_032["j1_1"]["temps (s)"] ) ):
	if _032_PREVAIL_LBM_032["j1_1"]["temps (s)"][ii] == 0:
		index_11__032 = ii 
		break
index_12__032 = 0
for ii in range(len(_032_PREVAIL_LBM_032["j1_2"]["temps (s)"] ) ):
	if _032_PREVAIL_LBM_032["j1_2"]["temps (s)"][ii] == 0:
		index_12__032 = ii 
		break
index_21__032 = 0
for ii in range(len(_032_PREVAIL_LBM_032["j2_1"]["temps (s)"] ) ):
	if _032_PREVAIL_LBM_032["j2_1"]["temps (s)"][ii] == 0:
		index_21__032 = ii 
		break
index_22__032 = 0
for ii in range(len(_032_PREVAIL_LBM_032["j2_2"]["temps (s)"] ) ):
	if _032_PREVAIL_LBM_032["j2_2"]["temps (s)"][ii] == 0:
		index_22__032 = ii 
		break		

length_list__032 = min(len(_032_PREVAIL_LBM_032["j1_1"]["temps (s)"][index_11__032:]) , len(_032_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12__032:]) , len(_032_PREVAIL_LBM_032["j2_1"]["temps (s)"][index_21__032:]) , len(_032_PREVAIL_LBM_032["j2_2"]["temps (s)"][index_22__032:] ) )

_032_Temperaturemean_signal=[]
_032_Temperaturemps = []
_032_Temperaturemms = []
_032_Temperaturestd_signal = []
_032_mean_signal=[]
_032_mps = []
_032_mms = []
_032_std_signal = []

for i in range(length_list__032):
	_032_mean_signal.append( np.mean( [ _032_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21__032+i],
		_032_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22__032+i] 
		] ) )
	# 
	_032_std_signal.append( np.std( [ _032_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21__032+i],
		_032_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22__032+i] 
		] ) )
	# 
	_032_mms.append( np.mean( [ _032_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21__032+i],
		_032_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22__032+i] 
		] )
	- 1.96* np.std( [ _032_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21__032+i],
		_032_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22__032+i] 
		] ) )
	# 
	_032_mps.append( np.mean( [ _032_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21__032+i],
		_032_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22__032+i] 
		] )
	+1.96*np.std( [ _032_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21__032+i],
		_032_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22__032+i] 
		] ) )
	# 
	_032_Temperaturemean_signal.append( np.mean( [ _032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__032+i],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__032+i] 
		] ) )
	# 
	_032_Temperaturestd_signal.append( np.std( [ _032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__032+i],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__032+i] 
		] ) )
	# 
	_032_Temperaturemms.append( np.mean( [ _032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__032+i],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__032+i] 
		] )
	-1.96*np.std( [ _032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__032+i],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__032+i] 
		] ) )
	# 
	_032_Temperaturemps.append( np.mean( [ _032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__032+i],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__032+i] 
		] )
	+1.96*np.std( [ _032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__032+i],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__032+i] 
		] ) )
	
_032_Temperaturemean_signal_raw=[]
_032_Temperaturemps_raw = []
_032_Temperaturemms_raw = []
_032_Temperaturestd_signal_raw = []
_032_mean_signal_raw=[]
_032_mps_raw = []
_032_mms_raw = []
_032_std_signal_raw = []

for i in range(length_list__032+110):
	_032_mean_signal_raw.append( np.mean( [ _032_PREVAIL_LBM_032["j1_1"]["PU"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["PU"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["PU"][index_21__032+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["PU"][index_22__032+i-110] 
		] ) )
	# 
	_032_std_signal_raw.append( np.std( [ _032_PREVAIL_LBM_032["j1_1"]["PU"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["PU"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["PU"][index_21__032+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["PU"][index_22__032+i-110] 
		] ) )
	# 
	_032_mms_raw.append( np.mean( [ _032_PREVAIL_LBM_032["j1_1"]["PU"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["PU"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["PU"][index_21__032+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["PU"][index_22__032+i-110] 
		] )
	-np.std( [ _032_PREVAIL_LBM_032["j1_1"]["PU"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["PU"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["PU"][index_21__032+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["PU"][index_22__032+i-110] 
		] ) )
	# 
	_032_mps_raw.append( np.mean( [ _032_PREVAIL_LBM_032["j1_1"]["PU"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["PU"][index_12__032+i-110] ,
		_032_PREVAIL_LBM_032["j2_1"]["PU"][index_21__032+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["PU"][index_22__032+i-110] 
		] )
	+np.std( [ _032_PREVAIL_LBM_032["j1_1"]["PU"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["PU"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["PU"][index_21__032+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["PU"][index_22__032+i-110] 
		] ) )
	# 
	_032_Temperaturemean_signal_raw.append( np.mean( [ _032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__032+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__032+i-110] 
		] ) )
	# 
	_032_Temperaturestd_signal_raw.append( np.std( [ _032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__032+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__032+i-110] 
		] ) )
	# 
	_032_Temperaturemms_raw.append( np.mean( [ _032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__032+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__032+i-110] 
		] )
	-np.std( [ _032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__032+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__032+i-110] 
		] ) )
	# 
	_032_Temperaturemps_raw.append( np.mean( [ _032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__032+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__032+i-110] 
		] )
	+np.std( [ _032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__032+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__032+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(_032_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12__032-110:index_12__032+length_list__032], _032_mean_signal_raw, linestyle='-', color='red', label="P$_8$", alpha=0.5)
plt.fill_between(_032_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12__032-110:index_12__032+length_list__032], _032_mms_raw, _032_mps_raw, color='red', alpha=0.2)
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
plt.plot(_032_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12__032-110:index_12__032+length_list__032], _032_Temperaturemean_signal_raw, linestyle='-', color='red', label="P$_8$", alpha=0.5)
plt.fill_between(_032_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12__032-110:index_12__032+length_list__032], _032_Temperaturemms_raw, _032_Temperaturemps_raw, color='red', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P8.jpg', bbox_inches='tight')









plt.close()


































# Load the Excel file
file_path = "033_PREVAIL_LBM_033.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
_033_PREVAIL_LBM_033 = {}
for sheet_name, data in all_sheets.items():
    _033_PREVAIL_LBM_033[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11__033 = 0
for ii in range(len(_033_PREVAIL_LBM_033["j1_1"]["temps (s)"] ) ):
	if _033_PREVAIL_LBM_033["j1_1"]["temps (s)"][ii] == 0:
		index_11__033 = ii 
		break
index_12__033 = 0
for ii in range(len(_033_PREVAIL_LBM_033["j1_2"]["temps (s)"] ) ):
	if _033_PREVAIL_LBM_033["j1_2"]["temps (s)"][ii] == 0:
		index_12__033 = ii 
		break
index_21__033 = 0
for ii in range(len(_033_PREVAIL_LBM_033["j2_1"]["temps (s)"] ) ):
	if _033_PREVAIL_LBM_033["j2_1"]["temps (s)"][ii] == 0:
		index_21__033 = ii 
		break
index_22__033 = 0
for ii in range(len(_033_PREVAIL_LBM_033["j2_2"]["temps (s)"] ) ):
	if _033_PREVAIL_LBM_033["j2_2"]["temps (s)"][ii] == 0:
		index_22__033 = ii 
		break		

length_list__033 = min(len(_033_PREVAIL_LBM_033["j1_1"]["temps (s)"][index_11__033:]) , len(_033_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12__033:]), len(_033_PREVAIL_LBM_033["j2_1"]["temps (s)"][index_21__033:]) , len(_033_PREVAIL_LBM_033["j2_2"]["temps (s)"][index_22__033:] ) )

_033_Temperaturemean_signal=[]
_033_Temperaturemps = []
_033_Temperaturemms = []
_033_Temperaturestd_signal = []
_033_mean_signal=[]
_033_mps = []
_033_mms = []
_033_std_signal = []

for i in range(length_list__033):
	_033_mean_signal.append( np.mean( [ _033_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22__033+i] 
		] ) )
	# 
	_033_std_signal.append( np.std( [ _033_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22__033+i] 
		] ) )
	# 
	_033_mms.append( np.mean( [ _033_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22__033+i] 
		] )
	- 1.96* np.std( [ _033_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22__033+i] 
		] ) )
	# 
	_033_mps.append( np.mean( [ _033_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22__033+i] 
		] )
	+1.96*np.std( [ _033_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22__033+i] 
		] ) )
	# 
	_033_Temperaturemean_signal.append( np.mean( [ _033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i] 
		] ) )
	# 
	_033_Temperaturestd_signal.append( np.std( [ _033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i] 
		] ) )
	# 
	_033_Temperaturemms.append( np.mean( [ _033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i] 
		] )
	-1.96*np.std( [ _033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i] 
		] ) )
	# 
	_033_Temperaturemps.append( np.mean( [ _033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i] 
		] )
	+1.96*np.std( [ _033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i] 
		] ) )
	
_033_Temperaturemean_signal_raw=[]
_033_Temperaturemps_raw = []
_033_Temperaturemms_raw = []
_033_Temperaturestd_signal_raw = []
_033_mean_signal_raw=[]
_033_mps_raw = []
_033_mms_raw = []
_033_std_signal_raw = []

for i in range(length_list__033+110):
	_033_mean_signal_raw.append( np.mean( [ _033_PREVAIL_LBM_033["j1_1"]["PU"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["PU"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["PU"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["PU"][index_22__033+i-110] 
		] ) )
	# 
	_033_std_signal_raw.append( np.std( [ _033_PREVAIL_LBM_033["j1_1"]["PU"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["PU"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["PU"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["PU"][index_22__033+i-110] 
		] ) )
	# 
	_033_mms_raw.append( np.mean( [ _033_PREVAIL_LBM_033["j1_1"]["PU"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["PU"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["PU"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["PU"][index_22__033+i-110] 
		] )
	-np.std( [ _033_PREVAIL_LBM_033["j1_1"]["PU"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["PU"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["PU"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["PU"][index_22__033+i-110] 
		] ) )
	# 
	_033_mps_raw.append( np.mean( [ _033_PREVAIL_LBM_033["j1_1"]["PU"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["PU"][index_12__033+i-110] ,
		_033_PREVAIL_LBM_033["j2_1"]["PU"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["PU"][index_22__033+i-110] 
		] )
	+np.std( [ _033_PREVAIL_LBM_033["j1_1"]["PU"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["PU"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["PU"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["PU"][index_22__033+i-110] 
		] ) )
	# 
	_033_Temperaturemean_signal_raw.append( np.mean( [ _033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i-110] 
		] ) )
	# 
	_033_Temperaturestd_signal_raw.append( np.std( [ _033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i-110] 
		] ) )
	# 
	_033_Temperaturemms_raw.append( np.mean( [ _033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i-110] 
		] )
	-np.std( [ _033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i-110] 
		] ) )
	# 
	_033_Temperaturemps_raw.append( np.mean( [ _033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i-110] 
		] )
	+np.std( [ _033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(_033_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12__033-110:index_12__033+length_list__033], _033_mean_signal_raw, linestyle='-', color='fuchsia', label="P$_9$", alpha=0.5)
plt.fill_between(_033_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12__033-110:index_12__033+length_list__033], _033_mms_raw, _033_mps_raw, color='fuchsia', alpha=0.2)
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
plt.plot(_033_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12__033-110:index_12__033+length_list__033], _033_Temperaturemean_signal_raw, linestyle='-', color='fuchsia', label="P$_9$", alpha=0.5)
plt.fill_between(_033_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12__033-110:index_12__033+length_list__033], _033_Temperaturemms_raw, _033_Temperaturemps_raw, color='fuchsia', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P9.jpg', bbox_inches='tight')




plt.close()



































# Load the Excel file
file_path = "034_PREVAIL_LBM_034.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
_034_PREVAIL_LBM_034 = {}
for sheet_name, data in all_sheets.items():
    _034_PREVAIL_LBM_034[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11__034 = 0
for ii in range(len(_034_PREVAIL_LBM_034["j1_1"]["temps (s)"] ) ):
	if _034_PREVAIL_LBM_034["j1_1"]["temps (s)"][ii] == 0:
		index_11__034 = ii 
		break
index_12__034 = 0
for ii in range(len(_034_PREVAIL_LBM_034["j1_2"]["temps (s)"] ) ):
	if _034_PREVAIL_LBM_034["j1_2"]["temps (s)"][ii] == 0:
		index_12__034 = ii 
		break
index_21__034 = 0
for ii in range(len(_034_PREVAIL_LBM_034["j2_1"]["temps (s)"] ) ):
	if _034_PREVAIL_LBM_034["j2_1"]["temps (s)"][ii] == 0:
		index_21__034 = ii 
		break
index_22__034 = 0
for ii in range(len(_034_PREVAIL_LBM_034["j2_2"]["temps (s)"] ) ):
	if _034_PREVAIL_LBM_034["j2_2"]["temps (s)"][ii] == 0:
		index_22__034 = ii 
		break		

length_list__034 = min(len(_034_PREVAIL_LBM_034["j1_1"]["temps (s)"][index_11__034:]) , len(_034_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12__034:]) , len(_034_PREVAIL_LBM_034["j2_1"]["temps (s)"][index_21__034:]) , len(_034_PREVAIL_LBM_034["j2_2"]["temps (s)"][index_22__034:] ) )

_034_Temperaturemean_signal=[]
_034_Temperaturemps = []
_034_Temperaturemms = []
_034_Temperaturestd_signal = []
_034_mean_signal=[]
_034_mps = []
_034_mms = []
_034_std_signal = []

for i in range(length_list__034):
	_034_mean_signal.append( np.mean( [ _034_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22__034+i] 
		] ) )
	# 
	_034_std_signal.append( np.std( [ _034_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22__034+i] 
		] ) )
	# 
	_034_mms.append( np.mean( [ _034_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22__034+i] 
		] )
	- 1.96* np.std( [ _034_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22__034+i] 
		] ) )
	# 
	_034_mps.append( np.mean( [ _034_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22__034+i] 
		] )
	+1.96*np.std( [ _034_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22__034+i] 
		] ) )
	# 
	_034_Temperaturemean_signal.append( np.mean( [ _034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i] 
		] ) )
	# 
	_034_Temperaturestd_signal.append( np.std( [ _034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i] 
		] ) )
	# 
	_034_Temperaturemms.append( np.mean( [ _034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i] 
		] )
	-1.96*np.std( [ _034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i] 
		] ) )
	# 
	_034_Temperaturemps.append( np.mean( [ _034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i] 
		] )
	+1.96*np.std( [ _034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i] 
		] ) )
	
_034_Temperaturemean_signal_raw=[]
_034_Temperaturemps_raw = []
_034_Temperaturemms_raw = []
_034_Temperaturestd_signal_raw = []
_034_mean_signal_raw=[]
_034_mps_raw = []
_034_mms_raw = []
_034_std_signal_raw = []

for i in range(length_list__034+110):
	_034_mean_signal_raw.append( np.mean( [ _034_PREVAIL_LBM_034["j1_1"]["PU"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["PU"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["PU"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["PU"][index_22__034+i-110] 
		] ) )
	# 
	_034_std_signal_raw.append( np.std( [ _034_PREVAIL_LBM_034["j1_1"]["PU"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["PU"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["PU"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["PU"][index_22__034+i-110] 
		] ) )
	# 
	_034_mms_raw.append( np.mean( [ _034_PREVAIL_LBM_034["j1_1"]["PU"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["PU"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["PU"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["PU"][index_22__034+i-110] 
		] )
	-np.std( [ _034_PREVAIL_LBM_034["j1_1"]["PU"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["PU"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["PU"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["PU"][index_22__034+i-110] 
		] ) )
	# 
	_034_mps_raw.append( np.mean( [ _034_PREVAIL_LBM_034["j1_1"]["PU"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["PU"][index_12__034+i-110] ,
		_034_PREVAIL_LBM_034["j2_1"]["PU"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["PU"][index_22__034+i-110] 
		] )
	+np.std( [ _034_PREVAIL_LBM_034["j1_1"]["PU"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["PU"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["PU"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["PU"][index_22__034+i-110] 
		] ) )
	# 
	_034_Temperaturemean_signal_raw.append( np.mean( [ _034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i-110] 
		] ) )
	# 
	_034_Temperaturestd_signal_raw.append( np.std( [ _034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i-110] 
		] ) )
	# 
	_034_Temperaturemms_raw.append( np.mean( [ _034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i-110] 
		] )
	-np.std( [ _034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i-110] 
		] ) )
	# 
	_034_Temperaturemps_raw.append( np.mean( [ _034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i-110] 
		] )
	+np.std( [ _034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(_034_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12__034-110:index_12__034+length_list__034], _034_mean_signal_raw, linestyle='-', color='darkgreen', label="P$_{10}$", alpha=0.5)
plt.fill_between(_034_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12__034-110:index_12__034+length_list__034], _034_mms_raw, _034_mps_raw, color='darkgreen', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
plt.ylim([0, 150])
plt.xlabel("Time [s]")
plt.ylabel("LDF [AU]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_raw_P10.jpg', bbox_inches='tight')

plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(_034_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12__034-110:index_12__034+length_list__034], _034_Temperaturemean_signal_raw, linestyle='-', color='darkgreen', label="P$_{10}$", alpha=0.5)
plt.fill_between(_034_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12__034-110:index_12__034+length_list__034], _034_Temperaturemms_raw, _034_Temperaturemps_raw, color='darkgreen', alpha=0.2)
plt.legend()
plt.xlim([-60, 780])
# plt.ylim([27, 33])
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.savefig('./Figures/superp_mean_std_temp_raw_P10.jpg', bbox_inches='tight')






























# Load the Excel file
file_path = "035_PREVAIL_LBM_035.xlsx"  # Change this to the path of your Excel file

# Read all sheets into a dictionary of DataFrames
all_sheets = pd.read_excel(file_path, sheet_name=None, engine="openpyxl")


# Optional: Extract specific columns for all sheets
column_names = ["temps (s)", "PU", "PU_pc", "Temperature", "baseline"]

# Extract data from each sheet into a structured dictionary
_035_PREVAIL_LBM_035 = {}
for sheet_name, data in all_sheets.items():
    _035_PREVAIL_LBM_035[sheet_name] = {col: data[col] for col in column_names if col in data.columns}

index_11__035 = 0
for ii in range(len(_035_PREVAIL_LBM_035["j1_1"]["temps (s)"] ) ):
	if _035_PREVAIL_LBM_035["j1_1"]["temps (s)"][ii] == 0:
		index_11__035 = ii 
		break
index_12__035 = 0
for ii in range(len(_035_PREVAIL_LBM_035["j1_2"]["temps (s)"] ) ):
	if _035_PREVAIL_LBM_035["j1_2"]["temps (s)"][ii] == 0:
		index_12__035 = ii 
		break
index_21__035 = 0
for ii in range(len(_035_PREVAIL_LBM_035["j2_1"]["temps (s)"] ) ):
	if _035_PREVAIL_LBM_035["j2_1"]["temps (s)"][ii] == 0:
		index_21__035 = ii 
		break
index_22__035 = 0
for ii in range(len(_035_PREVAIL_LBM_035["j2_2"]["temps (s)"] ) ):
	if _035_PREVAIL_LBM_035["j2_2"]["temps (s)"][ii] == 0:
		index_22__035 = ii 
		break		

length_list__035 = min(len(_035_PREVAIL_LBM_035["j1_1"]["temps (s)"][index_11__035:]) , len(_035_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12__035:]) , len(_035_PREVAIL_LBM_035["j2_1"]["temps (s)"][index_21__035:]) , len(_035_PREVAIL_LBM_035["j2_2"]["temps (s)"][index_22__035:] ) )

_035_Temperaturemean_signal=[]
_035_Temperaturemps = []
_035_Temperaturemms = []
_035_Temperaturestd_signal = []
_035_mean_signal=[]
_035_mps = []
_035_mms = []
_035_std_signal = []

for i in range(length_list__035):
	_035_mean_signal.append( np.mean( [ _035_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22__035+i] 
		] ) )
	# 
	_035_std_signal.append( np.std( [ _035_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22__035+i] 
		] ) )
	# 
	_035_mms.append( np.mean( [ _035_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22__035+i] 
		] )
	- 1.96* np.std( [ _035_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22__035+i] 
		] ) )
	# 
	_035_mps.append( np.mean( [ _035_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22__035+i] 
		] )
	+1.96*np.std( [ _035_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22__035+i] 
		] ) )
	# 
	_035_Temperaturemean_signal.append( np.mean( [ _035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i] 
		] ) )
	# 
	_035_Temperaturestd_signal.append( np.std( [ _035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i] 
		] ) )
	# 
	_035_Temperaturemms.append( np.mean( [ _035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i] 
		] )
	-1.96*np.std( [ _035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i] 
		] ) )
	# 
	_035_Temperaturemps.append( np.mean( [ _035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i] 
		] )
	+1.96*np.std( [ _035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i] 
		] ) )
	
_035_Temperaturemean_signal_raw=[]
_035_Temperaturemps_raw = []
_035_Temperaturemms_raw = []
_035_Temperaturestd_signal_raw = []
_035_mean_signal_raw=[]
_035_mps_raw = []
_035_mms_raw = []
_035_std_signal_raw = []

for i in range(length_list__035+110):
	_035_mean_signal_raw.append( np.mean( [ _035_PREVAIL_LBM_035["j1_1"]["PU"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["PU"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["PU"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["PU"][index_22__035+i-110] 
		] ) )
	# 
	_035_std_signal_raw.append( np.std( [ _035_PREVAIL_LBM_035["j1_1"]["PU"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["PU"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["PU"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["PU"][index_22__035+i-110] 
		] ) )
	# 
	_035_mms_raw.append( np.mean( [ _035_PREVAIL_LBM_035["j1_1"]["PU"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["PU"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["PU"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["PU"][index_22__035+i-110] 
		] )
	-np.std( [ _035_PREVAIL_LBM_035["j1_1"]["PU"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["PU"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["PU"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["PU"][index_22__035+i-110] 
		] ) )
	# 
	_035_mps_raw.append( np.mean( [ _035_PREVAIL_LBM_035["j1_1"]["PU"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["PU"][index_12__035+i-110] ,
		_035_PREVAIL_LBM_035["j2_1"]["PU"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["PU"][index_22__035+i-110] 
		] )
	+np.std( [ _035_PREVAIL_LBM_035["j1_1"]["PU"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["PU"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["PU"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["PU"][index_22__035+i-110] 
		] ) )
	# 
	_035_Temperaturemean_signal_raw.append( np.mean( [ _035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i-110] 
		] ) )
	# 
	_035_Temperaturestd_signal_raw.append( np.std( [ _035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i-110] 
		] ) )
	# 
	_035_Temperaturemms_raw.append( np.mean( [ _035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i-110] 
		] )
	-np.std( [ _035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i-110] 
		] ) )
	# 
	_035_Temperaturemps_raw.append( np.mean( [ _035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i-110] 
		] )
	+np.std( [ _035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i-110] 
		] ) )



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
# # 
plt.plot(_035_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12__035-110:index_12__035+length_list__035], _035_mean_signal_raw, linestyle='-', color='gold', label="P$_{11}$", alpha=0.5)
plt.fill_between(_035_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12__035-110:index_12__035+length_list__035], _035_mms_raw, _035_mps_raw, color='gold', alpha=0.2)
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
plt.plot(_035_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12__035-110:index_12__035+length_list__035], _035_Temperaturemean_signal_raw, linestyle='-', color='gold', label="P$_{11}$", alpha=0.5)
plt.fill_between(_035_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12__035-110:index_12__035+length_list__035], _035_Temperaturemms_raw, _035_Temperaturemps_raw, color='gold', alpha=0.2)
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
# for i in range(min(length_list__009+110,length_list+110,length_list__026+110,length_list__028+110,length_list__029+110,length_list__030+110,length_list__031+110,length_list__032+110,length_list__033+110,length_list__034+110,length_list__035+110)):
for i in range(min(length_list__009+110,length_list+110,length_list__026+110,length_list__028+110,length_list__029+110,length_list__030+110,length_list__031+110,length_list__032+110,length_list__033+110,length_list__034+110,length_list__035+110)):
	all_means.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["PU"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["PU"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["PU"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["PU"][index_22__009+i-110],
		# 
		_008_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
		_008_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
		_008_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
		_008_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110],
		# 
		_026_PREVAIL_LBM_026["j1_1"]["PU"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["PU"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["PU"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["PU"][index_22__026+i-110],
		# 
		_028_PREVAIL_LBM_028["j1_1"]["PU"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["PU"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["PU"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["PU"][index_22__028+i-110],  
		# 
		_029_PREVAIL_LBM_029["j1_1"]["PU"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["PU"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["PU"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["PU"][index_22__029+i-110], 
		# 
		_030_PREVAIL_LBM_030["j1_1"]["PU"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["PU"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["PU"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["PU"][index_22__030+i-110] ,
		# 
		_031_PREVAIL_LBM_031["j1_1"]["PU"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["PU"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["PU"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["PU"][index_22__031+i-110] , 
		# 
		_032_PREVAIL_LBM_032["j1_1"]["PU"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["PU"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["PU"][index_21__030+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["PU"][index_22__030+i-110] ,
		# #
		_033_PREVAIL_LBM_033["j1_1"]["PU"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["PU"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["PU"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["PU"][index_22__033+i-110] ,
		# # 
		_034_PREVAIL_LBM_034["j1_1"]["PU"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["PU"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["PU"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["PU"][index_22__034+i-110] ,
		# # 
		_035_PREVAIL_LBM_035["j1_1"]["PU"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["PU"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["PU"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["PU"][index_22__035+i-110] 
		# 
		] ) )
	# 
	all_std.append( np.std( [ _009_PREVAIL_LBM_009["j1_1"]["PU"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["PU"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["PU"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["PU"][index_22__009+i-110],
		# 
		_008_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
		_008_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
		_008_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
		_008_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110],
		# 
		_026_PREVAIL_LBM_026["j1_1"]["PU"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["PU"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["PU"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["PU"][index_22__026+i-110] ,
		# 
		_028_PREVAIL_LBM_028["j1_1"]["PU"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["PU"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["PU"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["PU"][index_22__028+i-110],  
		# 
		_029_PREVAIL_LBM_029["j1_1"]["PU"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["PU"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["PU"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["PU"][index_22__029+i-110], 
		# 
		_030_PREVAIL_LBM_030["j1_1"]["PU"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["PU"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["PU"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["PU"][index_22__030+i-110] ,
		# 
		_031_PREVAIL_LBM_031["j1_1"]["PU"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["PU"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["PU"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["PU"][index_22__031+i-110] ,
		# 
		# 
		_032_PREVAIL_LBM_032["j1_1"]["PU"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["PU"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["PU"][index_21__030+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["PU"][index_22__030+i-110] ,
		# #
		_033_PREVAIL_LBM_033["j1_1"]["PU"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["PU"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["PU"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["PU"][index_22__033+i-110] ,
		# # 
		_034_PREVAIL_LBM_034["j1_1"]["PU"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["PU"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["PU"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["PU"][index_22__034+i-110] ,
		# # 
		_035_PREVAIL_LBM_035["j1_1"]["PU"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["PU"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["PU"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["PU"][index_22__035+i-110] 
		] ) )
	# 
	all_means_temp.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i-110],
		# 
		_008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
		_008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
		_008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
		_008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110],
		# 
		_026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i-110],
		# 
		_028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i-110],  
		# 
		_029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i-110], 
		# 
		_030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i-110] ,
		# 
		_031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i-110] ,
		# 
		# 
		_032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__030+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__030+i-110] ,
		# #
		_033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i-110] ,
		# # 
		_034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i-110] ,
		# # 
		_035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i-110] 
		] ) )
	# 
	all_std_temp.append( np.std( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i-110],
		# 
		_008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
		_008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
		_008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
		_008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110],
		# 
		_026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i-110] ,
		# 
		_028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i-110],  
		# 
		_029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i-110], 
		# 
		_030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i-110] ,
		# 
		_031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i-110] ,
		# 
		_032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__030+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__030+i-110] ,
		# #
		_033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i-110] ,
		# # 
		_034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i-110] ,
		# # 
		_035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i-110] 
		# 
		] ) )

all_mms = [all_means[i]-all_std[i] for i in range(len(all_means))]
all_mps = [all_means[i]+all_std[i] for i in range(len(all_means))]

all_means_pc=[]
all_std_pc=[]
# Percentage
# for i in range(min(length_list__009,length_list,length_list__026,length_list__028,length_list__029,length_list__030,length_list__031,length_list__032,length_list__033,length_list__034,length_list__035)):
for i in range(min(length_list__009,length_list,length_list__026,length_list__028,length_list__029,length_list__030,length_list__031,length_list__032,length_list__033,length_list__035,length_list__034)):
	all_means_pc.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11__009+i],
		_009_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12__009+i],
		_009_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21__009+i],
		_009_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22__009+i],
		# 
		_008_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
		_008_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
		_008_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
		_008_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i],
		# 
		_026_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22__026+i]  ,
		# 
		_028_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22__028+i],  
		# 
		_029_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22__029+i], 
		# 
		_030_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22__030+i] ,
		# 
		_031_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22__031+i] ,
		# 
		_032_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21__030+i],
		_032_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22__030+i] ,
		# #
		_033_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22__033+i] ,
		# # 
		_034_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22__034+i] ,
		# # 
		_035_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22__035+i] 
		#                               
		] ) )
	all_std_pc.append( np.std( [ _009_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11__009+i],
		_009_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12__009+i],
		_009_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21__009+i],
		_009_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22__009+i],
		# 
		_008_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
		_008_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
		_008_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
		_008_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i],
		# 
		_026_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22__026+i] ,
		# 
		_028_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22__028+i],  
		# 
		_029_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22__029+i], 
		# 
		_030_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22__030+i] ,
		# 
		_031_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22__031+i] ,
		# 
		_032_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21__030+i],
		_032_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22__030+i] ,
		# #
		_033_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22__033+i] ,
		# # 
		_034_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22__034+i] ,
		# # 
		_035_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22__035+i] 
		#                                    
		 ] ) )

# 68% confidence
all_mms_pc = [all_means_pc[i] -  all_std_pc[i] for i in range(len(all_means_pc))]
all_mps_pc = [all_means_pc[i] +  all_std_pc[i] for i in range(len(all_means_pc))]

# # 95% Confidence
# all_mms_pc = [all_means_pc[i] - 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]
# all_mps_pc = [all_means_pc[i] + 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]




initial_temp = np.mean( all_means_temp[:110] )
std_initial_temp = np.mean( all_std_temp[:110] )
print('Initial Temperature °C:',initial_temp, '$\\pm$', std_initial_temp)


initial_baseline = np.mean( all_means[:110] )
std_initial = np.mean( all_std[:110] )
print('baseline LDF AU:',initial_baseline, '$\\pm$', std_initial)

first_occlusions = np.mean( [ all_means[60+110:160+110],all_means[400+110:500+110] ] )
std_fo = np.mean( [ np.mean(all_std[60+110:160+110]),np.mean(all_std[400+110:500+110]) ] )
print('First occlusion LDF AU:',first_occlusions, '$\\pm$', std_fo)

last_occlusions = np.mean( [ all_means[800+110:900+110],all_means[1200+110:1300+110] ] )
std_lo = np.mean( [ np.mean(all_std[800+110:900+110]),np.mean(all_std[1200+110:1300+110]) ] )
print('Last occlusion LDF AU:',last_occlusions, '$\\pm$', std_lo)


first_occlusions_pc = np.mean( [ all_means_pc[60:160],all_means_pc[400:500] ] )
std_fo_pc = np.mean( [ np.mean(all_std_pc[60:160]),np.mean(all_std_pc[400:500]) ] )
print('First occlusion LDF %:',first_occlusions_pc, '$\\pm$', std_fo_pc)

last_occlusions_pc = np.mean( [ all_means_pc[800:900],all_means_pc[1200:1300] ] )
std_lo_pc = np.mean( [ np.mean(all_std_pc[800:900]),np.mean(all_std_pc[1200:1300]) ] )
print('Last occlusion LDF %:',last_occlusions_pc, '$\\pm$', std_lo_pc)

# idex0 = 60+110
# idex = 160+110
# t = np.linspace(idex0,idex,idex-idex0)

# len_t = len(all_means)
# tcheck = np.linspace(0,_028_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12__028+length_list__028],len_t)

# plt.figure(figsize=(10, 6))
# plt.plot(tcheck, all_means, alpha=0.5)
# plt.plot(t, all_means[idex0:idex], label="Original Signal", alpha=0.5)
# plt.savefig('./Figures/check_time_range_occlusion_1.png')
# plt.close()

# idex0 = 400+110
# idex = 500+110
# t = np.linspace(idex0,idex,idex-idex0)

# plt.figure(figsize=(10, 6))
# plt.plot(tcheck, all_means, alpha=0.5)
# plt.plot(t, all_means[idex0:idex], label="Original Signal", alpha=0.5)
# plt.savefig('./Figures/check_time_range_occlusion_2.png')
# plt.close()



# idex0 = 800+110
# idex = 900+110
# t = np.linspace(idex0,idex,idex-idex0)

# plt.figure(figsize=(10, 6))
# plt.plot(tcheck, all_means, alpha=0.5)
# plt.plot(t, all_means[idex0:idex], label="Original Signal", alpha=0.5)
# plt.savefig('./Figures/check_time_range_occlusion_3.png')
# plt.close()



# idex0 = 1200+110
# idex = 1300+110
# t = np.linspace(idex0,idex,idex-idex0)

# plt.figure(figsize=(10, 6))
# plt.plot(tcheck, all_means, alpha=0.5)
# plt.plot(t, all_means[idex0:idex], label="Original Signal", alpha=0.5)
# plt.savefig('./Figures/check_time_range_occlusion_4.png')
# plt.close()



# len_t = len(all_means_pc)
# tcheck = np.linspace(0,_028_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12__028+length_list__028],len_t)


# idex0 = 60
# idex = 160
# t = np.linspace(idex0,idex,idex-idex0)

# plt.figure(figsize=(10, 6))
# plt.plot(tcheck, all_means_pc, alpha=0.5)
# plt.plot(t, all_means_pc[idex0:idex], label="Original Signal", alpha=0.5)
# plt.savefig('./Figures/check_time_range_occlusion_1pc.png')
# plt.close()


# idex0 = 400
# idex = 500
# t = np.linspace(idex0,idex,idex-idex0)

# plt.figure(figsize=(10, 6))
# plt.plot(tcheck, all_means_pc, alpha=0.5)
# plt.plot(t, all_means_pc[idex0:idex], label="Original Signal", alpha=0.5)
# plt.savefig('./Figures/check_time_range_occlusion_2pc.png')
# plt.close()



# idex0 = 800
# idex = 900
# t = np.linspace(idex0,idex,idex-idex0)

# plt.figure(figsize=(10, 6))
# plt.plot(tcheck, all_means_pc, alpha=0.5)
# plt.plot(t, all_means_pc[idex0:idex], label="Original Signal", alpha=0.5)
# plt.savefig('./Figures/check_time_range_occlusion_3pc.png')
# plt.close()



# idex0 = 1200
# idex = 1300
# t = np.linspace(idex0,idex,idex-idex0)

# plt.figure(figsize=(10, 6))
# plt.plot(tcheck, all_means_pc, alpha=0.5)
# plt.plot(t, all_means_pc[idex0:idex], label="Original Signal", alpha=0.5)
# plt.savefig('./Figures/check_time_range_occlusion_4pc.png')
# plt.close()



plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(_008_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-110:index_22+length_list], _008_mean_signal_raw, color='black', label="P$_1$", alpha=1, linewidth=3)
# 
plt.plot(_009_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22__009-110:index_22__009+length_list__009], _009_mean_signal_raw, color='blue', label="P$_2$", alpha=1, linewidth=3)
# 
plt.plot(_026_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12__026-110:index_12__026+length_list__026], _026_mean_signal_raw, color='red', label="P$_3$", alpha=1, linewidth=3)
plt.plot(_028_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12__029-110:index_12__029+length_list__029], _029_mean_signal_raw, color='darkgreen', label="P$_5$", alpha=1, linewidth=3)
plt.plot(_030_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12__030-110:index_12__030+length_list__030], _030_mean_signal_raw, color='gold', label="P$_6$", alpha=1, linewidth=3)
plt.plot(_033_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12__033-110:index_12__033+length_list__033], _033_mean_signal_raw, color='fuchsia', label="P$_9$", alpha=1, linewidth=3)
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
plt.plot(_028_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12__028-110:index_12__028+length_list__028], _028_mean_signal_raw, linestyle='-', color='black', label="P$_4$", alpha=1, linewidth=3)
plt.plot(_031_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12__031-110:index_12__031+length_list__031], _031_mean_signal_raw, linestyle='-', color='blue', label="P$_7$", alpha=1, linewidth=3)
plt.plot(_032_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12__032-110:index_12__032+length_list__032], _032_mean_signal_raw, linestyle='-', color='red', label="P$_8$", alpha=1, linewidth=3)
plt.plot(_034_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12__034-110:index_12__034+length_list__034], _034_mean_signal_raw, linestyle='-', color='darkgreen', label="P$_{10}$", alpha=1, linewidth=3)
plt.plot(_035_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12__035-110:index_12__035+length_list__035], _035_mean_signal_raw, linestyle='-', color='gold', label="P$_{11}$", alpha=1, linewidth=3)
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
plt.plot(_008_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22:index_22+length_list], _008_mean_signal, color='black', label="P$_1$", alpha=1, linewidth=3)
# 
plt.plot(_009_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22__009:index_22__009+length_list__009], _009_mean_signal, color='blue', label="P$_2$", alpha=1, linewidth=3)
# 
plt.plot(_026_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12__026:index_12__026+length_list__026], _026_mean_signal, color='red', label="P$_3$", alpha=1, linewidth=3)
plt.plot(_028_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12__029:index_12__029+length_list__029], _029_mean_signal, color='darkgreen', label="P$_5$", alpha=1, linewidth=3)
plt.plot(_030_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12__030:index_12__030+length_list__030], _030_mean_signal, color='gold', label="P$_6$", alpha=1, linewidth=3)
plt.plot(_033_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12__033:index_12__033+length_list__033], _033_mean_signal, color='fuchsia', label="P$_9$", alpha=1, linewidth=3)
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
plt.plot(_028_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12__028:index_12__028+length_list__028], _028_mean_signal, linestyle='-', color='black', label="P$_4$", alpha=1, linewidth=3)
plt.plot(_031_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12__031:index_12__031+length_list__031], _031_mean_signal, linestyle='-', color='blue', label="P$_7$", alpha=1, linewidth=3)
plt.plot(_032_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12__032:index_12__032+length_list__032], _032_mean_signal, linestyle='-', color='red', label="P$_8$", alpha=1, linewidth=3)
plt.plot(_034_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12__034:index_12__034+length_list__034], _034_mean_signal, linestyle='-', color='darkgreen', label="P$_{10}$", alpha=1, linewidth=3)
plt.plot(_035_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12__035:index_12__035+length_list__035], _035_mean_signal, linestyle='-', color='gold', label="P$_{11}$", alpha=1, linewidth=3)
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
plt.plot(_008_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22-110:index_22+length_list], _008_Temperaturemean_signal_raw, color='black', label="P$_1$", alpha=1, linewidth=3)
# 
plt.plot(_009_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22__009-110:index_22__009+length_list__009], _009_Temperaturemean_signal_raw, color='blue', label="P$_2$", alpha=1, linewidth=3)
# 
plt.plot(_026_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12__026-110:index_12__026+length_list__026], _026_Temperaturemean_signal_raw, color='red', label="P$_3$", alpha=1, linewidth=3)
plt.plot(_028_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12__029-110:index_12__029+length_list__029], _029_Temperaturemean_signal_raw, color='darkgreen', label="P$_5$", alpha=1, linewidth=3)
plt.plot(_030_PREVAIL_LBM_030["j1_2"]["temps (s)"][index_12__030-110:index_12__030+length_list__030], _030_Temperaturemean_signal_raw, color='gold', label="P$_6$", alpha=1, linewidth=3)
plt.plot(_033_PREVAIL_LBM_033["j1_2"]["temps (s)"][index_12__033-110:index_12__033+length_list__033], _033_Temperaturemean_signal_raw, color='fuchsia', label="P$_9$", alpha=0.5, linewidth=3)
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
plt.plot(_028_PREVAIL_LBM_028["j1_2"]["temps (s)"][index_12__028-110:index_12__028+length_list__028], _028_Temperaturemean_signal_raw, linestyle='-', color='black', label="P$_4$", alpha=1, linewidth=3)
plt.plot(_031_PREVAIL_LBM_031["j1_2"]["temps (s)"][index_12__031-110:index_12__031+length_list__031], _031_Temperaturemean_signal_raw, linestyle='-', color='blue', label="P$_7$", alpha=1, linewidth=3)
plt.plot(_032_PREVAIL_LBM_032["j1_2"]["temps (s)"][index_12__032-110:index_12__032+length_list__032], _032_Temperaturemean_signal_raw, linestyle='-', color='red', label="P$_8$", alpha=1, linewidth=3)
plt.plot(_034_PREVAIL_LBM_034["j1_2"]["temps (s)"][index_12__034-110:index_12__034+length_list__034], _034_Temperaturemean_signal_raw, linestyle='-', color='darkgreen', label="P$_{10}$", alpha=1, linewidth=3)
plt.plot(_035_PREVAIL_LBM_035["j1_2"]["temps (s)"][index_12__035-110:index_12__035+length_list__035], _035_Temperaturemean_signal_raw, linestyle='-', color='gold', label="P$_{11}$", alpha=1, linewidth=3)
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
column_names = ["dis_030cement_all",	"LDF_v_all",	"LDF_q_all",	"LDF_baseline_v",	"LDF_baseline_q",	"load_all",	"time_all"]


# Extract data from each sheet into a structured dictionary
model = {}
for sheet_name, data in all_sheets.items():
    model[sheet_name] = {col: data[col] for col in column_names if col in data.columns}
# 
plt.rcParams.update({'font.size': 25})
plt.figure(figsize=(20, 12))
plt.plot(_009_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22__009:index_22__009+length_list__009], _009_mean_signal, color='r', label="P$_2$", alpha=0.6)
plt.fill_between(_009_PREVAIL_LBM_009["j2_2"]["temps (s)"][index_22__009:index_22__009+length_list__009], _009_mms, _009_mps, color='r', alpha=0.2)
# 
plt.plot(_026_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12__026:index_12__026+length_list__026], _026_mean_signal , linestyle='-', color='g', label="P$_3$", alpha=0.6)
plt.fill_between(_026_PREVAIL_LBM_026["j1_2"]["temps (s)"][index_12__026:index_12__026+length_list__026], _026_mms, _026_mps , linestyle='-', color='g', alpha=0.2)
# 
# plt.plot(_029_PREVAIL_LBM_029["j1_2"]["temps (s)"][index_12__029:index_12__029+length_list__029], _029_mean_signal , linestyle='-', color='turquoise', label="P$_5$", alpha=0.6)
# plt.fill_between(_029_PREVAIL_LBM_029["j1_2"]["temps (s)"][index_12__029:index_12__029+length_list__029], _029_mms, _029_mps , linestyle='-', color='turquoise', alpha=0.2)
# plt.plot(_008_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22:index_22+length_list], all_means_pc, color='k', label="Exp", alpha=0.5, linewidth=3)
# plt.fill_between(_008_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22:index_22+length_list], all_mms_pc, all_mps_pc, color='k', alpha=0.2)
# plt.plot(_008_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22:index_22+length_list], all_means_pc, color='k', label="Exp", alpha=0.5, linewidth=3)
# plt.fill_between(_008_PREVAIL_LBM_008["j2_2"]["temps (s)"][index_22:index_22+length_list], _026_mms[:length_list], _009_mps[:length_list], color='k', alpha=0.2)

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
# for i in range(min(length_list__009+110,length_list+110,length_list__026+110,length_list__028+110,length_list__029+110,length_list__030+110,length_list__031+110,length_list__032+110,length_list__033+110,length_list__034+110,length_list__035+110)):
for i in range(min(length_list__009+110,length_list+110,length_list__026+110,length_list__029+110,length_list__030+110,length_list__033+110)):
	all_means.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["PU"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["PU"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["PU"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["PU"][index_22__009+i-110],
		# 
		_008_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
		_008_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
		_008_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
		_008_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110],
		# 
		_026_PREVAIL_LBM_026["j1_1"]["PU"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["PU"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["PU"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["PU"][index_22__026+i-110],
		# 
		# 
		_029_PREVAIL_LBM_029["j1_1"]["PU"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["PU"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["PU"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["PU"][index_22__029+i-110], 
		# 
		_030_PREVAIL_LBM_030["j1_1"]["PU"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["PU"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["PU"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["PU"][index_22__030+i-110] ,
		# 
		# 
		# #
		_033_PREVAIL_LBM_033["j1_1"]["PU"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["PU"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["PU"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["PU"][index_22__033+i-110] 
		# 
		] ) )
	# 
	all_std.append( np.std( [ _009_PREVAIL_LBM_009["j1_1"]["PU"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["PU"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["PU"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["PU"][index_22__009+i-110],
		# 
		_008_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110],
		_008_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110],
		_008_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110],
		_008_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110],
		# 
		_026_PREVAIL_LBM_026["j1_1"]["PU"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["PU"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["PU"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["PU"][index_22__026+i-110] ,
		# 
		_029_PREVAIL_LBM_029["j1_1"]["PU"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["PU"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["PU"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["PU"][index_22__029+i-110], 
		# 
		_030_PREVAIL_LBM_030["j1_1"]["PU"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["PU"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["PU"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["PU"][index_22__030+i-110] ,
		# 
		_031_PREVAIL_LBM_031["j1_1"]["PU"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["PU"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["PU"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["PU"][index_22__031+i-110] ,
		# 
		# #
		_033_PREVAIL_LBM_033["j1_1"]["PU"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["PU"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["PU"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["PU"][index_22__033+i-110] 
		] ) )
	# 
	all_means_temp.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i-110],
		# 
		_008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
		_008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
		_008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
		_008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110],
		# 
		_026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i-110],
		# 
		_029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i-110], 
		# 
		_030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i-110] ,
		# 
		# 
		# #
		_033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i-110] 
		] ) )
	# 
	all_std_temp.append( np.std( [ _009_PREVAIL_LBM_009["j1_1"]["Temperature"][index_11__009+i-110],
		_009_PREVAIL_LBM_009["j1_2"]["Temperature"][index_12__009+i-110],
		_009_PREVAIL_LBM_009["j2_1"]["Temperature"][index_21__009+i-110],
		_009_PREVAIL_LBM_009["j2_2"]["Temperature"][index_22__009+i-110],
		# 
		_008_PREVAIL_LBM_008["j1_1"]["Temperature"][index_11+i-110],
		_008_PREVAIL_LBM_008["j1_2"]["Temperature"][index_12+i-110],
		_008_PREVAIL_LBM_008["j2_1"]["Temperature"][index_21+i-110],
		_008_PREVAIL_LBM_008["j2_2"]["Temperature"][index_22+i-110],
		# 
		_026_PREVAIL_LBM_026["j1_1"]["Temperature"][index_11__026+i-110],
		_026_PREVAIL_LBM_026["j1_2"]["Temperature"][index_12__026+i-110],
		_026_PREVAIL_LBM_026["j2_1"]["Temperature"][index_21__026+i-110],
		_026_PREVAIL_LBM_026["j2_2"]["Temperature"][index_22__026+i-110] ,
		# 
		# 
		_029_PREVAIL_LBM_029["j1_1"]["Temperature"][index_11__029+i-110],
		_029_PREVAIL_LBM_029["j1_2"]["Temperature"][index_12__029+i-110],
		_029_PREVAIL_LBM_029["j2_1"]["Temperature"][index_21__029+i-110],
		_029_PREVAIL_LBM_029["j2_2"]["Temperature"][index_22__029+i-110], 
		# 
		_030_PREVAIL_LBM_030["j1_1"]["Temperature"][index_11__030+i-110],
		_030_PREVAIL_LBM_030["j1_2"]["Temperature"][index_12__030+i-110],
		_030_PREVAIL_LBM_030["j2_1"]["Temperature"][index_21__030+i-110],
		_030_PREVAIL_LBM_030["j2_2"]["Temperature"][index_22__030+i-110] ,
		# #
		_033_PREVAIL_LBM_033["j1_1"]["Temperature"][index_11__033+i-110],
		_033_PREVAIL_LBM_033["j1_2"]["Temperature"][index_12__033+i-110],
		_033_PREVAIL_LBM_033["j2_1"]["Temperature"][index_21__033+i-110],
		_033_PREVAIL_LBM_033["j2_2"]["Temperature"][index_22__033+i-110] 
		# 
		] ) )

all_mms = [all_means[i]-all_std[i] for i in range(len(all_means))]
all_mps = [all_means[i]+all_std[i] for i in range(len(all_means))]

all_means_pc=[]
all_std_pc=[]
# Percentage
# for i in range(min(length_list__009,length_list,length_list__026,length_list__028,length_list__029,length_list__030,length_list__031,length_list__032,length_list__033,length_list__034,length_list__035)):
for i in range(min(length_list__009,length_list,length_list__026,length_list__029,length_list__030,length_list__033)):
	all_means_pc.append( np.mean( [ _009_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11__009+i],
		_009_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12__009+i],
		_009_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21__009+i],
		_009_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22__009+i],
		# 
		_008_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
		_008_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
		_008_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
		_008_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i],
		# 
		_026_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22__026+i]  ,
		# 
		_029_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22__029+i], 
		# 
		_030_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22__030+i] ,
		# 
		# #
		_033_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22__033+i] 
		#                               
		] ) )
	all_std_pc.append( np.std( [ _009_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11__009+i],
		_009_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12__009+i],
		_009_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21__009+i],
		_009_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22__009+i],
		# 
		_008_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i],
		_008_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i],
		_008_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i],
		_008_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i],
		# 
		_026_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11__026+i],
		_026_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12__026+i],
		_026_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21__026+i],
		_026_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22__026+i] ,
		# 
		# 
		_029_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11__029+i],
		_029_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12__029+i],
		_029_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21__029+i],
		_029_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22__029+i], 
		# 
		_030_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11__030+i],
		_030_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12__030+i],
		_030_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21__030+i],
		_030_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22__030+i] ,
		# 
		# 
		# #
		_033_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11__033+i],
		_033_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12__033+i],
		_033_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21__033+i],
		_033_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22__033+i] 
		# # 
		#                                    
		 ] ) )

# 68% confidence
all_mms_pc = [all_means_pc[i] -  all_std_pc[i] for i in range(len(all_means_pc))]
all_mps_pc = [all_means_pc[i] +  all_std_pc[i] for i in range(len(all_means_pc))]

# # 95% Confidence
# all_mms_pc = [all_means_pc[i] - 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]
# all_mps_pc = [all_means_pc[i] + 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]



initial_temp = np.mean( all_means_temp[:110] )
std_initial_temp = np.mean( all_std_temp[:110] )
print('Initial Temperature °C: M',initial_temp, '$\\pm$', std_initial_temp)

initial_baseline = np.mean( all_means[:110] )
std_initial = np.mean( all_std[:110] )
print('baseline LDF AU: M',initial_baseline, '$\\pm$', std_initial)

first_occlusions = np.mean( [ all_means[60+110:160+110],all_means[400+110:500+110] ] )
std_fo = np.mean( [ np.mean(all_std[60+110:160+110]),np.mean(all_std[400+110:500+110]) ] )
print('First occlusion LDF AU: M',first_occlusions, '$\\pm$', std_fo)

last_occlusions = np.mean( [ all_means[800+110:900+110],all_means[1200+110:1300+110] ] )
std_lo = np.mean( [ np.mean(all_std[800+110:900+110]),np.mean(all_std[1200+110:1300+110]) ] )
print('Last occlusion LDF AU: M',last_occlusions, '$\\pm$', std_lo)


first_occlusions_pc = np.mean( [ all_means_pc[60:160],all_means_pc[400:500] ] )
std_fo_pc = np.mean( [ np.mean(all_std_pc[60:160]),np.mean(all_std_pc[400:500]) ] )
print('First occlusion LDF M %:',first_occlusions_pc, '$\\pm$', std_fo_pc)

last_occlusions_pc = np.mean( [ all_means_pc[800:900],all_means_pc[1200:1300] ] )
std_lo_pc = np.mean( [ np.mean(all_std_pc[800:900]),np.mean(all_std_pc[1200:1300]) ] )
print('Last occlusion LDF M %:',last_occlusions_pc, '$\\pm$', std_lo_pc)











# Raw

# A compléter

all_means=[]
all_means_temp = []
all_std=[]
all_std_temp =[]
# for i in range(min(length_list__009+110,length_list+110,length_list__026+110,length_list__028+110,length_list__029+110,length_list__030+110,length_list__031+110,length_list__032+110,length_list__033+110,length_list__034+110,length_list__035+110)):
for i in range(min(length_list__028+110,length_list__031+110,length_list__032+110,length_list__035+110,length_list__034+110)):
	all_means.append( np.mean( [ 
		_028_PREVAIL_LBM_028["j1_1"]["PU"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["PU"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["PU"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["PU"][index_22__028+i-110],  
		# 
		_031_PREVAIL_LBM_031["j1_1"]["PU"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["PU"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["PU"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["PU"][index_22__031+i-110] , 
		# 
		_032_PREVAIL_LBM_032["j1_1"]["PU"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["PU"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["PU"][index_21__030+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["PU"][index_22__030+i-110] ,
		# #
		# # 
		_034_PREVAIL_LBM_034["j1_1"]["PU"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["PU"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["PU"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["PU"][index_22__034+i-110] ,
		# # 
		_035_PREVAIL_LBM_035["j1_1"]["PU"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["PU"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["PU"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["PU"][index_22__035+i-110] 
		# 
		] ) )
	# 
	all_std.append( np.std( [ 
		# 
		_028_PREVAIL_LBM_028["j1_1"]["PU"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["PU"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["PU"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["PU"][index_22__028+i-110],  
		# 
		_031_PREVAIL_LBM_031["j1_1"]["PU"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["PU"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["PU"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["PU"][index_22__031+i-110] ,
		# 
		# 
		_032_PREVAIL_LBM_032["j1_1"]["PU"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["PU"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["PU"][index_21__030+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["PU"][index_22__030+i-110] ,
		# #
		# # 
		_034_PREVAIL_LBM_034["j1_1"]["PU"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["PU"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["PU"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["PU"][index_22__034+i-110] ,
		# # 
		_035_PREVAIL_LBM_035["j1_1"]["PU"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["PU"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["PU"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["PU"][index_22__035+i-110] 
		] ) )
	# 
	all_means_temp.append( np.mean( [ 
		# 
		_028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i-110],  
		# 
		# 
		_032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__030+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__030+i-110] ,
		# #
		# # 
		_034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i-110] ,
		# # 
		_035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i-110] 
		] ) )
	# 
	all_std_temp.append( np.std( [ 
		# 
		_028_PREVAIL_LBM_028["j1_1"]["Temperature"][index_11__028+i-110],
		_028_PREVAIL_LBM_028["j1_2"]["Temperature"][index_12__028+i-110],
		_028_PREVAIL_LBM_028["j2_1"]["Temperature"][index_21__028+i-110],
		_028_PREVAIL_LBM_028["j2_2"]["Temperature"][index_22__028+i-110],  
		# 
		# 
		_031_PREVAIL_LBM_031["j1_1"]["Temperature"][index_11__031+i-110],
		_031_PREVAIL_LBM_031["j1_2"]["Temperature"][index_12__031+i-110],
		_031_PREVAIL_LBM_031["j2_1"]["Temperature"][index_21__031+i-110],
		_031_PREVAIL_LBM_031["j2_2"]["Temperature"][index_22__031+i-110] ,
		# 
		_032_PREVAIL_LBM_032["j1_1"]["Temperature"][index_11__032+i-110],
		_032_PREVAIL_LBM_032["j1_2"]["Temperature"][index_12__032+i-110],
		_032_PREVAIL_LBM_032["j2_1"]["Temperature"][index_21__030+i-110],
		_032_PREVAIL_LBM_032["j2_2"]["Temperature"][index_22__030+i-110] ,
		# #
		# # 
		_034_PREVAIL_LBM_034["j1_1"]["Temperature"][index_11__034+i-110],
		_034_PREVAIL_LBM_034["j1_2"]["Temperature"][index_12__034+i-110],
		_034_PREVAIL_LBM_034["j2_1"]["Temperature"][index_21__034+i-110],
		_034_PREVAIL_LBM_034["j2_2"]["Temperature"][index_22__034+i-110] ,
		# # 
		_035_PREVAIL_LBM_035["j1_1"]["Temperature"][index_11__035+i-110],
		_035_PREVAIL_LBM_035["j1_2"]["Temperature"][index_12__035+i-110],
		_035_PREVAIL_LBM_035["j2_1"]["Temperature"][index_21__035+i-110],
		_035_PREVAIL_LBM_035["j2_2"]["Temperature"][index_22__035+i-110] 
		# 
		] ) )

all_mms = [all_means[i]-all_std[i] for i in range(len(all_means))]
all_mps = [all_means[i]+all_std[i] for i in range(len(all_means))]

all_means_pc=[]
all_std_pc=[]
# Percentage
# for i in range(min(length_list__009,length_list,length_list__026,length_list__028,length_list__029,length_list__030,length_list__031,length_list__032,length_list__033,length_list__034,length_list__035)):
for i in range(min(length_list__028,length_list__031,length_list__032,length_list__035,length_list__034)):
	all_means_pc.append( np.mean( [ 
		# 
		_028_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22__028+i],  
		# 
		_031_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22__031+i] ,
		# 
		_032_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21__030+i],
		_032_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22__030+i] ,
		# #
		# # 
		_034_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22__034+i] ,
		# # 
		_035_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22__035+i] 
		#                               
		] ) )
	all_std_pc.append( np.std( [ 
		# 
		# 
		# 
		_028_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11__028+i],
		_028_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12__028+i],
		_028_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21__028+i],
		_028_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22__028+i],  
		# 
		# 
		# 
		_031_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11__031+i],
		_031_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12__031+i],
		_031_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21__031+i],
		_031_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22__031+i] ,
		# 
		_032_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11__032+i],
		_032_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12__032+i],
		_032_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21__030+i],
		_032_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22__030+i] ,
		# #
		# # 
		_034_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11__034+i],
		_034_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12__034+i],
		_034_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21__034+i],
		_034_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22__034+i] ,
		# # 
		_035_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11__035+i],
		_035_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12__035+i],
		_035_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21__035+i],
		_035_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22__035+i] 
		#                                    
		 ] ) )

# 68% confidence
all_mms_pc = [all_means_pc[i] -  all_std_pc[i] for i in range(len(all_means_pc))]
all_mps_pc = [all_means_pc[i] +  all_std_pc[i] for i in range(len(all_means_pc))]

# # 95% Confidence
# all_mms_pc = [all_means_pc[i] - 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]
# all_mps_pc = [all_means_pc[i] + 1.96 * all_std_pc[i] for i in range(len(all_means_pc))]





initial_temp = np.mean( all_means_temp[:110] )
std_initial_temp = np.mean( all_std_temp[:110] )
print('Initial Temperature °C F:',initial_temp, '$\\pm$', std_initial_temp)

initial_baseline = np.mean( all_means[:110] )
std_initial = np.mean( all_std[:110] )
print('baseline LDF AU F:',initial_baseline, '$\\pm$', std_initial)

first_occlusions = np.mean( [ all_means[60+110:160+110],all_means[400+110:500+110] ] )
std_fo = np.mean( [ np.mean(all_std[60+110:160+110]),np.mean(all_std[400+110:500+110]) ] )
print('First occlusion LDF AU F:',first_occlusions, '$\\pm$', std_fo)

last_occlusions = np.mean( [ all_means[800+110:900+110],all_means[1200+110:1300+110] ] )
std_lo = np.mean( [ np.mean(all_std[800+110:900+110]),np.mean(all_std[1200+110:1300+110]) ] )
print('Last occlusion LDF AU F:',last_occlusions, '$\\pm$', std_lo)


first_occlusions_pc = np.mean( [ all_means_pc[60:160],all_means_pc[400:500] ] )
std_fo_pc = np.mean( [ np.mean(all_std_pc[60:160]),np.mean(all_std_pc[400:500]) ] )
print('First occlusion LDF F % :',first_occlusions_pc, '$\\pm$', std_fo_pc)

last_occlusions_pc = np.mean( [ all_means_pc[800:900],all_means_pc[1200:1300] ] )
std_lo_pc = np.mean( [ np.mean(all_std_pc[800:900]),np.mean(all_std_pc[1200:1300]) ] )
print('Last occlusion LDF F %:',last_occlusions_pc, '$\\pm$', std_lo_pc)



# Statistical test
# 
basal_value_male=[]
basal_value_female=[]
# 
first_ischaemia_male = []
first_ischaemia_female = []
second_ischaemia_male = []
second_ischaemia_female = []
third_ischaemia_male = []
third_ischaemia_female = []
fourth_ischaemia_male = []
fourth_ischaemia_female = []
first_ischaemia_male_pc = []
first_ischaemia_female_pc = []
second_ischaemia_male_pc = []
second_ischaemia_female_pc = []
third_ischaemia_male_pc = []
third_ischaemia_female_pc = []
fourth_ischaemia_male_pc = []
fourth_ischaemia_female_pc = []
# 
first_hyperaemia_male = []
first_hyperaemia_female = []
second_hyperaemia_male = []
second_hyperaemia_female = []
third_hyperaemia_male = []
third_hyperaemia_female = []
fourth_hyperaemia_male = []
fourth_hyperaemia_female = []
first_hyperaemia_male_pc = []
first_hyperaemia_female_pc = []
second_hyperaemia_male_pc = []
second_hyperaemia_female_pc = []
third_hyperaemia_male_pc = []
third_hyperaemia_female_pc = []
fourth_hyperaemia_male_pc = []
fourth_hyperaemia_female_pc = []
# 
# 

# Raw
# 1 male
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list+110):
	blank_a.append(_008_PREVAIL_LBM_008["j1_1"]["PU"][index_11+i-110])
	blank_b.append(_008_PREVAIL_LBM_008["j1_2"]["PU"][index_12+i-110])
	blank_c.append(_008_PREVAIL_LBM_008["j2_1"]["PU"][index_21+i-110])
	blank_d.append(_008_PREVAIL_LBM_008["j2_2"]["PU"][index_22+i-110])
first_ischaemia_male.append(np.mean(blank_a[60+110:160+110]))
first_ischaemia_male.append(np.mean(blank_b[60+110:160+110])) 
first_ischaemia_male.append(np.mean(blank_c[60+110:160+110]))
first_ischaemia_male.append(np.mean(blank_d[60+110:160+110]))
second_ischaemia_male.append(np.mean(blank_a[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_b[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_c[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_d[400+110:500+110]))
third_ischaemia_male.append(np.mean(blank_a[800+110:900+110]))
third_ischaemia_male.append(np.mean(blank_b[800+110:900+110]))
third_ischaemia_male.append(np.mean(blank_c[800+110:900+110])) 
third_ischaemia_male.append(np.mean(blank_d[800+110:900+110]))
fourth_ischaemia_male.append(np.mean(blank_a[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_b[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_c[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_d[1200+110:1300+110]))
first_hyperaemia_male.append(np.max(blank_a[160+110:400+110]))
first_hyperaemia_male.append(np.max(blank_b[160+110:400+110]))
first_hyperaemia_male.append(np.max(blank_c[160+110:400+110])) 
first_hyperaemia_male.append(np.max(blank_d[160+110:400+110]))
second_hyperaemia_male.append(np.max(blank_a[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_b[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_c[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_d[500+110:800+110]))
third_hyperaemia_male.append(np.max(blank_a[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_b[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_c[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_d[900+110:1200+110]))
fourth_hyperaemia_male.append(np.max(blank_a[1300+110:]))
fourth_hyperaemia_male.append(np.max(blank_b[1300+110:]))
fourth_hyperaemia_male.append( np.max(blank_c[1300+110:]))
fourth_hyperaemia_male.append(np.max(blank_d[1300+110:]))
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__009+110):
	blank_a.append(_009_PREVAIL_LBM_009["j1_1"]["PU"][index_11__009+i-110])
	blank_b.append(_009_PREVAIL_LBM_009["j1_2"]["PU"][index_12__009+i-110])
	blank_c.append(_009_PREVAIL_LBM_009["j2_1"]["PU"][index_21__009+i-110])
	blank_d.append(_009_PREVAIL_LBM_009["j2_2"]["PU"][index_22__009+i-110])
first_ischaemia_male.append(np.mean(blank_a[60+110:160+110]))
first_ischaemia_male.append(np.mean(blank_b[60+110:160+110])) 
first_ischaemia_male.append(np.mean(blank_c[60+110:160+110]))
first_ischaemia_male.append(np.mean(blank_d[60+110:160+110]))
second_ischaemia_male.append(np.mean(blank_a[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_b[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_c[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_d[400+110:500+110]))
third_ischaemia_male.append(np.mean(blank_a[800+110:900+110]))
third_ischaemia_male.append(np.mean(blank_b[800+110:900+110]))
third_ischaemia_male.append(np.mean(blank_c[800+110:900+110])) 
third_ischaemia_male.append(np.mean(blank_d[800+110:900+110]))
fourth_ischaemia_male.append(np.mean(blank_a[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_b[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_c[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_d[1200+110:1300+110]))
first_hyperaemia_male.append(np.max(blank_a[160+110:400+110]))
first_hyperaemia_male.append(np.max(blank_b[160+110:400+110]))
first_hyperaemia_male.append(np.max(blank_c[160+110:400+110])) 
first_hyperaemia_male.append(np.max(blank_d[160+110:400+110]))
second_hyperaemia_male.append(np.max(blank_a[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_b[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_c[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_d[500+110:800+110]))
third_hyperaemia_male.append(np.max(blank_a[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_b[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_c[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_d[900+110:1200+110]))
fourth_hyperaemia_male.append(np.max(blank_a[1300+110:]))
fourth_hyperaemia_male.append(np.max(blank_b[1300+110:]))
fourth_hyperaemia_male.append( np.max(blank_c[1300+110:]))
fourth_hyperaemia_male.append(np.max(blank_d[1300+110:]))
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__026+110):
	blank_a.append(_026_PREVAIL_LBM_026["j1_1"]["PU"][index_11__026+i-110])
	blank_b.append(_026_PREVAIL_LBM_026["j1_2"]["PU"][index_12__026+i-110])
	blank_c.append(_026_PREVAIL_LBM_026["j2_1"]["PU"][index_21__026+i-110])
	blank_d.append(_026_PREVAIL_LBM_026["j2_2"]["PU"][index_22__026+i-110])
first_ischaemia_male.append(np.mean(blank_a[60+110:160+110]))
first_ischaemia_male.append(np.mean(blank_b[60+110:160+110])) 
first_ischaemia_male.append(np.mean(blank_c[60+110:160+110]))
first_ischaemia_male.append(np.mean(blank_d[60+110:160+110]))
second_ischaemia_male.append(np.mean(blank_a[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_b[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_c[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_d[400+110:500+110]))
third_ischaemia_male.append(np.mean(blank_a[800+110:900+110]))
third_ischaemia_male.append(np.mean(blank_b[800+110:900+110]))
third_ischaemia_male.append(np.mean(blank_c[800+110:900+110])) 
third_ischaemia_male.append(np.mean(blank_d[800+110:900+110]))
fourth_ischaemia_male.append(np.mean(blank_a[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_b[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_c[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_d[1200+110:1300+110]))
first_hyperaemia_male.append(np.max(blank_a[160+110:400+110]))
first_hyperaemia_male.append(np.max(blank_b[160+110:400+110]))
first_hyperaemia_male.append(np.max(blank_c[160+110:400+110])) 
first_hyperaemia_male.append(np.max(blank_d[160+110:400+110]))
second_hyperaemia_male.append(np.max(blank_a[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_b[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_c[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_d[500+110:800+110]))
third_hyperaemia_male.append(np.max(blank_a[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_b[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_c[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_d[900+110:1200+110]))
fourth_hyperaemia_male.append(np.max(blank_a[1300+110:]))
fourth_hyperaemia_male.append(np.max(blank_b[1300+110:]))
fourth_hyperaemia_male.append( np.max(blank_c[1300+110:]))
fourth_hyperaemia_male.append(np.max(blank_d[1300+110:]))
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__029+110):
	blank_a.append(_029_PREVAIL_LBM_029["j1_1"]["PU"][index_11__029+i-110])
	blank_b.append(_029_PREVAIL_LBM_029["j1_2"]["PU"][index_12__029+i-110])
	blank_c.append(_029_PREVAIL_LBM_029["j2_1"]["PU"][index_21__029+i-110])
	blank_d.append(_029_PREVAIL_LBM_029["j2_2"]["PU"][index_22__029+i-110])
first_ischaemia_male.append(np.mean(blank_a[60+110:160+110]))
first_ischaemia_male.append(np.mean(blank_b[60+110:160+110])) 
first_ischaemia_male.append(np.mean(blank_c[60+110:160+110]))
first_ischaemia_male.append(np.mean(blank_d[60+110:160+110]))
second_ischaemia_male.append(np.mean(blank_a[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_b[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_c[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_d[400+110:500+110]))
third_ischaemia_male.append(np.mean(blank_a[800+110:900+110]))
third_ischaemia_male.append(np.mean(blank_b[800+110:900+110]))
third_ischaemia_male.append(np.mean(blank_c[800+110:900+110])) 
third_ischaemia_male.append(np.mean(blank_d[800+110:900+110]))
fourth_ischaemia_male.append(np.mean(blank_a[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_b[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_c[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_d[1200+110:1300+110]))
first_hyperaemia_male.append(np.max(blank_a[160+110:400+110]))
first_hyperaemia_male.append(np.max(blank_b[160+110:400+110]))
first_hyperaemia_male.append(np.max(blank_c[160+110:400+110])) 
first_hyperaemia_male.append(np.max(blank_d[160+110:400+110]))
second_hyperaemia_male.append(np.max(blank_a[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_b[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_c[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_d[500+110:800+110]))
third_hyperaemia_male.append(np.max(blank_a[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_b[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_c[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_d[900+110:1200+110]))
fourth_hyperaemia_male.append(np.max(blank_a[1300+110:]))
fourth_hyperaemia_male.append(np.max(blank_b[1300+110:]))
fourth_hyperaemia_male.append( np.max(blank_c[1300+110:]))
fourth_hyperaemia_male.append(np.max(blank_d[1300+110:]))
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__030+110):
	blank_a.append(_030_PREVAIL_LBM_030["j1_1"]["PU"][index_11__030+i-110])
	blank_b.append(_030_PREVAIL_LBM_030["j1_2"]["PU"][index_12__030+i-110])
	blank_c.append(_030_PREVAIL_LBM_030["j2_1"]["PU"][index_21__030+i-110])
	blank_d.append(_030_PREVAIL_LBM_030["j2_2"]["PU"][index_22__030+i-110])
first_ischaemia_male.append(np.mean(blank_a[60+110:160+110]))
first_ischaemia_male.append(np.mean(blank_b[60+110:160+110])) 
first_ischaemia_male.append(np.mean(blank_c[60+110:160+110]))
first_ischaemia_male.append(np.mean(blank_d[60+110:160+110]))
second_ischaemia_male.append(np.mean(blank_a[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_b[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_c[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_d[400+110:500+110]))
third_ischaemia_male.append(np.mean(blank_a[800+110:900+110]))
third_ischaemia_male.append(np.mean(blank_b[800+110:900+110]))
third_ischaemia_male.append(np.mean(blank_c[800+110:900+110])) 
third_ischaemia_male.append(np.mean(blank_d[800+110:900+110]))
fourth_ischaemia_male.append(np.mean(blank_a[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_b[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_c[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_d[1200+110:1300+110]))
first_hyperaemia_male.append(np.max(blank_a[160+110:400+110]))
first_hyperaemia_male.append(np.max(blank_b[160+110:400+110]))
first_hyperaemia_male.append(np.max(blank_c[160+110:400+110])) 
first_hyperaemia_male.append(np.max(blank_d[160+110:400+110]))
second_hyperaemia_male.append(np.max(blank_a[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_b[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_c[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_d[500+110:800+110]))
third_hyperaemia_male.append(np.max(blank_a[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_b[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_c[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_d[900+110:1200+110]))
fourth_hyperaemia_male.append(np.max(blank_a[1300+110:]))
fourth_hyperaemia_male.append(np.max(blank_b[1300+110:]))
fourth_hyperaemia_male.append( np.max(blank_c[1300+110:]))
fourth_hyperaemia_male.append(np.max(blank_d[1300+110:]))
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__033+110):
	blank_a.append(_033_PREVAIL_LBM_033["j1_1"]["PU"][index_11__033+i-110])
	blank_b.append(_033_PREVAIL_LBM_033["j1_2"]["PU"][index_12__033+i-110])
	blank_c.append(_033_PREVAIL_LBM_033["j2_1"]["PU"][index_21__033+i-110])
	blank_d.append(_033_PREVAIL_LBM_033["j2_2"]["PU"][index_22__033+i-110])
first_ischaemia_male.append(np.mean(blank_a[60+110:160+110]))
first_ischaemia_male.append(np.mean(blank_b[60+110:160+110])) 
first_ischaemia_male.append(np.mean(blank_c[60+110:160+110]))
first_ischaemia_male.append(np.mean(blank_d[60+110:160+110]))
second_ischaemia_male.append(np.mean(blank_a[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_b[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_c[400+110:500+110]))
second_ischaemia_male.append(np.mean(blank_d[400+110:500+110]))
third_ischaemia_male.append(np.mean(blank_a[800+110:900+110]))
third_ischaemia_male.append(np.mean(blank_b[800+110:900+110]))
third_ischaemia_male.append(np.mean(blank_c[800+110:900+110])) 
third_ischaemia_male.append(np.mean(blank_d[800+110:900+110]))
fourth_ischaemia_male.append(np.mean(blank_a[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_b[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_c[1200+110:1300+110]))
fourth_ischaemia_male.append(np.mean(blank_d[1200+110:1300+110]))
first_hyperaemia_male.append(np.max(blank_a[160+110:400+110]))
first_hyperaemia_male.append(np.max(blank_b[160+110:400+110]))
first_hyperaemia_male.append(np.max(blank_c[160+110:400+110])) 
first_hyperaemia_male.append(np.max(blank_d[160+110:400+110]))
second_hyperaemia_male.append(np.max(blank_a[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_b[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_c[500+110:800+110]))
second_hyperaemia_male.append(np.max(blank_d[500+110:800+110]))
third_hyperaemia_male.append(np.max(blank_a[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_b[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_c[900+110:1200+110]))
third_hyperaemia_male.append(np.max(blank_d[900+110:1200+110]))
fourth_hyperaemia_male.append(np.max(blank_a[1300+110:]))
fourth_hyperaemia_male.append(np.max(blank_b[1300+110:]))
fourth_hyperaemia_male.append( np.max(blank_c[1300+110:]))
fourth_hyperaemia_male.append(np.max(blank_d[1300+110:]))
# 


# 1 female 28 31 32 34 35
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__028+110):
	blank_a.append(_028_PREVAIL_LBM_028["j1_1"]["PU"][index_11__028+i-110])
	blank_b.append(_028_PREVAIL_LBM_028["j1_2"]["PU"][index_12__028+i-110])
	blank_c.append(_028_PREVAIL_LBM_028["j2_1"]["PU"][index_21__028+i-110])
	blank_d.append(_028_PREVAIL_LBM_028["j2_2"]["PU"][index_22__028+i-110])
first_ischaemia_female.append(np.mean(blank_a[60+110:160+110]))
first_ischaemia_female.append(np.mean(blank_b[60+110:160+110])) 
first_ischaemia_female.append(np.mean(blank_c[60+110:160+110]))
first_ischaemia_female.append(np.mean(blank_d[60+110:160+110]))
second_ischaemia_female.append(np.mean(blank_a[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_b[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_c[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_d[400+110:500+110]))
third_ischaemia_female.append(np.mean(blank_a[800+110:900+110]))
third_ischaemia_female.append(np.mean(blank_b[800+110:900+110]))
third_ischaemia_female.append(np.mean(blank_c[800+110:900+110])) 
third_ischaemia_female.append(np.mean(blank_d[800+110:900+110]))
fourth_ischaemia_female.append(np.mean(blank_a[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_b[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_c[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_d[1200+110:1300+110]))
first_hyperaemia_female.append(np.max(blank_a[160+110:400+110]))
first_hyperaemia_female.append(np.max(blank_b[160+110:400+110]))
first_hyperaemia_female.append(np.max(blank_c[160+110:400+110])) 
first_hyperaemia_female.append(np.max(blank_d[160+110:400+110]))
second_hyperaemia_female.append(np.max(blank_a[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_b[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_c[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_d[500+110:800+110]))
third_hyperaemia_female.append(np.max(blank_a[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_b[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_c[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_d[900+110:1200+110]))
fourth_hyperaemia_female.append(np.max(blank_a[1300+110:]))
fourth_hyperaemia_female.append(np.max(blank_b[1300+110:]))
fourth_hyperaemia_female.append( np.max(blank_c[1300+110:]))
fourth_hyperaemia_female.append(np.max(blank_d[1300+110:]))# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__031+110):
	blank_a.append(_031_PREVAIL_LBM_031["j1_1"]["PU"][index_11__031+i-110])
	blank_b.append(_031_PREVAIL_LBM_031["j1_2"]["PU"][index_12__031+i-110])
	blank_c.append(_031_PREVAIL_LBM_031["j2_1"]["PU"][index_21__031+i-110])
	blank_d.append(_031_PREVAIL_LBM_031["j2_2"]["PU"][index_22__031+i-110])
first_ischaemia_female.append(np.mean(blank_a[60+110:160+110]))
first_ischaemia_female.append(np.mean(blank_b[60+110:160+110])) 
first_ischaemia_female.append(np.mean(blank_c[60+110:160+110]))
first_ischaemia_female.append(np.mean(blank_d[60+110:160+110]))
second_ischaemia_female.append(np.mean(blank_a[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_b[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_c[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_d[400+110:500+110]))
third_ischaemia_female.append(np.mean(blank_a[800+110:900+110]))
third_ischaemia_female.append(np.mean(blank_b[800+110:900+110]))
third_ischaemia_female.append(np.mean(blank_c[800+110:900+110])) 
third_ischaemia_female.append(np.mean(blank_d[800+110:900+110]))
fourth_ischaemia_female.append(np.mean(blank_a[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_b[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_c[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_d[1200+110:1300+110]))
first_hyperaemia_female.append(np.max(blank_a[160+110:400+110]))
first_hyperaemia_female.append(np.max(blank_b[160+110:400+110]))
first_hyperaemia_female.append(np.max(blank_c[160+110:400+110])) 
first_hyperaemia_female.append(np.max(blank_d[160+110:400+110]))
second_hyperaemia_female.append(np.max(blank_a[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_b[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_c[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_d[500+110:800+110]))
third_hyperaemia_female.append(np.max(blank_a[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_b[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_c[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_d[900+110:1200+110]))
fourth_hyperaemia_female.append(np.max(blank_a[1300+110:]))
fourth_hyperaemia_female.append(np.max(blank_b[1300+110:]))
fourth_hyperaemia_female.append( np.max(blank_c[1300+110:]))
fourth_hyperaemia_female.append(np.max(blank_d[1300+110:]))# 
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__032+110):
	blank_a.append(_032_PREVAIL_LBM_032["j1_1"]["PU"][index_11__032+i-110])
	blank_b.append(_032_PREVAIL_LBM_032["j1_2"]["PU"][index_12__032+i-110])
	blank_c.append(_032_PREVAIL_LBM_032["j2_1"]["PU"][index_21__032+i-110])
	blank_d.append(_032_PREVAIL_LBM_032["j2_2"]["PU"][index_22__032+i-110])
first_ischaemia_female.append(np.mean(blank_a[60+110:160+110]))
first_ischaemia_female.append(np.mean(blank_b[60+110:160+110])) 
first_ischaemia_female.append(np.mean(blank_c[60+110:160+110]))
first_ischaemia_female.append(np.mean(blank_d[60+110:160+110]))
second_ischaemia_female.append(np.mean(blank_a[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_b[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_c[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_d[400+110:500+110]))
third_ischaemia_female.append(np.mean(blank_a[800+110:900+110]))
third_ischaemia_female.append(np.mean(blank_b[800+110:900+110]))
third_ischaemia_female.append(np.mean(blank_c[800+110:900+110])) 
third_ischaemia_female.append(np.mean(blank_d[800+110:900+110]))
fourth_ischaemia_female.append(np.mean(blank_a[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_b[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_c[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_d[1200+110:1300+110]))
first_hyperaemia_female.append(np.max(blank_a[160+110:400+110]))
first_hyperaemia_female.append(np.max(blank_b[160+110:400+110]))
first_hyperaemia_female.append(np.max(blank_c[160+110:400+110])) 
first_hyperaemia_female.append(np.max(blank_d[160+110:400+110]))
second_hyperaemia_female.append(np.max(blank_a[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_b[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_c[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_d[500+110:800+110]))
third_hyperaemia_female.append(np.max(blank_a[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_b[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_c[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_d[900+110:1200+110]))
fourth_hyperaemia_female.append(np.max(blank_a[1300+110:]))
fourth_hyperaemia_female.append(np.max(blank_b[1300+110:]))
fourth_hyperaemia_female.append( np.max(blank_c[1300+110:]))
fourth_hyperaemia_female.append(np.max(blank_d[1300+110:]))# 
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__034+110):
	blank_a.append(_034_PREVAIL_LBM_034["j1_1"]["PU"][index_11__034+i-110])
	blank_b.append(_034_PREVAIL_LBM_034["j1_2"]["PU"][index_12__034+i-110])
	blank_c.append(_034_PREVAIL_LBM_034["j2_1"]["PU"][index_21__034+i-110])
	blank_d.append(_034_PREVAIL_LBM_034["j2_2"]["PU"][index_22__034+i-110])
first_ischaemia_female.append(np.mean(blank_a[60+110:160+110]))
first_ischaemia_female.append(np.mean(blank_b[60+110:160+110])) 
first_ischaemia_female.append(np.mean(blank_c[60+110:160+110]))
first_ischaemia_female.append(np.mean(blank_d[60+110:160+110]))
second_ischaemia_female.append(np.mean(blank_a[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_b[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_c[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_d[400+110:500+110]))
third_ischaemia_female.append(np.mean(blank_a[800+110:900+110]))
third_ischaemia_female.append(np.mean(blank_b[800+110:900+110]))
third_ischaemia_female.append(np.mean(blank_c[800+110:900+110])) 
third_ischaemia_female.append(np.mean(blank_d[800+110:900+110]))
fourth_ischaemia_female.append(np.mean(blank_a[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_b[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_c[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_d[1200+110:1300+110]))
first_hyperaemia_female.append(np.max(blank_a[160+110:400+110]))
first_hyperaemia_female.append(np.max(blank_b[160+110:400+110]))
first_hyperaemia_female.append(np.max(blank_c[160+110:400+110])) 
first_hyperaemia_female.append(np.max(blank_d[160+110:400+110]))
second_hyperaemia_female.append(np.max(blank_a[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_b[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_c[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_d[500+110:800+110]))
third_hyperaemia_female.append(np.max(blank_a[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_b[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_c[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_d[900+110:1200+110]))
fourth_hyperaemia_female.append(np.max(blank_a[1300+110:]))
fourth_hyperaemia_female.append(np.max(blank_b[1300+110:]))
fourth_hyperaemia_female.append( np.max(blank_c[1300+110:]))
fourth_hyperaemia_female.append(np.max(blank_d[1300+110:]))# 
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__035+110):
	blank_a.append(_035_PREVAIL_LBM_035["j1_1"]["PU"][index_11__035+i-110])
	blank_b.append(_035_PREVAIL_LBM_035["j1_2"]["PU"][index_12__035+i-110])
	blank_c.append(_035_PREVAIL_LBM_035["j2_1"]["PU"][index_21__035+i-110])
	blank_d.append(_035_PREVAIL_LBM_035["j2_2"]["PU"][index_22__035+i-110])
first_ischaemia_female.append(np.mean(blank_a[60+110:160+110]))
first_ischaemia_female.append(np.mean(blank_b[60+110:160+110])) 
first_ischaemia_female.append(np.mean(blank_c[60+110:160+110]))
first_ischaemia_female.append(np.mean(blank_d[60+110:160+110]))
second_ischaemia_female.append(np.mean(blank_a[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_b[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_c[400+110:500+110]))
second_ischaemia_female.append(np.mean(blank_d[400+110:500+110]))
third_ischaemia_female.append(np.mean(blank_a[800+110:900+110]))
third_ischaemia_female.append(np.mean(blank_b[800+110:900+110]))
third_ischaemia_female.append(np.mean(blank_c[800+110:900+110])) 
third_ischaemia_female.append(np.mean(blank_d[800+110:900+110]))
fourth_ischaemia_female.append(np.mean(blank_a[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_b[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_c[1200+110:1300+110]))
fourth_ischaemia_female.append(np.mean(blank_d[1200+110:1300+110]))
first_hyperaemia_female.append(np.max(blank_a[160+110:400+110]))
first_hyperaemia_female.append(np.max(blank_b[160+110:400+110]))
first_hyperaemia_female.append(np.max(blank_c[160+110:400+110])) 
first_hyperaemia_female.append(np.max(blank_d[160+110:400+110]))
second_hyperaemia_female.append(np.max(blank_a[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_b[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_c[500+110:800+110]))
second_hyperaemia_female.append(np.max(blank_d[500+110:800+110]))
third_hyperaemia_female.append(np.max(blank_a[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_b[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_c[900+110:1200+110]))
third_hyperaemia_female.append(np.max(blank_d[900+110:1200+110]))
fourth_hyperaemia_female.append(np.max(blank_a[1300+110:]))
fourth_hyperaemia_female.append(np.max(blank_b[1300+110:]))
fourth_hyperaemia_female.append( np.max(blank_c[1300+110:]))
fourth_hyperaemia_female.append(np.max(blank_d[1300+110:]))# 
# 






















# pc
# 1 male
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list):
	blank_a.append(_008_PREVAIL_LBM_008["j1_1"]["PU_pc"][index_11+i])
	blank_b.append(_008_PREVAIL_LBM_008["j1_2"]["PU_pc"][index_12+i])
	blank_c.append(_008_PREVAIL_LBM_008["j2_1"]["PU_pc"][index_21+i])
	blank_d.append(_008_PREVAIL_LBM_008["j2_2"]["PU_pc"][index_22+i])
first_ischaemia_male_pc.append(np.mean(blank_a[60:160]))
first_ischaemia_male_pc.append(np.mean(blank_b[60:160])) 
first_ischaemia_male_pc.append(np.mean(blank_c[60:160]))
first_ischaemia_male_pc.append(np.mean(blank_d[60:160]))
second_ischaemia_male_pc.append(np.mean(blank_a[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_b[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_c[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_d[400:500]))
third_ischaemia_male_pc.append(np.mean(blank_a[800:900]))
third_ischaemia_male_pc.append(np.mean(blank_b[800:900]))
third_ischaemia_male_pc.append(np.mean(blank_c[800:900])) 
third_ischaemia_male_pc.append(np.mean(blank_d[800:900]))
fourth_ischaemia_male_pc.append(np.mean(blank_a[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_b[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_c[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_d[1200:1300]))
first_hyperaemia_male_pc.append(np.max(blank_a[160:400]))
first_hyperaemia_male_pc.append(np.max(blank_b[160:400]))
first_hyperaemia_male_pc.append(np.max(blank_c[160:400])) 
first_hyperaemia_male_pc.append(np.max(blank_d[160:400]))
second_hyperaemia_male_pc.append(np.max(blank_a[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_b[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_c[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_d[500:800]))
third_hyperaemia_male_pc.append(np.max(blank_a[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_b[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_c[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_d[900:1200]))
fourth_hyperaemia_male_pc.append(np.max(blank_a[1300:]))
fourth_hyperaemia_male_pc.append(np.max(blank_b[1300:]))
fourth_hyperaemia_male_pc.append( np.max(blank_c[1300:]))
fourth_hyperaemia_male_pc.append(np.max(blank_d[1300:]))
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__009):
	blank_a.append(_009_PREVAIL_LBM_009["j1_1"]["PU_pc"][index_11__009+i])
	blank_b.append(_009_PREVAIL_LBM_009["j1_2"]["PU_pc"][index_12__009+i])
	blank_c.append(_009_PREVAIL_LBM_009["j2_1"]["PU_pc"][index_21__009+i])
	blank_d.append(_009_PREVAIL_LBM_009["j2_2"]["PU_pc"][index_22__009+i])
first_ischaemia_male_pc.append(np.mean(blank_a[60:160]))
first_ischaemia_male_pc.append(np.mean(blank_b[60:160])) 
first_ischaemia_male_pc.append(np.mean(blank_c[60:160]))
first_ischaemia_male_pc.append(np.mean(blank_d[60:160]))
second_ischaemia_male_pc.append(np.mean(blank_a[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_b[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_c[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_d[400:500]))
third_ischaemia_male_pc.append(np.mean(blank_a[800:900]))
third_ischaemia_male_pc.append(np.mean(blank_b[800:900]))
third_ischaemia_male_pc.append(np.mean(blank_c[800:900])) 
third_ischaemia_male_pc.append(np.mean(blank_d[800:900]))
fourth_ischaemia_male_pc.append(np.mean(blank_a[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_b[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_c[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_d[1200:1300]))
first_hyperaemia_male_pc.append(np.max(blank_a[160:400]))
first_hyperaemia_male_pc.append(np.max(blank_b[160:400]))
first_hyperaemia_male_pc.append(np.max(blank_c[160:400])) 
first_hyperaemia_male_pc.append(np.max(blank_d[160:400]))
second_hyperaemia_male_pc.append(np.max(blank_a[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_b[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_c[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_d[500:800]))
third_hyperaemia_male_pc.append(np.max(blank_a[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_b[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_c[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_d[900:1200]))
fourth_hyperaemia_male_pc.append(np.max(blank_a[1300:]))
fourth_hyperaemia_male_pc.append(np.max(blank_b[1300:]))
fourth_hyperaemia_male_pc.append( np.max(blank_c[1300:]))
fourth_hyperaemia_male_pc.append(np.max(blank_d[1300:]))
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__026):
	blank_a.append(_026_PREVAIL_LBM_026["j1_1"]["PU_pc"][index_11__026+i])
	blank_b.append(_026_PREVAIL_LBM_026["j1_2"]["PU_pc"][index_12__026+i])
	blank_c.append(_026_PREVAIL_LBM_026["j2_1"]["PU_pc"][index_21__026+i])
	blank_d.append(_026_PREVAIL_LBM_026["j2_2"]["PU_pc"][index_22__026+i])
first_ischaemia_male_pc.append(np.mean(blank_a[60:160]))
first_ischaemia_male_pc.append(np.mean(blank_b[60:160])) 
first_ischaemia_male_pc.append(np.mean(blank_c[60:160]))
first_ischaemia_male_pc.append(np.mean(blank_d[60:160]))
second_ischaemia_male_pc.append(np.mean(blank_a[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_b[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_c[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_d[400:500]))
third_ischaemia_male_pc.append(np.mean(blank_a[800:900]))
third_ischaemia_male_pc.append(np.mean(blank_b[800:900]))
third_ischaemia_male_pc.append(np.mean(blank_c[800:900])) 
third_ischaemia_male_pc.append(np.mean(blank_d[800:900]))
fourth_ischaemia_male_pc.append(np.mean(blank_a[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_b[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_c[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_d[1200:1300]))
first_hyperaemia_male_pc.append(np.max(blank_a[160:400]))
first_hyperaemia_male_pc.append(np.max(blank_b[160:400]))
first_hyperaemia_male_pc.append(np.max(blank_c[160:400])) 
first_hyperaemia_male_pc.append(np.max(blank_d[160:400]))
second_hyperaemia_male_pc.append(np.max(blank_a[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_b[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_c[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_d[500:800]))
third_hyperaemia_male_pc.append(np.max(blank_a[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_b[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_c[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_d[900:1200]))
fourth_hyperaemia_male_pc.append(np.max(blank_a[1300:]))
fourth_hyperaemia_male_pc.append(np.max(blank_b[1300:]))
fourth_hyperaemia_male_pc.append( np.max(blank_c[1300:]))
fourth_hyperaemia_male_pc.append(np.max(blank_d[1300:]))
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__029):
	blank_a.append(_029_PREVAIL_LBM_029["j1_1"]["PU_pc"][index_11__029+i])
	blank_b.append(_029_PREVAIL_LBM_029["j1_2"]["PU_pc"][index_12__029+i])
	blank_c.append(_029_PREVAIL_LBM_029["j2_1"]["PU_pc"][index_21__029+i])
	blank_d.append(_029_PREVAIL_LBM_029["j2_2"]["PU_pc"][index_22__029+i])
first_ischaemia_male_pc.append(np.mean(blank_a[60:160]))
first_ischaemia_male_pc.append(np.mean(blank_b[60:160])) 
first_ischaemia_male_pc.append(np.mean(blank_c[60:160]))
first_ischaemia_male_pc.append(np.mean(blank_d[60:160]))
second_ischaemia_male_pc.append(np.mean(blank_a[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_b[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_c[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_d[400:500]))
third_ischaemia_male_pc.append(np.mean(blank_a[800:900]))
third_ischaemia_male_pc.append(np.mean(blank_b[800:900]))
third_ischaemia_male_pc.append(np.mean(blank_c[800:900])) 
third_ischaemia_male_pc.append(np.mean(blank_d[800:900]))
fourth_ischaemia_male_pc.append(np.mean(blank_a[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_b[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_c[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_d[1200:1300]))
first_hyperaemia_male_pc.append(np.max(blank_a[160:400]))
first_hyperaemia_male_pc.append(np.max(blank_b[160:400]))
first_hyperaemia_male_pc.append(np.max(blank_c[160:400])) 
first_hyperaemia_male_pc.append(np.max(blank_d[160:400]))
second_hyperaemia_male_pc.append(np.max(blank_a[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_b[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_c[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_d[500:800]))
third_hyperaemia_male_pc.append(np.max(blank_a[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_b[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_c[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_d[900:1200]))
fourth_hyperaemia_male_pc.append(np.max(blank_a[1300:]))
fourth_hyperaemia_male_pc.append(np.max(blank_b[1300:]))
fourth_hyperaemia_male_pc.append( np.max(blank_c[1300:]))
fourth_hyperaemia_male_pc.append(np.max(blank_d[1300:]))
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__030):
	blank_a.append(_030_PREVAIL_LBM_030["j1_1"]["PU_pc"][index_11__030+i])
	blank_b.append(_030_PREVAIL_LBM_030["j1_2"]["PU_pc"][index_12__030+i])
	blank_c.append(_030_PREVAIL_LBM_030["j2_1"]["PU_pc"][index_21__030+i])
	blank_d.append(_030_PREVAIL_LBM_030["j2_2"]["PU_pc"][index_22__030+i])
first_ischaemia_male_pc.append(np.mean(blank_a[60:160]))
first_ischaemia_male_pc.append(np.mean(blank_b[60:160])) 
first_ischaemia_male_pc.append(np.mean(blank_c[60:160]))
first_ischaemia_male_pc.append(np.mean(blank_d[60:160]))
second_ischaemia_male_pc.append(np.mean(blank_a[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_b[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_c[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_d[400:500]))
third_ischaemia_male_pc.append(np.mean(blank_a[800:900]))
third_ischaemia_male_pc.append(np.mean(blank_b[800:900]))
third_ischaemia_male_pc.append(np.mean(blank_c[800:900])) 
third_ischaemia_male_pc.append(np.mean(blank_d[800:900]))
fourth_ischaemia_male_pc.append(np.mean(blank_a[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_b[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_c[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_d[1200:1300]))
first_hyperaemia_male_pc.append(np.max(blank_a[160:400]))
first_hyperaemia_male_pc.append(np.max(blank_b[160:400]))
first_hyperaemia_male_pc.append(np.max(blank_c[160:400])) 
first_hyperaemia_male_pc.append(np.max(blank_d[160:400]))
second_hyperaemia_male_pc.append(np.max(blank_a[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_b[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_c[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_d[500:800]))
third_hyperaemia_male_pc.append(np.max(blank_a[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_b[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_c[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_d[900:1200]))
fourth_hyperaemia_male_pc.append(np.max(blank_a[1300:]))
fourth_hyperaemia_male_pc.append(np.max(blank_b[1300:]))
fourth_hyperaemia_male_pc.append( np.max(blank_c[1300:]))
fourth_hyperaemia_male_pc.append(np.max(blank_d[1300:]))
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__033):
	blank_a.append(_033_PREVAIL_LBM_033["j1_1"]["PU_pc"][index_11__033+i])
	blank_b.append(_033_PREVAIL_LBM_033["j1_2"]["PU_pc"][index_12__033+i])
	blank_c.append(_033_PREVAIL_LBM_033["j2_1"]["PU_pc"][index_21__033+i])
	blank_d.append(_033_PREVAIL_LBM_033["j2_2"]["PU_pc"][index_22__033+i])
first_ischaemia_male_pc.append(np.mean(blank_a[60:160]))
first_ischaemia_male_pc.append(np.mean(blank_b[60:160])) 
first_ischaemia_male_pc.append(np.mean(blank_c[60:160]))
first_ischaemia_male_pc.append(np.mean(blank_d[60:160]))
second_ischaemia_male_pc.append(np.mean(blank_a[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_b[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_c[400:500]))
second_ischaemia_male_pc.append(np.mean(blank_d[400:500]))
third_ischaemia_male_pc.append(np.mean(blank_a[800:900]))
third_ischaemia_male_pc.append(np.mean(blank_b[800:900]))
third_ischaemia_male_pc.append(np.mean(blank_c[800:900])) 
third_ischaemia_male_pc.append(np.mean(blank_d[800:900]))
fourth_ischaemia_male_pc.append(np.mean(blank_a[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_b[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_c[1200:1300]))
fourth_ischaemia_male_pc.append(np.mean(blank_d[1200:1300]))
first_hyperaemia_male_pc.append(np.max(blank_a[160:400]))
first_hyperaemia_male_pc.append(np.max(blank_b[160:400]))
first_hyperaemia_male_pc.append(np.max(blank_c[160:400])) 
first_hyperaemia_male_pc.append(np.max(blank_d[160:400]))
second_hyperaemia_male_pc.append(np.max(blank_a[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_b[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_c[500:800]))
second_hyperaemia_male_pc.append(np.max(blank_d[500:800]))
third_hyperaemia_male_pc.append(np.max(blank_a[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_b[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_c[900:1200]))
third_hyperaemia_male_pc.append(np.max(blank_d[900:1200]))
fourth_hyperaemia_male_pc.append(np.max(blank_a[1300:]))
fourth_hyperaemia_male_pc.append(np.max(blank_b[1300:]))
fourth_hyperaemia_male_pc.append( np.max(blank_c[1300:]))
fourth_hyperaemia_male_pc.append(np.max(blank_d[1300:]))
# 


# 1 female 28 31 32 34 35
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__028):
	blank_a.append(_028_PREVAIL_LBM_028["j1_1"]["PU_pc"][index_11__028+i])
	blank_b.append(_028_PREVAIL_LBM_028["j1_2"]["PU_pc"][index_12__028+i])
	blank_c.append(_028_PREVAIL_LBM_028["j2_1"]["PU_pc"][index_21__028+i])
	blank_d.append(_028_PREVAIL_LBM_028["j2_2"]["PU_pc"][index_22__028+i])
first_ischaemia_female_pc.append(np.mean(blank_a[60:160]))
first_ischaemia_female_pc.append(np.mean(blank_b[60:160])) 
first_ischaemia_female_pc.append(np.mean(blank_c[60:160]))
first_ischaemia_female_pc.append(np.mean(blank_d[60:160]))
second_ischaemia_female_pc.append(np.mean(blank_a[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_b[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_c[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_d[400:500]))
third_ischaemia_female_pc.append(np.mean(blank_a[800:900]))
third_ischaemia_female_pc.append(np.mean(blank_b[800:900]))
third_ischaemia_female_pc.append(np.mean(blank_c[800:900])) 
third_ischaemia_female_pc.append(np.mean(blank_d[800:900]))
fourth_ischaemia_female_pc.append(np.mean(blank_a[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_b[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_c[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_d[1200:1300]))
first_hyperaemia_female_pc.append(np.max(blank_a[160:400]))
first_hyperaemia_female_pc.append(np.max(blank_b[160:400]))
first_hyperaemia_female_pc.append(np.max(blank_c[160:400])) 
first_hyperaemia_female_pc.append(np.max(blank_d[160:400]))
second_hyperaemia_female_pc.append(np.max(blank_a[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_b[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_c[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_d[500:800]))
third_hyperaemia_female_pc.append(np.max(blank_a[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_b[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_c[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_d[900:1200]))
fourth_hyperaemia_female_pc.append(np.max(blank_a[1300:]))
fourth_hyperaemia_female_pc.append(np.max(blank_b[1300:]))
fourth_hyperaemia_female_pc.append( np.max(blank_c[1300:]))
fourth_hyperaemia_female_pc.append(np.max(blank_d[1300:]))# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__031):
	blank_a.append(_031_PREVAIL_LBM_031["j1_1"]["PU_pc"][index_11__031+i])
	blank_b.append(_031_PREVAIL_LBM_031["j1_2"]["PU_pc"][index_12__031+i])
	blank_c.append(_031_PREVAIL_LBM_031["j2_1"]["PU_pc"][index_21__031+i])
	blank_d.append(_031_PREVAIL_LBM_031["j2_2"]["PU_pc"][index_22__031+i])
first_ischaemia_female_pc.append(np.mean(blank_a[60:160]))
first_ischaemia_female_pc.append(np.mean(blank_b[60:160])) 
first_ischaemia_female_pc.append(np.mean(blank_c[60:160]))
first_ischaemia_female_pc.append(np.mean(blank_d[60:160]))
second_ischaemia_female_pc.append(np.mean(blank_a[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_b[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_c[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_d[400:500]))
third_ischaemia_female_pc.append(np.mean(blank_a[800:900]))
third_ischaemia_female_pc.append(np.mean(blank_b[800:900]))
third_ischaemia_female_pc.append(np.mean(blank_c[800:900])) 
third_ischaemia_female_pc.append(np.mean(blank_d[800:900]))
fourth_ischaemia_female_pc.append(np.mean(blank_a[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_b[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_c[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_d[1200:1300]))
first_hyperaemia_female_pc.append(np.max(blank_a[160:400]))
first_hyperaemia_female_pc.append(np.max(blank_b[160:400]))
first_hyperaemia_female_pc.append(np.max(blank_c[160:400])) 
first_hyperaemia_female_pc.append(np.max(blank_d[160:400]))
second_hyperaemia_female_pc.append(np.max(blank_a[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_b[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_c[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_d[500:800]))
third_hyperaemia_female_pc.append(np.max(blank_a[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_b[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_c[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_d[900:1200]))
fourth_hyperaemia_female_pc.append(np.max(blank_a[1300:]))
fourth_hyperaemia_female_pc.append(np.max(blank_b[1300:]))
fourth_hyperaemia_female_pc.append( np.max(blank_c[1300:]))
fourth_hyperaemia_female_pc.append(np.max(blank_d[1300:]))# 
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__032):
	blank_a.append(_032_PREVAIL_LBM_032["j1_1"]["PU_pc"][index_11__032+i])
	blank_b.append(_032_PREVAIL_LBM_032["j1_2"]["PU_pc"][index_12__032+i])
	blank_c.append(_032_PREVAIL_LBM_032["j2_1"]["PU_pc"][index_21__032+i])
	blank_d.append(_032_PREVAIL_LBM_032["j2_2"]["PU_pc"][index_22__032+i])
first_ischaemia_female_pc.append(np.mean(blank_a[60:160]))
first_ischaemia_female_pc.append(np.mean(blank_b[60:160])) 
first_ischaemia_female_pc.append(np.mean(blank_c[60:160]))
first_ischaemia_female_pc.append(np.mean(blank_d[60:160]))
second_ischaemia_female_pc.append(np.mean(blank_a[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_b[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_c[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_d[400:500]))
third_ischaemia_female_pc.append(np.mean(blank_a[800:900]))
third_ischaemia_female_pc.append(np.mean(blank_b[800:900]))
third_ischaemia_female_pc.append(np.mean(blank_c[800:900])) 
third_ischaemia_female_pc.append(np.mean(blank_d[800:900]))
fourth_ischaemia_female_pc.append(np.mean(blank_a[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_b[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_c[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_d[1200:1300]))
first_hyperaemia_female_pc.append(np.max(blank_a[160:400]))
first_hyperaemia_female_pc.append(np.max(blank_b[160:400]))
first_hyperaemia_female_pc.append(np.max(blank_c[160:400])) 
first_hyperaemia_female_pc.append(np.max(blank_d[160:400]))
second_hyperaemia_female_pc.append(np.max(blank_a[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_b[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_c[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_d[500:800]))
third_hyperaemia_female_pc.append(np.max(blank_a[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_b[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_c[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_d[900:1200]))
fourth_hyperaemia_female_pc.append(np.max(blank_a[1300:]))
fourth_hyperaemia_female_pc.append(np.max(blank_b[1300:]))
fourth_hyperaemia_female_pc.append( np.max(blank_c[1300:]))
fourth_hyperaemia_female_pc.append(np.max(blank_d[1300:]))# 
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__034):
	blank_a.append(_034_PREVAIL_LBM_034["j1_1"]["PU_pc"][index_11__034+i])
	blank_b.append(_034_PREVAIL_LBM_034["j1_2"]["PU_pc"][index_12__034+i])
	blank_c.append(_034_PREVAIL_LBM_034["j2_1"]["PU_pc"][index_21__034+i])
	blank_d.append(_034_PREVAIL_LBM_034["j2_2"]["PU_pc"][index_22__034+i])
first_ischaemia_female_pc.append(np.mean(blank_a[60:160]))
first_ischaemia_female_pc.append(np.mean(blank_b[60:160])) 
first_ischaemia_female_pc.append(np.mean(blank_c[60:160]))
first_ischaemia_female_pc.append(np.mean(blank_d[60:160]))
second_ischaemia_female_pc.append(np.mean(blank_a[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_b[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_c[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_d[400:500]))
third_ischaemia_female_pc.append(np.mean(blank_a[800:900]))
third_ischaemia_female_pc.append(np.mean(blank_b[800:900]))
third_ischaemia_female_pc.append(np.mean(blank_c[800:900])) 
third_ischaemia_female_pc.append(np.mean(blank_d[800:900]))
fourth_ischaemia_female_pc.append(np.mean(blank_a[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_b[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_c[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_d[1200:1300]))
first_hyperaemia_female_pc.append(np.max(blank_a[160:400]))
first_hyperaemia_female_pc.append(np.max(blank_b[160:400]))
first_hyperaemia_female_pc.append(np.max(blank_c[160:400])) 
first_hyperaemia_female_pc.append(np.max(blank_d[160:400]))
second_hyperaemia_female_pc.append(np.max(blank_a[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_b[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_c[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_d[500:800]))
third_hyperaemia_female_pc.append(np.max(blank_a[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_b[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_c[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_d[900:1200]))
fourth_hyperaemia_female_pc.append(np.max(blank_a[1300:]))
fourth_hyperaemia_female_pc.append(np.max(blank_b[1300:]))
fourth_hyperaemia_female_pc.append( np.max(blank_c[1300:]))
fourth_hyperaemia_female_pc.append(np.max(blank_d[1300:]))# 
# 
blank_a=[]
blank_b=[]
blank_c=[]
blank_d=[]
for i in range(length_list__035):
	blank_a.append(_035_PREVAIL_LBM_035["j1_1"]["PU_pc"][index_11__035+i])
	blank_b.append(_035_PREVAIL_LBM_035["j1_2"]["PU_pc"][index_12__035+i])
	blank_c.append(_035_PREVAIL_LBM_035["j2_1"]["PU_pc"][index_21__035+i])
	blank_d.append(_035_PREVAIL_LBM_035["j2_2"]["PU_pc"][index_22__035+i])
first_ischaemia_female_pc.append(np.mean(blank_a[60:160]))
first_ischaemia_female_pc.append(np.mean(blank_b[60:160])) 
first_ischaemia_female_pc.append(np.mean(blank_c[60:160]))
first_ischaemia_female_pc.append(np.mean(blank_d[60:160]))
second_ischaemia_female_pc.append(np.mean(blank_a[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_b[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_c[400:500]))
second_ischaemia_female_pc.append(np.mean(blank_d[400:500]))
third_ischaemia_female_pc.append(np.mean(blank_a[800:900]))
third_ischaemia_female_pc.append(np.mean(blank_b[800:900]))
third_ischaemia_female_pc.append(np.mean(blank_c[800:900])) 
third_ischaemia_female_pc.append(np.mean(blank_d[800:900]))
fourth_ischaemia_female_pc.append(np.mean(blank_a[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_b[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_c[1200:1300]))
fourth_ischaemia_female_pc.append(np.mean(blank_d[1200:1300]))
first_hyperaemia_female_pc.append(np.max(blank_a[160:400]))
first_hyperaemia_female_pc.append(np.max(blank_b[160:400]))
first_hyperaemia_female_pc.append(np.max(blank_c[160:400])) 
first_hyperaemia_female_pc.append(np.max(blank_d[160:400]))
second_hyperaemia_female_pc.append(np.max(blank_a[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_b[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_c[500:800]))
second_hyperaemia_female_pc.append(np.max(blank_d[500:800]))
third_hyperaemia_female_pc.append(np.max(blank_a[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_b[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_c[900:1200]))
third_hyperaemia_female_pc.append(np.max(blank_d[900:1200]))
fourth_hyperaemia_female_pc.append(np.max(blank_a[1300:]))
fourth_hyperaemia_female_pc.append(np.max(blank_b[1300:]))
fourth_hyperaemia_female_pc.append( np.max(blank_c[1300:]))
fourth_hyperaemia_female_pc.append(np.max(blank_d[1300:]))# 

print(fourth_hyperaemia_female_pc)

print(fourth_hyperaemia_male_pc)


from scipy import stats  
t_stat, p_val = stats.ttest_1samp(a=fourth_hyperaemia_female_pc, popmean = np.mean(fourth_hyperaemia_female_pc))
print("t-statistic = " + str(t_stat))  
print("p-value = " + str(p_val))

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.normaltest.html
t_stat, p_val = stats.normaltest(fourth_hyperaemia_female_pc)
print("t-statistic = " + str(t_stat))  
print("p-value = " + str(p_val))