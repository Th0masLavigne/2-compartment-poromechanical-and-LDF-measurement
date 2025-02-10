
def convert_to_seconds(time_str):
    days, time = time_str.split('-')
    hours, minutes, seconds = time.split(':')  
    days = int(days)
    hours = int(hours)
    minutes = int(minutes)
    seconds = int(seconds)
    total_seconds = days * 86400 + hours * 3600 + minutes * 60 + seconds
    return total_seconds

import pandas as pd
import matplotlib.pyplot as plt
data = pd.read_csv('speedup.dat', delimiter=" ")
data["Time"]=data['Time'].apply(convert_to_seconds)
data['Speedup'] = data['Time'].iloc[0] / data['Time']


plt.plot(data['CPU'], data['CPU'], label="ideal")
plt.plot(data['CPU'], data['Speedup'], 'o-', label="real")
plt.xlabel('Number of CPU')
plt.ylabel('speedup')
plt.legend()
plt.show()
