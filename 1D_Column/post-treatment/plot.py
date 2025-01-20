import pandas as pd
import matplotlib.pyplot as plt


# Import the CSV file as a Pandas DataFrame
data_table_1 = pd.read_csv('Results_2_comp_1.csv')
data_table_5 = pd.read_csv('Results_2_comp_5.csv')
data_table_20 = pd.read_csv('Results_2_comp_20.csv')
# avec epsilon 0.6
# data_table_epsilon0595= pd.read_csv('Results_2_comp_epsilon0595.csv')
# avec epsilon 0.5 => cell dans solid
data_table_mono_05= pd.read_csv('Results_2_comp_mono_05.csv')

index_ti,index_t2=15,90


# Access specific columns, for example:
load = data_table_1['load_all'][index_ti:index_t2]
load_ = data_table_1['load_all'][index_ti:]
time = data_table_1['time_all'][index_ti:index_t2]
time = [x-15 for x in time]
time2 = data_table_1['time_all'][index_ti:]
time2 = [x-15 for x in time2]

varepsilonb_1 = data_table_1['varepsilonb_all'][index_ti:]
varepsilonb_5 = data_table_5['varepsilonb_all'][index_ti:]
varepsilonb_20 = data_table_20['varepsilonb_all'][index_ti:]
varepsilonb_epsilon0595  = data_table_mono_05['varepsilonb_all'][index_ti:]

varepsilon_1 = data_table_1['varepsilon_all'][index_ti:]
varepsilon_5 = data_table_5['varepsilon_all'][index_ti:]
varepsilon_20 = data_table_20['varepsilon_all'][index_ti:]
varepsilon_epsilon0595  = data_table_mono_05['varepsilon_all'][index_ti:]


pressure_blood_1 = data_table_1['pressure_blood_all'][index_ti:index_t2]
pressure_blood_5 = data_table_5['pressure_blood_all'][index_ti:index_t2]
pressure_blood_20 = data_table_20['pressure_blood_all'][index_ti:index_t2]
pressure_blood_epsilon0595  = data_table_mono_05['pressure_blood_all'][index_ti:index_t2]


pressure_solid_1 = data_table_1['pressure_solid_all'][index_ti:index_t2]
pressure_solid_5 = data_table_5['pressure_solid_all'][index_ti:index_t2]
pressure_solid_20 = data_table_20['pressure_solid_all'][index_ti:index_t2]
pressure_solid_epsilon0595  = data_table_mono_05['pressure_solid_all'][index_ti:index_t2]

pressure_IF_1 = data_table_1['pressure_IF_all'][index_ti:index_t2]
pressure_IF_5 = data_table_5['pressure_IF_all'][index_ti:index_t2]
pressure_IF_20 = data_table_20['pressure_IF_all'][index_ti:index_t2]
pressure_IF_epsilon0595  = data_table_mono_05['pressure_IF_all'][index_ti:index_t2]

pressure_cell_1 = data_table_1['pressure_cell_all'][index_ti:index_t2]
pressure_cell_5 = data_table_5['pressure_cell_all'][index_ti:index_t2]
pressure_cell_20 = data_table_20['pressure_cell_all'][index_ti:index_t2]
pressure_cell_epsilon0595  = data_table_mono_05['pressure_cell_all'][index_ti:index_t2]




pressure_blood_1_ = data_table_1['pressure_blood_all'][index_ti:]
pressure_blood_5_ = data_table_5['pressure_blood_all'][index_ti:]
pressure_blood_20_ = data_table_20['pressure_blood_all'][index_ti:]
pressure_blood_epsilon0595_  = data_table_mono_05['pressure_blood_all'][index_ti:]


pressure_solid_1_ = data_table_1['pressure_solid_all'][index_ti:]
pressure_solid_5_ = data_table_5['pressure_solid_all'][index_ti:]
pressure_solid_20_ = data_table_20['pressure_solid_all'][index_ti:]
pressure_solid_epsilon0595_  = data_table_mono_05['pressure_solid_all'][index_ti:]

pressure_IF_1_ = data_table_1['pressure_IF_all'][index_ti:]
pressure_IF_5_ = data_table_5['pressure_IF_all'][index_ti:]
pressure_IF_20_ = data_table_20['pressure_IF_all'][index_ti:]
pressure_IF_epsilon0595_  = data_table_mono_05['pressure_IF_all'][index_ti:]

pressure_cell_1_ = data_table_1['pressure_cell_all'][index_ti:]
pressure_cell_5_ = data_table_5['pressure_cell_all'][index_ti:]
pressure_cell_20_ = data_table_20['pressure_cell_all'][index_ti:]
pressure_cell_epsilon0595_  = data_table_mono_05['pressure_cell_all'][index_ti:]




plt.rcParams.update({'font.size': 15})


fig1, ax1 = plt.subplots()
# ax1.plot(time,pressure_IF_all[index_ti:index_t2],linestyle='-',linewidth=3,color='navy', label='IF')
ax1.plot(time,load,linestyle='-',linewidth=3,color='black', label='LOAD')
ax1.plot(time,pressure_solid_epsilon0595,linestyle='--',linewidth=3,color='forestgreen', label='S, mono',alpha=1)
ax1.plot(time,pressure_blood_epsilon0595,linestyle='--',linewidth=3,color='forestgreen', label='B, mono')
ax1.plot(time,pressure_solid_20,linestyle=':',linewidth=3,color='cornflowerblue', label='S, $\\mu^c=20$Pa.s')
ax1.plot(time,pressure_blood_1,linestyle='-',linewidth=3,color='salmon', label='B, $\\mu^c=1$Pa.s')
ax1.plot(time,pressure_blood_5,linestyle='-.',linewidth=3,color='salmon', label='B, $\\mu^c=5$Pa.s')
ax1.plot(time,pressure_blood_20,linestyle=':',linewidth=3,color='salmon', label='B, $\\mu^c=20$Pa.s')
ax1.plot(time,pressure_solid_1,linestyle='-',linewidth=3,color='cornflowerblue', label='S, $\\mu^c=1$Pa.s')
ax1.plot(time,pressure_solid_5,linestyle='-.',linewidth=3,color='cornflowerblue', label='S, $\\mu^c=5$Pa.s')
# ax1.plot(ti2e_all,pressure_cell_all,linestyle='-.',linewidth=3,color='gray', label='C')
ax1.plot([-10, 30],[0, 0],linewidth=1,linestyle='-',color='darkgray')
ax1.plot([0, 0],[-20, 220],linewidth=1,linestyle='-',color='darkgray')
t=plt.text(20, 195, r'$Load$', fontsize = 12, color = 'black' )
t.set_bbox(dict(facecolor='white', alpha=0.7, linewidth=0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8,alpha=0.5)
plt.text(7, 105, r'$p^s~at~the~bottom~points$', fontsize = 12, color = 'cornflowerblue')
plt.text(7, 30, r'$p^b~at~the~bottom~points$', fontsize = 12, color = 'salmon')
# t=plt.text(5, 157, r'$4\%$', fontsize = 12, color = 'cornflowerblue')
# t.set_bbox(dict(facecolor='white', alpha=0.7, linewidth=0))
ax1.set_xlim([-10,30])
ax1.set_ylim([-20,220])
ax1.set_xlabel('time (s)')
ax1.set_ylabel('Pressure (Pa)')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
# ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
# ax1.legend()
fig1.tight_layout()
fig1.savefig('PT_pressures.jpg')

fig1, ax1 = plt.subplots()
ax1.plot(time,pressure_IF_1,linestyle='-',linewidth=3,color='tan', label='IF, $\\mu^c=1$Pa.s')
ax1.plot(time,pressure_IF_5,linestyle='-.',linewidth=3,color='tan', label='IF, $\\mu^c=5$Pa.s')
ax1.plot(time,pressure_IF_20,linestyle=':',linewidth=3,color='tan', label='IF, $\\mu^c=20$Pa.s')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('Intertitial Fluid Pressure (Pa)')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
# ax1.legend()
fig1.tight_layout()
fig1.savefig('PT_pressuresIF.jpg')


fig1, ax1 = plt.subplots()
ax1.plot(time,pressure_cell_1,linestyle='-',linewidth=3,color='thistle', label='C, $\\mu^c=1$Pa.s')
ax1.plot(time,pressure_cell_5,linestyle='-.',linewidth=3,color='thistle', label='C, $\\mu^c=5$Pa.s')
ax1.plot(time,pressure_cell_20,linestyle=':',linewidth=3,color='thistle', label='C, $\\mu^c=20$Pa.s')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('Cell Pressure (Pa)')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
# ax1.legend()
fig1.tight_layout()
fig1.savefig('PT_pressuresC.jpg')


fig1, ax1 = plt.subplots()
# ax1.plot(time,pressure_IF_all[index_ti:index_t2],linestyle='-',linewidth=3,color='navy', label='IF')
ax1.plot(time2,load_,linestyle='-',linewidth=3,color='black', label='LOAD')
ax1.plot(time2,pressure_solid_epsilon0595_,linestyle='--',linewidth=3,color='forestgreen', label='S, mono')
ax1.plot(time2,pressure_blood_epsilon0595_,linestyle='--',linewidth=3,color='forestgreen', label='B, mono')
ax1.plot(time2,pressure_solid_20_,linestyle=':',linewidth=3,color='cornflowerblue', label='S, $\\mu^c=20$Pa.s')
ax1.plot(time2,pressure_blood_1_,linestyle='-',linewidth=3,color='salmon', label='B, $\\mu^c=1$Pa.s')
ax1.plot(time2,pressure_blood_5_,linestyle='-.',linewidth=3,color='salmon', label='B, $\\mu^c=5$Pa.s')
ax1.plot(time2,pressure_blood_20_,linestyle=':',linewidth=3,color='salmon', label='B, $\\mu^c=20$Pa.s')
ax1.plot(time2,pressure_solid_1_,linestyle='-',linewidth=3,color='cornflowerblue', label='S, $\\mu^c=1$Pa.s')
ax1.plot(time2,pressure_solid_5_,linestyle='-.',linewidth=3,color='cornflowerblue', label='S, $\\mu^c=5$Pa.s')
# ax1.plot(ti2e_all,pressure_cell_all,linestyle='-.',linewidth=3,color='gray', label='C')
ax1.plot([-10, 120],[0, 0],linewidth=1,linestyle='-',color='darkgray')
ax1.plot([0, 0],[-20, 220],linewidth=1,linestyle='-',color='darkgray')
t=plt.text(20, 195, r'$Load$', fontsize = 12, color = 'black' )
t.set_bbox(dict(facecolor='white', alpha=0.7, linewidth=0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8,alpha=0.5)
plt.text(60, 105, r'$p^s~at~the~bottom~points$', fontsize = 12, color = 'cornflowerblue')
plt.text(7, 30, r'$p^b~at~the~bottom~points$', fontsize = 12, color = 'salmon')
# t=plt.text(5, 157, r'$4\%$', fontsize = 12, color = 'cornflowerblue')
# t.set_bbox(dict(facecolor='white', alpha=0.7, linewidth=0))
ax1.set_xlim([-10,120])
ax1.set_ylim([-20,220])
ax1.set_xlabel('time (s)')
ax1.set_ylabel('Pressure (Pa)')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
# ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
# ax1.legend()
fig1.tight_layout()
fig1.savefig('PT_pressures_.jpg')


fig1, ax1 = plt.subplots()
ax1.plot(time2,pressure_IF_1_,linestyle='-',linewidth=3,color='tan', label='IF, $\\mu^c=1$Pa.s')
ax1.plot(time2,pressure_IF_5_,linestyle='-.',linewidth=3,color='tan', label='IF, $\\mu^c=5$Pa.s')
ax1.plot(time2,pressure_IF_20_,linestyle=':',linewidth=3,color='tan', label='IF, $\\mu^c=20$Pa.s')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('Interstitial Fluid Pressure (Pa)')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
# ax1.legend()
ax1.plot([-10, 120],[370, 370],linewidth=1,linestyle='-',color='darkgray')
ax1.plot([0, 0],[350, 560],linewidth=1,linestyle='-',color='darkgray')
ax1.set_xlim([-10,120])
ax1.set_ylim([350,560])
fig1.tight_layout()
fig1.savefig('PT_pressuresIF_.jpg')


fig1, ax1 = plt.subplots()
ax1.plot(time2,pressure_cell_1_,linestyle='-',linewidth=3,color='thistle', label='C, $\\mu^c=1$Pa.s')
ax1.plot(time2,pressure_cell_5_,linestyle='-.',linewidth=3,color='thistle', label='C, $\\mu^c=5$Pa.s')
ax1.plot(time2,pressure_cell_20_,linestyle=':',linewidth=3,color='thistle', label='C, $\\mu^c=20$Pa.s')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('Cell Pressure (Pa)')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
# ax1.legend()
ax1.plot([-10, 120],[-1870,-1870],linewidth=1,linestyle='-',color='darkgray')
ax1.plot([0, 0],[-2000, -1600],linewidth=1,linestyle='-',color='darkgray')
ax1.set_xlim([-10,120])
ax1.set_ylim([-1875,-1675])
fig1.tight_layout()
fig1.savefig('PT_pressuresC_.jpg')

fig1, ax1 = plt.subplots()
ax1.plot(time2,varepsilon_epsilon0595,linestyle='--',linewidth=3,color='forestgreen')
ax1.plot(time2,varepsilon_1,linestyle='-',linewidth=3,color='tan', label='IF, $\\mu^c=1$Pa.s')
ax1.plot(time2,varepsilon_5,linestyle='-.',linewidth=3,color='tan', label='IF, $\\mu^c=5$Pa.s')
ax1.plot(time2,varepsilon_20,linestyle=':',linewidth=3,color='tan', label='IF, $\\mu^c=20$Pa.s')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('$\\varepsilon$ (-) at bottom points')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
# ax1.legend()
ax1.plot([-10, 120],[0.5, 0.5],linewidth=1,linestyle='-',color='darkgray')
ax1.plot([0, 0],[0, 0.6],linewidth=1,linestyle='-',color='darkgray')
ax1.set_xlim([-10,120])
ax1.set_ylim([0.485,0.507])
# plt.text(0, 0.2, r'$Extra-vascular~porosity~at~the~bottom~points$', fontsize = 12, color = 'cornflowerblue')
fig1.tight_layout()
fig1.savefig('PT_epsilon.jpg')

fig1, ax1 = plt.subplots()
ax1.plot(time2,varepsilonb_epsilon0595,linestyle='--',linewidth=3,color='forestgreen')
ax1.plot(time2,varepsilonb_1,linestyle='-',linewidth=3,color='salmon', label='B, $\\mu^c=1$Pa.s')
ax1.plot(time2,varepsilonb_5,linestyle='-.',linewidth=3,color='salmon', label='B, $\\mu^c=5$Pa.s')
ax1.plot(time2,varepsilonb_20,linestyle=':',linewidth=3,color='salmon', label='B, $\\mu^c=20$Pa.s')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('$\\varepsilon^b$ (-) at bottom points')
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.grid(color='lightgray', linestyle=':', linewidth=0.8)
# plt.text(30, 0.01, r'$Vascular~porosity~at~the~bottom~points$', fontsize = 12, color = 'salmon')
# ax1.legend()
ax1.plot([-10, 120],[0.04, 0.04],linewidth=1,linestyle='-',color='darkgray')
ax1.plot([0, 0],[0, 0.045],linewidth=1,linestyle='-',color='darkgray')
ax1.set_xlim([-10,120])
ax1.set_ylim([0.032,0.045])
fig1.tight_layout()
fig1.savefig('PT_epsilonb.jpg')

x=[None, None]
y=[None, None]
fig1d, ax1d = plt.subplots()
ax1d.plot(x,y,linestyle='--',linewidth=1.5,label='Sciumè et al. (2021)',color='forestgreen')
ax1d.plot(x,y,linestyle='-',linewidth=1.5,label='$\\mu^c=1$Pa.s',color='forestgreen')
ax1d.plot(x,y,linestyle='-.',linewidth=1.5,label='$\\mu^c=5$Pa.s',color='forestgreen')
ax1d.plot(x,y,linestyle=':',linewidth=1.5,label='$\\mu^c=20$Pa.s',color='forestgreen')
ax1d.axis('off')
fig1d.legend(loc='center', prop={'size': 12})
fig1d.tight_layout()
fig1d.savefig('PT_Figure_legend.jpg')