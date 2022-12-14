import yt
import fnmatch
import os
import pandas as pd

# Select the time
istart = 0 # Start index
iskip  = 1 # Skip index

path = os.getcwd()

list_output=fnmatch.filter(os.listdir(path), 'plt*')
list_output.sort()
print("Found output   : {}".format(list_output))

list_process=list_output[istart:len(list_output):iskip]
# list_process = ["plt01119_scaldis"]
print("To be processed: {} \n".format(list_process))

d = {}
d['time'] = []
d['temperature'] = []
d['rho_O'] = []
d['rho_adv_0'] = []
d['adv_0'] = []
for file in list_process:
    ds = yt.load("{}/{}".format(path,file))
    pt = ds.point([0.0,0.0,0.0])

    d['time'].append(float(ds.current_time))
    d['temperature'].append(float(pt['Temp']))
    d['rho_adv_0'].append(float(pt['rho_adv_0']))
    d['adv_0'].append(float(pt['adv_0']))
    d['rho_O'].append(float(pt['Y(O)']*pt['density']))


df = pd.DataFrame(d)
df = df.sort_values(by='time', ascending=True)
df.to_csv('data_center_discharge.csv')
