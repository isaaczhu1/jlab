# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 11:28:06 2024

@author: evanb
"""
from readfile import read_data
import numpy as np
import matplotlib.pyplot as plt

# Extract file name
filename = "2DayRun1-marin" 



start_time = (0,12,28) # 12*3600+28*60 # days, hours, minutes
end_time = (2,8,40) # 2*24*3600+8*3600+40*60

def extract_times(clock):
    """
    

    Parameters
    ----------
    clock : string
        the clock time of the data.

    Returns
    -------
    integer
        seconds after midnight that the detection was made.

    """
    time_splits = clock.split(":")
    time = (3600*int(time_splits[0])+60*int(time_splits[1])+int(time_splits[2]))+3600*start_time[1]+60*start_time[2]
    return time % (24*3600)

def get_hourly_rate(filename, start_time, end_time):
    data = read_data(filename)
    # Get times after midnight for every single coincidence
    #times = [] 
    counts = [0]*24
    for datum in data:
        time = extract_times(datum[1])
        counts[time // 3600] += 1
        # Normalize based on how many hours we segmented over.
    counts = np.array(counts)
    num_hours = [0]*24

    for i in range(24):
        num_hours[i] = end_time[0] + 1
        if i > end_time[1]:
            num_hours[i] -= 1
        if i < start_time[1]:
            num_hours[i] -= 1
        if i == end_time[1]:
            num_hours[i] -= (60-end_time[2])/60
        if i == start_time[1]:
            num_hours[i] -= (start_time[2])/60
    
    num_hours = np.array(num_hours)
    errors = np.sqrt(counts)
    return counts, errors, num_hours
   
def normalize(counts, hours):
    for i in range(len(counts)):
        counts[i] /= hours[i]
    return counts


#start_time = (0,12,5)
#end_time = (5,9,27)
counts = np.zeros(24)
errors = np.zeros(24)
num_hours = np.zeros(24)
#counts, errors, num_hours =  get_hourly_rate(filename,start_time,end_time) 

filenames = ["2DayRun1-marin"]
for file in filenames:
    tcount, terr, thour = get_hourly_rate(file, start_time, end_time)
    counts += tcount
    errors += terr
    num_hours += thour


counts = normalize(counts, num_hours)
errors = normalize(errors, num_hours)
average = sum(counts)/24 # Average value

chi2 = 0
for i in range(24):
    chi2 += (counts[i]-average)**2/(errors[i])**2
    
print(f"chi-squared is {chi2} // 23")

plt.errorbar([i+0.5 for i in range(24)], counts, yerr=errors, xerr=0.5,linestyle='none',label="Data")
plt.plot([0,24], [average, average],label="Average")

plt.xlabel("Hours since midnight")
plt.xlim(0,24)
plt.ylabel("Count rate per hour")
plt.legend()


night = (sum(counts[0:6])+sum(counts[17:]))/(6+24-17)
night_err = np.sqrt((sum(errors[0:6]**2)+sum(errors[17:]**2)))/(6+24-17)
day = sum(counts[7:16])/(16-7)
day_err = np.sqrt(sum(errors[7:16]**2))/(16-7)

print(f"Day-time average: {day} +- {day_err} counts per hour")
print(f"Night-time average: {night} +- {night_err} counts per hour")


#plt.hist(times) # See if there is some correlation in time
    


