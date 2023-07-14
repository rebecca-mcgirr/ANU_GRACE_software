#Import required programs, load Data file and variables
import sys
import numpy as np
import time
from time import gmtime
from datetime import datetime, timedelta
print('Y     M  D KBR R  Start End')
yr = str(sys.argv[1])
month = str(sys.argv[2])
day = str(sys.argv[3])

f = open('KBR1B_{yr}-{month}-{day}_X_03.asc'.format(yr = yr, month = month, day = day), 'r')
data = f.readlines()

split = np.array_split(data,len(data))

KBR = np.zeros(len(data)-48, dtype='int')
i = 48
DATA_START = data[12]

split_start = np.array_split(data,len(data))
string_start = str(split_start[11])
res_start = int(''.join(map(str, string_start[35:44])))

x_start = str(split_start[11])

x = (''.join(map(str, x_start[64:72])))

a = time.strptime(x, "%H:%M:%S")

DIFF = timedelta(hours=a.tm_hour, minutes=a.tm_min, seconds=a.tm_sec).seconds


while 48 <= i < len(data):
	split = np.array_split(data,len(data))
	string = str(split[i])
	res = int(''.join(map(str, string[3:13])))
	KBR[i-48] = int(res)
	i = i + 1
f.close()



def LastDigit(num):
    return num % 10

def find_missing(lst):
        return [i for x, y in zip(lst, lst[1:]) 
            for i in range(x + 1, y) if y - x >= 1]
    
array_int = np.array(KBR, dtype='int')
lst = array_int
missing = np.array(find_missing(lst), dtype='int')
missed_measurements = np.zeros((int(len(missing))))
for i in range(len(missing)):   
            number = missing[i]
            if LastDigit(number) == 0:
                missed_measurements[i] = number
            elif LastDigit(number) == 5:
                missed_measurements[i] = number

res = missed_measurements[missed_measurements != 0]
missed_val_int = np.array(res, dtype='int')


def consecutive(data, stepsize=5):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

DECIMATE_SECONDS = np.zeros((len(consecutive(missed_val_int)[0:]),2))
for i in range(len(KBR)):
    if i < len(consecutive(missed_val_int)[0:]):
        minimum = np.min(consecutive(missed_val_int)[i])
        maximum = np.max(consecutive(missed_val_int)[i])
        DECIMATE_SECONDS[i,0] = minimum
        DECIMATE_SECONDS[i,1] = maximum
        
        
        
        
#GRACE time start of 2000-01-01 12:00:00:00
grace_start_time = 946728000

###### FIRST TIME OF THE DAY
grace_time = res_start
#print(grace_time)

grace_date = time.gmtime(grace_time + grace_start_time)


if a != '00:00:00':
	print('{yr} {month} {day} KBR RA {start} {end}'.format(yr = yr, month = month, day = day, start = int(1), end = int(DIFF/5)))
	print('{yr} {month} {day} KBR RR {start} {end}'.format(yr = yr, month = month, day = day, start = int(1), end = int(DIFF/5)))

for i in range(len(DECIMATE_SECONDS)):
    yr = int(grace_date.tm_year)
    month = int(grace_date.tm_mon)
    day = int(grace_date.tm_mday)


    ###### START AND END TIME #####
    start_sec = DECIMATE_SECONDS[i,0]
    end_sec = DECIMATE_SECONDS[i,1]


    epoch_start = ((start_sec - grace_time))/5
    epoch_end = ((end_sec - grace_time))/5
   
    
    print('{yr} {month} {day} KBR RA {start} {end}'.format(yr = yr, month = month, day = day, start = int(epoch_start), end = int(epoch_end)))
    print('{yr} {month} {day} KBR RR {start} {end}'.format(yr = yr, month = month, day = day, start = int(epoch_start), end = int(epoch_end)))

