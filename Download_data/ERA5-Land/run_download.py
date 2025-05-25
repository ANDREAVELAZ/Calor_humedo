import os
import sys,glob
import numpy as np 
month_list=['01','02','03','04','05','06','07','08','09','10','11','12']
start_year=1979
end_year=2023
for year in np.arange(start_year,end_year+1):
	for month in month_list:
		os.system('python download_ERA5-Land.py '+str(year)+' '+month)


