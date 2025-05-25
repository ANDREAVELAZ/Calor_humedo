#Importando paqueterías
import requests
import numpy as np
import os
import pandas as pd


#Cargando el id de las estaciones HadISD
df=pd.read_fwf("idstation.txt",header=None) #Leyendo el archivo de texto
#df

#Selección de id de las estaciones en México
#Las estaciones de México abarcan desde tijuana hasta Tapachula, el número 76 al inicio indica la zona de México.
ID=df[(df[0]>="760011-99999")&(df[0]<="769043-99999")] 


ID=ID[0].str.split(pat="-", expand=True) #Separando los datos

#Descargando los archivos pertenecientes a los datos HadISD de México
for i,ii in enumerate(ID[0]):
    url = 'https://data.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadISD/subdaily/HadISDTable/r1/v3-4-0-2023f/station_data'
    filename=(f'{ii}99999_HadISD_HadOBS_19310101-20240101_v3-4-0-2023f.nc')
    #print(filename)
    os.system('wget '+url+'/'+filename)