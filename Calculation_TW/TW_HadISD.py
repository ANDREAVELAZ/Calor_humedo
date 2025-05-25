import requests
import xarray as xr
import numpy as np
import metpy
import metpy.calc as mpycalc
from metpy.units import units
import matplotlib.pyplot as plt 
from wetbulb_dj08_spedup import WetBulb
import os
import time
import pandas as pd

#Función para seleccionar las variables
def selec_variables(filename,ii):
    datos=xr.open_dataset(filename)
    if datos.stnlp.isnull().all()==False: #Para encontrar si todos los valores son nulos
        datos=datos[["tds","stnlp","tas"]]  
    else:
        datos=datos[["tds","psl","tas"]]
    name=(f"{ii}99999.nc") 
    datos.to_netcdf(name)
    
    
#Función para limpieza de los datos
def limpieza(ii):
    name=(f"{ii}99999.nc")
    data=xr.open_dataset(name)
    va=list(data.variables.keys()) #Lista de las variables 
   
    #Quitando datos erróneos
    dew_point=data[va[0]].where(data[va[0]]!=data[va[0]].flagged_value)
    presion=data[va[1]].where(data[va[1]]!=data[va[1]].flagged_value)
    temp=data[va[2]].where(data[va[2]]!=data[va[2]].flagged_value)
    
    #Quitando datos faltantes
    dew_point=dew_point.dropna("time")
    presion=presion.dropna("time")
    temp=temp.dropna("time")
    
    #Encontrando datos similares
    dew_point=dew_point.where((dew_point.time==presion.time)&(dew_point.time==temp.time),drop=True)
    presion=presion.where((dew_point.time==presion.time)&(dew_point.time==temp.time),drop=True)
    temp=temp.where((dew_point.time==presion.time)&(dew_point.time==temp.time),drop=True)

    #Condición en caso de que no existan datos similares

    if dew_point.isnull().all()==True:
        pass
    else:

        #Seleccionando años completos y con más de dos mediciones al día
        anios=[dew_point.time.dt.year.min()]
        a=dew_point.time.dt.year.min()
        for i in range(0,(np.array(dew_point.time.dt.year.max())+1)-(np.array(dew_point.time.dt.year.min())+1)):
            a=a+1
            anios.append(a)
        anio_med=[]
        for i in anios:
            time1=dew_point.where(dew_point.time.dt.year==i,drop=True)
            if len(time1)>730:
                anio_med.append(i)
            
        #Cortando los datos con los años completos
        dew_point=dew_point.where(dew_point.time.dt.year.isin(anio_med), drop=True)
        presion=presion.where(presion.time.dt.year.isin(anio_med), drop=True)
        temp=temp.where(temp.time.dt.year.isin(anio_med), drop=True)
        
    return dew_point,temp,data,name

def TW(dew_point,datos):
    time_vec=dew_point.time
    HS=np.zeros(time_vec.shape)
    Twb=np.zeros(time_vec.shape)
    va=list(datos.variables.keys())
    for it,ti in enumerate(time_vec):
       
        slice_data=datos.sel(time=ti) #Recortando los datos
        
        # convertimos usando metpy
        hs=mpycalc.specific_humidity_from_dewpoint(slice_data[va[1]],slice_data[va[0]])
        SpecificHumidity=float(hs.values)
        HS[it]=SpecificHumidity #Para que guarde la humedad específica para cada punto
        
        
        
        Temperature=float(slice_data[va[2]])
        
        Pressure=float(slice_data[va[1]].values)*100
        #utilizando el codigo en wetbulb_dj08_spedup
        Twb[it]=WetBulb(Temperature,Pressure,SpecificHumidity,0)[0]+273.15
        
    return Twb,HS,time_vec

#Función que hace los cálculos de TWmax,TWmean,Tmax
def calculos (Tw,Hs,time_vec,temp):
    arr = xr.DataArray(Tw, dims=("time"),   #Haciando un xarray.DataArray para poder encontrar de forma más rápida 
                                                                    #los datos de TWmax,TWmean,Tmax usando resample
                   coords={"time": time_vec},
                   attrs={"units": "K"})
    arr2 = xr.DataArray(Hs, dims=("time"),   #Haciando un xarray.DataArray para poder encontrar de forma más rápida 
                                                                    #los datos de TWmax,TWmean,Tmax usando resample
                   coords={"time": time_vec},
                   attrs={"units": "kg/kg"})
    #Temperatura min de bulbo húmedo
    Twmin=arr.resample(time="D").min() 
    Twmin.attrs=attrs={"units": "K", "description": "Temperatura mínima diaria de bulbo húmedo"} #Cambiando los atributos 
    Twmin.name="TWmin"
    #Temperatura max de bulbo húmedo
    Twmax=arr.resample(time="D").max() 
    Twmax.attrs=attrs={"units": "K", "description": "Temperatura máxima diaria de bulbo húmedo"} #Cambiando los atributos 
    Twmax.name="TWmax" 
    #Temperatura promedio de bulbo húmedo
    Twmean=arr.resample(time="D").mean()
    Twmean.attrs=attrs={"units": "K", "description": "Temperatura promedio diaria de bulbo húmedo promedio"}
    Twmean.name="TWmean"  
     #Temperatura min 
    Tmin=temp.resample(time="D").min() 
    Tmin.attrs=attrs={"units": "K", "description": "Temperatura mínima diaria"} #Cambiando los atributos 
    Tmin.name="Tmin"
    #Temperatura máxima
    Tmax=temp.resample(time="D").max()
    Tmax.attrs=attrs={"units": "K", "description": "Temperatura máxima diaria"}
    Tmax.name="Tmax"
    #Temperatura promedio
    Tmean=temp.resample(time="D").mean()
    Tmean.attrs=attrs={"units": "K", "description": "Temperatura promedio diaria"}
    Tmean.name="Tmean"  
     #Temperatura min de bulbo húmedo
    SHmin=arr2.resample(time="D").min() 
    SHmin.attrs=attrs={"units": "kg/kg", "description": "Humedad específica mínima diaria"} #Cambiando los atributos 
    SHmin.name="SHmin"
    #Humedad específica máxima
    SHmax=arr2.resample(time="D").max()
    SHmax.attrs=attrs={"units": "kg/kg", "description": "Humedad específica máxima diaria"}
    SHmax.name="SHmax"
   
    #Humedad específica promedio
    SHmean=arr2.resample(time="D").mean()
    SHmean.attrs=attrs={"units": "kg/kg", "description": "Humedad específica promedio diaria"}
    SHmean.name="SHmean" 
    return Twmin,Twmax,Twmean,Tmin,Tmax,Tmean,SHmin,SHmax,SHmean,arr,arr2


def archivo (arr,arr2,temp,Twmin,Twmax,Twmean,Tmin,Tmax,Tmean,SHmin,SHmax,SHmean,time_vec,name,ii):
    #Creando un nuevo archivo NetCDF 
# Crear un xr.Dataset
    esta=xr.Dataset(data_vars={"Twmin":Twmin,"Twmax":Twmax,"Twmean":Twmean,"Tmin":Tmin,"Tmax":Tmax,"Tmean":Tmean,"SHmin":SHmin,"SHmax":SHmax,"SHmean":SHmean},
                 #Para cambiar el formato de la fechas
                         
                   coords={"time": np.array(Twmax.time).astype("<M8[ns]")},
                   attrs={"description": "Estadísticos diarios:Twmax,Twmean,Tmax,Tmean,SHmax,SHmean"})
    newdata = xr.Dataset(data_vars={"Tw":arr,"T":temp,"SH":arr2},
                 #Para cambiar el formato de la fechad
                   coords={"time": np.array(arr.time).astype("<M8[ns]")},
                   attrs={"description": "Variables subdiarias"})
    newdata.to_netcdf(name)
    name1=(f"{ii}99999_estadisticos.nc")
    esta.to_netcdf(name1)

#Leyendo el archivo de texto con la identificación de las estaciones
df=pd.read_fwf("idstation.txt",header=None) 

#Recortando la lista a sólo los ID de las estaciones en México (Tijuana-Tapachula=107 estaciones)
ID=df[(df[0]>="760011-99999")&(df[0]<="769043-99999")] #el número 76 es el inicio de las estaciones dentro de México

ID=ID[0].str.split(pat="-", expand=True) #Separando los datos


#Descarga y cálculo de TW
txt= open("ID.txt","w")
txt.close()
for i,ii in enumerate(ID[0]):
    url = 'https://data.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadISD/subdaily/HadISDTable/r1/v3-4-0-2023f/station_data'
    filename=(f'{ii}99999_HadISD_HadOBS_19310101-20240101_v3-4-0-2023f.nc')
    #print(filename)
    os.system('wget '+url+'/'+filename)

    #Seleccionando variables a ocupar
    selec_variables(filename,ii)

    #Limpieza de datos
    dew_point,temp,datos,name=limpieza(ii)

    #Condición para el caso de que la estación no coincida en ninguna fecha, cambie a la siguiente estación
    if dew_point.isnull().all()==True:
        continue
       
    else:
         #Cálculo de Tw y humedad específica (Hs)
        Tw,Hs,time_vec=TW(dew_point,datos=datos)
    
        #Cálculo de estadísticas diarias
        Twmin,Twmax,Twmean,Tmin,Tmax,Tmean,SHmin,SHmax,SHmean,arr,arr2=calculos(Tw,Hs,time_vec,temp)
    
        #Guardando los datos
        archivo(arr,arr2,temp,Twmin,Twmax,Twmean,Tmin,Tmax,Tmean,SHmin,SHmax,SHmean,time_vec,name,ii)
        
    