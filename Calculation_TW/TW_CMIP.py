import xarray as xr
import numpy as np
import metpy
import metpy.calc as mpycalc
from metpy.units import units
import matplotlib.pyplot as plt
from wetbulb_dj08_spedup import WetBulb
import os
import cftime
import pandas as pd

#Este script se corrió dos veces, una vez por modelo

#Función para unir todos los archivos en uno sólo y así, posteriormente, calcular TW.

def model(escenario,periodo_i):
 df=pd.read_excel("names.xlsx")
 var=["huss_day","tas_day","psl_day"]

 df=df.loc[(df["scenario"]==escenario) & (df["period"]>periodo_i)]
 huss=df.loc[df["var"]=="huss_day"].reset_index()
 psl=df.loc[df["var"]=="psl_day"].reset_index()
 tas=df.loc[df["var"]=="tas_day"].reset_index()

#Datos concatenados
 for ij,jj in enumerate(huss.index):
    print(ij)
    
    #Abriendo los archivos
    data_huss=xr.open_dataset(f"{huss.iloc[ij,1]}_{huss.iloc[ij,2]}_{huss.iloc[ij,3]}_{huss.iloc[ij,4]}_{huss.iloc[ij,5]}.nc")
    data_psl=xr.open_dataset(f"{psl.iloc[ij,1]}_{psl.iloc[ij,2]}_{psl.iloc[ij,3]}_{psl.iloc[ij,4]}_{psl.iloc[ij,5]}.nc")
    data_tas=xr.open_dataset(f"{tas.iloc[ij,1]}_{tas.iloc[ij,2]}_{tas.iloc[ij,3]}_{tas.iloc[ij,4]}_{tas.iloc[ij,5]}.nc")

    #Recortando los datos
    data_huss=data_huss["huss"].sel(lat=slice(5,37), lon=slice(360-125,360-75))
    data_psl=data_psl["psl"].sel(lat=slice(5,37), lon=slice(360-125,360-75))
    data_tas=data_tas["tas"].sel(lat=slice(5,37), lon=slice(360-125,360-75))

    # Crear un xr.Dataset
    if ij==0:
        newdata = xr.Dataset(data_vars={"huss":data_huss,"psl":data_psl,"tas":data_tas},
                     #Para cambiar el formato de la fecha
                       coords={"time": data_huss.time, "latitude": data_huss.lat, "longitude": data_huss.lon},
                       attrs={"description": f"Datos escenario {escenario}"})
    else:
        newdata1 = xr.Dataset(data_vars={"huss":data_huss,"psl":data_psl,"tas":data_tas},
                     #Para cambiar el formato de la fecha
                       coords={"time": np.array(data_huss.time), "latitude": data_huss.lat, "longitude": data_huss.lon},
                       attrs={"description": f"Datos escenario {escenario}"})
        newdata_actualizado = xr.concat([newdata, newdata1], dim="time")
        newdata=newdata_actualizado
        
    if escenario=="historical":

    newdata_actualizado.to_netcdf("historic_HADGEM.nc") #Nombre del archivo concatenado 

    else:
    newdata_actualizado.to_netcdf("ssp5852_HADGEM.nc")

    #Corriendo la función para cada escenario
model("historical","1949")
model("ssp585","2013")

#Función para calcular TW

def TW(lats,lons,time_vec,datos):
    HS=np.zeros(datos["tas"].shape)
    Twb=np.zeros(datos["tas"].shape) #Haciendo un array para rellenar los datos, dado que la rutina dos es más rápida, se ocupará ese método para calcular TW
    
    for ilat,lat in enumerate(lats):
        #print(ilat)
        for ilon,lon in enumerate(lons):
            #print(ilon)
            slice_data=datos.isel(lat=ilat,lon=ilon)
    		# checamos si esta coordenada tiene datos o no, si son nulos, no calculamos nada
            if all(slice_data.psl.isnull()):
                continue
            else:
                hs= slice_data.huss.values
 
                for it,ti in enumerate(time_vec):
                    time_data=slice_data.isel(time=it)
                    SpecificHumidity=hs[it]
                    Temperature=time_data.tas.values-273.15
                   # print(Temperature)
                    Pressure=time_data.psl.values
    				#print(Temperature,Humidity,Pressure)
                    Twb[it,ilat,ilon]=WetBulb(Temperature,Pressure,SpecificHumidity,0)[0]+273.15
                    
    return Twb,hs,Temperature

datos=xr.open_dataset("historic_HADGEM.nc")
datos1=xr.open_dataset("ssp5852_HADGEM.nc")

lats=datos.huss.lat
lons=datos.huss.lon
lats1=datos1.huss.lat
lons1=datos1.huss.lon

time_vec=datos.huss.time
time_vec1=datos1.huss.time

TWhistorical,SHhistorical,Temperaturehistorical=TW(lats,lons,time_vec,datos)
arr = xr.DataArray(TWhistorical, dims=("time", "lat","lon"),   #Haciando un xarray.DataArray para poder encontrar de forma más rápida 
                                                                    #los datos de TWmax,TWmean,Tmax usando resample
                   coords={"time": time_vec, "lat": lats,
                "lon": lons},
                   attrs={"units": "K"})
newdata=xr.Dataset(data_vars={"Tw":arr,"T":datos.tas,"SH":datos.huss},
                 #Para cambiar el formato de la fecha
                   coords={"time": time_vec, "lat": lats, "lon": lons},
                   attrs={"description": "Valores calculados historical"})
TWssp585,SHssp585,Temperaturessp585=TW(lats1,lons1,time_vec1,datos1)
arr1 = xr.DataArray(TWssp585, dims=("time", "lat","lon"),   #Haciando un xarray.DataArray para poder encontrar de forma más rápida 
                                                                    #los datos de TWmax,TWmean,Tmax usando resample
                   coords={"time": time_vec1, "lat": lats1,
                "lon": lons1},
                   attrs={"units": "K"})
newdata1 = xr.Dataset(data_vars={"Tw":arr1,"T":datos1.tas,"SH":datos1.huss},
                 #Para cambiar el formato de la fecha
                   coords={"time":time_vec1, "lat":lats1, "lon":lons1},
                   attrs={"description": "Valores calculados ssp585"})


#Guardándolo como archivo tipo netCDF
newdata.to_netcdf("cal_historical_HADGEM.nc") #Nombre del archivo con los cálculos de TW
newdata1.to_netcdf("cal_ssp585_HADGEM.nc")
 

