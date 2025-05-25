import xarray as xr
import metpy.calc as mpycalc
from wetbulb_dj08_spedup import WetBulb  # Asegúrate de que este es tu módulo correcto
import numpy as np
from metpy.units import units
import multiprocessing as mp
import time
import matplotlib.pyplot as plt
import os 
cpu_affinity=os.sched_getaffinity(0)

#print(cpu_affinity)

#print(len(cpu_affinity))
# Cargar los datos y la máscara de tierra-mar
anio=1979
while anio<=2023:
    from metpy.units import units
    print(anio)
    filename = f"{anio}*.nc"
    datos = xr.open_mfdataset(filename)
    # Definir los límites de latitud y longitud para una porción de México
    lat_min = 12.0  # Latitud mínima (ej. sur de México)
    lat_max = 35.0  # Latitud máxima (ej. norte de México)
    lon_min = -118.0  # Longitud mínima (ej. oeste de México)
    lon_max = -86.0  # Longitud máxima (ej. este de México)
   # print(datos)
    # Hacer un subset de los datos en esa región geográfica
    datos = datos.sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))
    datos=datos.rename({'valid_time':'time'})
    mask_filename = 'lsm_1279l4_0.1x0.1.grb_v4_unpack.nc'
    mask = xr.open_dataset(mask_filename)
    mask.coords['longitude'] = (mask.coords['longitude'] + 180) % 360 - 180
    mask = mask.sortby('longitude')
    lsm = mask['lsm'].squeeze()
    
    # Interpolación de la máscara al grid de los datos
    lsm_interp = lsm.interp(longitude=datos.longitude, latitude=datos.latitude)
    
    # Filtrar las coordenadas de latitud y longitud válidas (donde la máscara es 1)
    valid_coords = np.column_stack(np.where(lsm_interp != 0))
    latitudes = datos.latitude.values[valid_coords[:, 0]]
    longitudes = datos.longitude.values[valid_coords[:, 1]]
    
    # Crear una lista de pares (latitud, longitud) con buena data
    lat_lon_pairs = list(zip(latitudes, longitudes))
   # print('grid points to analyse',len(lat_lon_pairs))
    
    # Convertir los datos a las unidades adecuadas
 #   import pint

# Crear un objeto UnitRegistry
   # ureg = pint.UnitRegistry()
    p_with_units = datos.sp * units.Pa
    print(p_with_units) 
    dewpoint_with_units = datos.d2m * units.kelvin
    temperature_c = datos.t2m - 273.15  # Convertir a Celsius
    temperature_c.load()
    p_with_units.load()
    def save_output(da, name, units, path,long_name):
        da.name=name
        da.attrs={'units':units,'long_name':long_name}
        return da
    #    da.to_netcdf(path+name+'.nc')
    #    return
    
    def create_array(data, lats, lons, time_vec, units, operation, name):
        arr = xr.DataArray(
            data, dims=("time", "latitude", "longitude"),
            coords={"time": time_vec, "latitude": lats, "longitude": lons},
            name=name, attrs={"units": units}
        )
        arr = arr.where(arr != 0)
        if operation == 'daily':
            return arr.resample(time='1D').mean(), arr.resample(time='1D').max()
        else:
            return arr
    
    qs = mpycalc.specific_humidity_from_dewpoint(p_with_units, dewpoint_with_units).values
    HS = create_array(qs, datos.latitude, datos.longitude, datos.time, 'kg kg-1', None, 'specific_humidity')
    
    # Función que calcula la temperatura de bulbo húmedo para un punto (lat, lon) en todos los tiempos
    def calculate_wet_bulb_for_point(lat, lon):
    #    print(lat,lon)
        timei=time.time()
    #    try:
        temp = np.asarray(temperature_c.sel(latitude=lat, longitude=lon).values)
    
        pressure = np.asarray(p_with_units.sel(latitude=lat, longitude=lon).values)
    
        humidity = np.asarray(HS.sel(latitude=lat,longitude=lon).values)
    
           # mpycalc.specific_humidity_from_dewpoint(
           #     pressure, 
           #     dewpoint_with_units.sel(latitude=lat, longitude=lon).values
           # )
            
            # Asegurarse de que no haya valores inválidos
        if np.isnan(temp).any() or np.isnan(pressure).any() or np.isnan(humidity).any():
            return None
        lenh=len(humidity)
        wet_bulb_temp=np.zeros(lenh) 
        for ij in range(lenh):
        # Calcular la temperatura de bulbo húmedo en todos los tiempos
            wet_bulb_temp[ij] = WetBulb(temp[ij], pressure[ij], humidity[ij], 0)[0]
    #    print(lon,lat,time.time()-timei,np.max(wet_bulb_temp))
        return lat, lon, wet_bulb_temp  # Devolver lat, lon y los valores calculados
    #except ZeroDivisionError:
    #    print(f"División por cero en el punto ({lat}, {lon})")
    #    return None
    def resample_max_mean(arr):
        return arr.resample(time='1D').mean(), arr.resample(time='1D').max()
    # Función para paralelizar el cálculo usando multiprocessing
    def parallel_wet_bulb_calculation(lat_lon_pairs):
        with mp.Pool(processes=10) as pool:
            results = pool.starmap(calculate_wet_bulb_for_point, lat_lon_pairs)
        return results
    #test code
    #print(mp.cpu_count())
    #for ipp,pair in enumerate(lat_lon_pairs):
    #    print(calculate_wet_bulb_for_point(pair[0],pair[1]))
    #    if ipp>20:
    #        break
    time_vec = datos.time.values
    n_time = len(time_vec)
    # Ejecutar la paralelización
    results = parallel_wet_bulb_calculation(lat_lon_pairs)#[0:8000])
    twb_full_array = np.full((n_time, len(datos.latitude), len(datos.longitude)), np.nan)
    
    # Colocar los resultados en las posiciones correctas
    for result in results:
        if result is not None:
            lat, lon, twb = result
    
            # Encontrar el índice de la latitud y longitud en 'datos'
            lat_idx = np.where(datos.latitude == lat)[0][0]
            lon_idx = np.where(datos.longitude == lon)[0][0]
            
            # Asignar los valores de wet bulb temperature (twb) al array correspondiente
            twb_full_array[:, lat_idx, lon_idx] = twb
    #print(np.nanmax(twb_full_array))
    # Convertir la matriz final en un DataArray usando los mismos ejes que 'datos'
    twb = xr.DataArray(
        twb_full_array,
        dims=("time", "latitude", "longitude"),
            coords={"time": time_vec, "latitude": datos.latitude.values, "longitude": datos.longitude.values},
        name="wet_bulb_temperature",
        attrs={"units": "Celsius"}
    )
   # print(twb.max())
    
    #twb_mean,twb_max=resample_max_mean(twb)
    #T_mean,T_max=resample_max_mean(temperature_c)
    #qs_mean,qs_max=resample_max_mean(HS)
   # print(qs_mean)
    array_list=[twb,temperature_c,HS]
    names=['TWB_d','T_d','q_d',]
    long_names=['Wet bulb temperature daily ','Dry bulb temperature daily ','Daily mean specific humidity']
    units=['degrees Celsius','degrees Celsius','g kg-1']
    da_list={}
    for ii,array in enumerate(array_list):
            da=save_output(array, names[ii],units[ii],'',long_names[ii])
            da_list[names[ii]]=da
    #def save_output(da, name, units, path,long_name):
   # if anio==1979:
    #    ds=xr.Dataset(da_list)
        #print(ds)
       # quit()
    #else:
    nuevos_datos=xr.Dataset(da_list)
   #     newdata_actualizado = xr.concat([ds, nuevos_datos], dim="time")
    #    ds=newdata_actualizado
    nuevos_datos.to_netcdf(f"ERA_{anio}_hourly.nc") 
   # print(ds)

    #if anio==2023:
     #   ds.to_netcdf("analysis_ERA2.nc")
 #   el#se:
  #      pass        
    anio=anio+1




    # Crear el DataArray final usando la función create_array
   # twbmean,twbmax = create_array(wet_bulb_values_array, np.array(latitudes), np.array(longitudes), time_vec, "Celsius", "daily","wet_bulb_temperature")
   # Tmean,Tmax = create_array(wet_bulb_values_array, np.array(latitudes), np.array(longitudes), time_vec, "Celsius", "daily","wet_bulb_temperature")
   # shmean,shmax = create_array(wet_bulb_values_array, np.array(latitudes), np.array(longitudes), time_vec, "Celsius", "daily","wet_bulb_temperature")
   # print(twb)
   # newdata = xr.Dataset(da_list)#ta_vars={"Twmax":TWMAX,"Twmean":TWMEAN,"Tmax":TMAX,"Tmean":TMEAN,"SHmax":SHMAX,"SHmean":SHMEAN},
                     #Para cambiar el formato de la fecha
    #                   coords={"time": np.array(time_vec).astype("<M8[ns]"), "latitude": TWMAX.latitude, "longitude": TWMAX.longitude},
    #                   attrs={"description": "Valores calculados Twmax,Twmean,Tmax,Tmean, SHmax,SHmean"})
  #  print(newdata)


   # anio=anio+1
ds = xr.open_mfdataset("ERA_hourly*.nc", combine="by_coords")
ds.to_netcdf("ERA_hourly.nc")
