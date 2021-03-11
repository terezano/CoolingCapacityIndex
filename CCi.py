#!/usr/bin/env python
# -*- coding: utf-8 -*-
try: 
    import csv
    import datetime
    import json
    import os
    import pprint
    import sys

    import earthpy.spatial as es
    import geopandas as gpd
    import numpy as np
    import pandas as pd
    #import pycrs
    import rasterio
    import tifffile as tf
    from fiona.crs import from_epsg
    from rasterio.mask import mask
    from scipy import ndimage
    from shapely.geometry import box

    import funkce_cci
    print('I solemnly swear that I am up to no good.')
except:
    print('Import FAIL')

############################################################################################
############################################################################################

"""

SCRIPT_DIR      -- 
IN_FOLDER       -- Folder with L8 satellite images in TIF
OUT_FOLDER      -- Folder for result TIF
CLIPPED         -- Folder for clipped original L8 img TIF 
DEM             -- TIF with elevation in meters
df_METEODATA    -- Panda Data Frame with meteorolofical data in CSV, delimited with ','

"""
SCRIPT_DIR = os.path.dirname((os.path.abspath(__file__)))
#print(SCRIPT_DIR)
IN_FOLDER = (os.path.join(SCRIPT_DIR, r'input_L8'))
OUT_FOLDER = (os.path.join(SCRIPT_DIR, r'output_L8'))
CLIPPED = (os.path.join(SCRIPT_DIR, r'output_L8'))
DEM = (os.path.join(SCRIPT_DIR, r'input_L8\DEM_Olomouc.TIF'))


df_METEODATA = pd.read_csv((os.path.join(SCRIPT_DIR, 'KGI_libelium_T_RH.csv')), delimiter =',')
df_METEODATA["datetime_1"] = pd.to_datetime(df_METEODATA['time'])
df_METEODATA.drop(["time", "id"], 1, inplace=True)


df_RadShortIn = pd.read_csv((os.path.join(SCRIPT_DIR, 'ShortInRad.csv')), delimiter =',')

"""
EMIS_ATM          -- Emmisivity of atmosphere
"""

EMIS_ATM     = 0.8 #[unitless]

"""
BBOX                    -- Boudning box in WGS coordinates
EPSG_CODE               -- EPSG of desired coordinate system   
SATELLITE_SENSING_START -- Start of approximate Satellite sensing [time can befound in metadata]
SATELLITE_SENSING_END   -- End of Satellite Approximate sensing [time can be found in metadata]
L8_dict                 -- Dictionary {sensing_date : imgs}

"""
#Olomouc
BBOX = box(17.146144, 49.523127, 17.412798, 49.672685) #(minx, miny, maxx, maxy)
#BBOX = box(17.057542, 49.488845, 17.455584, 49.692628)
EPSG_CODE = 32633 
SATELLITE_SENSING_START = '09:30:00' #[HH:MM:SS]
SATELLITE_SENSING_END   = '09:55:00' #[HH:MM:SS]
L8_dict = {}
############################################################################################
############################################################################################

#Find all paths with valid sufix (TIF for L8)
def find_path(input_folder, file_name):
    
    path_list_folder = []
    for root, dirs, files in os.walk(input_folder, topdown=False):
        for name in files:
            #loops through all files and dirs
            if name.endswith(file_name) :
                #if suffix in name:
                    #list of paths of all files in condition
                path_list_folder.append(os.path.join(root, name))
        #pprint.pprint(path_list_folder)
    return path_list_folder

#List of strings paths to Landsat Images
landsat_tif = find_path(IN_FOLDER, ".TIF")
#pprint.pprint(landsat_tif)

#List of Landsat paths with valid name
for landsat_path in landsat_tif:
    L8_metadata = {}
    landsat_name = os.path.basename(landsat_path)
    
    #Split the valid files with '_', then make them into a key in a dictionary
    if 'LC08_L1TP' in landsat_name:
    
        #List of each part between '_' -> ['LC08', 'L1TP', '189026', '20200801', '20200807', '01', 'T1', 'B1.TIF']
        #Make keys for dictionary, consists of unique sensing dates
        date = landsat_name.split('_')[3]
        #Create dictionary with date (key) and path string to each band (value) -> filtering data by sensing date
        #print(landsat_path)
        if date in L8_dict:
            L8_dict[date].append(landsat_path)
        else:
            L8_dict[date] = [landsat_path]
#pprint.pprint(L8_dict)
        
for dates, list_of_paths in L8_dict.items():
    list_of_paths_clipped = []
    txt_path = ''
    df_air_temp_average = -99
    df_average_humidity = -99
    
    ShortInRad_value = df_RadShortIn[df_RadShortIn['date'] == int(dates)]['value']
    
    #print(ShortInRad_value)

    if ShortInRad_value.empty:
        ShortInRad_value = df_RadShortIn["value"].mean()
    else:
        ShortInRad_value = ShortInRad_value.to_numpy()[0]
    #print('datum', dates, '=> hodnota:', ShortInRad_value)

    #new data frames for humidity and temperature from CSV meteorological
    df_air_temp = df_METEODATA[df_METEODATA['sensor'] == 'TCB']
    df_humidity = df_METEODATA[df_METEODATA['sensor'] == 'HUMB']
    #print(df_air_temp)  
    #Converting SATELLITE_SENSING_START and SATELLITE_SENSING_END to datetime format
    date_time_start = datetime.datetime.strptime(f'{dates} {SATELLITE_SENSING_START}', '%Y%m%d %H:%M:%S')
    date_time_end = datetime.datetime.strptime(f'{dates} {SATELLITE_SENSING_END}', '%Y%m%d %H:%M:%S')
    
    #print(date_time_start.strftime("%Y-%m-%d %H:%M:%S %Z%z"))
    #print(date_time_end.strftime("%Y-%m-%d %H:%M:%S %Z%z"))

    #Extracting the meteorological information from CSV measured during the satellite sensing period
    df_air_temp_extract = df_air_temp[(df_air_temp['datetime_1'] > date_time_start) & (df_air_temp['datetime_1'] <= date_time_end)]
    df_humidity_extract = df_humidity[(df_humidity['datetime_1'] > date_time_start) & (df_humidity['datetime_1'] <= date_time_end)]
    #print(df_air_temp_extract)
    
    #print(df_humidity_extract)
    
    #If there is more than one number measured by meteo sensors in CSV, calculate average and round to 2 dec. places
    if len(df_air_temp_extract) > 1:
        df_air_temp_average = round((np.average(df_air_temp_extract["value"])),2)
        #print(df_air_temp_average)

    #If there is exactly one number measured by meteo sensors in CSV, use the value
    elif len(df_air_temp_extract) == 1:
        df_air_temp_average = df_air_temp_extract["value"].item()
    
    #If there is anything else, use predefined number
    else:
        df_air_temp_average = 10
    
    #If there is more than one number measured by meteo sensors in CSV, calculate average and round to 2 dec. places
    if len(df_humidity_extract) > 1:
        df_average_humidity = round((np.average(df_humidity_extract["value"])),2)
    
    #If there is exactly one number measured by meteo sensors in CSV, use the value
    elif len(df_humidity_extract) == 1:
        df_average_humidity = df_humidity_extract["value"].item()
    
    #If there is anything else, use predefined number
    else:
        df_average_humidity = 70
    
    if df_air_temp_average == -99 or df_average_humidity == -99:
        raise Exception('TEMP or HUMIDITY value error')
    

    
    #print()
    #print(df_average_humidity)
    #print(df_humidity_extract)
    #print()
    
    # read  L8 metadata text file and store values for one day in dict
    for path in list_of_paths:
        if 'B1.TIF' in path:
            txt_path = path.replace('B1.TIF', 'MTL.txt')
            #print(txt_path)
            metadata_file = open(txt_path, 'r')
            for line in metadata_file:
                #print(line)
                # last value is just END without any value
                if '=' in line:
                    key, value = line.strip().split(' = ')
                    #print(key,value)
                    L8_metadata[key] = value                                 
            metadata_file.close()
            # after find B1.TIF exit the loop
            break
          
    for path in list_of_paths:
            #print(txt_path)
        image_name = os.path.basename(path).replace('.TIF', '')
        
        #print(list_of_paths)
        out_path = CLIPPED + '\\' + image_name + '_Clipped.TIF'
        #print(out_path)
        funkce_cci.fClipSatelliteImage(path, out_path, BBOX, EPSG_CODE)
        
        list_of_paths_clipped.append(out_path)
    #pprint.pprint(list_of_paths_clipped)
    
    #Extract constants from metadata L8
    for path in list_of_paths_clipped:
        #print(path)  
        if 'B1_' in path:
            b1 = path
            reflectance_MULT_B1 = float((L8_metadata['REFLECTANCE_MULT_BAND_1']))
            reflectance_ADD_B1 = float((L8_metadata['REFLECTANCE_ADD_BAND_1']))
        if 'B2' in path:
            b2 = path
        if 'B3' in path:
            b3 = path
            reflectance_MULT_B3 = float((L8_metadata['REFLECTANCE_MULT_BAND_3']))
            reflectance_ADD_B3 = float((L8_metadata['REFLECTANCE_ADD_BAND_3']))
        if 'B4' in path:
            b4 = path
            reflectance_MULT_B4 = float((L8_metadata['REFLECTANCE_MULT_BAND_4']))
            reflectance_ADD_B4 = float((L8_metadata['REFLECTANCE_ADD_BAND_4']))
        if 'B5' in path:
            b5 = path
            reflectance_MULT_B5 = float((L8_metadata['REFLECTANCE_MULT_BAND_5']))
            reflectance_ADD_B5 = float((L8_metadata['REFLECTANCE_ADD_BAND_5']))
        if 'B6' in path:
            b6 = path
        if 'B7' in path:
            b7 = path
            reflectance_MULT_B7 = float((L8_metadata['REFLECTANCE_MULT_BAND_7']))
            reflectance_ADD_B7 = float((L8_metadata['REFLECTANCE_ADD_BAND_7']))
        if 'B10' in path:
            b10 = path
            radiance_MULT_B10 = float((L8_metadata['RADIANCE_MULT_BAND_10']))
            radiance_ADD_B10 = float((L8_metadata['RADIANCE_ADD_BAND_10']))
            K1_CONSTANT_BAND_10 = float((L8_metadata['K1_CONSTANT_BAND_10']))
            K2_CONSTANT_BAND_10 = float((L8_metadata['K2_CONSTANT_BAND_10']))
            
    #print(radiance_MULT_B10)
    #print(radiance_ADD_B10)
     
    
    ############################################################################################
    ############################################################################################
    #Computing fuctions for Cooling Capacity Index
    
    shade_path = funkce_cci.fshade(DEM, 315, 45, OUT_FOLDER, name = 'shade')

    toa_path_b10 = funkce_cci.fTOA_Rad(b10, radiance_MULT_B10, radiance_ADD_B10, OUT_FOLDER,'TOARad_b10_' + dates)  
    
    BriTemp_path = funkce_cci.fBT(toa_path_b10, K1_CONSTANT_BAND_10, K2_CONSTANT_BAND_10, OUT_FOLDER, 'BriTemp_' + dates)
    ndvi_path = funkce_cci.fndvi(b4, b5, OUT_FOLDER, 'NDVI_' + dates)
    VegCov_path = funkce_cci.fVc(ndvi_path, OUT_FOLDER, 'VegCov_' + dates)
    SurfEmis_path = funkce_cci.fEmis(VegCov_path, ndvi_path, b4, OUT_FOLDER, 'SurfEmis_' + dates)
    LST_path = funkce_cci.fLST(BriTemp_path, SurfEmis_path, OUT_FOLDER, 'LST_' + dates)
       
   
    toa_refl_b1 = funkce_cci.fTOA_Refl(b1, reflectance_ADD_B1, reflectance_MULT_B1, OUT_FOLDER, 'TOARef_b1_' + dates)
    toa_refl_b3 = funkce_cci.fTOA_Refl(b3, reflectance_ADD_B3, reflectance_MULT_B3, OUT_FOLDER, 'TOARef_b3_' + dates)
    toa_refl_b4 = funkce_cci.fTOA_Refl(b4, reflectance_ADD_B4, reflectance_MULT_B4, OUT_FOLDER, 'TOARef_b4_' + dates) 
    toa_relf_b5 = funkce_cci.fTOA_Refl(b5, reflectance_ADD_B5, reflectance_MULT_B5, OUT_FOLDER, 'TOARef_b5_' + dates) 
    toa_refl_b7 = funkce_cci.fTOA_Refl(b7, reflectance_ADD_B7, reflectance_MULT_B7, OUT_FOLDER, 'TOARef_b7_' + dates)
    albedo_path = funkce_cci.fAlbedo(toa_refl_b1, toa_refl_b3, toa_refl_b4, toa_relf_b5, toa_refl_b7, OUT_FOLDER, 'Albedo_' + dates)
    

    RadLongIn = funkce_cci.fRadLongIn(EMIS_ATM, df_air_temp_average, name = "RadTLongIn")
    RadLongOut_path = funkce_cci.fRadLongOut(LST_path, SurfEmis_path, OUT_FOLDER, 'RadLongOut_' + dates)
    
    
    RadShortOut_path = funkce_cci.fRadShortOut(ShortInRad_value, albedo_path, OUT_FOLDER, 'RadShortOut_' + dates)
    
    
    TotalRadiation_path = funkce_cci.fTotalRad(ShortInRad_value, RadShortOut_path, RadLongIn, RadLongOut_path, OUT_FOLDER, 'TotalRad_' + dates)
    GroundHeatFlux_path = funkce_cci.fGroundHeatFl(albedo_path, ndvi_path, LST_path, TotalRadiation_path, OUT_FOLDER, 'GroundHeatFlux_' + dates)
    
    AvailableEvapo_path = funkce_cci.fAvailableEvapo(TotalRadiation_path, GroundHeatFlux_path, OUT_FOLDER, 'AvailableEvapo_' + dates)
    EvapoFraction_path = funkce_cci.fEvapoFraction(LST_path, df_air_temp_average, OUT_FOLDER, 'EvapoFraction_' + dates)
    LatentHeatFlux_path = funkce_cci.fLatentHeatFlux(EvapoFraction_path, AvailableEvapo_path, OUT_FOLDER, "LatentHeatFlux" +dates)  
    
    SensibleHeatFlux_path = funkce_cci.fSensibHeatFl(LatentHeatFlux_path, AvailableEvapo_path, OUT_FOLDER, 'SensinbleHeatFlux_' + dates)
    
    SaturVapPress = funkce_cci.fSaturVapPress(df_air_temp_average)
    #print(SaturVapPress)
    
    PsychCons = funkce_cci.fPsychCons(101.325, 2.45, 0.001013, 0.622)
    #print(PsychCons)
    
    SlopeVapPress = funkce_cci.fSlopeSaturVapPress(df_air_temp_average)
    #print(SlopeVapPress)

    ET0_path = funkce_cci.fET0(1.6, SlopeVapPress, 4.5, SaturVapPress, GroundHeatFlux_path, PsychCons, TotalRadiation_path, df_air_temp_average, OUT_FOLDER, 'ET0_'+dates)  
    
    Kc_path = funkce_cci.fKc(b4, b5, OUT_FOLDER,'Kc_' + dates)
 
    ETI_path = funkce_cci.fETI(Kc_path, ET0_path, OUT_FOLDER, 'ETI_' + dates)
    
    CCi_path = funkce_cci.fCCi(albedo_path, ETI_path, shade_path, OUT_FOLDER, 'CCI_' + dates )  
    
pprint.pprint("Mischief managed.")

