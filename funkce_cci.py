
try: 
    import json
    import os
    import sys

    import earthpy as et
    import earthpy.spatial as es
    import geopandas as gpd
    import imagecodecs
    import numpy as np
    import rasterio
    import tifffile as tf
    from fiona.crs import from_string, from_epsg
    from osgeo import gdal, osr
    from rasterio.mask import mask
    from rasterio.plot import show, show_hist
    from scipy import ndimage
    from shapely.geometry import box
    import pprint

    #print('Import OK.')
    
except ImportError as IE:
    print('Import not OK - ', IE )
    sys.exit('exit - bad import')


"""
input_array         --  Input TIF photo in numpy array format
src_dataset_path    --  path to input array
output_path         --  Output TIF path
"""

#https://gist.github.com/jkatagi/a1207eee32463efd06fb57676dcf86c8
def fGeoref(input_array, src_dataset_path, output_path):
        cols = input_array.shape[1]
        rows = input_array.shape[0]
        
        dataset = gdal.Open(src_dataset_path, gdal.GA_ReadOnly)
        originX, pixelWidth, b, originY, d, pixelHeight = dataset.GetGeoTransform() 
        driver = gdal.GetDriverByName('GTiff')
        band_num = 1
        GDT_dtype = gdal.GDT_Float32
        outRaster = driver.Create(output_path, cols, rows, band_num, GDT_dtype)
        outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
        outband = outRaster.GetRasterBand(band_num)
        outband.WriteArray(input_array)
        prj=dataset.GetProjection()
        outRasterSRS = osr.SpatialReference(wkt=prj)
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        
        return True

"""
input_tif       -- Input TIF photo 
output_tif      -- Output TIF photo
bbox            -- Borders of bounding box (minx, miny, maxx, maxy) in WGS84 decimal degrees
epsg_code       -- EPSG of output TIF
"""
#https://automating-gis-processes.github.io/CSC/notebooks/L5/clipping-raster.html
def fClipSatelliteImage (input_tif, output_tif, bbox, epsg_code):
    data = rasterio.open(input_tif)
    #print(data.meta)
    #Creating GeoPanda dataframe for bbox
    #old way
    #geo = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=from_epsg(4326))
    geo = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs="EPSG:4326")

    # Project the Polygon into same CRS as the grid
    #geo = geo.to_crs(crs=data.crs.data)
    geo = geo.to_crs(data.crs)
    geo.crs
    
    
    #pprint.pprint(data.crs)
    #Function to parse features from GeoDataFrame in such a manner that rasterio wants them
    def getFeatures(gdf):
        return [json.loads(gdf.to_json())['features'][0]['geometry']]
    
    #Get Geometry coordinates
    coords = getFeatures(geo)
    
    # Clip the raster with Polygon
    out_img, out_transform = mask(dataset=data, shapes=coords, crop=True)
    #Copy metadata and update
    out_meta = data.meta.copy()
    out_meta.update({"driver": "GTiff",
                     "height": out_img.shape[1],
                     "width": out_img.shape[2],
                     "transform": out_transform,
                     "crs": epsg_code}
                    )
    #Save the clipped raster
    with rasterio.open(output_tif, "w", **out_meta) as dest:
            dest.write(out_img)
    return out_img

"""
red_path        -- String path to RED band (band4 in Landsat8)
nir_path        -- String path to NIR band (band5 in Landsat8) 
out_folder      -- Output folder
"""
def fndvi(red_path, nir_path, out_folder, name = "_NDVI"):

    red= tf.imread(red_path)
    nir = tf.imread(nir_path)

    zero_exception = np.seterr(all = "ignore")
    red_array = np.array(red).astype(np.float32)
    nir_array = np.array(nir).astype(np.float32)
    
    result_ndvi = np.divide((nir_array-red_array),(red_array+nir_array))
      
    result_ndvi_path = os.path.join(out_folder, name + ".TIF") 


    fGeoref(result_ndvi, red_path, result_ndvi_path)
    
    return result_ndvi_path

"""
band_path   -- String path band
Ap          -- Constant from L8 metadata REFLECTANCE_ADD_BAND_X
Mp          -- Constant from L8 metadata REFLECTANCE_MULT_BAND_X
out_folder  -- Output folder
"""
def fTOA_Refl(band_path, Ap, Mp, out_folder, name = "_TOA_Ref"):   
    
    band= tf.imread(band_path)
    
    result_toa = Mp * band + Ap
    result_toa_path = os.path.join(out_folder, name + ".TIF")
    
    #tf.imsave(result_toa_path , result_toa)
    
    fGeoref(result_toa, band_path, result_toa_path)
    
    return result_toa_path

"""
band_path   -- String path band
ML          -- Constant from L8 metadata RADIANCE_MULT_BAND_X, constant oloaded from dict
AL          -- Constant from L8 metadata RADIANCE_ADD_BAND_X, constant oloaded from dict
out_folder  -- Output folder
"""
def fTOA_Rad(thermo_path, ML, AL, out_folder, name = "TOA_Rad"):   
        
    thermo= tf.imread(thermo_path)
    
    result_toa_rad = (ML  * thermo) + AL
    result_toa_rad_path = os.path.join(out_folder, name + ".TIF")
    
    #tf.imsave(result_toa_rad_path , result_toa_rad)
    
    fGeoref(result_toa_rad, thermo_path, result_toa_rad_path)
    
    return result_toa_rad_path

"""
input_path_toa   -- String path to TOA TIF
K1               -- Constant from L8 metadata K1_CONSTANT_BAND_10, constant oloaded from dict
K2               -- Constant from L8 metadata K2_CONSTANT_BAND_10, constant oloaded from dict
out_folder       -- Output folder
"""
#vyzarena teplota bez korekce pro emisivitu

def fBT(input_path_toa, K1, K2, out_folder, name = "BT"):
    
    #K1 = 774.8853 #Kelvin constant
    #K2 = 1321.0789 #Kelvin constant
    
    toa_array = np.array(tf.imread(input_path_toa))
    
    zero_exception = np.seterr(all = "ignore")
        
    result_BT = K2/np.log(((K1/toa_array+1))) - 273.15
    #result_BT[result_BT == np.inf] = 0    
    
    result_bt_path = os.path.join(out_folder, name + ".TIF") 
    
    fGeoref(result_BT, input_path_toa, result_bt_path)
    
    #tf.imsave(result_bt_path, result_BT)
    
    return result_bt_path

"""
ndvi        -- String path to NDVI TIF
out_folder  -- Output folder
"""
def fVc(ndvi, out_folder, name = "Vc"):
    ndvi_array = np.array(tf.imread(ndvi))
    result_Vc = np.divide(np.power(ndvi_array,2),0.3)

    result_Vc_path = os.path.join(out_folder, name +".TIF")

    #tf.imsave(result_Vc_path, result_Vc)

    fGeoref(result_Vc, ndvi, result_Vc_path)

    return result_Vc_path

"""
vc          -- String path to vegetation cover TIF
ndvi        -- String path to ndvi TIF
red         -- String path to red (band4 in L8) TIF
out_folder  -- Output folder
"""
def fEmis(vc, ndvi, red, out_folder, name = "Emis"):
        
    ndvi_array = np.array(tf.imread(ndvi))
    vc_array = np.array(tf.imread(vc))
    red = np.array(tf.imread(red))

    result_e = (0.004 * vc_array) + 0.986
    #print(result_e)
    result_e = np.where(ndvi_array < 0.2, 1 - red, result_e) 
    result_e = np.where(ndvi_array > 0.5, 0.99, result_e)              
    result_e[result_e > 1] = 0.99                           
    result_e[result_e < 0.8] = 0.8
    #print(result_e)

    result_e_path = os.path.join(out_folder, name +".TIF")

    fGeoref(result_e, vc, result_e_path)

    #tf.imsave(result_e_path, result_e)
    return result_e_path

"""
bt          -- String path to brightness temperature TIF
surf_emis   -- String path to surface emissivity TIF
out_folder  -- Output folder
"""
def fLST(bt, surf_emis, out_folder, name = "LST"):
        bt_array = np.array(tf.imread(bt))
        e_array = np.array(tf.imread(surf_emis))
        
        result_LST = (bt_array / (1 + ((0.0015 * bt_array)/1.4488) * np.log(e_array)))
        
        result_LST_path = os.path.join(out_folder, name + ".TIF")
        
        fGeoref(result_LST, bt, result_LST_path)
        
        #tf.imsave(result_LST_path, result_LST)
        return result_LST_path


"""
lst         -- String path to land surface temperature TIF
surf_emis   -- String path to surface emissivity TIF
out_folder  -- Output folder
"""
def fRadLongOut(lst, surf_emis, out_folder, name = "RadLongOut"):
    lst_array = np.array(tf.imread(lst))
    e_array = np.array(tf.imread(surf_emis))
    
    result_long_out = e_array * 5.67 * 10 ** (-8.0) * (lst_array + 273.15)**4
        
    long_out_path = os.path.join(out_folder, name +".TIF")
    
    fGeoref(result_long_out, lst, long_out_path)
    
    #tf.imsave(long_out_path, result_long_out)
    #print(result_long_out)
    return long_out_path

"""
toa_band_1  -- String path to TOA reflectance of band1 TIF
toa_band_3  -- String path to TOA reflectance of band3 TIF
toa_band_4  -- String path to TOA reflectance of band4 TIF
toa_band_5  -- String path to TOA reflectance of band5 TIF
toa_band_7  -- String path to TOA reflectance of band7 TIF
out_folder  -- Output folder
"""
def fAlbedo(toa_band_1, toa_band_3, toa_band_4, toa_band_5, toa_band_7, out_folder, name = 'albedo'):
    band1_array = np.array(tf.imread(toa_band_1))
    band3_array = np.array(tf.imread(toa_band_3))
    band4_array = np.array(tf.imread(toa_band_4))
    band5_array = np.array(tf.imread(toa_band_5))
    band7_array = np.array(tf.imread(toa_band_7))
    
    result_albedo = ((0.356 * band1_array) + (0.130 * band3_array) + (0.373 * band4_array) + (0.085 * band5_array) + (0.072 * band7_array))-0.0018/1.016
    result_albedo_path = os.path.join(out_folder, name + ".TIF")
    
    fGeoref(result_albedo, toa_band_4, result_albedo_path)
      
    #tf.imsave(result_albedo_path, result_albedo)
    
    return result_albedo_path

"""
albedo      -- String path to albedo TIF
out_folder  -- Output folder
"""
def fRadShortOut(RadShIn, albedo, out_folder, name ="RadShortOut"):
    
    albedo_array = np.array(tf.imread(albedo))

    result_RadShortOut = albedo_array * RadShIn
    
    result_short_out_path = os.path.join(out_folder, name +".TIF")
    
    fGeoref(result_RadShortOut, albedo, result_short_out_path)
    
    #tf.imsave(result_short_out_path, result_RadShortOut)
    
    return result_short_out_path

"""
emis_atmo     -- String path to Short Incoming Radiation TIF
air_temp      -- Air Temperature in [°C], constant loaded from CSV
out_folder    -- Output folder
"""
def fRadLongIn(emis_atmo, air_temp, name = "RadTLongIn"):
    
    result_long_in = emis_atmo * 5.6703 * 10.0 ** (-8.0) * (air_temp + 273.15) **4
    
    return result_long_in

"""
RadShortIn    -- String path to Short Incoming Radiation TIF
RadShortOut   -- Air Temperature in [°C], constant loaded from CSV
RadLongOut    -- String path to Long Outcoming Radiation TIF
RadLongIn     -- Constant, Long Incoming Radiation
out_folder    -- Output folder
"""
def fTotalRad(RadShortIn, RadShortOut, RadLongIn, RadLongOut, out_folder, name = "TotalRad"):
    
    RadShortOut_array = np.array(tf.imread(RadShortOut))
    RadLongOut_array = np.array(tf.imread(RadLongOut))
    
    result_TotalRad = (RadShortIn - RadShortOut_array) + (RadLongIn - RadLongOut_array)
    
    result_TotalRad_path = (os.path.join(out_folder, name + ".TIF"))
    
    fGeoref(result_TotalRad, RadLongOut, result_TotalRad_path)
    
    #tf.imsave(result_TotalRad_path, result_TotalRad)
    
    return result_TotalRad_path

"""
albedo       -- String path to Albedo TIF
ndvi         -- String path to NDVI TIF
lst          -- String path to Land Surface Temperature TIF
TotalRad     -- String path to Total Radiation (Net Radiation) TIF
out_folder   -- Output folder
"""
def fGroundHeatFl(albedo, ndvi, lst, TotalRad, out_folder, name = "GroundHeatFlux"):
    albedo_array = np.array(tf.imread(albedo))
    ndvi_array = np.array(tf.imread(ndvi))
    lst_array = np.array(tf.imread(lst))
    TotalRad_array = np.array(tf.imread(TotalRad))
    
    result_ghf_a = (0.0038 * albedo_array) + (0.0074 * albedo_array)
    result_ghf_b = (1 - 0.98) * ndvi_array
    result_ghf_c = (np.power(result_ghf_a,2) * np.power(result_ghf_b,4)) * TotalRad_array
    result_ghf_d = albedo_array * result_ghf_c
    result_ghf = np.divide(lst_array, result_ghf_d)
 
    result_ghf = (lst_array/albedo_array) * ((0.0038 * albedo_array) + 0.0074 * albedo_array)**2 * (1 - (0.98 * ndvi_array**4)) * TotalRad_array
    result_ghf_path = (os.path.join(out_folder, name + ".TIF"))
    
    fGeoref(result_ghf, ndvi, result_ghf_path)
    #tf.imsave(result_ghf_path, result_ghf)
    
    return result_ghf_path

"""
lst          -- String path to Land Surface Temperature TIF
air_temp     -- Air Temperature in [°C], constant loaded from CSV
out_folder   -- Output folder
"""
def fEvapoFraction(lst, air_temp, out_folder, name = 'EvapoFraction'):
    lst_array = np.array(tf.imread(lst))
    
    lst_median = ndimage.median_filter(lst_array, 5)
    lst_median = lst_median[~np.isnan(lst_median)]
    Tmax = np.nanmax(lst_array)
    evapo_frac = np.divide(Tmax - lst_array, Tmax - air_temp)
    
    evapo_frac_path = (os.path.join(out_folder, name + ".TIF"))
    fGeoref(evapo_frac, lst, evapo_frac_path)
    
    #tf.imsave(evapo_frac_path, evapo_frac)
    
    return evapo_frac_path

"""
TotalRad            -- String path to Total Radiation (Net Radiation) TIF
ground_heat_flux    -- String path to Ground Heat Flux (Net Radiation) TIF
out_folder          -- Output folder
"""
def fAvailableEvapo(TotalRad, ground_heat_flux, out_folder, name = 'AvailableEvapo'):
    TotalRad_array = np.array(tf.imread(TotalRad))
    ground_heat_flux_array= np.array(tf.imread(ground_heat_flux))
    
    avialable_evapo = TotalRad_array - ground_heat_flux_array #available energy for evaporation

    avialable_evapo_path = (os.path.join(out_folder,  name + ".TIF"))  
    fGeoref(avialable_evapo, TotalRad, avialable_evapo_path)
    
    #tf.imsave(avialable_evapo_path, avialable_evapo)
    
    return avialable_evapo_path

"""
TotalRad            -- String path to Total Radiation (Net Radiation) TIF
ground_heat_flux    -- String path to Ground Heat Flux (Net Radiation) TIF
out_folder          -- Output folder
"""
def fLatentHeatFlux(evapo_frac, available_evapo, out_folder, name = "LatentHeatFlux" ):    
    evapo_frac_array = np.array(tf.imread(evapo_frac))
    avialable_evapo_array = np.array(tf.imread(available_evapo))
    
    latent_heat_flux = evapo_frac_array * avialable_evapo_array
    
    latent_heat_flux_path = (os.path.join(out_folder,  name + ".TIF"))
    fGeoref(latent_heat_flux, evapo_frac, latent_heat_flux_path)
    
    #tf.imsave(latent_heat_flux_path, latent_heat_flux)
    
    return latent_heat_flux_path

"""
latent_heat_flux    -- String path to Latent Heat Flux TIF
available_evapo     -- String path to Available Evaporation (Net Radiation) TIF
out_folder          -- Output folder
"""
def fSensibHeatFl(latent_heat_flux, available_evapo, out_folder, name = "SensibleHeatFlux" ):
    latent_heat_flux_array = np.array(tf.imread(latent_heat_flux))
    available_evapo_array = np.array(tf.imread(available_evapo))

    result_sensible_heat_flux = available_evapo_array - latent_heat_flux_array

    result_sensible_heat_flux_path = (os.path.join(out_folder, name + ".TIF"))

    
    fGeoref(result_sensible_heat_flux, available_evapo, result_sensible_heat_flux_path)
   
    
    #tf.imsave(result_sensible_heat_flux_path, result_sensible_heat_flux)
    
    
    return result_sensible_heat_flux_path
  

"""
air_temp    -- Air Temperature in [°C], constant loaded from CSV
"""
def fSaturVapPress(air_temp): 
    SatVapPress =  0.6108 * 2.7183 ** ((17.27 * air_temp)/(air_temp+237.3))
    return SatVapPress

"""
AtmPress    -- atmospheric pressure [kPa], constant loaded from CSV
latent_heat -- String path to Latent Heat Flux TIF
cp          -- specific heat at constant pressure, 1.013 10-3 [MJ kg-1 °C-1]
epsilon     -- ratio molecular weight of water vapour/dry air = 0.622
"""
def fPsychCons(AtmPress, latent_heat, cp, epsilon):
    Y = (cp * AtmPress) / (epsilon * latent_heat)
    return Y


"""
air_temp    -- Air Temperature in [°C], constant loaded from CSV
"""
def fSlopeSaturVapPress(air_temp):
    SlopeSatVapPress = round((4098 * 0.6108 * 2.7183 ** 
    ((17.27 * air_temp)/(air_temp+237.3)) / (air_temp + 237.3) ** 2),3)
    return SlopeSatVapPress


"""
ActVapPress     -- http://www.fao.org/3/X0490E/x0490e07.htm#atmospheric%20pressure%20(p), mean saturation vap pres
SatVapCurve     -- Saturation Vapour by Curve, constant
WindSp          -- Wind Speed [m/s], constant
SaturVapPress   -- Saturated Vapour Pressure, constant
GroundHeatFlux  -- String Path to Ground Heat Flux TIF
psych           -- Psychometric constant, fPsychCons
TotalRad        -- String path to Total Radiation TIF
air_temp        -- air temperature [°C], from CSV
out_folder      -- Output Folder
"""
#ActVapPress -> http://www.fao.org/3/X0490E/x0490e07.htm#atmospheric%20pressure%20(p), mean saturation vap press
def fET0(ActVapPress, SatVapCurve, WindSp, SaturVapPress, GroundHeatFlux, psych, TotalRad, air_temp, out_folder, name= "ET0"):
    ground_heat_flux_array = np.array(tf.imread(GroundHeatFlux))
    total_rad_array = np.array(tf.imread(TotalRad))

    ET0 = ((0.408 * SatVapCurve * (total_rad_array - ground_heat_flux_array) + psych * (900/(air_temp + 273.15)) * WindSp * (SaturVapPress - ActVapPress))/(SatVapCurve + psych * (1+0.34*WindSp)))/10
    
    ET0_path = (os.path.join(out_folder, name + '.TIF'))     
    
    fGeoref(ET0, TotalRad, ET0_path)
    #tf.imsave(ET0_path, ET0)

    return ET0_path

"""
red_path        -- String path to RED band (band4 in Landsat8)
nir_path        -- String path to NIR band (band5 in Landsat8) 
out_folder      -- Output folder
"""
def fKc(red, nir, out_folder, name= "Kc"):
    red_array = np.array(tf.imread(red)).astype(np.float32)
    nir_array = np.array(tf.imread(nir)).astype(np.float32)

    np.seterr(all = "ignore")

    RVI_result = nir_array / red_array
    Kc_result = 1.1 * (1- np.exp(-1.5 * RVI_result))
    Kc_path = (os.path.join(out_folder, name + '.TIF'))

    fGeoref(Kc_result, nir, Kc_path)

    return Kc_path


"""
Kc_path         -- String path to crop coefficient TIF
ET0_path        -- String path to availale evapotranspiration TIF
out_folder      -- Output folder
"""
def fETI(Kc_path, ET0_path, out_folder, name = "ETI"):
    Kc_array = np.array(tf.imread(Kc_path))
    ET0_array = np.array(tf.imread(ET0_path))

    #nanmax -> ignore nan values
    ETI_result = (Kc_array * ET0_array) /  np.nanmax(ET0_array)
    ETI_path = (os.path.join(out_folder, name + ".TIF"))

    fGeoref(ETI_result, Kc_path, ETI_path)
    #tf.imsave(ETI_path, ETI_result)
    return ETI_path

"""
DEM_paht    -- String path to Digital Elevation Model in meters
azimuth     -- azimuth of the sun [default = 315°]
altitude    -- altitude of the sun above horizon [deafult = 45°]
out_folder  -- Output folder
"""
def fshade(DEM_path, azimuth, altitude, out_folder, name = "shade"):
    elevation = tf.imread(DEM_path)
    elevation_array = np.array(elevation).astype(np.float32)
    #Divided by 255 to standardize from hillshade scale (0-255) to 0-1 scale
    result_shade = es.hillshade(elevation_array, azimuth, altitude)/255
    
    result_shade_path = os.path.join(out_folder, name + ".TIF")
    
    fGeoref(result_shade, DEM_path, result_shade_path)
    
    #tf.imsave(result_shade_path, result_shade)
    return result_shade_path

"""
albedo      -- String path to albedo TIF
eti         -- String path to evapotranspiration index TIF
shade       -- String path to shade TIF
out_folder  -- Output Folder
"""
def fCCi(albedo, eti, shade, out_folder, name = 'CCi'):
    shade_array = np.array(tf.imread(shade))
    albedo_array = np.array(tf.imread(albedo))
    ETI_array = np.array(tf.imread(eti))
    
    CCi = (0.6 * shade_array) + (0.2 * albedo_array) + (0.2 * ETI_array)
    CCi_path = (os.path.join(out_folder, name + ".TIF"))
    
    fGeoref(CCi, albedo, CCi_path)

    #tf.imsave(CCi_path, CCi)
    return CCi_path


