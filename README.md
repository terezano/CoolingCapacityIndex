# CoolingCapacityIndex

A short script generating Cooling Capacity Index made for Landsat 8 images as main part of master thesis. Based on InVEST <a href = https://naturalcapitalproject.stanford.edu/software/invest> Natural Capital project by Standford University<a/>.

# Specifications

- Python 3 (originally made for Python 3.8.8 using Miniconda3)
- Input folder with subfolders containing unzipped images with all bands
- Output folder
- Digital Elevation model in meters for your region of interest
- Sensor data for air temperature and humidity in CSV
- Short Wave Incoming Radiation in CSV


# Obtaining satellite images
Landsat 8 images can be obtained from <a href = https://earthexplorer.usgs.gov>EarthExplorer by USGS </a>, where you need to register and then download image bundle, unzip it and locate in input folder:

# Sensor Data in CSV
Your sensor data should look a bit like this: 

![image](https://user-images.githubusercontent.com/60270092/112305646-9c060780-8ca7-11eb-8b80-034a91609ab9.png)

- Sensor -- name of the sensor
- Value  -- actual value of temperature in °C/humidity in %
- Time   -- Timestamp in RRRR-MM-DD HH:MM:SS
- Id     -- Id of each value (can be ommited)

# Short Wave Incoming Radiation CSV
![image](https://user-images.githubusercontent.com/60270092/112306032-074fd980-8ca8-11eb-8b13-93734ad4085c.png)

- Date   -- Date in format of RRRRMMDD
- Value  -- actual value of radiation in watts per square meter


# Running with Miniconda3
All kinds of env activation with miniconda are available <a href= https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>here</a>.
If your conda is connected to another env you should deactivate it first using:

```
conda deactivate
```

and then create new env using conda create env from file:

```
conda env create -f environment.yml
```

Conda asks to procees, press Y

```
proceed ([y]/n)?
```

If you are using Miniconda3 in MicrosofT Visual Code Studio, <a href = https://medium.com/@udiyosovzon/how-to-activate-conda-environment-in-vs-code-ce599497f20d>Conda may have not been recognized</a>, you need to connect it manually in JSON settings:
```
{
    ... # any other settings you have already added (remove this line)

    "terminal.integrated.shell.windows": "<path-to-conda-cmd.exe>\\cmd.exe",
    “terminal.integrated.shellArgs.windows”: [“/K”, “C:\\<path-to-conda-installation>\\Scripts\\activate.bat C:\\<path-to-conda-installation>”]
    "python.condaPath": "C:\\<path-to-conda.exe>conda.exe"
}
```

# Running the script
Conda should be able to start using:
```
python CCi.py
```
