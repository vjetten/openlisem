# Contents of the quick start

* [Setup & Installation](https://github.com/vjetten/openlisem/wiki/Getting-started#setup--installation)
* [The user interface](https://github.com/vjetten/openlisem/wiki/Getting-started#the-user-interface)
* [Required input](https://github.com/vjetten/openlisem/wiki/Getting-started#required-input)
* [Output & results](https://github.com/vjetten/openlisem/wiki/Getting-started#output--results)

A more extensive guide to using OpenLISEM, and detailed explanation of all parts of the model can be found in the [manual](). Throughout the quick start we link to sections of the manual were applicable. 

## Setup & Installation

**NOTE: this setup works for windows systems, for linux systems users are referred to the manual**

To use OpenLISEM, the model has to be installed on your system, and some other software is required for a full use of the model. The OpenLISEM model is available as an executable program for Windows (64-bit), and can be [compiled]() both under Windows and Linux systems. The OpenLISEM model uses the [PCRaster](https://pcraster.geo.uu.nl/) GIS environment for input and output maps. To work with this PCRaster and user interface, Nutshell, is available. The easiest way to use PCRaster is by installing it through a conda environment on your system. Besides this minimal installation, a dedicated GIS program is usefull for database preparation, for example the opensource [QGIS](https://qgis.org/en/site/).

Follow the steps below to setup your PC:

1. Install miniconda
1. Install PCRaster software with miniconda
1. download and install Nutshell
1. download and install OpenLISEM

### 1. Install miniconda
The PCRaster software we use is installed in **Conda**. Conda is an installation environment independent of programming language and platform (windows, unix). On Windows machines, first download and install the "minimal conda environment" **miniconda**.

[Download](https://docs.conda.io/en/latest/miniconda.html) and install the latest 64bit version.

Anaconda/miniconda is a repository that works with so called 'environments'. These are separated folders under which you can install programs and libraries. Each environment is independent from other environments, if something goes wrong you can delete the entire environment without damaging other installs. We will create an environment called “lisem” and install everything in there.  To open miniconda go to your program list and look under ’A’ for Anaconda3 (64-bit) and under that folder click on Anaconda Prompt (miniconda3). 

This opens a command window. Note that the prompt says: 

> (base) PS C:\path\to\username>

Now miniconda is installed, we can proceed by installing PCRaster in a separate environment.

### 2. Install PCRaster with miniconda

In the miniconda shell, we are now in the 'base' environment. For OpenLISEM we want to create a new environment, we suggest to call this 'lisem'.

> conda create --name lisem

Now activate the 'lisem' environment. Whatever you install now will be done in this workspace/environment. In that way you can separate different projects:

> conda activate lisem

Install the following packages (answer “y” to prompts):

> conda install -c conda-forge pcraster gdal

This takes a while, but this should be all. Check if it worked:  

> pcrcalc

Should return a help syntax. If you already have installed python on your laptop, it may cause interference. Please remove any python references from your windows path.

### 3. Download and install Nutshell

PCRaster does not have a menu of interface, it can be operated from the command line. For ease of operations this was created separately and is called [Nutshell](https://github.com/vjetten/NutShell/releases/tag/NutShell). This program is not part of the conda system and is downloaded as a zip file, which you can unzip on your system. 

**TIP: store the software in a dedicated programs folder on your system, so you can use it across multiple prohjects, e.g. '~/user/programs/'.**

[Download Nutshell](https://github.com/vjetten/NutShell/releases/download/NutShell/nutshell-v5.14.1-mapedit-v3.3_220923.zip) and unzip it in your programs folder. Open the NutShell.exe. The program checks all environments and looks for PCRaster and gives a warning if it finds an environment that does not have the proper libs. In the interface, click on Files and select options. This opens a window which lets you select the PCRaster environment location. It already points to the Miniconda environment folder. If you have multiple environments, select one. Click on the folder icon and select the “lisem” environment:

<p align="center">
  <img width="569" height="346" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/NutShell_interface.png">
</p>

To test if everything works:

> pcrcalc

This should again show the help syntax. PCRaster and NutShell are now ready for use! An introduction to using PCRaster and Nutshell can be found [here](https://github.com/vjetten/openlisem/wiki/Introduction-PCRaster-&-Nutshell).

### 4. Download and install OpenLISEM

Download the latest release from [OpenLISEM](https://github.com/vjetten/openlisem/releases/). Unzip the files in your programs folder. Find the Lisem.exe in the folder and open the application. This is the basic setup of OpenLISEM, and your system is now ready to start using the model!

## The user interface

The OpenLISEM interface has 3 pages, see the figures below:  
1. Has the main categories for running the model. On the left side are the INPUT definitions and options, to the right the OUTPUT definitions and options.  
1. A tree structure that shows the names of all input maps. You can click on each name and change names to your own names, but it is far more convenient to use default names and make sure that different scenarios are stored in different folders.  
1. Model results interface that shows summarized data, hydrographs and spatial output while running.  

The main interface of the software shows the model settings. Here, the model controls, process settings and output settings are visible. 

<p align="center">
  <img width="850" height="540" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/OpenLISEM_input_interface.png">
</p>

When opening the map database, an overview of required input maps is available. Double-clicking a map name allows you to select a file, although it is usually easier to use the default filenames. The lighter text indicates an input map is not required for the current simulation, as the relevant process has not been activated. Input maps are usually constructed using the PCRaster + Nutshell software -or the new LISEM database generator-. All input maps must be of identical rows and columns, but might be of any format supported by GDAL (.map, GeoTiff, ascii, etc..).

<p align="center">
  <img width="818" height="597" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/OpenLISEM_mapDB.png">
</p>

When the simulation is running, real-time output can be visualized in the simulation tab. Here, hydrological statistics of the simulation are provided in catchment-average mm equivalent. Additionally, a spatial view of many model variables can be shown, such as infiltration, flow heights, flow velocities and soil loss. For specified points in the simulation, such as channel outlets or bridges, a discharge curve (hydrograph) can be viewed during the simulation.

<p align="center">
  <img width="782" height="503" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/OpenLISEM_output_interface.png">
</p>

## Required input

OpenLISEM needs a minimum of 24 raster maps depending on the input options selected in the interface such as special surfaces, sediment modelling, and 2-layer sediment transport. In theory these maps can be derived from 5 base sources:

1. Rainfall
1. Topography
1. Land use/cover
1. Soil type
1. Infrastructure

The number of input maps needed for a OpenLISEM run depends on the selected options. Most of these options can be set from the user interface. Configuration of a OpenLISEM model run can be saved and loaded in a .run file. This is an easy to edit ASCII file containing a line for every option, with the option name and its value. Default values are used to fill up any missing data in a .run file. The options can also be used as command line arguments when OpenLISEM is used without user interface. For more information about this, see the section about running OpenLISEM.

Besides options about functionality, many options about technical parts of the model are available to users. Minimum time steps, courant factors and approximation methods are some of these options. In general, when not sure about the influence of an option on the model, it is best to use the default value for that option. The default values have been chosen in such a way that they should ensure optimal performance in a wide range of problems. Changing these values to outside their proper range can generate errors within the model, and lead to unrealistic behavior.

The raster format causes openLISEM to need a large data set. For instance a river channel is characterized by its network, shape, width, depth, bed slope, resistance, strength and infiltration rate or stationary baseflow (8 maps). Do not be scared by this!

For information about the preparation of all specific input maps and data we refer to the [manual](https://github.com/vjetten/openlisem/wiki/Manual). There we discuss the preparation of the following input data:

1. Rainfall [🔨 ](https://github.com/vjetten/openlisem/wiki/Prepare-Rainfall) [🔣 ](https://github.com/vjetten/openlisem/wiki/Meteorological-input)
1. Catchment [🔨](https://github.com/vjetten/openlisem/wiki/Preparing-Topography)
1. Land use [🔨](https://github.com/vjetten/openlisem/wiki/Prepare-Land-Use)
1. Soil properties & Infiltration [🔨](https://github.com/vjetten/openlisem/wiki/Prepare-Soil-Properties)
1. Channels [🔨](https://github.com/vjetten/openlisem/wiki/Prepare-Channels) [🔣 ](https://github.com/vjetten/openlisem/wiki/General-process-options#channels--rivers)
1. Infrastructure [🔨](https://github.com/vjetten/openlisem/wiki/Prepare-Infrastructure) [🔣 ](https://github.com/vjetten/openlisem/wiki/General-process-options#infrastructure) 

## Output & results

When running openLISEM, on the simulation tab, the real time values for many variables can be displayed over the catchment as well as at the outlet. The graph in the simulation tab will show the discharge, precipitation and if sediment dynamics are modeled also the suspended sediment load and concentration. These values are displayed for the outlet, but additional points can be added to the database which also will be monitored by OpenLISEM. More information on output interpretation is available in the [manual]().

<p align="center">
  <img width="682" height="384" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/Hydrograph_output.png">
</p>

Besides the discharge at specific locations, the OpenLISEM model can create insight in spatial dynamics within the modeled catchment. On the map display the dynamics for water or sediment related variables can be viewed during the model run.

<p align="center">
  <img width="782" height="503" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/Infiltration_spatial_output.png">
</p>

For further analysis the display of OpenLISEM is not very suitable. However the model stores the hydrograph values in the selected 'results' folder. Moreover it is possible to store map series of variables of interest, this can be selected at the right half of the input tab. Map series can be analyzed and viewed after the model run is finished, for example with NutShell.