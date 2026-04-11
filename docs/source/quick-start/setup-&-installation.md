A more extensive guide to using OpenLISEM, and detailed explanation of all parts of the model can be found in the [manual](../using-openlisem/manual/index.md). Throughout the quick start we link to sections of the manual were applicable. 

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