# Run OpenLISEM with the user interface

🚧 add information here

# OpenLISEM from the command line

In some cases it is required to run OpenLISEM multiple times. For example for sensitivity or uncertainty analyses or for auto-calibration. This can be automated by running OpenLISEM from the command line. Two option exist; a batch mode including the graphical user interface (GUI) and a version without GUI. The last option is designed to work on systems without a graphical interface, for example high performance computing systems.
All settings which can be set in the GUI are stored in a runfile. To make sure the command line execution of OpenLISEM works smoothly we advise to first run the model from the GUI. If this works and the results are within expectations the runfile can be loaded from the command line.

In the default batch mode, OpenLISEM will automatically open the user interface for each run. When the `-ni` option is used, only minor run information is shown in the command line application.

To load a `.run` file from a command line application, use the following command:
```
Lisem [-ni] -r runfilename.run
```
`-r`    required     provide the runfile name (including full path and extension)  
`-ni`   optional     use OpenLISEM without GUI

## Solving errors when running from the command line

When running OpenLISEM from the command line, the runfile must be the latest version, if this is not the case an error message will be returned. Besides that OpenLISEM can currently only run with absolute paths in the runfile.

### ERROR - Filename not found for map
This error occurs because the runfile has absolute pathnames, which do not point to the location where the map database currently is located. Adjust the pathnames for:

`Map Directory`, `Result Directory`, `Rainfall Directory` and `Rainfall Map Directory` to point to the correct absolute paths.

### ERROR - Map or Variable not found! You could be using an old runfile, or a map is not present

This error can occur because a map the is required is not present in the map database. If the database however does run from the User Interface, the most likely reason is that the runfile is from an older version and some maps or variable options are not present in the runfile. To solve this error follow these steps:

1. Open the runfile in the user interface of OpenLISEM. Many warning messages will be shown, accept them all.
2. In the Input section set the Database options to the correct folder, on the Rainfall tab, also adjust the directories.
3. (Optional) Check if all other settings are as you expect for your run.
4. Save these settings to a **new** runfile, with the 'save as' option. This will show a message that the runfile will be updated to the newest version.

This new runfile should run without errors from the command line.

## Container to run OpenLISEM

To run OpenLISEM on a high performance computing cluster (HPC), it can be convenient to have the model within a container which can easily be distributed, including all dependencies etc. This is available at: [apptainer openlisem](https://git.wur.nl/comme002/singularity-openlisem)
