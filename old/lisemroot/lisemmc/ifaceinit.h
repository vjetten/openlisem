//---------------------------------------------------------------------------

#ifndef ifaceinitH
#define ifaceinitH
//---------------------------------------------------------------------------

#define hintlisemstart "Goto startscreen LISEM: undo all changes, reset the interface and unload runfiles."
#define hintfileopen "Open one or more runfiles. They are added to a drop down list at the bottom. Each file is executed in sequence, when stopped the next runfile is stared. Press the stop-all button to end all runs."
#define hintfilesave "Save the current run file, all changes made in the interface will be stored."
#define hintfilesaveas "Save the current run file under a different name. This is the easiest way to make a new runfile."
#define hintE_MapDir       "Directory of the input maps. This name will we added to all mapnames. If this is changed all input map pathnames will be changed too."
#define hintE_RainfallName "Rainfall file in LISEM format (not the same as PCRaster format). The first column is\
 the time in minutes. Other columns have rainfall intensity in mm/h (one column per station). The intensity is\
 assigned to the interval as follows: if the first line is \"0.0 0.0\" and the second is \"10.0 15.0\", LISEM assumes 15 mm/h from minute 0 to 10."
#define hintE_SnowmeltName "Snowmelt file in LISEM format. The format is the same as the rainfall format."
#define hintE_begintime    ""
#define hintE_Endtime      ""
#define hintE_Timestep     ""




#endif
