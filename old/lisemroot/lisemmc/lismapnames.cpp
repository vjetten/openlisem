//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "lismapnames.h"
#include "iface.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)

//---------------------------------------------------------------------------
void __fastcall TLisIFace::SizeMapNames(TStringGrid *N, int w1, int w2,int w3)
{
     N->ColWidths[0]=w1;
     N->ColWidths[1]=w2;
     N->ColWidths[2]=w3;
     N->ColWidths[3]=0;//PageControl->Width-(w1+w2+w3)-48;
     N->ColWidths[4]=0;
     N->ColWidths[5]=0;

     N->Cells[0][0] = "Variable";
     N->Cells[1][0] = "Mapname";
     N->Cells[2][0] = "Description";
     N->Cells[3][0] = "Path";
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::SizeMapNamesL(TStringGrid *N, int w1, int w2,int w3, int w4)
{
     N->ColWidths[0]=w1;
     N->ColWidths[1]=w2;
     N->ColWidths[2]=w3;
     N->ColWidths[3]=0;//w4;//PageControl->Width-(w1+w2+w3)-32-N->Left;
     N->ColWidths[4]=0;
     N->ColWidths[5]=0;

     N->Cells[0][0] = "Variable";
     N->Cells[1][0] = "Mapname";
     N->Cells[2][0] = "Description";
     N->Cells[3][0] = "Path";
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::FillMapNames(TStringGrid *N, int i,
           AnsiString Varname, AnsiString Filename, AnsiString Desc, AnsiString ID)
{
     if (!CheckAdjustMapDirectoryName)
     {
          N->Cells[0][i] = Varname;
          N->Cells[1][i] = Filename;
          N->Cells[2][i] = Desc;
     }

     if (E_MapDir->Text.IsEmpty())
       N->Cells[4][i] = Filename;
     else
     if (N->Cells[3][i] == OldMapDir || ForceMapdir)
     {
       N->Cells[3][i] = E_MapDir->Text;
       if (E_MapDir->Text.IsPathDelimiter(E_MapDir->Text.Length()))
          N->Cells[4][i] = E_MapDir->Text+Filename;
       else
          N->Cells[4][i] = E_MapDir->Text+"\\"+Filename;
     }

     if (!CheckAdjustMapDirectoryName)
     {
     N->Cells[5][i] = ID;
     }
     N->Width = 5+3 + N->ColWidths[0] + N->ColWidths[1] + N->ColWidths[2] + N->ColWidths[3];
     if (i==1)
        N->Height = N->RowHeights[0]+2+3;
     N->Height += N->RowHeights[i]+1;
     N->RowCount =i+1;

}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::FillMapNamesO(TStringGrid *N, int i,
           AnsiString Varname, AnsiString Filename, AnsiString Desc, AnsiString ID)
{
     N->Cells[0][i] = Varname;
     N->Cells[1][i] = Filename;
     N->Cells[2][i] = Desc;
     N->Cells[3][i] = E_ResultDir->Text;
     if (E_ResultDir->Text.IsEmpty())
       N->Cells[4][i] = Filename;
     else
     {
       if (E_ResultDir->Text.IsPathDelimiter(E_ResultDir->Text.Length()))
          N->Cells[4][i] = E_ResultDir->Text+Filename;
       else
          N->Cells[4][i] = E_ResultDir->Text+"\\"+Filename;
     }
     N->Cells[5][i] = ID.UpperCase();

     N->Width = 5+3 + N->ColWidths[0] + N->ColWidths[1] + N->ColWidths[2] + N->ColWidths[3];
     if (i==1)
        N->Height = N->RowHeights[0]+2+3;
     N->Height += N->RowHeights[i]+1;
     N->RowCount =i+1;

}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::InitMapNames()
{
//VJ 050812: soildepth text indicated cm should be mm 
    int i = 1;
    int n = 0;
    int w1 = 80, w2=140, w3=380;
    SizeMapNames(MapsCatchment,w1,w2,w3);
 //   FillMapNames(MapsCatchment,i,"Area","area.map","Catchment area mask, boundary check for all other maps","area");i++;
    FillMapNames(MapsCatchment,i,"Gradient","grad.map","Slope gradient in direction of flow","grad");i++;
    FillMapNames(MapsCatchment,i,"LDD","ldd.map","Local Drain Direction network, steepest slope","ldd");i++;
//    FillMapNames(MapsCatchment,i,"LDDtill","lddtill.map","Local Drain Direction network, tillage","lddtill");i++;
    FillMapNames(MapsCatchment,i,"Outlet","outlet.map","Main catchment outlet, corresponding to LDD map","outlet");i++;
    FillMapNames(MapsCatchment,i,"ID","id.map","Raingauge zone ID numbers","ID");i++;
    FillMapNames(MapsCatchment,i,"OutPoint","outpoint.map","Reporting points (cells > 0) for runoff, main outlet must have value 1","outpoint");i++;
    n = i-1;
    i = 1;
    SizeMapNames(MapsLanduse,w1,w2,w3);
    FillMapNames(MapsLanduse,i,"Cover","per.map","Fraction surface cover by vegetation and residue","cover");i++;
    FillMapNames(MapsLanduse,i,"LAI","lai.map","Leaf area index of the plant cover in a gridcell (m2/m2)","lai");i++;
    FillMapNames(MapsLanduse,i,"Height","ch.map","Plant height (m)","ch");i++;
    FillMapNames(MapsLanduse,i,"Smax","smax.map","Plant canopy storage (mm)","smax");i++;
    FillMapNames(MapsLanduse,i,"Road width","roadwidt.map","Width of impermeable roads (m)","road");i++;
    FillMapNames(MapsLanduse,i,"Grass strips","grasswid.map","Width of grass strips (m)","grasswidth");i++;
    n = n+i-1;
// VJ 040514 include buffers
    i = 1;
    SizeMapNames(MapsBuffers,w1,w2,w3);
    FillMapNames(MapsBuffers,i,"Buffer ID nr","bufferid.map","ID number for each buffer starting with 1 (0 is outside area)","bufferID");i++;
    FillMapNames(MapsBuffers,i,"Buffer volume","buffervol.map","Buffer volumes at the locations of the buffers (m3)","bufferVolume");i++;
 //   FillMapNames(MapsBuffers,i,"Buffer area","bufferarea.map","Buffer area at locations of the buffers (m2)","bufferarea");i++;
//    FillMapNames(MapsBuffers,i,"Buffer outflow","bufferQ.map","Buffer gage outflow discharge (m3/s)","bufferdischarge");i++;
    n = n+i-1;
//VJ 080423 Snowmelt
    i = 1;
    SizeMapNames(MapsSnowmelt,w1,w2,w3);
    FillMapNames(MapsSnowmelt,i,"Snowmelt ID","snowid.map","Snowmelt zone ID number for snowmelt file starting with 1 (0 is non-snow area)","SnowID");i++;
    n = n+i-1;

 //   FillMapNames(MapsLanduse,i,"Manning Grass","grassman.map","mannings n of grass strips (-)","grassman");i++;

    i = 1;
    SizeMapNames(MapsSurface,w1,w2,w3);
    FillMapNames(MapsSurface,i,"RR","rr.map","Random Roughness (here standard deviation of heights) (cm)","rr");i++;
    FillMapNames(MapsSurface,i,"n","n.map","Manning's n (-)","manning");i++;
    FillMapNames(MapsSurface,i,"Crust","crustfrc.map","Fraction of gridcell covered with Crust (-)","crustfrc");i++;
    FillMapNames(MapsSurface,i,"Compacted","compfrc.map","Fraction of gridcell compacted (-)","compfrc");i++;
    FillMapNames(MapsSurface,i,"Stoniness","stonefrc.map","Fraction of gridcell covered by stones (-)","stonefrc");i++;
    FillMapNames(MapsSurface,i,"Hard Surface","hardsurf.map","no interc/infil/detach (1); normal surface (0), do not use for roads","hardsurf");i++;
    n = n+i-1;
    i = 1;
    SizeMapNames(MapsErosion,w1,w2,w3);
    FillMapNames(MapsErosion,i,"Cohesion","coh.map","Cohesion (kPa)","coh");i++;
    FillMapNames(MapsErosion,i,"Cohesion","cohadd.map","Extra cohesion factor by e.g. plant roots (kPa)","cohadd");i++;
    FillMapNames(MapsErosion,i,"Aggregates","aggrstab.map","Aggregate stability for splash erosion (-)","aggrstab");i++;
    FillMapNames(MapsErosion,i,"D50","d50.map","Median of the texture of the suspendeed matter (mu)","d50");i++;
    n = n+i-1;
    i = 1;
    SizeMapNames(MapsTexture,w1,w2,200);
    FillMapNames(MapsTexture,i,"Class 1","mu0.map","Clay fraction (MUST BE CLAY < 2mu) ","fractionmu0");i++;
    FillMapNames(MapsTexture,i,"Class 2","mu1.map","Soil texture fraction for class 1 (-) ","fractionmu1");i++;
    FillMapNames(MapsTexture,i,"Class 3","mu2.map","Soil texture fraction for class 2 (-) ","fractionmu2");i++;
    FillMapNames(MapsTexture,i,"Class 4","mu3.map","Soil texture fraction for class 3 (-) ","fractionmu3");i++;
    FillMapNames(MapsTexture,i,"Class 5","mu4.map","Soil texture fraction for class 4 (-) ","fractionmu4");i++;
    FillMapNames(MapsTexture,i,"Class 6","mu5.map","Soil texture fraction for class 5 (-) ","fractionmu5");i++;
    n = n+i-1;
    i = 1;
    SizeMapNames(MapsInfilSwatre,w1,w2,w3);
    FillMapNames(MapsInfilSwatre,i,"Profile table","profile.inp","Text file with land unit descriptions and links to tables","profinp");i++;
    FillMapNames(MapsInfilSwatre,i,"Profile ","profile.map","ID numbers corresponding to land units in profile table","profmap");i++;
    FillMapNames(MapsInfilSwatre,i,"Prof. Crust","profcrst.map","ID numbers of crusted soils (using also profile table)","profcrst");i++;
    FillMapNames(MapsInfilSwatre,i,"Prof. Wheel","profwltr.map","ID numbers of compacted wheel tracks (using also profile table)","profwltr");i++;
    FillMapNames(MapsInfilSwatre,i,"Prof. Grass","profgras.map","ID numbers of grasstrips (using also profile table)","profgras");i++;
    FillMapNames(MapsInfilSwatre,i,"Init. suction","inithead","Initial matrix potential of layers 001 to nnn, give only filename, not extention (cm)","inithead");i++;
    FillMapNames(MapsInfilSwatre,i,"Output head","headout.map","Locations to write tables of the matrix potential","headout");i++;
    n = n+i-1;
    i = 1;
    SizeMapNames(MapsInfilGA1,w1,w2,w3);
    FillMapNames(MapsInfilGA1,i,"Ksat1","ksat1.map","Layer 1: Saturated Hydraulic Conductivity (mm/h)","ksat1");i++;
    FillMapNames(MapsInfilGA1,i,"Psi1","psi1.map","Layer 1: Average suction at the wetting front (cm) ","psi1");i++;
    FillMapNames(MapsInfilGA1,i,"Thetas1","thetas1.map","Layer 1: Porosity (-)","thetas1");i++;
    FillMapNames(MapsInfilGA1,i,"Thetai1","thetai1.map","Layer 1: Initial moisture content (-)","thetai1");i++;
    FillMapNames(MapsInfilGA1,i,"Depth1","soildep1.map","Layer 1: Depth (mm) to bottom of layer 1","soildep1");i++;
    n = n+i-1;
    i = 1;
    SizeMapNames(MapsInfilGA2,w1,w2,w3);
    FillMapNames(MapsInfilGA2,i,"Ksat2","ksat2.map","Layer 2: Saturated Hydraulic Conductivity (mm/h)","ksat2");i++;
    FillMapNames(MapsInfilGA2,i,"Psi2","psi2.map","Layer 2: Average suction at the wetting front (cm) ","psi2");i++;
    FillMapNames(MapsInfilGA2,i,"Thetas2","thetas2.map","Layer 2: Porosity (-)","thetas2");i++;
    FillMapNames(MapsInfilGA2,i,"Thetai2","thetai2.map","Layer 2: Initial moisture content (-)","thetai2");i++;
    FillMapNames(MapsInfilGA2,i,"Depth2","soildep2.map","Layer 2: Depth (mm) to bottom of layer 2","soildep2");i++;
    n = n+i-1;
    i = 1;
    SizeMapNames(MapsInfilMorel,w1,w2,w3);
    FillMapNames(MapsInfilMorel,i,"Ksat1","ksat1.map","Layer 1: Saturated Hydraulic Conductivity (mm/h)","ksat1");i++;
    FillMapNames(MapsInfilMorel,i,"G1","G1.map","Layer 1: Driving force ","psi1");i++;
    FillMapNames(MapsInfilMorel,i,"Thetas1","thetas1.map","Layer 1: Porosity (-)","thetas1");i++;
    FillMapNames(MapsInfilMorel,i,"Thetai1","thetai1.map","Layer 1: Initial moisture content (-)","thetai1");i++;
    FillMapNames(MapsInfilMorel,i,"Depth1","soildep1.map","Layer 1: Depth (mm)","soildep1");i++;
    n = n+i-1;
    i = 1;
    SizeMapNames(MapsInfilSmith,w1,w2,w3);
    FillMapNames(MapsInfilSmith,i,"Ksat1","ksat1.map","Layer 1: Saturated Hydraulic Conductivity (mm/h)","ksat1");i++;
    FillMapNames(MapsInfilSmith,i,"G1","G1.map","Layer 1: Driving force ","psi1");i++;
    FillMapNames(MapsInfilSmith,i,"Thetas1","thetas1.map","Layer 1: Porosity (-)","thetas1");i++;
    FillMapNames(MapsInfilSmith,i,"Thetai1","thetai1.map","Layer 1: Initial moisture content (-)","thetai1");i++;
    FillMapNames(MapsInfilSmith,i,"Depth1","soildep1.map","Layer 1: Depth (mm)","soildep1");i++;
    n = n+i-1;
    i = 1;
    SizeMapNames(MapsInfilExtra,w1,w2,w3);
    FillMapNames(MapsInfilExtra,i,"Ksat Crust","ksatcrst.map","Ksat of crusts (all models except SWATRE) (mm/h)","ksatcrst");i++;
    FillMapNames(MapsInfilExtra,i,"Ksat Compact","ksatcomp.map","Ksat of compacted areas (all models except SWATRE) (mm/h)","ksatcomp");i++;
    FillMapNames(MapsInfilExtra,i,"Ksat Grass","ksatgras.map","Ksat of grassstrips (all models except SWATRE) (mm/h)","ksatgras");i++;
 //   FillMapNames(MapsInfilExtra,i,"Ksat Wheel","ksatwt.map","Ksat of wheel tracks (all models except SWATRE) (mm/h)","ksatwt");i++;
    i = 1;
//VJ 050812 drainage with GA
    SizeMapNames(MapsInfilDrainage,w1,w2,w3);
    FillMapNames(MapsInfilDrainage,i,"Drain. fact.","drfactor.map","Drainage exponent in k=ks*(moist/pore)^d (-)","drfactor");
    i = 1;
//VJ 080521 drainage with GA
    SizeMapNames(MapsInfilKsat,w1,w2,w3);
    FillMapNames(MapsInfilKsat,i,"Ksat1","ksat1.map","Saturated Hydraulic Conductivity (mm/h)","ksat1");
    i = 1;
    SizeMapNames(MapsInfilHoltan,w1,w2,w3);
    FillMapNames(MapsInfilHoltan,i,"A","A","A","A");i++;
    FillMapNames(MapsInfilHoltan,i,"FP","FP","FP","FP");i++;
    FillMapNames(MapsInfilHoltan,i,"P","P","P","P");i++;

    i = 1;
    SizeMapNames(MapsChannels,w1,w2,w3);
    FillMapNames(MapsChannels,i,"LDD ","lddchan.map","LDD of main channel (must be 1 branch connected to the outlet)","lddchan");i++;
    FillMapNames(MapsChannels,i,"Width","chanwidt.map","Channel width (m)","chanwidth");i++;
    FillMapNames(MapsChannels,i,"Side angle","chanside.map","Channel side angle (between channel side and surface: 0 is rectangular, 1 = 45o)","chanside");i++;
    FillMapNames(MapsChannels,i,"Gradient","changrad.map","Slope gradient of channel bed (-)","changrad");i++;
    FillMapNames(MapsChannels,i,"N","chanman.map","Mannings n of channel bed (-)","chanman");i++;
    FillMapNames(MapsChannels,i,"Cohesion","chancoh.map","Cohesion of channel bed (kPa)","chancoh");i++;
    i = 1;
    SizeMapNames(MapsChannelinfil,w1,w2,w3);
    FillMapNames(MapsChannelinfil,i,"Ksat","chanksat.map","Infiltration rate of channel bed (mm/h)","chanksat");i++;
//    FillMapNames(MapsChannelinfil,i,"Storage","chanstor.map","Storage capacity of channel (mm)","chanstor");i++;

//VJ 080217 add baseflow variables
    i = 1;
    SizeMapNames(MapsChannelBaseflow,w1,w2,w3);
    FillMapNames(MapsChannelBaseflow,i,"Inflow flux","chanbaseflux.map","Incoming flux into channel from the two sides (m3/s)","chanbaseflux");i++;
    FillMapNames(MapsChannelBaseflow,i,"Increase in baseflow","chanincrease.map","Increase in basevolume during rainstorm (-)","chanincrease");i++;
    FillMapNames(MapsChannelBaseflow,i,"Initial volume","chanvini.map","Initial baseflow water volume in channel (m3)","chanvolini");i++;

    i = 1;
    SizeMapNames(MapsWheeltrack,w1,w2,w3);
    FillMapNames(MapsWheeltrack,i,"LDD ","lddwheel.map","LDD of wheeltrack network (can be separate branches with pits)","lddwheel");i++;
    FillMapNames(MapsWheeltrack,i,"Number","wheelnbr.map","Number of wheeltrack channels in a gridcell (-)","wheelnbr");i++;
    FillMapNames(MapsWheeltrack,i,"Width","wheelwid.map","Sum of widths of wheeltracks in a gridcell (m)","wheelwidth");i++;
    FillMapNames(MapsWheeltrack,i,"Depth","wheeldep.map","Wheel track overflow depth (cm)","wheeldepth");i++;
    FillMapNames(MapsWheeltrack,i,"Gradient","wheelgrd.map","Slope gradient of wheel tracks (-)","wheelgradient");i++;
    FillMapNames(MapsWheeltrack,i,"N","wheelman.map","Mannings n of Wheel tracks (-)","wheelman");i++;
    FillMapNames(MapsWheeltrack,i,"Cohesion","wheelcoh.map","Cohesion of wheel tracks (kPa)","wheelcohesion");i++;
    FillMapNames(MapsWheeltrack,i,"Ksat","ksatwt.map","Saturated hydraulic conductivity of wheel tracks (mm/h)","ksatwt");i++;
    i = 1;
    SizeMapNames(MapsNutsBD,w1,w2,w3);
    FillMapNames(MapsNutsBD,i,"Bulk Dens.","bulkdens.map","Bulk density of the topsoil (kg/m3)","bulk");i++;
    i = 1;
    SizeMapNames(MapsNutsP,w1,w2,w3);
    FillMapNames(MapsNutsP,i,"Content","pcont.map","Phosphate (P) content of the soil (kg/kg)","pcont");i++;
    FillMapNames(MapsNutsP,i,"Solute","psolut.map","Initial solute store P in surface layer (kg/m2)","psolute");i++;
    FillMapNames(MapsNutsP,i,"Efficiency","peff.map","Extraction efficiency (s-1)","pefficiency");i++;
    FillMapNames(MapsNutsP,i,"Sorption","Psorp.map","Sorption isotherm kd (m3/kg)","psorp");i++;
    FillMapNames(MapsNutsP,i,"Conversion","Pconv.map","Conversion P from soil content to clay content(-)","pconv");i++;
    i = 1;
    SizeMapNames(MapsNutsNH4,w1,w2,w3);
    FillMapNames(MapsNutsNH4,i,"Content","nh4cont.map","Ammonium (NH4+) content of the soil (kg/kg)","nh4cont");i++;
    FillMapNames(MapsNutsNH4,i,"Solute","nh4solut.map","Initial solute store NH4 in surface layer (kg/m2)","nh4solute");i++;
    FillMapNames(MapsNutsNH4,i,"Efficiency","nh4eff.map","Extraction efficiency (s-1)","nh4efficiency");i++;
    FillMapNames(MapsNutsNH4,i,"Sorption","NH4sorp.map","Sorption isotherm kd (m3/kg)","nh4sorp");i++;
    FillMapNames(MapsNutsNH4,i,"Conversion","NH4conv.map","Conversion NH4 from soil content to clay content(-)","nh4conv");i++;
    i = 1;
    SizeMapNames(MapsNutsNO3,w1,w2,w3);
    FillMapNames(MapsNutsNO3,i,"Content","NO3cont.map","Nitrate (NO3-) content of the soil (kg/kg)","no3cont");i++;
    FillMapNames(MapsNutsNO3,i,"Solute","NO3solut.map","Initial solute store NO3 in surface layer (kg/m2)","no3solute");i++;
    FillMapNames(MapsNutsNO3,i,"Efficiency","NO3eff.map","Extraction efficiency (s-1)","no3efficiency");i++;
    FillMapNames(MapsNutsNO3,i,"Sorption","NO3sorp.map","Sorption isotherm kd (m3/kg)","no3sorp");i++;
    FillMapNames(MapsNutsNO3,i,"Conversion","NO3conv.map","Conversion NO3 from soil content to clay content(-)","no3conv");i++;
    i = 1;
    SizeMapNames(MapsGully,w1,w2,w3);
    FillMapNames(MapsGully,i,"DEM","dem.map","Digital elevation model (m)","dem");i++;
//    FillMapNames(MapsGully,i,"Initial  Width","gulwidth.map","Initial gully width (m)","gullyinit");i++;
    FillMapNames(MapsGully,i,"mannings N","gullyman.map","manning's n gully bottom (-)","gullyn");i++;
    FillMapNames(MapsGully,i,"BulkDensity","bulkdens.map","Bulkdensity of topsoil (kg/m3)","bulkdens1");i++;
    FillMapNames(MapsGully,i,"Ksat","gulksat1.map","Ksat of topsoil for gully infil (mm/h)","gulksat1");i++;
    FillMapNames(MapsGully,i,"Depth layer 2","soilDep2.map","Depth to subsoil (cm)","gullydep");i++;
    FillMapNames(MapsGully,i,"Cohesion 2","coh2.map","Cohesion of subsoil (kPa)","gullycoh");i++;
    FillMapNames(MapsGully,i,"BulkDensity 2","bulkden2.map","Bulkdensity of subsoil (kg/m3)","bulkdens2");i++;
    FillMapNames(MapsGully,i,"Ksat 2","gulksat2.map","Ksat of subsoil for gully infil (mm/h)","gulksat2");i++;
    FillMapNames(MapsGully,i,"Exclude","noncrit.map","areas to be excluded (roads etc.) = 1, rest = 0","nonfcrit");i++;
//VJ 040331 Included init gully dimensions
    i = 1;
    SizeMapNames(MapsGullyInit,w1,w2,w3);
    FillMapNames(MapsGullyInit,i,"Gully Width","gulwinit.map","Initial gully width (m)","gulwinit");i++;
    FillMapNames(MapsGullyInit,i,"Gully Depth","guldinit.map","Initial gully depth (m)","guldinit");i++;
    InitOutMapNames();
    if (!CheckAdjustMapDirectoryName)
    {
      TextureClass->Cells[0][0] = "Class (mu)";
      TextureClass->Cells[0][1]= "2";
      TextureClass->Cells[0][2]= "16";
      TextureClass->Cells[0][3]= "32";
      TextureClass->Cells[0][4]= "53";
      TextureClass->Cells[0][5]= "75";
      TextureClass->Cells[0][6]= "105";
    }

}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::InitOutMapNames()
{
    int i = 1;
    int w1 = 86, w2=60, w3= 230, w4 = 250;
    SizeMapNamesL(MapsOutputBASIC,w1,w2,w3,w4);
    FillMapNamesO(MapsOutputBASIC,i,"Runoff","ro"," runoff (l/s)","outrunoff");i++;
    FillMapNamesO(MapsOutputBASIC,i,"Concentration","conc"," sediment concentration (g/l)","outconc");i++;
    FillMapNamesO(MapsOutputBASIC,i,"Water height","wh"," water height on surface (mm)","outwh");i++;
    FillMapNamesO(MapsOutputBASIC,i,"Cumulative WH","whc","Cumulative runoff water height on surface (mm)","outrwh");i++;
    FillMapNamesO(MapsOutputBASIC,i,"Transport cap.","tc"," transport capacity (excl. channels)(g/l)","outtc");i++;
    FillMapNamesO(MapsOutputBASIC,i,"Detachment","det"," detachment (kg)","outeros");i++;
    FillMapNamesO(MapsOutputBASIC,i,"Deposition","depo"," deposition (kg)","outdepo");i++;
    FillMapNamesO(MapsOutputBASIC,i,"Velocity","velo"," velocity (m/s)","outvelo");i++;
    FillMapNamesO(MapsOutputBASIC,i,"Infiltration","inf","Cumulative infil (mm)","outinf");i++;
    FillMapNamesO(MapsOutputBASIC,i,"Surface Storage","sstor"," Surface storage (mm)","outss");i++;
    FillMapNamesO(MapsOutputBASIC,i,"Channel volume","chanvol"," Channel water volume (m3)","outchvol");i++;
    i = 1;
    SizeMapNamesL(MapsOutputNut,w1,w2,w3,w4);
    FillMapNamesO(MapsOutputNut,i,"Nuts P solut.","NPsol"," P in solution (kg)","outpsolut");i++;
    FillMapNamesO(MapsOutputNut,i,"Nuts P susp.","NPsus"," P in suspension (kg)","outpsus");i++;
    FillMapNamesO(MapsOutputNut,i,"Nuts P inf.","NPinf"," P in solution (kg)","outpinf");i++;
    FillMapNamesO(MapsOutputNut,i,"Nuts NH4 solut.","NNH4sol"," NH4 in solution (kg)","outnh4solut");i++;
    FillMapNamesO(MapsOutputNut,i,"Nuts NH4 susp.","NNH4sus"," NH4 in suspension (kg)","outnh4sus");i++;
    FillMapNamesO(MapsOutputNut,i,"Nuts NH4 inf.","NNH4inf"," NH4 in solution (kg)","outnh4inf");i++;
    FillMapNamesO(MapsOutputNut,i,"Nuts NO3 solut.","NNO3sol"," NO3 in solution (kg)","outNO3solut");i++;
    FillMapNamesO(MapsOutputNut,i,"Nuts NO3 susp.","NNO3sus"," NO3 in suspension (kg)","outNO3sus");i++;
    FillMapNamesO(MapsOutputNut,i,"Nuts NO3 inf.","NNO3inf"," NO3 in solution (kg)","outno3inf");i++;
    //VJ 030415 changed suspension to solution in infiltration text
    i = 1;
    SizeMapNamesL(MapsOutputNutErosDep,w1,w2,w3,w4);
    FillMapNamesO(MapsOutputNutErosDep,i,"P Dep.","NPdep.map"," P in clay deposition (kg/m2)","outpdep");i++;
    FillMapNamesO(MapsOutputNutErosDep,i,"NH4 Dep.","NNH4dep.map"," NH4 in clay deposition (kg/m2)","outnh4dep");i++;
    FillMapNamesO(MapsOutputNutErosDep,i,"NO3 Dep.","NNO3dep.map"," NO3 in clay deposition (kg/m2)","outno3dep");i++;
    FillMapNamesO(MapsOutputNutErosDep,i,"P Detach.","NPdet.map"," P in clay detachment (kg/m2)","outpdet");i++;
    FillMapNamesO(MapsOutputNutErosDep,i,"NH4 Detach.","NNH4det.map"," NH4 in clay detachment (kg/m2)","outnh4det");i++;
    FillMapNamesO(MapsOutputNutErosDep,i,"NO3 Detach.","NNO3det.map"," NO3 in clay detachment (kg/m2)","outno3det");i++;
    i = 1;
    SizeMapNamesL(MapsOutputGul,w1,w2,w3,w4);
    FillMapNamesO(MapsOutputGul,i,"Gully depth","guld","gully depth (m)","outguld");i++;
    FillMapNamesO(MapsOutputGul,i,"Gully width","gulw","gully width (m)","outgulw");i++;
    FillMapNamesO(MapsOutputGul,i,"Gully area","gula","gully cross section area (m2)","outgula");i++;
    FillMapNamesO(MapsOutputGul,i,"Gully fraction","gulf","cross section, fraction of foot^2 (-)","outgulf");i++;
    FillMapNamesO(MapsOutputGul,i,"Gully DEM","guldem","digital elevation model (m)","outguldem");i++;
    i = 1;
    SizeMapNamesL(MapsOutputMC,w1,w2,w3,w4);
    FillMapNamesO(MapsOutputMC,i,"Sed. Conc mu0","smu0","suspended sediment flux class 0 (mu)","outmu0");i++;
    FillMapNamesO(MapsOutputMC,i,"Sed. Conc mu1","smu1","suspended sediment flux class 1 (mu)","outmu1");i++;
    FillMapNamesO(MapsOutputMC,i,"Sed. Conc mu2","smu2","suspended sediment flux class 2 (mu)","outmu2");i++;    FillMapNamesO(MapsOutputMC,i,"Sed. Conc mu3","smu3","suspended sediment flux class 3 (mu)","outmu3");i++;
    FillMapNamesO(MapsOutputMC,i,"Sed. Conc mu4","smu4","suspended sediment flux class 4 (mu)","outmu4");i++;
    FillMapNamesO(MapsOutputMC,i,"Sed. Conc mu5","smu5","suspended sediment flux class 5 (mu)","outmu5");i++;
    FillMapNamesO(MapsOutputMC,i,"D50","D50s","D50 suspended sediment (mu)","outD50susp");i++;

}
