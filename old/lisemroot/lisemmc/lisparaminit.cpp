// -----------------------------------------------------------------------------
// ******** calculate general variables needed *********************************
// -----------------------------------------------------------------------------

     //------------depression storage-------------------------------------------


/* OBSOLETE =====> moved to surfstor3.cpp

     _spatial(REAL4, DeprStoreMaxH);
     calc(" DeprStoreMaxH = max(10*(0.112*RR + 0.031*sqr(RR) - "
          "           0.012*RR*(Gradient*100)),0) ");
        // Onstad (1984); maximum storage in micro-depressions
        // RR (Random roughness) in cm; DeprStoreMaxH now in mm
        // DeprStoreMaxH cannot be less than 0

     _spatial(REAL4, DeprStoreRainH);
     calc(" DeprStoreRainH = max(10*(0.329*RR + 0.073*sqr(RR) - "
          "            0.018*RR*(Gradient*100)),0) ");
        // Onstad (1984); net rainfall needed to fill all micro-depressions
        // RR (Random roughness) in cm; DeprStoreRainH now in mm
        // DeprStoreRainH cannot be less than 0

    _spatial(REAL4, DeprStoreStartH);
     calc(" DeprStoreStartH = max(DeprStoreRainH*(0.0527*RR - 0.0049*(Gradient*100)),0) ");
        // threshold net rainfall after which runoff starts
        // runoff starts before all micro-depressions are filled
        // LISEM project (1993); DeprStoreStartH in mm
        // DeprStoreStartH cannot be less than    0

     _spatial(REAL4, PondAreaFractMax);
     calc(" PondAreaFractMax = max(0.152*RR - 0.008*sqr(RR) - 0.008*RR*(Gradient*100),0) ");
        // Onstad (1984); maximum fraction of surface covered with water
        // PondAreaFractMax is assumed not to be less than 0.10; lower values give
        // mass balance error for steep slopes!)
  //   _spatial(REAL4, PAFcoeff);
//     calc(" PAFcoeff = 1.4057*(RR*10)**-0.9422 ");
     //RR in MM in this equation!
//VJ 040224 MAJOR BUG hier stond RR/10 ipv maal 10, terwijl dit in mm moet zijn
//VJ 050302 mag weg, is vervangen door de formule
//      calc(" PondAreaFract = 1-exp(-1.88*WH/(RR*10)) ");

====> OBSOLETE */

// RR is given in cm in the maps
     _spatial(REAL4, MDS);
     calc(" MDS = max(0, 0.243*RR*10 + 0.010*RR*RR*100 - 0.012*RR*10*Gradient*100) ");
        // new threshold var according to kamphorst et al.  in MM !!!!!!
//VJ 050302 simplified this
     _spatial(REAL4, SDS);
//     calc(" SDS = -log10(0.9)/PAFcoeff ");
//     calc(" SDS = min(SDS, MDS*0.9)");
     calc(" SDS = 0.1*MDS ");
      // MDS can be 0 because of steep slopes, SDS should be smaller than MDS always
     // start of runoff depression
     // based on fraction ponded area: fpa = 1-exp(-aWH)
     // runoff starts if fpa is at 10% -> 0.1 = 1-exp(-aWH)
     // -0.9 = -exp(-aWH) -> WH = ln(0.9)/-a, with a = PAFcoeff

//     calc(" SDS = max(MDS*(0.0527*RR - 0.0049*(Gradient*100)),0) ");
       // based on previous, in MM !!!!!!
       

     //------------overland flow------------------------------------------------

     _spatial(REAL4, infilsurplus);
     calc("infilsurplus = 0");
     // infiltration surplus for kinematic wavr

     _nonspatial(REAL4, QOutflow );
     QOutflow = 0.0;

     _spatial(REAL4, WH);
     calc(" WH = 0.0");
        // WH is the waterheight, in mm

     _spatial(REAL4, WHRoad);
     calc(" WHRoad = 0.0");
        // WHRoad is the waterheight on roads, in mm

     _spatial(REAL4, WHWheelTrack);
     calc(" WHWheelTrack = 0.0");
        // WHWheelTrack is the waterheight on wheeltracks, in mm

     _spatial(REAL4, WHGrass);
     calc(" WHGrass = 0.0");
        // grasstrips waterheight

     _spatial(REAL4, WHCrust);
     calc(" WHCrust = 0.0");

     _spatial(REAL4, WHCompact);
     calc(" WHCompact = 0.0");

// VJ 040514 include buffers
     _spatial(REAL4, TotalBufferVolume);
     calc(" TotalBufferVolume = 0.0 ");

     _spatial(REAL4, BufferVolumeInit);
     _spatial(REAL4, BufferVolumeInitChannel);
     _spatial(REAL4, BufferVolumeCurrent);
     calc(" BufferVolumeCurrent = 0.0 ");
     _spatial(REAL4, BufferVolumeCurrentChannel);
     calc(" BufferVolumeCurrentChannel = 0.0 ");

     _spatial(REAL4, BufferSedVolumeCurrent);
     calc(" BufferSedVolumeCurrent = 0.0 ");
     _spatial(REAL4, BufferSedVolumeCurrentChannel);
     calc(" BufferSedVolumeCurrentChannel = 0.0 ");

     if (SwitchBuffers)
     {
       if (SwitchIncludeChannel)
       {
          calc(" BufferVolumeCurrentChannel = cover(mif(ChannelLDD, BufferVolume, 0),0) ");
          calc(" BufferVolumeInitChannel = BufferVolumeCurrentChannel ");
       }
       calc(" BufferVolumeCurrent = mif(BufferVolumeCurrentChannel eq 0, BufferVolume, 0) ");
       calc(" BufferVolumeInit = BufferVolumeCurrent ");

       // adapt parameters for buffer, avoid errors in kin wave
       if (!SwitchSedtrap)
       {
          calc(" Gradient = mif(BufferID gt 0, 0.001, Gradient) ");
          calc(" MDS = mif(BufferID gt 0, 0.1, MDS) ");
          calc(" SDS = mif(BufferID gt 0, 0, SDS) ");
          calc(" RR = mif(BufferID gt 0, 0.05, RR) ");
       }

       if (SwitchIncludeChannel)
       {
         calc(" ChannelGradient = mif(BufferID gt 0, 0.001, ChannelGradient) ");
         calc(" ChannelN = mif(BufferID gt 0, 0.5, ChannelN) ");
       }  
     }
     _nonspatial(REAL4, buffersize);
     buffersize = 0;
     _nonspatial(REAL4, sedbuffersize);
     sedbuffersize = 0;

     // needed for kin wave that is not using buffers
     _spatial(REAL4, Bufferdummy);
     calc(" Bufferdummy = 0 ");

     //-------------------------channel-----------------------------------------

     _spatial(REAL4, ChannelVolout);
     calc(" ChannelVolout = 0 ");
     //VJ 080217 moved from lisPchannel.cpp because of output

     _spatial(REAL4, Outlet1);
     calc(" Outlet1 = mif(Outlet gt 0 and Outlet lt 2, 1.0, 0.0) ");
        // remaining part of the map becomes 0
     _spatial(REAL4, Outlet2);
     calc(" Outlet2 = mif(Outlet gt 1 and Outlet lt 3, 1.0, 0.0) ");
        // remaining part of the map becomes MV
     _spatial(REAL4, Outlet3);
     calc(" Outlet3 = mif(Outlet gt 2 and Outlet lt 4, 1.0, 0.0) ");
        // remaining part of the map becomes MV

     //-------------------------erosion-----------------------------------------
     _nonspatial(REAL4, SedOutflow );
     SedOutflow = 0.0;

     _spatial(REAL4, SedConcentration);
     calc(" SedConcentration = 0 ");

     _spatial(REAL4, Sedin);
     calc(" Sedin = 0 ");

     _spatial(REAL4, TransportCapOUT);
     calc(" TransportCapOUT = 0 ");
     //needed for timeseries output



    if (!SwitchNoErosion)
    {

//VJ 031211 removed < 100 check for cohesion
     _spatial(REAL4, Y);
     calc(" Y = 1.0/(0.89+0.56*CohesionTotal)");
//     calc(" Y = mif(CohesionTotal <100, 1.0/(0.89+0.56*CohesionTotal), 0) ");
//     calc(" Y = 0.79*exp(-0.85*CohesionTotal) ");
     calc(" Y = max(0, min(Y,1.0)) ");
        // Y is the flow detachment efficiency coefficient
        // Ugmin = 1.0 cm/s (Rauws & Govers, 1988)
        // Ugcrit = 0.89+0.56*(CohesionSoil+CohesionRoot)
        // CohesionSoil is the soil cohesion (kPa), which is 9.806*Torvane value (kg/cm2)
        // CohesionRoot is the extra cohesion provided by plant roots
        // enter high CohesionSoil and CohesionTotal values (9999) for non-erodible surfaces
        // thus, Y becomes 0

     calc(" Y = mif(hardsurface gt 0, 0, Y) ");
     //VJ 080613 include non-erodible surfaces

     _spatial(REAL4, SettlingVelocity);
     calc(" SettlingVelocity = 2*(2650-1000)*9.80*sqr(D50/2000000)/(9*0.001) ");
         // settling velocity of sediment, in m/s!, at 20 C
         // 2650 = particle density kgm-3
         // 1000 = densityof water kgm-3
         // 9.80 = acceleration of gravity
         // 0.001 = viscosity of water
         // Settling Velocity = 0.00080903 m/s for D50 = 30 mu at 20 C
         // Stokes' Law) from "Soil Physics" (Marshall & Holmes, '79; p.24/5)
//VJ new from EUROSEM manual v2
     _spatial(REAL4, CGovers);
     calc(" CGovers = ((D50+5)/0.32)**-0.6 ");
        // calc(" CGovers = 0.015061+exp(-2.33860-0.014059*D50) ");
        // experimentally derived coefficient
        // depending on D50 (Govers, 1990; EUROSEM, 1992)
        // CSS fit from EUROSEM manual data: r2=0.9979
        // CGovers and DGovers are calculated from D50, and thus now spatial

//VJ new from EUROSEM manual v2
     _spatial(REAL4, DGovers);
     calc(" DGovers = ((D50+5)/300)**0.25 ");
        // calc(" DGovers = log10(2.431200+0.027716*D50) ");
        // experimentally derived coefficient
        // depending on D50 (Govers, 1990; EUROSEM, 1992)
        // CSS fit from EUROSEM manual data: r2=0.979


//     _nonspatial(REAL4, VCourant);
//      VCourant = DX/DTSEC;
//      if (VCourant < 0.5)
//         LisemWarning("DX/DT is smaller than 0.5 m/s! EFFECT: possible mass balance errors ADVISE: choose timestep DT (sec) <= pixel size DX (m)");
//ONZIN de courant vergelijking is: snelheid * dx/dt < 1
    /*
     _spatial(REAL4, VCount);
     calc("VCount = 0.0 ");
     _spatial(REAL4, VMax);
     calc(" VMax = 0.0 ");
        // check how many times flow rate is larger than DX/DTSEC
        // and constant variables, kinematic wave check
    */

}


