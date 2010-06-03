//VJ 090520
//setting up sub-pixel gridcell surfaces

      _spatial(REAL4, WheelWidthDX);
      _spatial(REAL4, StoneWidthDX);
      _spatial(REAL4, SoilWidthDX);

     // LISEM discrimates the following surfaces within DX:
     // - ChannelWidthUpDX: upper width of channels in m
     // - RoadWidthDX: width of roads in m
     // - WheelWidthDX: width of wheeltracks in m (based on wheelwid.map)
     // - SoilWidthDX: width of soil surface in m (based on stonefrc.map)
     // - StoneWidthDX: width of stone surface in m.
     // furthermore: LISEM deals with:
     // - crusted surfaces (based on crustfrc.map), for infiltration only

//VJ 030701 changed definition: avoid confusion: wheelwidth, gullywidth, channelwidth and roadwidth
// are all zero when not chosen or present !!

      calc(" WheelWidthDX = 0");
      if (SwitchWheelPresent) //wheeltracks as channels, no gullies!
      {
         calc(" WheelWidthDX = cover(max(DX-RoadWidthDX-ChannelWidthUpDX,0)/DX*WheelWidth,0) ");
      }
      calc(" StoneWidthDX = 0");// max(DX-ChannelWidthUpDX-GullyWidthDX-RoadWidthDX-WheelWidthDX,0)*StoneFraction ");
      calc(" SoilWidthDX = max(DX-ChannelWidthUpDX-GullyWidthDX-RoadWidthDX-WheelWidthDX-StoneWidthDX,0) ");


            /* old junk, MOGELIJK WEL EEN GOED IDEE!
           // Ingrid: Indien deel van gridcel bestaat uit Road of Channel wordt het
           //         aantal Wheeltracks berekend uitgaande van het 'interrill'-deel van het oppervlak
           calc(" WheelWidthDX = WheelWidth/WheelNumber "); // WheelWidth voor 1 wheeltrack
           calc(" WheelNumber = roundoff(cover(max(DX-ChannelWidthUpDX-RoadWidthDX,0)*WheelNumber/DX,0)) ");
           // Totale wheelwidth gecorrigeerd voor cellen met roads en channels
           calc(" WheelWidthDX = WheelNumber*WheelWidthDX ");
           */


 