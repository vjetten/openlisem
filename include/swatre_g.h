#ifndef SWATRE_G_H
#define SWATRE_G_H


//---------------------------------------------------------------------------
/// SWATRE structure for soil profile information and swatre options

/** SWATRE: structure for soil profile information and swatre options.\n
This model is declared for normal soils, crustd soils, compacted soils and grass strips. \n
*/

typedef struct SOIL_MODEL {
   struct PIXEL_INFO  *pixel; //defined in swatre_p.h
  // double precision;
   double minDt;
  // double calibrationfactor;
  // bool geometric;
  // bool swatreBottomClosed;
  // long nrCells;
} SOIL_MODEL;
//---------------------------------------------------------------------------

#endif // SWATRE_G_H
