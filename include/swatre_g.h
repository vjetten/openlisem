#ifndef SWATRE_G_H
#define SWATRE_G_H


//---------------------------------------------------------------------------
/** SWATRE: structure containing the 3D soil model. This model is declared for normal soils,\n
crustd soils, compacted soils and grass strips. PIXEL_INFO is defined in
  swatre_p.h (the private declarations)
*/
typedef struct SOIL_MODEL {
   struct PIXEL_INFO  *pixel;
   double precision;
   double minDt;
   double calibrationfactor;
   bool geometric;
   bool swatreBottomClosed;
   long nrCells;
} SOIL_MODEL;
//---------------------------------------------------------------------------

#endif // SWATRE_G_H
