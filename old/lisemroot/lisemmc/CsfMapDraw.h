//---------------------------------------------------------------------------
#ifndef CsfMapDrawH
#define CsfMapDrawH
//---------------------------------------------------------------------------

#include <SysUtils.hpp>
#include <Controls.hpp>
#include <Classes.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>

#include "CsfMap.h"

//---------------------------------------------------------------------------

#define MAXCLASS 32
#define MAXSTEP 8

#define CT_Lin 0
//#define CT_Exp 1
#define CT_Log 1
#define CT_Log2 2

//typedef struct _Legend {
//  int Left, Top, Width, Height;
//  bool show;
//  Graphics::TBitmap *bmp;
//} Legend;

typedef struct _CInfo {
    int r, c;
    double X, Y;
    double dx, dy;
    REAL4 V;
    bool MV;
    } CInfo;

typedef struct _RGB { int r,g,b;}RGB;

class TMapDraw : public TMap
{
protected:
    void GetMinMax();
    void GetMinMaxMAP();
    void MakeClassTable();
    void GetMinMaxI();
    void GetMinMaxC();
    void OrdinalClassify();
    void ScalarClassify();

    float Class[65];
    bool isScalar;

//    INT4 **CData;
    int **CData;

    double Dx, Dy;
    double *coordx;
    double *coordy;
    void Coordinates();
    void MakePalette();
    void CreateDrawMap();
    RGB pal[5][34];
public:
    void KillDrawMap();
    bool digitize;
    int NrClasses;
    Graphics::TBitmap *bmp;
    Graphics::TBitmap *bmpleg;
    int PaletteNr;
    char ClassType;
    bool getmin, getmax;
    float MinV, MaxV;
    bool ShowGrid;
    bool ShowLegend;
    int c1, r1, c2, r2; //display rows and cols
    AnsiString Title;
    bool DoUserMax;
    double UserMax;


    CInfo CellInfo;
    TColor backgroundcolor;
    void GetCellInfo(int X, int Y);
    void Classify();
    void SetClassifyType(char type = 0, int nrcl = 32)
         {ClassType = type; NrClasses = nrcl;}
    void ReversePalette();

    void DrawMapLegend();
    void DrawMap();
    void DrawMap1();
    void LoadFromFile(MEM_HANDLE *M, bool start);
    void TMapDraw::RrowCol2Coords(
		    double row,  /* Row number (relates to y position). */
		    double col,  /* Column number (relates to x position). */
		    double *x,   /* write-only. x co-ordinate */
		    double *y);   /* write-only. y co-ordinate */


    __fastcall TMapDraw(TComponent* Owner);
    __fastcall ~TMapDraw();
};

#endif
