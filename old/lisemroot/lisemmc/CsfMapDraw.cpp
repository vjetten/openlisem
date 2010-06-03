//---------------------------------------------------------------------------
#include <vcl.h>
#include <stdlib.h>
#include <math.h>
//---------------------------------------------------------------------------
#pragma hdrstop

#include "CsfMapDraw.h"
//#include "CsfMapLeg.h"
//#include "CsfMapEdit.h"

//---------------------------------------------------------------------
__fastcall TMapDraw::TMapDraw(TComponent* Owner)
    : TMap()
{
    coordx = NULL;
    coordy = NULL;
    CData = NULL;
    bmp = NULL;
    bmpleg = NULL;

    getmin = true;
    getmax = true;
    SetClassifyType();
    MakePalette();

//    ReversePalette();
    //!!!!!!!!!!1
    PaletteNr = 1;
    ShowGrid = false;
    ShowLegend = false;
    digitize = false;
    c1 = 0; r1 = 0; c2 = 0; r2 = 0;
    NrClasses = MAXCLASS;
    Title = "";
    UserMax = 0;
    DoUserMax = false;
}
//---------------------------------------------------------------------------
__fastcall TMapDraw::~TMapDraw()
{
     KillDrawMap();
}
//---------------------------------------------------------------------------
void TMapDraw::ReversePalette()
{
   RGB paltmp[33];
   int c, i;
   for (i = 0; i < 4; i++)
   {
     for (c=0; c<33; c++)
        paltmp[c] = pal[i][c];
     for (c=0; c<33; c++)
     {
        pal[i][c].r = paltmp[32-c].r;
        pal[i][c].g = paltmp[32-c].g;
        pal[i][c].b = paltmp[32-c].b;
     }
   }
}
//---------------------------------------------------------------------------
void TMapDraw::MakePalette()
{
    int i, c;
    RGB paltmp[33];

   /* red-yellow-blue */
   for (c=0; c<8; c++)
   {
      pal[0][c].r = 0;
      pal[0][c].g = 35*(1-exp(-pow(2.84/16*c,2)))+17+c;
      pal[0][c].b = max(63-pow((-0.22*(4*c-12)),2),0.0);
   }

   for (c=0; c<8; c++)
   {
      pal[0][c+8].r = min(8*c,63);
      pal[0][c+8].g = 35*(1-exp(-pow(2.84/16*(c+8),2)))+27;
      pal[0][c+8].b = max(63-pow((-0.22*(4*(c+8)-12)),2),0.0);
   }

   for (c=0; c<8; c++)
   {
      pal[0][c+16].r = 63;
      pal[0][c+16].g = min(67-4*c,63);
      pal[0][c+16].b = 0;
   }

   for (c=0; c<8; c++)
   {
      pal[0][c+24].r = 63-c;
      pal[0][c+24].g = 36-4*c;
      pal[0][c+24].b = 0;
   }

   /* blue to yellow, functions fitted with statistica
      through handmade points !!!*/
   for (c=0; c<32; c++)
   {
//      pal[1][c].r = min(50*(1-exp(-pow(2.0/32*c,5)))+9,63.0);
//      pal[1][c].g = min(50*(1-exp(-pow(3.0/32*c,2)))+17,63.0);
//      pal[1][c].b = max(60-pow((-0.17*(2.0*c-12)),2),0.0);

      pal[1][c].r = min(50*(1-exp(-pow(1.77/32*c,5)))+9,63.0);
      pal[1][c].g = min(45*(1-exp(-pow(2.84/32*c,2)))+17,63.0);
      pal[1][c].b = max(50-pow((-0.22*(2*c-12)),2),0.0);

//      pal[1][c].r = (int) 64 * (float)c/32.0;
//      pal[1][c].g = (int) 64 * (float)c/32.0;
//      pal[1][c].b = (int) 40 * (1-(float)c/32.0);
   }

   for (c=0; c<4; c++)
   {
      pal[1][(3-c)].r -= c*2;
      pal[1][(3-c)].g -= c*2;
      pal[1][(3-c)].b -= c*2;
   }
   for (c=24; c<32; c++)
   {
 //     pal[1][c].r -= (c-23);
 //     pal[1][c].g -= (c-23);
   }
     pal[1][33] = pal[1][32] = pal[1][31];
     for (c=0; c<33; c++)
        paltmp[c] = pal[1][c];
     for (c=0; c<33; c++)
     {
        pal[1][c].r = paltmp[32-c].r;
        pal[1][c].g = paltmp[32-c].g;
        pal[1][c].b = paltmp[32-c].b;
     }



   /* red to green */
   for (c=0; c<16; c++)
   {
      pal[2][c].r = 4*c;
      pal[2][c].g = 2*c+32;
      pal[2][c].b = 0;
   }
   for (c=0; c<16; c++)
   {
      pal[2][c+16].r = 63-c*0.5;
      pal[2][c+16].g = 62-3*c;
      pal[2][c+16].b = 0;
   }
   /*grey*/
   for (c=0; c<32; c++)
   {
      pal[3][c].r = 1.8*c+7;
      pal[3][c].g = 1.8*c+7;
      pal[3][c].b = 1.8*c+7;
   }

   for (i = 0; i < 4; i++)
   {
     for (c=0; c<33; c++)
     {
        pal[i][c].r *= 255.0/63.0;
        pal[i][c].g *= 255.0/63.0;
        pal[i][c].b *= 255.0/63.0;
     }
     pal[i][33] = pal[i][32] = pal[i][31];
   }

}
//---------------------------------------------------------------------------
void TMapDraw::GetMinMax()
{
    int r, c;

    MaxV = -9e19;
    MinV = 9e19;
    for(r = 0; r < nrRows; r++)
      for(c = 0; c < nrCols; c++)
        if (!IS_MV_REAL4(&(REAL4)Data[r][c]))
        {
            if (getmax)
               MaxV = max(MaxV, Data[r][c]);
            if (getmin)
               MinV = min(MinV, Data[r][c]);
        }
}
//---------------------------------------------------------------------------
void TMapDraw::GetMinMaxMAP()
{
    MAP *m = Mopen(MapName.c_str(), M_READ);
    REAL4 V;

    if(getmax)
    {
        RgetMaxVal(m, &V);
        MaxV = V;
    }

    if(getmin && MRH.valueScale != VS_DIRECTION)
    {
        RgetMinVal(m, &V);
        MinV = V;
    }

    if (m)
       Mclose(m);

}
//---------------------------------------------------------------------------
void TMapDraw::MakeClassTable()
{
    int i;
    float minm = 0.01;

    GetMinMax();
    if (DoUserMax)
       MaxV = UserMax;

 //   if (MinV==MaxV)
   //    NrClasses = 1;

    if(MRH.valueScale == VS_DIRECTION)
    {
        MinV = 0;
        MaxV = 360;
    }


    if (ClassType == CT_Lin)
    {
      for (i = 0; i < NrClasses; i++)
         Class[i] = MinV + i*(MaxV-MinV)/(float)NrClasses;
      Class[NrClasses] = MaxV;
    }
    else
    if (ClassType == CT_Log && MinV >= 0 && MaxV >= 0)
    {
       Class[0] = minm;
       for (i = 1; i < NrClasses; i++)
         Class[i] = pow(10,log10(minm) + i*(log10(MaxV)-log10(minm))/(float)NrClasses);
       Class[NrClasses] = MaxV;
    }
    else
    if (ClassType == CT_Log && MinV < 0 && MaxV > 0)
    {
       int Nrc2 = NrClasses/2;
       Class[0] = MinV;
       for (i = 1; i < Nrc2; i++)
         Class[i] = -pow(10,log10(-MinV) - i*(log10(-MinV)-log10(minm))/(float)Nrc2);
       Class[Nrc2] = 0;
//       for (i = Nrc2 + 1; i < NrClasses; i++)
       for (i = 1; i < Nrc2; i++)
         Class[i+Nrc2] = pow(10,log10(minm) + (i)*(log10(MaxV)-log10(minm))/(float)Nrc2);
       Class[NrClasses] = MaxV;
    }
    else
    if (ClassType == CT_Log && MinV < 0 && MaxV <= 0)
    {
       Class[0] = MinV;
       for (i = 1; i < NrClasses; i++)
         Class[i] = -pow(10,log10(-MinV) + i*(log10(minm)-log10(-MinV))/(float)NrClasses);
       Class[NrClasses] = -minm;
    }

}
//---------------------------------------------------------------------------
void TMapDraw::ScalarClassify()
{
    int r, c, i;

   // NrClasses = MAXCLASS;
    MakeClassTable();

    for(r = 0; r < nrRows; r++)
      for(c = 0; c < nrCols; c++)
        if (!IS_MV_REAL4(&Data[r][c]))
        {
          for (i = 0; i <= NrClasses; i++)
             if (Class[i] >= Data[r][c])
                break;
          CData[r][c] = i;
        }
        else
          CData[r][c] = -1;


    if (MRH.valueScale == VS_DIRECTION)
       for(r = 0; r < nrRows; r++)
         for(c = 0; c < nrCols; c++)
           if (!IS_MV_REAL4(&Data[r][c]) && Data[r][c] < 0)
             CData[r][c] = 0;
}
//---------------------------------------------------------------------------
void TMapDraw::GetMinMaxI()
{
    int r, c;

    for(r = 0; r < nrRows; r++)
      for(c = 0; c < nrCols; c++)
      {
        if (CData[r][c] != MV_INT4)
        {
            if (getmax)
               MaxV = max(MaxV, (float)CData[r][c]);
            if (getmin) MinV =
               min(MinV, (float)CData[r][c]);
        }
      }
}
//---------------------------------------------------------------------------
void ClassSort(int n, float *arr)
{
    int i,j;
    float a;

    for(j = 1; j <= n; j++)
    {
        a = arr[j];
        i = j-1;
        while(i > 0 && arr[i] > a)
        {
            arr[i+1] = arr[i];
            i--;
        }
        arr[i+1] = a;
    }
}
//---------------------------------------------------------------------------
void TMapDraw::GetMinMaxC()
{
    int r, c;

    if (getmax) MaxV = -9e19;
    if (getmin) MinV = 9e19;

    for(r = 0; r < nrRows; r++)
      for(c = 0; c < nrCols; c++)
        if (!IS_MV_REAL4(&Data[r][c]))
        {
            if (getmax)// && Data[r][c] > 0)
               MaxV = max(MaxV, Data[r][c]);
            if (getmin)// && Data[r][c] > 0)
               MinV = min(MinV, Data[r][c]);   //vj 080915 min was max in version 1.7 
        }

}
//---------------------------------------------------------------------------
void TMapDraw::OrdinalClassify()
{
    int r, c, i=0, j;

    GetMinMaxC();
    Class[0] = MinV;

    /* fill Class with all seperate classes */
    for(r = 0; r < nrRows; r++)
     for(c = 0; c < nrCols; c++)
      if (!IS_MV_REAL4(&Data[r][c]))
      {
          bool found = false;

          for(j = 0; j <= i; j++)
              if ((long)Data[r][c] == (long)Class[j]) found = true;

          if(!found && i < MAXCLASS)
          {
              i++;
              Class[i] = Data[r][c];
          }
      }

     if (i == 0) i++;
    NrClasses = min(i, MAXCLASS);
    if (NrClasses == 0)
       NrClasses = 1;

    /* sort classes */
    ClassSort(NrClasses, Class);

    /* reclass map for colors */
     for(r = 0; r < nrRows; r++)
      for(c = 0; c < nrCols; c++)
        if (!IS_MV_REAL4(&Data[r][c]))
        {
//            C[r][c] = 1+(((long)Data[r][c] % MAXCLASSC)+ MAXCLASSC ) % MAXCLASSC;
//           CData[r][c] = 1+(((long)Data[r][c] % NrClasses)+ NrClasses) % NrClasses;

          // for (i = 1; i <= NrClasses; i++)
            // if (Data[r][c] == Class[i])
               CData[r][c] = (int)Data[r][c] ;//i*(float)MAXCLASS/(float)NrClasses;//+1;//(int)Data[r][c];

        }       
        else
           CData[r][c] = -1;

}
//---------------------------------------------------------------------------
void TMapDraw::Classify()
{
   if (!Created)
      return;

//   if (MRH.valueScale == VS_LDD)
  //     NrClasses=8;

   if (MRH.valueScale == VS_SCALAR)
      ScalarClassify();

 //  if (MRH.valueScale == VS_NOMINAL || MRH.valueScale == VS_ORDINAL)
   //   OrdinalClassify();
      /*
 * NOTE new VS_* types must be different from CR_* values
 *      VS_BOOLEAN       0xE0 => 224
 *      VS_NOMINAL       0xE2 => 226
 *      VS_ORDINAL       0xF2 => 242
 *      VS_SCALAR        0xEB => 235
 *      VS_DIRECTION     0xFB => 251
 *      VS_LDD           0xF0 => 240
 */
}
//---------------------------------------------------------------------------
void TMapDraw::Coordinates()
{
    int r, c;
    int nrR = max(10, r2-r1);
    int nrC = max(10, c2-c1);

    //Dy, Dx =s number of screen pixels per gridcell

    //first look at the shape of bmp window and see if it can fit nrR or nrC
    if (bmp->Width > bmp->Height)
    {
       Dy = (double)bmp->Height/(double)nrR;
       Dx = Dy;
    }
    else
    {
       Dx = (double)bmp->Width/(double)nrC;
       Dy = Dx;
    }

    //now correct
    if (int(Dy*nrR) > bmp->Height)
    {
       Dy = (double)bmp->Height/(double)nrR;
       Dx = Dy;
    }
    if (int(Dx*nrC) > bmp->Width)
    {
       Dx = (double)bmp->Width/(double)nrC;
       Dy = Dx;
    }

     for (c = 0; c <= nrCols+MAXSTEP-1; c++)
        coordx[c] = c * Dx;
     for (r = 0; r <= nrRows+MAXSTEP-1; r++)
        coordy[r] = r * Dy;
}
//---------------------------------------------------------------------------
void TMapDraw::RrowCol2Coords(
		    double row,  /* Row number (relates to y position). */
		    double col,  /* Column number (relates to x position). */
		    double *x,   /* write-only. x co-ordinate */
		    double *y)   /* write-only. y co-ordinate */
{
	double cs     = RgiveCellSizeX();
//	double c      = m->angleCos;
//	double s      = m->angleSin;  // ignored
	double yRow   = cs * row;
	double xCol   = cs * col;
	double xCol_t = xCol;
	double yRow_t = yRow;

	*x = RgiveX0() + xCol_t;
	if (MgiveProjection() == PT_YINCT2B)
		*y = RgiveY0() + yRow_t;
	else  /* all other projections */
		*y = RgiveY0() - yRow_t;

}
//---------------------------------------------------------------------------
void TMapDraw::GetCellInfo(int X, int Y)   //X and Y are screen coordinates on the bitmap
{
   // get real coordinates from row and col
   double ra=r1+(double)Y/Dy, ca = c1+(double)X/Dx;
   double xa, ya;

   RrowCol2Coords(ra, ca, &xa, &ya);

   CellInfo.c = (int) ca;
   CellInfo.r = (int) ra;
   CellInfo.X = xa; //map coords
   CellInfo.Y = ya; //map coords

   CellInfo.MV = (CellInfo.r > nrRows-1 || CellInfo.c > nrCols-1);
   if (!CellInfo.MV)
      CellInfo.MV = IS_MV_REAL4(&Data[CellInfo.r][CellInfo.c]);
   if (!CellInfo.MV)
      CellInfo.V = Data[CellInfo.r][CellInfo.c];
   CellInfo.dx = Dx; //???????  necessary?
   CellInfo.dy = Dy; //???????
}
//---------------------------------------------------------------------------
void TMapDraw::DrawMap()
{
   int r, c, col;
   bmp->Canvas->Brush->Style = bsSolid;
   bmp->Canvas->Brush->Color = backgroundcolor;
   bmp->Canvas->FillRect(Rect(0,0,bmp->Width, bmp->Height));

   Coordinates();

   for (r = 0; r < nrRows; r++)
    for (c = 0; c < nrCols; c++)
    {
       if (CData[r][c] < 0)
          bmp->Canvas->Brush->Color = backgroundcolor;//clWhite;
      else
      {
          int col = (int)CData[r][c]*MAXCLASS/NrClasses;
           bmp->Canvas->Brush->Color =
                static_cast<TColor>RGB(pal[PaletteNr][col].r,pal[PaletteNr][col].g,pal[PaletteNr][col].b);
           bmp->Canvas->FillRect(Rect((c-c1)*Dx,(r-r1)*Dy,(c-c1+1)*Dx,(r-r1+1)*Dy));
       }
        if (ShowGrid)
        {
          bmp->Canvas->Pen->Width = 1;
          bmp->Canvas->Pen->Color = clBlack;
          bmp->Canvas->Pen->Style = psDot;
          bmp->Canvas->MoveTo((c-c1)*Dx,(r-r1)*Dy);
          bmp->Canvas->LineTo((c-c1)*Dx+1,(r-r1)*Dy+1);
        }

    }

}

//---------------------------------------------------------------------------
void TMapDraw::DrawMap1()
{
   int r, c;

   bmp->Canvas->Brush->Style = bsSolid;

   //wipe bmp
   bmp->Canvas->Brush->Color = clBlue;//backgroundcolor;
   bmp->Canvas->FillRect(Rect(0,0,bmp->Width, bmp->Height));

   //get correct Dx, Dy
   Coordinates();

   //adjust r2 and c2 to fill screen
   r2 = min(nrRows, bmp->Height/Dy + r1);
   c2 = min(nrCols, bmp->Width/Dx + c1);

   for (r = r1; r < r2; r++)
    for (c = c1; c < c2; c++)
    {
       if (CData[r][c] < 0)
          bmp->Canvas->Brush->Color = backgroundcolor;//clWhite;
      else
      {
        if (MRH.valueScale != VS_LDD)
        {
              if (MRH.valueScale == VS_SCALAR)
                    bmp->Canvas->Brush->Color = static_cast<TColor>
                     RGB(pal[PaletteNr][(int)CData[r][c]*MAXCLASS/(NrClasses)].r,
                         pal[PaletteNr][(int)CData[r][c]*MAXCLASS/(NrClasses)].g,
                         pal[PaletteNr][(int)CData[r][c]*MAXCLASS/(NrClasses)].b);
              bmp->Canvas->FillRect(Rect((c-c1)*Dx,(r-r1)*Dy,(c-c1+1)*Dx,(r-r1+1)*Dy));

              if (ShowGrid)
              {
                bmp->Canvas->Pen->Width = 1;
                bmp->Canvas->Pen->Color = clBlack;
                bmp->Canvas->Pen->Style = psDot;
                bmp->Canvas->MoveTo((c-c1)*Dx,(r-r1)*Dy);
                bmp->Canvas->LineTo((c-c1)*Dx+1,(r-r1)*Dy+1);
              }
        }
        else   //LDD
        {
            int i = 0, j = 0;

            bmp->Canvas->Brush->Color = backgroundcolor;
            bmp->Canvas->Pen->Width = 1;
            bmp->Canvas->Pen->Color = clBlack;
            bmp->Canvas->Pen->Style = psSolid;
            bmp->Canvas->MoveTo((c-c1+0.5)*Dx,(r-r1+0.5)*Dy);
            if (Data[r][c] == 1){ i = -1; j = +1;}
            if (Data[r][c] == 2){ i =  0; j = +1;}
            if (Data[r][c] == 3){ i = +1; j = +1;}
            if (Data[r][c] == 4){ i = -1; j =  0;}
            if (Data[r][c] == 5){ i =  0; j =  0;}
            if (Data[r][c] == 6){ i = +1; j =  0;}
            if (Data[r][c] == 7){ i = -1; j = -1;}
            if (Data[r][c] == 8){ i =  0; j = -1;}
            if (Data[r][c] == 9){ i = +1; j = -1;}
            bmp->Canvas->LineTo((c-c1+i+0.5)*Dx,(r-r1+j+0.5)*Dy);

            if (Data[r][c] == 5)
            {
               bmp->Canvas->Brush->Color = clBlack;
               bmp->Canvas->FillRect(Rect((c-c1+0.3)*Dx,(r-r1+0.3)*Dy,(c-c1+0.7)*Dx+2,(r-r1+0.7)*Dy));
            }

              if (ShowGrid)
              {
                bmp->Canvas->Pen->Width = 1;
                bmp->Canvas->Pen->Color = clBlack;
                bmp->Canvas->Pen->Style = psDot;
                bmp->Canvas->MoveTo((c-c1)*Dx,(r-r1)*Dy);
                bmp->Canvas->LineTo((c-c1)*Dx+1,(r-r1)*Dy+1);
              }
        } //LDD
      }// background color       
    }

    Screen->Cursor = crArrow;
}

//---------------------------------------------------------------------------
void TMapDraw::CreateDrawMap()
{
//   if (start)//Created)
  // {
     // bitmap to blit
     bmp = new Graphics::TBitmap();
     bmpleg = new Graphics::TBitmap();

     //data structure for display cell classes
     CData = new int*[nrRows];
     for(int r=0; r < nrRows; r++)
        CData[r] = new int[nrCols];
     //initialize view rows and cols
     c1 = 0;
     r1 = 0;
     c2 = nrCols;
     r2 = nrRows;

     coordx = new double[nrCols+MAXSTEP];
     coordy = new double[nrRows+MAXSTEP];
//   }
}
//---------------------------------------------------------------------------
void TMapDraw::KillDrawMap()
{
   if (Created)
   {
     if (bmp)
     {
        delete bmp;
        bmp = NULL;
     }
     if (bmpleg)
     {
        delete bmpleg;
        bmpleg = NULL;
     }
    if (coordx)
       delete[] coordx;
    if (coordy)
       delete[] coordy;
     coordx = NULL;
     coordy = NULL;
     if (CData)
     {
        for(int r=0; r < nrRows; r++)
           delete[] CData[r];
        delete[] CData;
        CData = NULL;
     }
   }
}
//---------------------------------------------------------------------------
void TMapDraw::DrawMapLegend()
{
   int i;
   int h = 24, w = 8;
   String S = "";

   bmpleg->Canvas->Brush->Style = bsSolid;
   bmpleg->Canvas->Font->Name = "Consolas";
   bmpleg->Canvas->Font->Size = 8;
   bmpleg->Canvas->Brush->Color = clBtnFace;
   bmpleg->Canvas->FillRect(Rect(0,0,bmpleg->Width,bmpleg->Height));

  // S = (String)ExtractFileName(MapName);
  // bmpleg->Canvas->TextOut(w,2,S);

   if (MRH.valueScale == VS_LDD)
   {
     h = 64;
     w = 48;
     int i, j, k = 0, kk;
     for (int j = -1; j <= 1; j++)
      for (int i = -1; i <= 1; i++)
      {
        int x = 0, y = 0, z = 15;
        k++;
        if (k % 2 == 1) z = 13;
        bmpleg->Canvas->MoveTo(w+i*z,h+j*z);
        bmpleg->Canvas->LineTo(w+i*4,h+j*4);
        kk = 10 - k;
        if (kk == 7) kk = 9;
        else
        if (kk == 9) kk = 7;
        else
        if (kk == 1) kk = 3;
        else
        if (kk == 3) kk = 1;
        else
        if (kk == 4) kk = 6;
        else
        if (kk == 6) kk = 4;
        S = kk;
        if (k == 2) y = -3;
        if (k == 8) y =  3;
        if (k == 4) x = -3;
        if (k == 6) x =  3;
        bmpleg->Canvas->TextOut(w-3+i*20+x,h-6+j*22+y,S);

        if (k == 1)
        {
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z);
          bmpleg->Canvas->LineTo(w+i*z+4,h+j*z);
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z);
          bmpleg->Canvas->LineTo(w+i*z,h+j*z+4);
        }
        if (k == 2)
        {
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z-2);
          bmpleg->Canvas->LineTo(w+i*z-3,h+j*z-2+3);
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z-2);
          bmpleg->Canvas->LineTo(w+i*z+3,h+j*z-2+3);
        }
        if (k == 3)
        {
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z);
          bmpleg->Canvas->LineTo(w+i*z-4,h+j*z);
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z);
          bmpleg->Canvas->LineTo(w+i*z,h+j*z+4);
        }
        if (k == 4)
        {
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z);
          bmpleg->Canvas->LineTo(w+i*z+3,h+j*z-3);
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z);
          bmpleg->Canvas->LineTo(w+i*z+3,h+j*z+3);
        }
        if (k == 6)
        {
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z);
          bmpleg->Canvas->LineTo(w+i*z-3,h+j*z-3);
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z);
          bmpleg->Canvas->LineTo(w+i*z-3,h+j*z+3);
        }
        if (k == 7)
        {
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z);
          bmpleg->Canvas->LineTo(w+i*z,h+j*z-3);
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z);
          bmpleg->Canvas->LineTo(w+i*z+3,h+j*z);
        }
        if (k == 8)
        {
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z+2);
          bmpleg->Canvas->LineTo(w+i*z-3,h+j*z-3+2);
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z+2);
          bmpleg->Canvas->LineTo(w+i*z+3,h+j*z-3+2);
        }
        if (k == 9)
        {
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z);
          bmpleg->Canvas->LineTo(w+i*z,h+j*z-3);
          bmpleg->Canvas->MoveTo(w+i*z,h+j*z);
          bmpleg->Canvas->LineTo(w+i*z-3,h+j*z);
        }
      }
   }
   else
   {
       char buff[255];
       int step = 4, xtra = 0;

      if (MRH.valueScale == VS_NOMINAL || MRH.valueScale == VS_ORDINAL)
      {
         step = 12;
         xtra = 1;
      }

       h = 32;
       int j = NrClasses;
       for (i = 0; i <= NrClasses; i++)
       {
           if (i < NrClasses+xtra)
           {
             bmpleg->Canvas->Brush->Color = static_cast<TColor>
                   RGB(pal[PaletteNr][(int)j*MAXCLASS/NrClasses].r,
                       pal[PaletteNr][(int)j*MAXCLASS/NrClasses].g,
                       pal[PaletteNr][(int)j*MAXCLASS/NrClasses].b);
             bmpleg->Canvas->FillRect(Rect(w,h+i*(step+1),32,h+(i+1)*(step+1)));
           }

           if (MRH.valueScale == VS_NOMINAL || MRH.valueScale == VS_ORDINAL)
           {
                  sprintf(buff," %g",Class[i]);
                  S = "";
                  S = buff;
                  bmpleg->Canvas->MoveTo(32, h-1+i*(step+1));
                  bmpleg->Canvas->LineTo(36, h-1+i*(step+1));
                  bmpleg->Canvas->Brush->Color = clBtnFace;
                  bmpleg->Canvas->TextOut(32+w,h+i*(step+1),S);
           }
           else
           if (NrClasses > 1)
             if (i % step == 0)
             {
                sprintf(buff," %.4g",Class[i]);
                S = "";
                S = buff;
                bmpleg->Canvas->MoveTo(32, h-1+j*(step+1));
                bmpleg->Canvas->LineTo(36, h-1+j*(step+1));
                bmpleg->Canvas->Brush->Color = clBtnFace;
                bmpleg->Canvas->TextOut(32+w,h-9+j*(step+1),S);
             }
           j--;
         }
         bmpleg->Canvas->TextOut(w,h-24,Title);
         bmpleg->Canvas->Brush->Style = bsClear;
         bmpleg->Canvas->Rectangle(w,h-1,32,h+(NrClasses+xtra)*(step+1));

    }
}
//---------------------------------------------------------------------------
//void TMapDraw::ReadDrawMap(AnsiString Name)
void TMapDraw::LoadFromFile(MEM_HANDLE *M, bool start)
{
    // if old maps exist kill them first
   if (start)
      KillDrawMap();

   TMap::LoadFromFile(M, start);

   if (start)
      CreateDrawMap();
   Created = true;

}
//---------------------------------------------------------------------------


