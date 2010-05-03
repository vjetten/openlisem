//---------------------------------------------------------------------------
#include <stdlib.h>
#include <math.h>
#include <qfile.h>
#include <qdir.h>
//---------------------------------------------------------------------------
#include "CsfMap.h"
#include "error.h"


//---------------------------------------------------------------------------
cTMap::cTMap()
{
    Data = NULL;
    Created = false;
    nrRows = 0;
    nrCols = 0;
    MapName = "";
    PathName = "";  
}
//---------------------------------------------------------------------------
cTMap::~cTMap()
{
    KillMap();
}
//---------------------------------------------------------------------------
void cTMap::KillMap()
{
    if (Data && nrRows > 0)
    {
       for (int r=0; r < nrRows; r++)
            delete[] Data[r];
       delete[] Data;
       Data = NULL;
    }
    Created = false;
}
//---------------------------------------------------------------------------
void cTMap::GetMapHeader(QString Name)
{
   MAP *m = Mopen(Name.toLatin1(), M_READ);
   MH = m->raster;

   // easier to separate these
   nrRows = (int)MH.nrRows;
   nrCols = (int)MH.nrCols;

   Mclose(m);
}
//---------------------------------------------------------------------------
void cTMap::CreateMap(QString Name)
{
     // if exist delete main map data structure and set Created to false
     if (Created)
        KillMap();

     // now get the header for nrrows and nrcols
     GetMapHeader(Name);

//     if (MH.cellRepr == CR_REAL4)
//     {
        Data = new REAL4*[nrRows];
        for(int r=0; r < nrRows; r++)
           Data[r] = new REAL4[nrCols];
//     }

     if (Data == NULL)
        return;
     Created = true;
}
//---------------------------------------------------------------------------
bool cTMap::LoadFromFile()
{
    MAP *m;
    //QFile fff(PathName);
    QFileInfo fi(PathName);

    if (!fi.exists())
      return(false);
    //fff.close();

    // make map structure
    CreateMap(PathName);

    if (!Created)
       return(false);

    MapName = PathName;

    m = Mopen(MapName.toLatin1(), M_READ);
    if (!m)
       return(false);
    RuseAs(m, CR_REAL4); //RgetCellRepr(m));//CR_REAL4);
    for(int r=0; r < nrRows; r++)
       RgetSomeCells(m, (UINT4)r*nrCols, (UINT4)nrCols, Data[r]);

    if (RgetCellRepr(m) == CR_REAL4)
       ResetMinMax();

    Mclose(m);
    return(true);
}
//---------------------------------------------------------------------------
void cTMap::_MakeMap(cTMap *dup, REAL4 value)
{
     if (dup == NULL)
        return;
     //MH = dup->MH;
     MH.minVal = value;
     MH.maxVal = value;
     nrRows = dup->nrRows;
     nrCols = dup->nrCols;
     MH.nrRows = nrRows;
     MH.nrCols = nrCols;
     MH.valueScale = dup->MH.valueScale;
     MH.cellRepr = dup->MH.cellRepr;
     MH.xUL = dup->MH.xUL;
	  MH.yUL = dup->MH.yUL;
     MH.cellSizeX  = dup->MH.cellSizeX;
     MH.cellSizeY  = dup->MH.cellSizeX;
     MH.angle = dup->MH.angle;

     Data = new REAL4*[nrRows];
     for(int r=0; r < nrRows; r++)
        Data[r] = new REAL4[nrCols];

     if (Data == NULL)
        return;

     for(int r = 0; r < nrRows; r++)
	 SetMemMV(Data[r],nrCols,CR_REAL4);

     for(int r=0; r < nrRows; r++)
      for(int c=0; c < nrCols; c++)
      if (!IS_MV_REAL4(&dup->Data[r][c]))
      {
          Data[r][c] = value;
      }

      Created = true;
}
//---------------------------------------------------------------------------
void cTMap::ResetMinMax(void)
{
     REAL8 minv = 1e30, maxv = -1e30;

     if (!Created)
        return;

      for(int r=0; r < nrRows; r++)
        for(int c=0; c < nrCols; c++)
        if (!IS_MV_REAL4(&Data[r][c]))
        {
           if (maxv < Data[r][c]) maxv = Data[r][c];
           if (minv > Data[r][c]) minv = Data[r][c];
        }

     MH.minVal = minv;
     MH.maxVal = maxv;
}
//---------------------------------------------------------------------------
void cTMap::WriteMap(QString Name)
{
    MAP *out;
    int r;//, c;
    REAL4 *buf;
    REAL8 *buf8;
    UINT1 *bufi1;
    INT4  *bufi4;

    if (!Created)
      {
        return;
      }
    MH.cellRepr = CR_REAL4; 
    if (MH.cellRepr == CR_REAL4)    
    ResetMinMax();

   out = Rcreate(Name.toLatin1(),nrRows, nrCols, (CSF_CR)MH.cellRepr, VS_SCALAR,
                 PT_YINCT2B, MH.xUL, MH.yUL, MH.angle, MH.cellSizeX);
    RputAllMV(out);

    if (out == NULL)
    {
      ErrorString = "Cannot write file: " + Name;
      throw 4;
    }
    //return;

    switch (MH.cellRepr)
    {
      case CR_REAL4 :
        buf = (REAL4 *)Rmalloc(out, nrCols);
        for(r=0; r < nrRows; r++)
        {
           memcpy(buf, Data[r], sizeof(REAL4)*nrCols);
           RputRow(out, r, buf);
        }
        free(buf);
        break;
      case CR_REAL8 :
        buf8 = (REAL8 *)Rmalloc(out, nrCols);
        for(r=0; r < nrRows; r++)
        {
           memcpy(buf8, Data[r], sizeof(REAL8)*nrCols);
           RputRow(out, r, buf8);
        }
        free(buf8);
        break;
      case CR_INT4 :
        bufi4 = (INT4 *)Rmalloc(out, nrCols);
        for(r=0; r < nrRows; r++)
        {
           memcpy(bufi4, Data[r], sizeof(INT4)*nrCols);
           RputRow(out, r, bufi4);
        }
        free(bufi4);
        break;
      case CR_UINT1 :
        bufi1 = (UINT1 *)Rmalloc(out, nrCols);
        for(r=0; r < nrRows; r++)
        {
           memcpy(bufi1, Data[r], sizeof(UINT1)*nrCols);
           RputRow(out, r, bufi1);
        }
        free(bufi1);
        break;
      }

     Mclose(out);
}
//---------------------------------------------------------------------------
void cTMap::WriteMapSeries(QString Dir, QString Name, int count)
{
    QString path;
    QString nam, dig;
    QFileInfo fi(Name);

    if (Name.indexOf(".") < 0)
    {

    	nam = Name + "00000000";

    	//    nam = fi.baseName();
    	nam.remove(8, 80);
    	dig.setNum(count);

    	if (count > 999)
    	{
    		nam.remove(8,1);
    		nam = nam + dig;
    		nam.insert(9,QString("."));
    	}
    	else
    		if (count > 99)
    			nam = nam + "." + dig;
    		else
    			if (count > 9)
    				nam = nam + "." + "0" + dig;
    			else
    				nam = nam + "." + "00" + dig;
    }
    else
       nam = Name;

    path = Dir + nam;
    WriteMap(path);
}


