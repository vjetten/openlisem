#include "io.h"
#include <QFileInfo>
#include "CsfMap.h"
#include "error.h"


/* for info:
typedef struct CSF_RASTER_HEADER
{
	 UINT2    valueScale;
	 UINT2    cellRepr;
	 CSF_VAR_TYPE minVal;
	 CSF_VAR_TYPE maxVal;
	 REAL8    xUL;
	 REAL8    yUL;
	 UINT4    nrRows;
	 UINT4    nrCols;
	 REAL8    cellSizeX;
	 REAL8    cellSizeY;
	 REAL8    angle;
} CSF_RASTER_HEADER;
*/


bool LoadFromFile(
    cTMap& raster)
{
    if(!QFileInfo(raster.PathName()).exists()) {
        Error(QString("Map %1 does not exist.").arg(raster.PathName()));
    }

    raster.Data = MaskedRaster<double>();

    // make map structure
    MapPtr m{Mopen(raster.PathName().toAscii().constData(), M_READ),
        close_csf_map};
    if(m == nullptr) {
        Error(QString("Map %1 cannot be opened.").arg(raster.PathName()));
    }

    auto header = m->raster;

    raster.Data = MaskedRaster<REAL8>(header.nrRows, header.nrCols,
        header.yUL, header.xUL, header.cellSize);

    if(!raster.created()) {
        return(false);
    }

    raster.setMapName(raster.PathName());

    RuseAs(m.get(), CR_REAL8); //RgetCellRepr(m));
    RgetSomeCells(m.get(), 0, raster.Data.nr_cells(), raster.Data[0]);

    /// KDJ: not needed ?! if (RgetCellRepr(m) == CR_REAL8)
    /// KDJ: not needed ?!   ResetMinMax();

    return true;
}


/// write a map to disk
void WriteMap(
    cTMap const& raster,
    QString Name)
{
    if(!raster.created()) {
        return;
    }

    if(Name.isEmpty()) {
        ErrorString = "Cannot write file, file name empty";
        throw 1;
    }

    // KDJ: Not needed?!: raster.ResetMinMax();

    // Make an array for output.
    std::unique_ptr<REAL4[]> Dt{new REAL4[raster.nrCols()]};

    /// KDJ: _MH.cellRepr = CR_REAL4;
    auto deleter = [](MAP* map) { Mclose(map); };
    std::unique_ptr<MAP, decltype(deleter)> out{Rcreate(
        Name.toAscii().constData(), raster.nrRows(),
        raster.nrCols(), CR_REAL4, VS_SCALAR, PT_YDECT2B, raster.west(),
        raster.north(), 0.0, raster.cellSize()), deleter};
    RuseAs(out.get(), CR_REAL4);

    for(long r=0; r < raster.nrRows(); ++r) {
        for(long c=0; c < raster.nrCols(); c++) {
            Dt[c] = static_cast<REAL4>(raster.Data[r][c]);
        }

        if(RputRow(out.get(), r, Dt.get()) != static_cast<UINT4>(raster.nrCols())) {
            ErrorString = "rputrow write error with" + Name;
            throw 1;
        }
    }
}


/// makes mapname if (name.map) or mapseries (name0000.001 to name0009.999)
void WriteMapSeries(
    cTMap const& raster,
    QString Dir,
    QString Name,
    int count)
{
    QString path;
    QFileInfo fi(Name);

    if(Name.indexOf(".") < 0) {
        QString nam, dig;

        nam = Name + "00000000";

        nam.remove(7, 10);
        dig = QString("%1").arg(count, 4, 10, QLatin1Char('0'));
        dig.insert(1, ".");
        Name = nam + dig;
    }

    path = Dir + Name;
    WriteMap(raster, path);
}
