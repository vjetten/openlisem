#include "io.h"
#include <cassert>
#include <cstring>
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


cTMap readRaster(
    QString const& pathName)
{
    if(!QFileInfo(pathName).exists()) {
        Error(QString("Map %1 does not exist.").arg(pathName));
    }

    MapPtr csfMap{Mopen(pathName.toAscii().constData(), M_READ), close_csf_map};
    if(!csfMap) {
        Error(QString("Map %1 cannot be opened.").arg(pathName));
    }
    RuseAs(csfMap.get(), CR_REAL8);

    auto header = csfMap->raster;
    MaskedRaster<double> raster_data(header.nrRows, header.nrCols,
        header.yUL, header.xUL, header.cellSize);
    RgetSomeCells(csfMap.get(), 0, raster_data.nr_cells(), raster_data[0]);

    return cTMap(std::move(raster_data), pathName);
}


void writeRaster(
    cTMap const& raster,
    QString pathName)
{
    if(raster.nrRows() == 0 || raster.nrCols() == 0) {
        return;
    }

    if(pathName.isEmpty()) {
        ErrorString = "Cannot write file, file name empty";
        throw 1;
    }

    // Create and configure CSF map.
    MapPtr csfMap{Rcreate(pathName.toAscii().constData(), raster.nrRows(),
        raster.nrCols(), CR_REAL4, VS_SCALAR, PT_YDECT2B, raster.west(),
        raster.north(), 0.0, raster.cellSize()), close_csf_map};
    RuseAs(csfMap.get(), CR_REAL8);

    // Copy cells to write to new buffer.
    auto const& raster_data(raster.Data);
    std::unique_ptr<double[]> buffer{new double[raster_data.nr_cells()]};
    std::memcpy(buffer.get(), raster_data[0], 8 * raster_data.nr_cells());

    // Write cells from buffer to file.
    size_t nr_cells_written = RputSomeCells(csfMap.get(), 0,
        raster_data.nr_cells(), buffer.get());

    if(!nr_cells_written == raster_data.nr_cells()) {
        Error("rputsomecells write error with " + pathName);
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
    writeRaster(raster, path);
}
