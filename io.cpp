#include "io.h"
#include <cassert>
#include <cstring>
#include <iostream>
#include <QFileInfo>
#include <gdal_priv.h>
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


auto close_gdal_dataset = [](GDALDataset* dataset) { GDALClose(dataset); };
using GDALDatasetPtr = std::unique_ptr<GDALDataset, decltype(
    close_gdal_dataset)>;


cTMap readRaster(
    QString const& pathName)
{
    // Open raster dataset and obtain some properties.
    GDALDatasetPtr dataset(static_cast<GDALDataset*>(GDALOpen(
        pathName.toAscii().constData(), GA_ReadOnly)), close_gdal_dataset);
    if(!dataset) {
        Error(QString("Map %1 cannot be opened.").arg(pathName));
    }

    int nr_bands = dataset->GetRasterCount();
    if(nr_bands == 0) {
        Error(QString("Map %1 does not contain any bands.").arg(pathName));
    }

    double transformation[6];
    dataset->GetGeoTransform(transformation);


    // Read the first raster band.
    GDALRasterBand* band{dataset->GetRasterBand(1)};
    assert(band);

    int const nr_rows{band->GetYSize()};
    int const nr_cols{band->GetXSize()};
    double const west{transformation[0]};
    double const north{transformation[3]};
    double const cell_size{transformation[1]};

    MaskedRaster<double> raster_data(nr_rows, nr_cols, north, west, cell_size);

    // All raster values are read into doubles. PCRaster value scales are not
    // taken into account.
    if(band->RasterIO(GF_Read, 0, 0, nr_cols, nr_rows, raster_data[0],
            nr_cols, nr_rows, GDT_Float64, 0, 0) != CE_None) {
        Error(QString("Raster band %1 cannot be read.").arg(pathName));
    }

    int hasNoDataValue{false};
    double noDataValue{band->GetNoDataValue(&hasNoDataValue)};
    if(hasNoDataValue) {
        raster_data.replace_with_mv(noDataValue);
    }

    return cTMap(std::move(raster_data), pathName);
}


void writePCRasterRaster(
    cTMap const& raster,
    QString pathName)
{
    // Create and configure CSF map.
    MapPtr csfMap{Rcreate(pathName.toAscii().constData(), raster.nrRows(),
        raster.nrCols(), CR_REAL4, VS_SCALAR, PT_YDECT2B, raster.west(),
        raster.north(), 0.0, raster.cellSize()), close_csf_map};

    if(!csfMap) {
        Error(QString("Dataset %1 cannot be created.").arg(pathName));
    }

    RuseAs(csfMap.get(), CR_REAL8);

    // Copy cells to write to new buffer.
    auto const& raster_data(raster.Data);
    std::unique_ptr<double[]> buffer{new double[raster_data.nr_cells()]};
    std::memcpy(buffer.get(), raster_data[0], sizeof(double) *
        raster_data.nr_cells());

    // Write cells from buffer to file.
    size_t nr_cells_written = RputSomeCells(csfMap.get(), 0,
        raster_data.nr_cells(), buffer.get());

    if(!nr_cells_written == raster_data.nr_cells()) {
        Error("rputsomecells write error with " + pathName);
    }
}


void writeGDALRaster(
    cTMap const& raster,
    QString const& pathName,
    GDALDriver& driver)
{
    // Create new dataset.
    int const nrRows{raster.nrRows()};
    int const nrCols{raster.nrCols()};
    int const nrBands{1};
    GDALDatasetPtr dataset{driver.Create(pathName.toAscii().constData(),
        nrCols, nrRows, nrBands, GDT_Float32, nullptr), close_gdal_dataset};

    if(!dataset) {
        Error(QString("Dataset %1 cannot be created.").arg(pathName));
    }

    MaskedRaster<double> const& raster_data{raster.Data};

    // Set some metadata.
    double transformation[]{
        raster_data.west(),
        raster_data.cell_size(),
        0.0,
        raster_data.north(),
        0.0,
        raster_data.cell_size()};
    dataset->SetGeoTransform(transformation);

    // PCRaster supports value scales, but other formats don't. We set the
    // value scale as a meta data item in the raster. If the format supports
    // setting meta data items, this allows for round tripping values scale
    // information back to the PCRaster format, in case the raster is
    // translated to PCRaster format later.
    dataset->SetMetadataItem("PCRASTER_VALUESCALE", "VS_SCALAR");

    // Write values to the raster band.
    auto band = dataset->GetRasterBand(1);

    band->SetNoDataValue(-FLT_MAX);

    if(band->RasterIO(GF_Write, 0, 0, nrCols, nrRows,
            const_cast<double*>(&raster_data.cell(0)),
            nrCols, nrRows, GDT_Float64, 0, 0) != CE_None) {
        Error(QString("Raster band %1 cannot be written.").arg(pathName));
    }
}


void writeRaster(
    cTMap const& raster,
    QString const& pathName,
    QString const& format)
{
    if(raster.nrRows() == 0 || raster.nrCols() == 0) {
        return;
    }

    if(pathName.isEmpty()) {
        ErrorString = "Cannot write file, file name empty";
        throw 1;
    }


    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName(
        format.toAscii().constData());

    if(!driver) {
        Error(QString("Format driver %1 not available.").arg(
            format.toAscii().constData()));
    }

    char** metadata{driver->GetMetadata()};
    bool driverSupportsCreate{CSLFetchBoolean(metadata, GDAL_DCAP_CREATE,
        FALSE) != FALSE};

    if(!driverSupportsCreate && format != "PCRaster") {
        Error(QString(
            "Format driver %1 cannot be used to create datasets.").arg(
                format.toAscii().constData()));
    }


    if(format == "PCRaster") {
        writePCRasterRaster(raster, pathName);
    }
    else {
        writeGDALRaster(raster, pathName, *driver);
    }
}


/// makes mapname if (name.map) or mapseries (name0000.001 to name0009.999)
void WriteMapSeries(
    cTMap const& raster,
    QString const& Dir,
    QString Name,
    int count,
    QString const& format)
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
    writeRaster(raster, path, format);
}
