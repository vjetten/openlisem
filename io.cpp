/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/
#include "io.h"
#include <cassert>
#include <cstring>
#include <iostream>
#include <memory>
#include <QFileInfo>
#include <gdal_priv.h>

#include "csf.h"
#include "CsfMap.h"
#include "CsfRGBMap.h"
#include "error.h"
#include "QDebug"


//! Function to close a CSF MAP.
auto close_csf_map = [](MAP* map) { Mclose(map); };

//! Auto-ptr type for CSF MAPs.
using MapPtr = std::unique_ptr<MAP, decltype(close_csf_map)>;

//! Function to close a GDAL GDALDataset.
auto close_gdal_dataset = [](GDALDataset* dataset) { GDALClose(dataset); };

//! Auto-ptr type for GDAL GDALDatasets.
using GDALDatasetPtr = std::unique_ptr<GDALDataset, decltype(
    close_gdal_dataset)>;


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


/*!
    @brief      Return whether raster @a pathName can be opened for reading.
*/
bool rasterCanBeOpenedForReading(
    QString const& pathName)
{
    GDALDatasetPtr dataset(static_cast<GDALDataset*>(GDALOpen(
        pathName.toLatin1().constData(), GA_ReadOnly)), close_gdal_dataset);
    bool result{dataset};

    return result;
}


/*!
    @brief      Read raster @a pathName and return the result.
*/
cTMap readRaster(
    QString const& pathName)
{
    // Open raster dataset and obtain some properties.
    GDALDatasetPtr dataset(static_cast<GDALDataset*>(GDALOpen(
        pathName.toLatin1().constData(), GA_ReadOnly)), close_gdal_dataset);
    if(!dataset) {
        Error(QString("Map %1 cannot be opened.").arg(pathName));
    }

    int nr_bands = dataset->GetRasterCount();
    if(nr_bands == 0) {
        Error(QString("Map %1 does not contain any bands.").arg(pathName));
    }

    double transformation[6];
    dataset->GetGeoTransform(transformation);

    QString projection{dataset->GetProjectionRef()};


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

    return cTMap(std::move(raster_data), projection, pathName);
}

/*!
    @brief      Read raster @a pathName and return the result.
*/
cTRGBMap *readRasterImage(
    QString const& pathName)
{
    // Open raster dataset and obtain some properties.
    GDALDatasetPtr dataset(static_cast<GDALDataset*>(GDALOpen(
        pathName.toLatin1().constData(), GA_ReadOnly)), close_gdal_dataset);
    if(!dataset) {
        Error(QString("Map %1 cannot be opened.").arg(pathName));
    }

    int nr_bands = dataset->GetRasterCount();
    if(nr_bands == 0) {
        Error(QString("Map %1 does not contain any bands.").arg(pathName));
    }

    double transformation[6];
    dataset->GetGeoTransform(transformation);

    QString projection{dataset->GetProjectionRef()};


    // Read the first raster band.
    GDALRasterBand* band{dataset->GetRasterBand(1)};
    assert(band);

    int const nr_rows{band->GetYSize()};
    int const nr_cols{band->GetXSize()};
    double const west{transformation[0]};
    double const north{transformation[3]};
    double const cell_size{transformation[1]};

    //qDebug() << "loading image " << pathName<< west << north << cell_size << nr_rows << nr_cols << nr_bands << band->GetRasterDataType();


    if(nr_bands == 1 || nr_bands == 2)
    {
        MaskedRaster<float> raster_data(nr_rows, nr_cols, north, west, cell_size);

        if(band->RasterIO(GF_Read, 0, 0, nr_cols, nr_rows, raster_data[0],
                nr_cols, nr_rows, GDT_Float32, 0, 0) != CE_None) {
            Error(QString("Raster band %1 cannot be read.").arg(pathName));
        }
        int s = 0;
        int s2 = 0;
        double min = band->GetMinimum(&s);
        double max = band->GetMaximum(&s2);

        int hasNoDataValue{false};
        double noDataValue{band->GetNoDataValue(&hasNoDataValue)};
        if(hasNoDataValue) {
            raster_data.replace_with_mv(noDataValue);
        }


        double av = 0;
        int n = 0;
        for(int i =0; i < nr_rows; i++)
        {
            for(int j = 0; j < nr_cols; j++)
            {
                if(!pcr::isMV(raster_data[i][j]))
                {
                     av += raster_data[i][j];
                     n++;
                }
            }
        }
        av = av/n;

        MaskedRaster<char> raster_data_int(nr_rows, nr_cols, north, west, cell_size);

        for(int i =0; i < nr_rows; i++)
        {
            for(int j = 0; j < nr_cols; j++)
            {
                raster_data_int[i][j] = char(int(std::min(255.0,std::max(1.0,(255.0*(raster_data[i][j] - min)/(std::min(2.0*av,max)-min))))));
            }
        }

        return new cTRGBMap(std::move(raster_data_int),projection,pathName);

    }else if(nr_bands > 2)
    {

        MaskedRaster<char> raster_data_int1(nr_rows, nr_cols, north, west, cell_size);
        MaskedRaster<char> raster_data_int2(nr_rows, nr_cols, north, west, cell_size);
        MaskedRaster<char> raster_data_int3(nr_rows, nr_cols, north, west, cell_size);


        {
            MaskedRaster<float> raster_data(nr_rows, nr_cols, north, west, cell_size);

            if(band->RasterIO(GF_Read, 0, 0, nr_cols, nr_rows, raster_data[0],
                    nr_cols, nr_rows, GDT_Float32, 0, 0) != CE_None) {
                Error(QString("Raster band %1 cannot be read.").arg(pathName));
            }
            int s = 0;
            int s2 = 0;


            double min = band->GetMinimum(&s);
            double max = band->GetMaximum(&s2);

            int hasNoDataValue{false};
            double noDataValue{band->GetNoDataValue(&hasNoDataValue)};
            if(hasNoDataValue) {
                raster_data.replace_with_mv(noDataValue);
            }

            double av = 0;
            int n = 0;
            for(int i =0; i < nr_rows; i++)
            {
                for(int j = 0; j < nr_cols; j++)
                {
                    if(!pcr::isMV(raster_data[i][j]))
                    {
                         av += raster_data[i][j];
                         n++;
                    }
                }
            }
            av = av/n;

            //qDebug() << av << min << max;
            for(int i =0; i < nr_rows; i++)
            {
                for(int j = 0; j < nr_cols; j++)
                {
                    raster_data_int1[i][j] = char(int(std::min(255.0,std::max(1.0,255.0*(raster_data[i][j] - min)/(std::min(2.0*av,max)-min)))));
                }
            }
        }

        band = dataset->GetRasterBand(2);
        assert(band);

        if(!(band->GetYSize() == nr_rows && band->GetXSize() == nr_cols))
        {
            Error(QString("band 2 is not identical in size to band 1.").arg(pathName));
        }

        {
            MaskedRaster<float> raster_data(nr_rows, nr_cols, north, west, cell_size);

            if(band->RasterIO(GF_Read, 0, 0, nr_cols, nr_rows, raster_data[0],
                    nr_cols, nr_rows, GDT_Float32, 0, 0) != CE_None) {
                Error(QString("Raster band %1 cannot be read.").arg(pathName));
            }
            int s = 0;
            int s2 = 0;
            double min = band->GetMinimum(&s);
            double max = band->GetMaximum(&s2);

            int hasNoDataValue{false};
            double noDataValue{band->GetNoDataValue(&hasNoDataValue)};
            if(hasNoDataValue) {
                raster_data.replace_with_mv(noDataValue);
            }

            double av = 0;
            int n = 0;
            for(int i =0; i < nr_rows; i++)
            {
                for(int j = 0; j < nr_cols; j++)
                {
                    if(!pcr::isMV(raster_data[i][j]))
                    {
                         av += raster_data[i][j];
                         n++;
                    }
                }
            }
            av = av/n;

            for(int i =0; i < nr_rows; i++)
            {
                for(int j = 0; j < nr_cols; j++)
                {
                    raster_data_int2[i][j] = char(int(std::min(255.0,(255.0*(raster_data[i][j] - min)/(std::min(2.0*av,max)-min)))));
                }
            }
        }

        band = dataset->GetRasterBand(3);
        assert(band);

        if(!(band->GetYSize() == nr_rows && band->GetXSize() == nr_cols))
        {
            Error(QString("band 3 is not identical in size to band 1.").arg(pathName));
        }

        {
            MaskedRaster<float> raster_data(nr_rows, nr_cols, north, west, cell_size);

            if(band->RasterIO(GF_Read, 0, 0, nr_cols, nr_rows, raster_data[0],
                    nr_cols, nr_rows, GDT_Float32, 0, 0) != CE_None) {
                Error(QString("Raster band %1 cannot be read.").arg(pathName));
            }
            int s = 0;
            int s2 = 0;
            double min = band->GetMinimum(&s);
            double max = band->GetMaximum(&s2);


            int hasNoDataValue{false};
            double noDataValue{band->GetNoDataValue(&hasNoDataValue)};
            if(hasNoDataValue) {
                raster_data.replace_with_mv(noDataValue);
            }

            double av = 0;
            int n = 0;
            for(int i =0; i < nr_rows; i++)
            {
                for(int j = 0; j < nr_cols; j++)
                {
                    if(!pcr::isMV(raster_data[i][j]))
                    {
                         av += raster_data[i][j];
                         n++;
                    }
                }
            }
            av = av/n;


            for(int i =0; i < nr_rows; i++)
            {
                for(int j = 0; j < nr_cols; j++)
                {
                    raster_data_int3[i][j] = char(int((std::min(255.0,255.0*(raster_data[i][j] - min)/(std::min(2.0*av,max)-min)))));
                }
            }
        }

        return new cTRGBMap(std::move(raster_data_int1),std::move(raster_data_int2),std::move(raster_data_int3),projection,pathName);
    }
    return nullptr;
}


void writePCRasterRaster(
    cTMap const& raster,
    QString pathName)
{
    // Create and configure CSF map.
    MapPtr csfMap{Rcreate(pathName.toLatin1().constData(), raster.nrRows(),
        raster.nrCols(), CR_REAL4, VS_SCALAR, PT_YDECT2B, raster.west(),
        raster.north(), 0.0, raster.cellSize()), close_csf_map};

    if(!csfMap) {
        Error(QString("Dataset %1 cannot be created.").arg(pathName));
    }

    RuseAs(csfMap.get(), CR_REAL8);

    // Copy cells to write to new buffer.
    auto const& raster_data(raster.data);
    std::unique_ptr<double[]> buffer{new double[raster_data.nr_cells()]};
    std::memcpy(buffer.get(), raster_data[0], sizeof(double) *
        raster_data.nr_cells());

    // Write cells from buffer to file.
    size_t nr_cells_written = RputSomeCells(csfMap.get(), 0,
        raster_data.nr_cells(), buffer.get());

    if(nr_cells_written != raster_data.nr_cells()) {
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
    GDALDatasetPtr dataset{driver.Create(pathName.toLatin1().constData(),
        nrCols, nrRows, nrBands, GDT_Float32, nullptr), close_gdal_dataset};

    if(!dataset) {
        Error(QString("Dataset %1 cannot be created.").arg(pathName));
    }

    MaskedRaster<double> const& raster_data{raster.data};

    // Set some metadata.
    double transformation[]{
        raster_data.west(),
        raster_data.cell_size(),
        0.0,
        raster_data.north(),
        0.0,
        raster_data.cell_size()};
    dataset->SetGeoTransform(transformation);

    dataset->SetProjection(raster.projection().toLatin1().constData());

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


/*!
    @brief      Write raster @a raster to @a pathName using format driver
                @a format.
    @param      raster Raster to write.
    @param      pathName Name of dataset to write.
    @param      format Name of driver to use for writing. Not that only drivers
                that implement the Create() method can be used to write
                rasters. (In case of the PCRaster driver an exception is made.
                Its driver doesn't implement Create() yet, but we handle
                saving to PCRaster raster format ourselves.)
*/
void writeRaster(
    cTMap const& raster,
    QString const& pathName,
    QString const& format)
{
//    if(raster.nrRows() == 0 || raster.nrCols() == 0) {
//        return;
//    }

//    if(pathName.isEmpty()) {
//        ErrorString = "Cannot write file, file name empty";
//        throw 1;
//    }

    if (format == "PCRaster") {
        // OK, until PCRaster supports Create(), we'll handle writing to
        // PCRaster format ourselves. Work is underway to add support
        // for Create() to the GDAL PCRaster driver.
        writePCRasterRaster(raster, pathName);
    } else {
        GDALDriver* driver = GetGDALDriverManager()->GetDriverByName(
            format.toLatin1().constData());

        if(!driver) {
            Error(QString("Format driver %1 not available.").arg(
                format.toLatin1().constData()));
        }

        char** metadata{driver->GetMetadata()};
        bool driverSupportsCreate{CSLFetchBoolean(metadata, GDAL_DCAP_CREATE, FALSE) != FALSE};
        if(driverSupportsCreate) {
            // All is well, write using GDAL.
            writeGDALRaster(raster, pathName, *driver);
        } else {
            Error(QString(
                "Format driver %1 cannot be used to create datasets.").arg(
                    format.toLatin1().constData()));
        }
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
