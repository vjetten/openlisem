#pragma once
#include <QString>
#include <masked_raster.h>

class cTMap;

bool               rasterCanBeOpenedForReading(
                                        QString const& pathName);

cTMap              readRaster          (QString const& pathName);

void               writeRaster         (cTMap const& raster,
                                        QString const& Name,
                                        QString const& format="PCRaster");

void               WriteMapSeries      (cTMap const& raster,
                                        QString const& Dir,
                                        QString Name,
                                        int count,
                                        QString const& format="PCRaster");

class cTRGBMap;

cTRGBMap  *readRasterImage        (QString const& pathName);
