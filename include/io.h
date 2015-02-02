#pragma once
#include <memory>
#include <QString>
#include "csf.h"


class cTMap;

//! Function to close a CSF MAP.
auto close_csf_map = [](MAP* map) { Mclose(map); };

//! Auto-ptr type for CSF MAPs.
using MapPtr = std::unique_ptr<MAP, decltype(close_csf_map)>;

bool               LoadFromFile        (cTMap& raster);

void               WriteMap            (cTMap const& raster,
                                        QString Name);

void               WriteMapSeries      (cTMap const& raster,
                                        QString Dir,
                                        QString Name,
                                        int count);
