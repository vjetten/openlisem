#pragma once
#include <cassert>
#include "fern/core/argument_traits.h"
#include "CsfMap.h"


namespace fern {

template<>
struct ArgumentTraits<cTMap>
{

    using value_type = double;

    using const_reference = double const&;

    using argument_category= array_2d_tag;

};


size_t size(
    cTMap const& raster,
    size_t index)
{
    assert(index == 0 || index == 1);
    return index == 0 ? raster.nrRows : raster.nrCols;
}


double const& get(
    cTMap const& raster,
    size_t row,
    size_t col)
{
    assert(row < size(raster, 0));
    assert(col < size(raster, 1));
    return raster.Data[row][col];
}


double& get(
    cTMap& raster,
    size_t row,
    size_t col)
{
    assert(row < size(raster, 0));
    assert(col < size(raster, 1));
    return raster.Data[row][col];
}

} // namespace fern
