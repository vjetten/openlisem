#include <cstdlib>
#include "error.h"
#include "raster_traits.h"
#include "fern/algorithm/algebra/elementary/add.h"


// Define global static variable and function. There is no LISEM lib, so we
// need to define them here again.
QString ErrorString;

void Error(
    QString string)
{
    ErrorString = string;
    throw 1;
}


void init_raster(
    cTMap& raster,
    int nr_rows,
    int nr_cols)
{
    raster.nrRows = nr_rows;
    raster.nrCols = nr_cols;

    raster.Data = new REAL8*[raster.nrRows];
    for(int row = 0; row < raster.nrRows; ++row) {
        raster.Data[row] = new REAL8[raster.nrCols];
    }
}


int main(
    int /* argc */,
    char** /* argv */)
{
    // Create rasters.
    int const nr_rows{6000};
    int const nr_cols{4000};

    cTMap raster1;
    init_raster(raster1, nr_rows, nr_cols);
    cTMap raster2;
    init_raster(raster2, nr_rows, nr_cols);
    cTMap raster3;
    init_raster(raster3, nr_rows, nr_cols);

    // Fill input rasters.
    // [0, 1, 2, 3, ...]
    for(int row = 0; row < nr_rows; ++row) {
        std::iota(raster1.Data[row], raster1.Data[row] + nr_cols,
            0 + row * nr_cols);
    }

    // [5, 5, 5, ...]
    for(int row = 0; row < nr_rows; ++row) {
        std::fill(raster2.Data[row], raster2.Data[row] + nr_cols, 5);
    }

    // Add input rasters.
    // [5, 6, 7, 8, ...]
    fern::algebra::add(raster1, raster2, raster3);

    // Verify result raster.
    assert(raster3.Data[0][100] == 100 + 5);
    assert(raster3.Data[100][100] == (100 * nr_cols + 100) + 5);

    return EXIT_SUCCESS;
}
