#include "fixture.h"
#include <gdal_priv.h>


Fixture::Fixture()
{
    // Register all GDAL drivers.
    GDALAllRegister();

    // GDAL mustn't throw in case of an error.
    CPLSetErrorHandler(CPLQuietErrorHandler);
}


Fixture::~Fixture()
{
}
