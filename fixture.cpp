#include "fixture.h"
#include <gdal_priv.h>


/*!
    @brief      Initializes the Lisem runtime environment.

    All GDAL I/O drivers are registered and GDAL is configured not to throw
    exceptions, but to return error codes.
*/
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
