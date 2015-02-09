#pragma once

/*!
    @brief      A fixture performs necessary setup and cleanup for Lisem to run.

    Upon construction a fixture initializes the runtime environment for Lisem.
    At the moment this involves configuring the GDAL library.

    This class can also be used as a test fixture in unit tests that depend
    on the initializations.
*/
class Fixture
{

public:

                   Fixture             ();

                   ~Fixture            ();

private:

};
