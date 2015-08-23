#pragma once

#include <QString>


#define ADD 0
#define SUB 1
#define MUL 2
#define DIV 3
#define POW 4
#define MIN 5  //VJ 041120 added this functionality
#define MAX 6
#define LARGER 7
#define SMALLER 8
#define LARGEREQUAL 9
#define SMALLEREQUAL 10
#define HIGHER 11
#define LOWER 12


class cTMap;

void               copy                (cTMap& raster,
                                        cTMap const& other);

int                countUnits          (cTMap const& raster);

void               fill                (cTMap& raster,
                                        double value);

double             mapTotal            (cTMap const& raster);

double             mapAverage          (cTMap const& raster);

double             mapMinimum          (cTMap const& raster);

double             mapMaximum          (cTMap const& raster);

double             getWindowAverage    (cTMap const& raster,
                                        int r,
                                        int c,
                                        bool center);

void               cover               (cTMap& raster,
                                        cTMap const& value1,
                                        double value2);

void               calcValue           (cTMap& raster,
                                        double value,
                                        int oper);

void               calcMap             (cTMap& raster,
                                        cTMap const& value,
                                        int oper);

void               calc2Maps           (cTMap& raster,
                                        cTMap const& value1,
                                        cTMap const& value2,
                                        int oper);

void               calcMapValue        (cTMap& raster,
                                        cTMap const& value1,
                                        double value2,
                                        int oper);

void               checkMap            (cTMap const& raster,
                                        int oper,
                                        double value,
                                        QString SS);
