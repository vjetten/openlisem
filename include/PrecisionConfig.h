#ifndef PRECISIONCONFIG_H
#define PRECISIONCONFIG_H

// PrecisionConfig.h
#pragma once

#ifdef USE_FLOAT
using Real = float;
//#define GDALfloat 6
//using RREAL = REAL4;
#else
using Real = double;
//#define GDALfloat 7
//using RREAL = REAL8;
#endif

#endif // PRECISIONCONFIG_H
