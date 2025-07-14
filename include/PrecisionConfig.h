#ifndef PRECISIONCONFIG_H
#define PRECISIONCONFIG_H

// PrecisionConfig.h
#pragma once

#ifdef USE_FLOAT
using Real = float;
#else
using Real = double;
#endif

#endif // PRECISIONCONFIG_H
