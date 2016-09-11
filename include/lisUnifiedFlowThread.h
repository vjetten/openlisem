/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Author: Victor Jetten
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

#ifndef LisemThreadH
#define LisemThreadH

/*!
 \file lisUnifiedFlowThread.h
 \brief

functions: \n
*/
#include <QtGui>
#include "CsfMap.h"
#include <algorithm>
#include "model.h"
#include "operation.h"
#include <thread>
#include <functional>
#include "lisUnifiedFlowThreadPool.h"


class LisemThreadPool;
class LisemThread;

class ThreadFunction
{
public:
    int maskindex;

    enum ThreadType
    {
        THREAD_SIMPLE,
        THREAD_SPATIAL,
        THREAD_SPATIAL_ORDERED
    };
    ThreadType type;
    std::function<void (LisemThread*)> f;
    int waitid;
    bool initial;
    bool final;
};

class TWorld;

class LisemThread : public QObject
{
    Q_OBJECT

public:

    int threadindex;
    bool active = false;
    LisemThreadPool *ThreadPool;

    std::thread threadobject;

    cTMap * mask;
    cTMap * cellR;
    cTMap * cellC;

    cTMap * UF_t1;
    cTMap * UF_t2;
    cTMap * UF_t3;
    cTMap * UF_t4;
    cTMap * UF_t5;
    cTMap * UF_t6;
    cTMap * UF_t7;
    cTMap * UF_t8;
    cTMap * UF_t9;
    cTMap * UF_t10;
    cTMap * UF_t11;

    ThreadFunction *functionreference;

    void CreateResources(TWorld * world, LisemThreadPool * Pool, int id);
    void Start();
    void Start_intern();
    void Quit();
    void CloseResources();

};

#endif
