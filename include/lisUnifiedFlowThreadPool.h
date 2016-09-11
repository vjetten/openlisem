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


/*!
 \file lisUnifiedFlowThreadPool.h
 \brief

functions: \n
*/


#ifndef LisemThreadPoolH
#define LisemThreadPoolH

#include "CsfMap.h"
#include <algorithm>
#include "model.h"
#include "operation.h"
#include <functional>
#include <thread>
#include "lisUnifiedFlowThread.h"
#include <QtGui>
#include <mutex>

//Main way to use functionaly of lisUnifiedFlowThreadPool
//use these macros around a function in the UnifiedFlowCode
//Running threaded means that each thread takes a part of the map as calculation domain
//running on a thread means one thread will calculate the entire function!
#define RUN_THREADED_ORDENED(x,c,self,...) ((ThreadPool->RunFunctionThreadedOrdered(std::bind(&(c::x),self,std::placeholders::_1 ,__VA_ARGS__ ))))
#define RUN_THREADED(x,c,self,...) ((ThreadPool->RunFunctionThreaded(std::bind(&(c::x),self,std::placeholders::_1 ,__VA_ARGS__ ))))
#define RUN_ON_THREAD(x,c,self,...) ((ThreadPool->RunFuctionOnThread(std::bind(&(c::x),self,std::placeholders::_1 ,__VA_ARGS__ ))))

class TWorld;
class LisemThread;
class ThreadFunction;

class LisemThreadPool: public QObject
{
    Q_OBJECT

public:

    int TP_NumberOfCores;

    bool initialized = false;

    std::mutex ThreadListMutex;
    QList<LisemThread*> ThreadList;
    QList<ThreadFunction*> ThreadCalls;
    QList<cTMap *> MaskList;
    QList<cTMap *> CellRList;
    QList<cTMap *> CellCList;
    QList<cTMap *> MaskListOrdered;
    QList<cTMap *> CellRListOrdered;
    QList<cTMap *> CellCListOrdered;
    cTMap * CoreMask;
    cTMap * InverseDTMask;
    cTMap * Mask;
    cTMap * DTMask;
    cTMap * cellR;
    cTMap * cellC;
    cTMap * minusone;
    cTMap * emptymask;
    int Priv_Iterator = 0;

    int _nrRows;
    int _nrCols;

    void InitThreads(TWorld * world);
    void SetMaskInitial(cTMap * _demmask);
    void SetMask(cTMap * _demmask, cTMap * _dtmask, cTMap * _cellR, cTMap * _cellC, bool create_subdivides = true);
    bool StartThread();
    void ThreadDone(LisemThread *thread);
    bool Allow_Work(int id);
    bool Allow_Final(int id);
    void RunFuctionOnThread(std::function<void (LisemThread*)> f);
    void RunFunctionThreaded(std::function<void (LisemThread*)> f);
    void RunFunctionThreadedOrdered(std::function<void (LisemThread*)> f);
    void WaitForAll();
    void clear();

};


/*
 * Getting the number of cores with c++ versions older than c++11 is more difficult
 * Example below is almost OS independently working (thanks to stackoverflow)
 *


#ifdef _WIN32
#include <windows.h>
#elif MACOS
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

int getNumCores() {
#ifdef WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#elif MACOS
    int nm[2];
    size_t len = 4;
    uint32_t count;

    nm[0] = CTL_HW; nm[1] = HW_AVAILCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);

    if(count < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &count, &len, NULL, 0);
        if(count < 1) { count = 1; }
    }
    return count;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}
*/

/*LisemThread * tr = 0;
std::function<void(LisemThread*)> f = std::bind(&(TWorld::UF2D_FluidSource),this,std::placeholders::_1 ,dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fn);
std::thread thread1 = std::thread(f,tr);
std::function<void(LisemThread*)> f2 = std::bind(&(TWorld::UF2D_SolidSource),this,std::placeholders::_1 ,dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_sn);
std::thread thread2 = std::thread(f2,tr);
thread1.join();
thread2.join();*/


#endif
