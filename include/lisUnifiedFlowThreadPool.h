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


    QList<cTMap *> CellRDerListOrdered2d;
    QList<cTMap *> CellCDerListOrdered2d;
    QList<cTMap *> CellRMaskListOrdered2d;
    QList<cTMap *> CellCMaskListOrdered2d;
    QList<cTMap *> CellRListOrdered2d;
    QList<cTMap *> CellCListOrdered2d;

    QList<cTMap *> CellRDerListOrdered1d;
    QList<cTMap *> CellCDerListOrdered1d;
    QList<cTMap *> CellRListOrdered1d;
    QList<cTMap *> CellCListOrdered1d;
    QList<cTMap *> CellRMaskListOrdered1d;
    QList<cTMap *> CellCMaskListOrdered1d;

    QList<cTMap *>  UF_t1;
    QList<cTMap *>  UF_t2;
    QList<cTMap *>  UF_t3;
    QList<cTMap *>  UF_t4;
    QList<cTMap *>  UF_t5;
    QList<cTMap *>  UF_t6;
    QList<cTMap *>  UF_t7;
    QList<cTMap *>  UF_t8;
    QList<cTMap *>  UF_t9;
    QList<cTMap *>  UF_t10;
    QList<cTMap *>  UF_t11;

    cTMap * CoreMask1d;
    cTMap * InverseDTMask1d;
    cTMap * Mask1d;
    cTMap * DTMask1d;
    cTMap * cellR1d;
    cTMap * cellC1d;

    cTMap * CoreMask2d;
    cTMap * InverseDTMask2d;
    cTMap * Mask2d;
    cTMap * DTMask2d;
    cTMap * cellR2d;
    cTMap * cellC2d;

    cTMap * minusone;
    cTMap * emptymask;

    int Priv_Iterator = 0;

    int _nrRows;
    int _nrCols;

    cTMap * tm;
    cTMap * tma;
    cTMap * tmb;
    cTMap * tmc;

    LisemThread * reportthread;

    std::chrono::high_resolution_clock::time_point time_since_call;

    void StartReportThread(TWorld * world);
    void RunReportFunction(std::function<void (int)> f);
    void WaitForReportThread();
    void Close();

    void InitThreads(TWorld * world);
    void SetMaskInitial(cTMap * _demmask, cTMap * _ldd);
    void SetMask(cTMap * _demmask, cTMap * _dtmask2d, cTMap * _cellR2d, cTMap * _cellC2d,cTMap * _ldd,cTMap * _dtmask1d, cTMap * _cellR1d, cTMap * _cellC1d);
    bool StartThread();
    void WaitForThreadSafe(LisemThread *thread);
    void ThreadDone(LisemThread *thread);
    bool Allow_Work(int id);
    bool Allow_Final(int id);
    void RunDynamicCompute(std::function<void (int)> f);
    void RunCellCompute(std::function<void (int)> f);
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


#endif
