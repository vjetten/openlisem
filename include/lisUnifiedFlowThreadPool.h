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

class TWorld;
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
    int core;
    std::function<void (int)> f;
    int waitid;
    bool initial;
    bool final;
};

class LisemThreadPool: public QObject
{
    Q_OBJECT

public:

    int TP_NumberOfCores;

    bool initialized = false;

    std::mutex ThreadListMutex;
    QList<LisemThread*> ThreadList;
    QList<ThreadFunction*> ThreadCalls;


    QList<cTMap *> CellRDerListOrdered;
    QList<cTMap *> CellCDerListOrdered;
    QList<cTMap *> CellRMaskListOrdered;
    QList<cTMap *> CellCMaskListOrdered;
    QList<cTMap *> CellRListOrdered;
    QList<cTMap *> CellCListOrdered;

    QList<cTMap *>  UF_t1;

    QList<double>  Double_Out1;

    cTMap * CoreMask;
    cTMap * CoreMask2dfull;
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
    void SetMaskInitial(cTMap * _demmask);
    void SetMask(cTMap * _demmask, cTMap * _dtmask, cTMap * _cellR, cTMap * _cellC);
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
 * Example below is almost OS independently working (source: stackoverflow)
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
    sysctl(nm, 2, &count, &len, nullptr, 0);

    if(count < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &count, &len, nullptr, 0);
        if(count < 1) { count = 1; }
    }
    return count;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}
*/


#endif
