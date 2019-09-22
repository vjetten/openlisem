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
#include "operation.h"
#include <thread>
#include <functional>
#include "lisUnifiedFlowThreadPool.h"
#include <mutex>
#include <condition_variable>
#include <atomic>

class LisemThreadPool;
class ThreadFunction;

class LisemThread
{

public:

    int threadindex;

    std::atomic<bool> active;
    std::atomic<bool> done;

    std::unique_lock<std::mutex> tpl;


    LisemThreadPool *ThreadPool;

    int is_active();

    std::thread threadobject;

    ThreadFunction *functionreference;
    std::mutex mutex_fr;

    std::condition_variable cv;
    std::atomic<bool> quit;



    std::chrono::high_resolution_clock::time_point time_used_in_last_function;

    void Start();
    void Start_intern();
    void Quit();

};

#endif
