
#include "lisUnifiedFlowThread.h"
#include "lisUnifiedFlowThreadPool.h"
#include <chrono>

void LisemThreadPool::StartReportThread(TWorld * world)
{
//    tm = world->NewMap(0.0);
//    tma = world->NewMap(0.0);
//    tmb = world->NewMap(0.0);
//    tmc = world->NewMap(0.0);

    reportthread = new LisemThread();

    reportthread->tpl = std::unique_lock<std::mutex>(reportthread->mutex_fr);
    reportthread->active.store(false);
    reportthread->quit.store(false);
    reportthread->Start();
    reportthread->cv.notify_all();

}

void LisemThreadPool::RunReportFunction(std::function<void (int)> f)
{

    ThreadFunction *tf = new ThreadFunction();
    tf->type = ThreadFunction::THREAD_SIMPLE;
    tf->core = -1;
    tf->f = f;
    tf->final = false;
    LisemThread*tr =reportthread;

    std::unique_lock<std::mutex> lock(tr->mutex_fr);

    reportthread->active.store(true);
    reportthread->done.store(true);

    reportthread->functionreference = tf;

    if(lock.owns_lock())
    {
        lock.unlock();
    }
    tr->cv.notify_all();
}

void LisemThreadPool::WaitForReportThread()
{

    if(reportthread->active.load() || !reportthread->active.load() )
    {
        while(!reportthread->done.load())
        {

        }

        LisemThread*tr =reportthread;
        std::unique_lock<std::mutex> lock(tr->mutex_fr);

        if(lock.owns_lock())
        {
            lock.unlock();
        }
        tr->cv.notify_all();
    }
}

void LisemThreadPool::Close()
{

    LisemThread*tr =reportthread;

    std::unique_lock<std::mutex> rlock(reportthread->mutex_fr);
    reportthread->quit.store(true);
    reportthread->functionreference = 0;

    if(rlock.owns_lock())
    {
        rlock.unlock();
    }
    reportthread->cv.notify_all();
    for(int i = 0; i < ThreadList.length(); i++)
    {
        LisemThread *thread = ThreadList.at(i);
        std::unique_lock<std::mutex> lock(thread->mutex_fr);

        thread->quit.store(true);
        thread->functionreference = 0;
        if(lock.owns_lock())
        {
            lock.unlock();
        }
        thread->cv.notify_all();
    }

    qDebug() << "closed threads";
    //delete reportthread;
    //ThreadList.clear();
}

void LisemThreadPool::InitThreads(TWorld * world)
{
    //first get the number of cores
    TP_NumberOfCores = std::thread::hardware_concurrency();

    //minimum equals one, else nothing would happen and we wouldn't be here!
    TP_NumberOfCores = std::max(TP_NumberOfCores,1);
    if((!(TP_NumberOfCores % 2 == 0)) && TP_NumberOfCores != 1)
    {
        TP_NumberOfCores -= 1;
    }

    if (world->SwitchUserCores)
    {
        TP_NumberOfCores = std::min(world->userCores, TP_NumberOfCores);
    }
    //VJ
//qDebug() << "cores" << TP_NumberOfCores;
    minusone = world->NewMap(0.0);
    emptymask = world->NewMap(0.0);

    CoreMask = world->NewMap(-1.0);
    CoreMask2dfull = world->NewMap(-1.0);
    InverseDTMask = world->NewMap(0.0);
    Mask = world->NewMap(0.0);
    DTMask = world->NewMap(0.0);
    cellR = world->NewMap(0.0);
    cellC = world->NewMap(0.0);

    UF_t1.clear();

    Double_Out1.clear();

    _nrRows = world->_nrRows;
    _nrCols = world->_nrCols;

    clear();

    Priv_Iterator = 0;

    for(int i = 0; i < TP_NumberOfCores; i++)
    {
        LisemThread *thread = new LisemThread();

        thread->ThreadPool = this;
        thread->threadindex = i;
        ThreadList.append(thread);

        UF_t1.append(world->NewMap(0.0));

        Double_Out1.append(0.0);

        CellRDerListOrdered.append(world->NewMap(0.0));
        CellCDerListOrdered.append(world->NewMap(0.0));
        CellRMaskListOrdered.append(world->NewMap(0.0));
        CellCMaskListOrdered.append(world->NewMap(0.0));
        CellRListOrdered.append(world->NewMap(0.0));
        CellCListOrdered.append(world->NewMap(0.0));

        CellRDerListOrdered.at(i)->setAllMV();
        CellCDerListOrdered.at(i)->setAllMV();
        CellRMaskListOrdered.at(i)->setAllMV();
        CellCMaskListOrdered.at(i)->setAllMV();
        CellRListOrdered.at(i)->setAllMV();
        CellCListOrdered.at(i)->setAllMV();
    }

    LisemThread *thread = new LisemThread();

    thread->ThreadPool = this;
    thread->threadindex = TP_NumberOfCores;
    ThreadList.append(thread);

    CellRDerListOrdered.append(world->NewMap(0.0));
    CellCDerListOrdered.append(world->NewMap(0.0));
    CellRMaskListOrdered.append(world->NewMap(0.0));
    CellCMaskListOrdered.append(world->NewMap(0.0));
    CellRListOrdered.append(world->NewMap(0.0));
    CellCListOrdered.append(world->NewMap(0.0));

    UF_t1.append(world->NewMap(0.0));

    CellRDerListOrdered.at(TP_NumberOfCores)->setAllMV();
    CellCDerListOrdered.at(TP_NumberOfCores)->setAllMV();
    CellRMaskListOrdered.at(TP_NumberOfCores)->setAllMV();
    CellCMaskListOrdered.at(TP_NumberOfCores)->setAllMV();
    CellRListOrdered.at(TP_NumberOfCores)->setAllMV();
    CellCListOrdered.at(TP_NumberOfCores)->setAllMV();

    for(int i = 0; i < TP_NumberOfCores+1; i++)
    {
        LisemThread * tr = ThreadList.at(i);
        tr->tpl = std::unique_lock<std::mutex>(tr->mutex_fr);
        tr->active.store(false);
        tr->quit.store(false);
        tr->Start();
        tr->cv.notify_all();
    }

}

void LisemThreadPool::SetMaskInitial(cTMap * _demmask)
{

    bool allinend = false;
    int minr = _nrRows;
    int maxr = 0;
    int minc = _nrCols;
    int maxc = 0;
    emptymask->setAllMV();

    QList<int> trcl;
    QList<int> tccl;

    for(int i = 0; i < TP_NumberOfCores + 1; i ++)
    {

        trcl.append(0);
        tccl.append(0);
    }

    //double count = 0;
    for(int r = 0; r < _nrRows; r++)
    {
        for (int c = 0; c < _nrCols; c++)
        {
            minusone->Drc = -1;

            if(!pcr::isMV(_demmask->data[r][c]) )
            {
                 minr = std::min(r,minr);
                 maxr = std::max(r,maxr);
                 minc = std::min(c,minc);
                 maxc = std::max(c,maxc);
            }
        }
    }
    int dr, dc,c2,rsize,csplit;

    if(!(TP_NumberOfCores == 1))
    {
        dr = (maxr - minr);
        dc = (maxc - minc);
        if((dr < TP_NumberOfCores * 4) || (dc < TP_NumberOfCores * 4))
        {
            allinend = true;
        }
        c2 = TP_NumberOfCores/2;
        rsize = std::max(0,(dr)/c2);

        csplit = minc + std::floor((maxc-minc)/2.0);
    }else
    {
        allinend = true;
    }

    for(int rc = minr; rc < maxr+1; rc++)
    {
        for (int cc = minc; cc < maxc+1; cc++)
        {
         //   if(!pcr::isMV(_demmask->data[rc][cc]) )
        //    {
            bool add = false;
            if(!pcr::isMV(_demmask->data[rc][cc]) )
            {
                add = true;
            }

            bool inend = allinend;
            if(!allinend)
            {
                int rend = (rc - minr) % rsize;
                if(((cc == csplit) || (cc == csplit + 1) || (rend == 0) || (rend == 1)) )
                {
                    if(((rend == 0) || (rend == 1)) && !((cc == csplit) || (cc == csplit + 1)))
                    {
                        if(((rc - minr > 2) && ((rc - minr < dr - 3)) ))
                        {
                            inend = true;
                        }
                    }else
                    {
                        inend = true;
                    }
                }
            }

            if((!allinend) && (!inend))
            {


                int core = std::floor((rc - minr)/rsize);
                core = std::min( (TP_NumberOfCores/2) -1,core);

                if(cc >(csplit + 1) )
                {
                    core += c2;
                }
                core = std::max(0,std::min( TP_NumberOfCores - 1,core));
                CoreMask2dfull->data[rc][cc] = core;

                if(add)
                {
                    CoreMask->data[rc][cc] = core;

                    int trc = trcl.at(core);
                    int tcc = tccl.at(core);

                    CellRMaskListOrdered.at(core)->data[trc][tcc] = rc;
                    CellCMaskListOrdered.at(core)->data[trc][tcc] = cc;

                    tcc ++;
                    if(cc == _nrCols)
                    {
                        trc ++;
                        tcc = 0;
                    }
                    if(INSIDE(trc,tcc))
                    {
                        CellRMaskListOrdered.at(core)->data[trc][tcc] = -1;
                        CellCMaskListOrdered.at(core)->data[trc][tcc] = -1;
                    }
                    trcl.replace(core,trc);
                    tccl.replace(core,tcc);
                }
            } else
                {
                    int core = TP_NumberOfCores;
                    if(TP_NumberOfCores == 1)
                    {
                        core= 0;
                    }
                    CoreMask2dfull->data[rc][cc] = core;
                    if(add)
                    {

                        CoreMask->data[rc][cc] = core;

                        int trc = trcl.at(core);
                        int tcc = tccl.at(core);

                        CellRMaskListOrdered.at(core)->data[trc][tcc] = rc;
                        CellCMaskListOrdered.at(core)->data[trc][tcc] = cc;

                        tcc ++;
                        if(tcc == _nrCols)
                        {
                            trc ++;
                            tcc = 0;
                        }
                        if(INSIDE(trc,tcc))
                        {
                            CellRMaskListOrdered.at(core)->data[trc][tcc] = -1;
                            CellCMaskListOrdered.at(core)->data[trc][tcc] = -1;
                        }
                        trcl.replace(core,trc);
                        tccl.replace(core,tcc);
                    }


                }
          //  }
        }
    }

    trcl.clear();
    tccl.clear();

    for(int i = 0; i < TP_NumberOfCores + 1; i ++)
    {
        trcl.append(0);
        tccl.append(0);
    }

    for(int i = 0; i < TP_NumberOfCores +1 ; i++)
    {

        for(int r = 0; r < _nrRows; r++)
        {
            for (int c = 0; c < _nrCols; c++)
            {
                bool found = false;
                if(CoreMask2dfull->Drc == i)
                {
                    found = true;
                }
                if(INSIDE(r,c-1)){if(CoreMask2dfull->data[r][c-1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r,c+1)){if(CoreMask2dfull->data[r][c+1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r-1,c)){if(CoreMask2dfull->data[r-1][c] == i)
                {
                    found = true;
                }}
                if(INSIDE(r+1,c)){if(CoreMask2dfull->data[r+1][c] == i)
                {
                    found = true;
                }}
                if(INSIDE(r+1,c-1)){if(CoreMask2dfull->data[r+1][c-1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r-1,c+1)){if(CoreMask2dfull->data[r-1][c+1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r-1,c-1)){if(CoreMask2dfull->data[r-1][c-1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r+1,c+1)){if(CoreMask2dfull->data[r+1][c+1] == i)
                {
                    found = true;
                }}
                UF_t1.at(0)->Drc = found? 1.0:0.0;

            }
        }


        for(int r = 0; r < _nrRows; r++)
        {
            for (int c = 0; c < _nrCols; c++)
            {
                if(!pcr::isMV(_demmask->data[r][c]) )
                {
                    if(UF_t1.at(0)->Drc == 1)
                    {
                        int trc = trcl.at(i);
                        int tcc = tccl.at(i);

                        CellRDerListOrdered.at(i)->data[trc][tcc] = r;
                        CellCDerListOrdered.at(i)->data[trc][tcc] = c;

                        tcc ++;
                        if(tcc == _nrCols)
                        {
                            trc ++;
                            tcc = 0;
                        }
                        if(INSIDE(trc,tcc))
                        {
                            CellRDerListOrdered.at(i)->data[trc][tcc] = -1;
                            CellCDerListOrdered.at(i)->data[trc][tcc] = -1;
                        }
                        trcl.replace(i,trc);
                        tccl.replace(i,tcc);
                    }

                }
            }
        }
    }

}

void LisemThreadPool::SetMask(cTMap * _demmask, cTMap * _dtmask, cTMap * _cellR, cTMap * _cellC)
{
    for(int r = 0; r < _nrRows; r++)
    {
        for (int c = 0; c < _nrCols; c++)
        {

            cellC->Drc = _cellC->Drc;
            cellR->Drc = _cellR->Drc;
            Mask->Drc =  _demmask->Drc;
            DTMask->Drc =  _dtmask->Drc;



            for(int i = 0; i < TP_NumberOfCores + 1; i ++)
            {

                CellRListOrdered.at(i)->data[r][c] = -1;
                CellCListOrdered.at(i)->data[r][c] = -1;
            }
    }}

    QList<int> trcl;
    QList<int> tccl;

    for(int i = 0; i < TP_NumberOfCores + 1; i ++)
    {
        trcl.append(0);
        tccl.append(0);

    }
    double count = 0;
    for(int rc = 0; rc < _nrRows; rc++)
    {
        bool out = false;
        for (int cc = 0; cc < _nrCols; cc++)
        {
            int r = (int) (cellR->data[rc][cc]);
            int c = (int) (cellC->data[rc][cc]);

            if(!INSIDE(r,c)){out = true; break;}
            if(pcr::isMV(_demmask->Drc)){out = true; break;}

            int core = CoreMask->Drc;

            int trc = trcl.at(core);
            int tcc = tccl.at(core);

            CellRListOrdered.at(core)->data[trc][tcc] = r;
            CellCListOrdered.at(core)->data[trc][tcc] = c;

            count ++;

            tcc ++;
            if(tcc == _nrCols)
            {
                trc ++;
                tcc = 0;
            }

            if(INSIDE(trc,tcc))
            {
                CellRListOrdered.at(core)->data[trc][tcc] = -1;
                CellCListOrdered.at(core)->data[trc][tcc] = -1;
            }
            trcl.replace(core,trc);
            tccl.replace(core,tcc);
        }
        //if(out){break;}
    }
}

bool LisemThreadPool::StartThread()
{

    bool started = false;
    bool CouldStart = false;

    for(int j =  0; j < ThreadCalls.length();j++)
    {
        CouldStart = false;
        ThreadFunction * tf = ThreadCalls.at(j);
        if(tf->type == ThreadFunction::THREAD_SIMPLE)
        {

            if(tf->final)
            {
                if(ThreadCalls.length() == 1)
                {
                    CouldStart = true;
                    for(int i = 0; i <  ThreadList.length(); i++)
                    {
                        if((ThreadList.at(i)->active.load()))
                        {
                            CouldStart = false;
                        }
                    }

                }else
                {
                    CouldStart = false;
                }
            }else
            {
                CouldStart = true;
            }
        }

        if(CouldStart)
        {
            CouldStart = false;
            for(int i = 0; i <  TP_NumberOfCores; i++)
            {
                if(!(ThreadList.at(i)->active.load()) || tf->final)
                {
                        if(tf->final)
                        {
                            i =TP_NumberOfCores;
                        }

                        //We now know an unactive thread exist, and a threading call is available
                        //Start the thread!
                        LisemThread * thread = ThreadList.at(i);

                        std::unique_lock<std::mutex> lock(thread->mutex_fr);

                        thread->functionreference = tf;

                        ThreadCalls.removeAt(j);

                        thread->active.store(true);
                        thread->done.store(false);

                        if(lock.owns_lock())
                        {
                            lock.unlock();
                        }
                        thread->cv.notify_all();

                        started = true;
                        CouldStart = true;
                        j--;
                        break;
                }
            }
        }

    }

    return started;
}



void LisemThreadPool::WaitForThreadSafe(LisemThread *thread)
{


}

void LisemThreadPool::ThreadDone(LisemThread *thread)
{

}

void LisemThreadPool::RunCellCompute(std::function<void (int)> f)
{
    for(int i = 0; i < TP_NumberOfCores; i++)
    {


        ThreadFunction *tf = new ThreadFunction();
        tf->type = ThreadFunction::THREAD_SIMPLE;
        tf->core = i;
        tf->f = f;
        tf->final = false;
        ThreadCalls.append(tf);
        StartThread();

    }

}

void LisemThreadPool::RunDynamicCompute(std::function<void (int)> f)
{

    time_since_call = std::chrono::high_resolution_clock::now();

    for(int i = 0; i < TP_NumberOfCores; i++)
    {


        ThreadFunction *tf = new ThreadFunction();
        tf->type = ThreadFunction::THREAD_SIMPLE;
        tf->core = i;
        tf->f = f;
        tf->final = false;
        ThreadCalls.append(tf);
    StartThread();

    }

    if(TP_NumberOfCores > 1)
    {

        ThreadFunction *tf = new ThreadFunction();
        tf->type = ThreadFunction::THREAD_SIMPLE;
        tf->core = TP_NumberOfCores;
        tf->final = true;
        tf->f = f;
        ThreadCalls.append(tf);


    }
    StartThread();

}

void LisemThreadPool::WaitForAll()
{

    bool wait = true;


    while(wait == true)
    {
        wait = false;
        LisemThread *tr;
        int index = 0;
        for(int i = 0; i < ThreadList.length(); i++)
        {
            if(ThreadList.at(i)->active.load())
            {
                index = i;
                tr = (ThreadList.at(i));

                wait = true;
                break;
            }
        }

        if(wait)
        {
            while(!(tr->done.load()))
            {

            }

            std::unique_lock<std::mutex> lock(tr->mutex_fr);

            tr->active.store(false);
            if(lock.owns_lock())
            {
                lock.unlock();
            }
            tr->cv.notify_all();

            StartThread();

        }else
        {
            int l = ThreadCalls.length();


            if(l > 0)
            {
                if(StartThread())
                {
                    wait = true;
                }else
                {
                    break;
                }
            }else
            {
                break;
            }
        }
    }


  //  std::chrono::duration<double> timespan = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now()-time_since_call);
    //qDebug() << "TOTAL TIME:   " << timespan.count();
}


void LisemThreadPool::clear()
{

}
