
#include "lisUnifiedFlowThread.h"
#include "lisUnifiedFlowThreadPool.h"

void LisemThreadPool::StartReportThread(TWorld * world)
{
    tm = world->NewMap(0.0);
    tma = world->NewMap(0.0);
    tmb = world->NewMap(0.0);
    tmc = world->NewMap(0.0);

    reportthread = new LisemThread();
    reportthread->Start();

}

void LisemThreadPool::RunReportFunction(std::function<void (int)> f)
{

    ThreadFunction *tf = new ThreadFunction();
    tf->type = ThreadFunction::THREAD_SIMPLE;
    tf->core = -1;
    tf->f = f;
    tf->final = false;
    LisemThread*tr =reportthread;

    //std::unique_lock<std::mutex> lock(tr->mutex_fr);
    //tr->cv.wait(lock, [tr]{return !(tr->active);});
    std::unique_lock<std::mutex> lock = WaitForThreadSafe(tr);

    reportthread->active = true;

    reportthread->functionreference = tf;

    lock.unlock();
    reportthread->cv.notify_one();

}

void LisemThreadPool::WaitForReportThread()
{

    if(reportthread->active)
    {
        LisemThread*tr =reportthread;
        //std::unique_lock<std::mutex> lock(tr->mutex_fr);
        //tr->cv.wait(lock, [tr]{return !(tr->active);});
        std::unique_lock<std::mutex> lock = WaitForThreadSafe(tr);

        lock.unlock();
    }
}

void LisemThreadPool::Close()
{
    std::unique_lock<std::mutex> lock(reportthread->mutex_fr);

    reportthread->quit = true;

    reportthread->functionreference = 0;

    lock.unlock();
    reportthread->cv.notify_all();



    for(int i = 0; i < ThreadList.length(); i++)
    {
        LisemThread *thread = ThreadList.at(i);

        std::unique_lock<std::mutex> lock(thread->mutex_fr);

        thread->quit = true;

        thread->functionreference = 0;

        lock.unlock();
        thread->cv.notify_all();
    }

    for(int i = 0; i < ThreadList.length(); i++)
    {
        LisemThread *thread = ThreadList.at(i);
        if(thread->threadobject.joinable())
        {
            thread->threadobject.join();
        }
    }

    if(reportthread->threadobject.joinable())
    {
        //reportthread->threadobject.join();
    }

    for(int i = ThreadList.length()-1; i > -1; i--)
    {
        delete ThreadList.at(i);
        ThreadList.removeAt(i);
    }
    //delete reportthread;
    ThreadList.clear();
}

void LisemThreadPool::InitThreads(TWorld * world)
{
    //first get the number of cores
    TP_NumberOfCores = 1;//std::thread::hardware_concurrency();
    //make this 1, instead of the actual cores


    //minimum equals one, else nothing would happen and we wouldn't be here!
    TP_NumberOfCores = std::max(TP_NumberOfCores,1);
    if((!(TP_NumberOfCores % 2 == 0)) && TP_NumberOfCores != 1)
    {
        TP_NumberOfCores -= 1;
    }

    minusone = world->NewMap(0.0);
    emptymask = world->NewMap(0.0);

    CellRDerListOrdered2d.clear();
    CellCDerListOrdered2d.clear();
    CellRListOrdered2d.clear();
    CellCListOrdered2d.clear();
    CellRListOrdered2d.clear();
    CellCListOrdered2d.clear();
    CellRDerListOrdered1d.clear();
    CellCDerListOrdered1d.clear();
    CellRListOrdered1d.clear();
    CellCListOrdered1d.clear();
    CellRListOrdered1d.clear();
    CellCListOrdered1d.clear();

    CoreMask1d = world->NewMap(-1.0);
    InverseDTMask1d = world->NewMap(0.0);
    Mask1d = world->NewMap(0.0);
    DTMask1d = world->NewMap(0.0);
    cellR1d = world->NewMap(0.0);
    cellC1d = world->NewMap(0.0);

    CoreMask2d = world->NewMap(-1.0);
    InverseDTMask2d = world->NewMap(0.0);
    Mask2d = world->NewMap(0.0);
    DTMask2d = world->NewMap(0.0);
    cellR2d = world->NewMap(0.0);
    cellC2d = world->NewMap(0.0);

    UF_t1.clear();
    UF_t2.clear();
    UF_t3.clear();
    UF_t4.clear();
    UF_t5.clear();
    UF_t6.clear();
    UF_t7.clear();
    UF_t8.clear();
    UF_t9.clear();
    UF_t10.clear();
    UF_t11.clear();

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
        UF_t2.append(world->NewMap(0.0));
        UF_t3.append(world->NewMap(0.0));
        UF_t4.append(world->NewMap(0.0));
        UF_t5.append(world->NewMap(0.0));
        UF_t6.append(world->NewMap(0.0));
        UF_t7.append(world->NewMap(0.0));
        UF_t8.append(world->NewMap(0.0));
        UF_t9.append(world->NewMap(0.0));
        UF_t10.append(world->NewMap(0.0));
        UF_t11.append(world->NewMap(0.0));


        CellRDerListOrdered1d.append(world->NewMap(0.0));
        CellCDerListOrdered1d.append(world->NewMap(0.0));
        CellRMaskListOrdered1d.append(world->NewMap(0.0));
        CellCMaskListOrdered1d.append(world->NewMap(0.0));
        CellRListOrdered1d.append(world->NewMap(0.0));
        CellCListOrdered1d.append(world->NewMap(0.0));

        CellRDerListOrdered2d.append(world->NewMap(0.0));
        CellCDerListOrdered2d.append(world->NewMap(0.0));
        CellRMaskListOrdered2d.append(world->NewMap(0.0));
        CellCMaskListOrdered2d.append(world->NewMap(0.0));
        CellRListOrdered2d.append(world->NewMap(0.0));
        CellCListOrdered2d.append(world->NewMap(0.0));



        CellRDerListOrdered1d.at(i)->setAllMV();
        CellCDerListOrdered1d.at(i)->setAllMV();
        CellRMaskListOrdered1d.at(i)->setAllMV();
        CellCMaskListOrdered1d.at(i)->setAllMV();
        CellRListOrdered1d.at(i)->setAllMV();
        CellCListOrdered1d.at(i)->setAllMV();


        CellRDerListOrdered2d.at(i)->setAllMV();
        CellCDerListOrdered2d.at(i)->setAllMV();
        CellRMaskListOrdered2d.at(i)->setAllMV();
        CellCMaskListOrdered2d.at(i)->setAllMV();
        CellRListOrdered2d.at(i)->setAllMV();
        CellCListOrdered2d.at(i)->setAllMV();
    }

    CellRDerListOrdered1d.append(world->NewMap(0.0));
    CellCDerListOrdered1d.append(world->NewMap(0.0));
    CellRMaskListOrdered1d.append(world->NewMap(0.0));
    CellCMaskListOrdered1d.append(world->NewMap(0.0));
    CellRListOrdered1d.append(world->NewMap(0.0));
    CellCListOrdered1d.append(world->NewMap(0.0));

    CellRDerListOrdered2d.append(world->NewMap(0.0));
    CellCDerListOrdered2d.append(world->NewMap(0.0));
    CellRMaskListOrdered2d.append(world->NewMap(0.0));
    CellCMaskListOrdered2d.append(world->NewMap(0.0));
    CellRListOrdered2d.append(world->NewMap(0.0));
    CellCListOrdered2d.append(world->NewMap(0.0));

    UF_t1.append(world->NewMap(0.0));
    UF_t2.append(world->NewMap(0.0));
    UF_t3.append(world->NewMap(0.0));
    UF_t4.append(world->NewMap(0.0));
    UF_t5.append(world->NewMap(0.0));
    UF_t6.append(world->NewMap(0.0));
    UF_t7.append(world->NewMap(0.0));
    UF_t8.append(world->NewMap(0.0));
    UF_t9.append(world->NewMap(0.0));
    UF_t10.append(world->NewMap(0.0));
    UF_t11.append(world->NewMap(0.0));


    CellRDerListOrdered1d.at(TP_NumberOfCores)->setAllMV();
    CellCMaskListOrdered1d.at(TP_NumberOfCores)->setAllMV();
    CellRMaskListOrdered1d.at(TP_NumberOfCores)->setAllMV();
    CellCMaskListOrdered1d.at(TP_NumberOfCores)->setAllMV();
    CellRListOrdered1d.at(TP_NumberOfCores)->setAllMV();
    CellCListOrdered1d.at(TP_NumberOfCores)->setAllMV();

    CellRDerListOrdered2d.at(TP_NumberOfCores)->setAllMV();
    CellCDerListOrdered2d.at(TP_NumberOfCores)->setAllMV();
    CellRMaskListOrdered2d.at(TP_NumberOfCores)->setAllMV();
    CellCMaskListOrdered2d.at(TP_NumberOfCores)->setAllMV();
    CellRListOrdered2d.at(TP_NumberOfCores)->setAllMV();
    CellCListOrdered2d.at(TP_NumberOfCores)->setAllMV();

    for(int i = 0; i < TP_NumberOfCores; i++)
    {
        LisemThread * tr = ThreadList.at(i);
        tr->Start();
    }

}

void LisemThreadPool::SetMaskInitial(cTMap * _demmask,cTMap *_ldd)
{

    bool allinend = false;
    int minr = _nrRows;
    int maxr = 0;
    int minc = _nrCols;
    int maxc = 0;
    emptymask->setAllMV();

    QList<int> trcl2d;
    QList<int> tccl2d;
    QList<int> trcl1d;
    QList<int> tccl1d;

    for(int i = 0; i < TP_NumberOfCores + 1; i ++)
    {
        trcl1d.append(0);
        tccl1d.append(0);
        trcl2d.append(0);
        tccl2d.append(0);
    }

    double count = 0;
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
    int dr,dc,c2,rsize,csplit;

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
            if(!pcr::isMV(_demmask->data[rc][cc]) )
            {
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
                    CoreMask2d->data[rc][cc] = core;

                    int trc = trcl2d.at(core);
                    int tcc = tccl2d.at(core);

                    CellRMaskListOrdered2d.at(core)->data[trc][tcc] = rc;
                    CellCMaskListOrdered2d.at(core)->data[trc][tcc] = cc;

                    tcc ++;
                    if(cc == _nrCols)
                    {
                        trc ++;
                        tcc = 0;
                    }
                    if(INSIDE(trc,tcc))
                    {
                        CellRMaskListOrdered2d.at(core)->data[trc][tcc] = -1;
                        CellCMaskListOrdered2d.at(core)->data[trc][tcc] = -1;
                    }
                    trcl2d.replace(core,trc);
                    tccl2d.replace(core,tcc);

                    if(_ldd != 0)
                    {
                        if(!pcr::isMV(_ldd->data[rc][cc]) )
                        {
                            CoreMask1d->data[rc][cc] = core;

                            trc = trcl1d.at(core);
                            tcc = tccl1d.at(core);

                            CellRMaskListOrdered1d.at(core)->data[trc][tcc] = rc;
                            CellCMaskListOrdered1d.at(core)->data[trc][tcc] = cc;

                            tcc ++;
                            if(tcc == _nrCols)
                            {
                                trc ++;
                                tcc = 0;
                            }
                            if(INSIDE(trc,tcc))
                            {
                                CellRMaskListOrdered1d.at(core)->data[trc][tcc] = -1;
                                CellCMaskListOrdered1d.at(core)->data[trc][tcc] = -1;
                            }
                            trcl1d.replace(core,trc);
                            tccl1d.replace(core,tcc);

                        }
                    }
                }else
                {
                    int core = TP_NumberOfCores;
                    if(TP_NumberOfCores == 1)
                    {
                        core= 0;
                    }
                    CoreMask2d->data[rc][cc] = core;

                    int trc = trcl2d.at(core);
                    int tcc = tccl2d.at(core);

                    CellRMaskListOrdered2d.at(core)->data[trc][tcc] = rc;
                    CellCMaskListOrdered2d.at(core)->data[trc][tcc] = cc;

                    tcc ++;
                    if(tcc == _nrCols)
                    {
                        trc ++;
                        tcc = 0;
                    }
                    if(INSIDE(trc,tcc))
                    {
                        CellRMaskListOrdered2d.at(core)->data[trc][tcc] = -1;
                        CellCMaskListOrdered2d.at(core)->data[trc][tcc] = -1;
                    }
                    trcl2d.replace(core,trc);
                    tccl2d.replace(core,tcc);

                    if(_ldd != 0)
                    {
                        if(!pcr::isMV(_ldd->data[rc][cc]))
                        {

                            CoreMask1d->data[rc][cc] = core;

                            trc = trcl1d.at(core);
                            tcc = tccl1d.at(core);

                            CellRMaskListOrdered1d.at(core)->data[trc][tcc] = rc;
                            CellCMaskListOrdered1d.at(core)->data[trc][tcc] = cc;

                            tcc ++;
                            if(tcc == _nrCols)
                            {
                                trc ++;
                                tcc = 0;
                            }
                            if(INSIDE(trc,tcc))
                            {
                                CellRMaskListOrdered1d.at(core)->data[trc][tcc] = -1;
                                CellCMaskListOrdered1d.at(core)->data[trc][tcc] = -1;
                            }
                            trcl1d.replace(core,trc);
                            tccl1d.replace(core,tcc);

                        }
                    }
                }
            }
        }
    }

    trcl2d.clear();
    tccl2d.clear();
    trcl1d.clear();
    tccl1d.clear();

    for(int i = 0; i < TP_NumberOfCores + 1; i ++)
    {
        trcl1d.append(0);
        tccl1d.append(0);
        trcl2d.append(0);
        tccl2d.append(0);
    }

    for(int i = 0; i < TP_NumberOfCores +1 ; i++)
    {

        for(int r = 0; r < _nrRows; r++)
        {
            for (int c = 0; c < _nrCols; c++)
            {
                bool found = false;
                if(CoreMask2d->Drc == i)
                {
                    found = true;
                }
                if(INSIDE(r,c-1)){if(CoreMask2d->data[r][c-1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r,c+1)){if(CoreMask2d->data[r][c+1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r-1,c)){if(CoreMask2d->data[r-1][c] == i)
                {
                    found = true;
                }}
                if(INSIDE(r+1,c)){if(CoreMask2d->data[r+1][c] == i)
                {
                    found = true;
                }}
                if(INSIDE(r+1,c-1)){if(CoreMask2d->data[r+1][c-1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r-1,c+1)){if(CoreMask2d->data[r-1][c+1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r-1,c-1)){if(CoreMask2d->data[r-1][c-1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r+1,c+1)){if(CoreMask2d->data[r+1][c+1] == i)
                {
                    found = true;
                }}
                UF_t1.at(0)->Drc = found? 1.0:0.0;
                found = false;
                if(CoreMask1d->Drc == i)
                {
                    found = true;
                }
                if(INSIDE(r,c-1)){if(CoreMask1d->data[r][c-1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r,c+1)){if(CoreMask1d->data[r][c+1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r-1,c)){if(CoreMask1d->data[r-1][c] == i)
                {
                    found = true;
                }}
                if(INSIDE(r+1,c)){if(CoreMask1d->data[r+1][c] == i)
                {
                    found = true;
                }}
                if(INSIDE(r-1,c-1)){if(CoreMask1d->data[r-1][c-1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r-1,c+1)){if(CoreMask1d->data[r-1][c+1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r+1,c+1)){if(CoreMask1d->data[r+1][c+1] == i)
                {
                    found = true;
                }}
                if(INSIDE(r+1,c-1)){if(CoreMask1d->data[r+1][c-1] == i)
                {
                    found = true;
                }}
                UF_t2.at(0)->Drc = (found? 1.0:0.0);
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
                        int trc = trcl2d.at(i);
                        int tcc = tccl2d.at(i);

                        CellRDerListOrdered2d.at(i)->data[trc][tcc] = r;
                        CellCDerListOrdered2d.at(i)->data[trc][tcc] = c;

                        tcc ++;
                        if(tcc == _nrCols)
                        {
                            trc ++;
                            tcc = 0;
                        }
                        if(INSIDE(trc,tcc))
                        {
                            CellRDerListOrdered2d.at(i)->data[trc][tcc] = -1;
                            CellCDerListOrdered2d.at(i)->data[trc][tcc] = -1;
                        }
                        trcl2d.replace(i,trc);
                        tccl2d.replace(i,tcc);
                    }

                }

                if(_ldd != 0)
                {
                    if(!pcr::isMV(_ldd->data[r][c]) && !pcr::isMV(_demmask->data[r][c]))
                    {
                        if(UF_t2.at(0)->Drc == 1)
                        {
                            int trc = trcl1d.at(i);
                            int tcc = tccl1d.at(i);

                            CellRDerListOrdered1d.at(i)->data[trc][tcc] = r;
                            CellCDerListOrdered1d.at(i)->data[trc][tcc] = c;

                            tcc ++;
                            if(tcc == _nrCols)
                            {
                                trc ++;
                                tcc = 0;
                            }
                            if(INSIDE(trc,tcc))
                            {
                                CellRDerListOrdered1d.at(i)->data[trc][tcc] = -1;
                                CellCDerListOrdered1d.at(i)->data[trc][tcc] = -1;
                            }
                            trcl1d.replace(i,trc);
                            tccl1d.replace(i,tcc);
                        }

                    }
                }
            }
        }
    }

}

void LisemThreadPool::SetMask(cTMap * _demmask, cTMap * _dtmask2d, cTMap * _cellR2d, cTMap * _cellC2d,cTMap * _ldd,cTMap * _dtmask1d, cTMap * _cellR1d, cTMap * _cellC1d)
{
    for(int r = 0; r < _nrRows; r++)
    {
        for (int c = 0; c < _nrCols; c++)
        {

            cellC2d->Drc = _cellC2d->Drc;
            cellR2d->Drc = _cellR2d->Drc;
            Mask2d->Drc =  _demmask->Drc;
            DTMask2d->Drc =  _dtmask2d->Drc;

            cellC1d->Drc = _cellC1d->Drc;
            cellR1d->Drc = _cellR1d->Drc;
            if(_ldd != 0)
            {
                Mask1d->Drc = _ldd->Drc;
            }
            DTMask1d->Drc =  _dtmask1d->Drc;

            for(int i = 0; i < TP_NumberOfCores + 1; i ++)
            {
                CellRListOrdered2d.at(i)->data[r][c] = -1;
                CellCListOrdered2d.at(i)->data[r][c] = -1;
                CellRListOrdered1d.at(i)->data[r][c] = -1;
                CellCListOrdered1d.at(i)->data[r][c] = -1;
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
            for (int cc = 0; cc < _nrCols; cc++)
            {
                int r = (int) (cellR2d->data[rc][cc]);
                int c = (int) (cellC2d->data[rc][cc]);

                if(!INSIDE(r,c)){break;}

                int core = CoreMask2d->Drc;

                int trc = trcl.at(core);
                int tcc = tccl.at(core);

                CellRListOrdered2d.at(core)->data[trc][tcc] = r;
                CellCListOrdered2d.at(core)->data[trc][tcc] = c;

                count ++;

                tcc ++;
                if(tcc == _nrCols)
                {
                    trc ++;
                    tcc = 0;
                }

                if(INSIDE(trc,tcc))
                {
                    CellRListOrdered2d.at(core)->data[trc][tcc] = -1;
                    CellCListOrdered2d.at(core)->data[trc][tcc] = -1;
                }
                trcl.replace(core,trc);
                tccl.replace(core,tcc);
            }
        }


    if(_ldd != 0)
    {

        QList<int> trcl;
        QList<int> tccl;

        for(int i = 0; i < TP_NumberOfCores + 1; i ++)
        {
            trcl.append(0);
            tccl.append(0);
            CellRListOrdered1d.at(i)->data[0][0] = -1;
            CellCListOrdered1d.at(i)->data[0][0] = -1;
        }

        for(int rc = 0; rc < _nrRows; rc++)
        {
            for (int cc = 0; cc < _nrCols; cc++)
            {


                    int r = (int) (cellR1d->data[rc][cc]);
                    int c = (int) (cellC1d->data[rc][cc]);
                    if(!INSIDE(r,c) ){break;}

                    if(!pcr::isMV(Mask1d->data[r][c]) && !pcr::isMV(Mask2d->data[r][c]))
                    {

                        int core = CoreMask1d->Drc;
                        count ++;
                        core = core;

                        int trc = trcl.at(core);
                        int tcc = tccl.at(core);

                        CellRListOrdered1d.at(core)->data[trc][tcc] = r;
                        CellCListOrdered1d.at(core)->data[trc][tcc] = c;

                        tcc ++;
                        if(tcc == _nrCols)
                        {
                            trc ++;
                            tcc = 0;
                        }
                        if(INSIDE(trc,tcc))
                        {
                            CellRListOrdered1d.at(core)->data[trc][tcc] = -1;
                            CellCListOrdered1d.at(core)->data[trc][tcc] = -1;
                        }

                        trcl.replace(core,trc);
                        tccl.replace(core,tcc);
                    }

            }
        }
    }

}

bool LisemThreadPool::StartThread()
{
    bool started = false;
    //ThreadListMutex.lock();
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
            for(int i = 0; i <  ThreadList.length(); i++)
            {
                if(!(ThreadList.at(i)->active))
                {
                        //qDebug() << "starting thread at" << i << " for core " << tf->core << "list length"  << ThreadCalls.length();

                        //We now know an unactive thread exist, and a threading call is available
                        //Start the thread!
                        LisemThread * thread = ThreadList.at(i);


                        //std::unique_lock<std::mutex> lock(thread->mutex_fr);
                        //thread->cv.wait(lock, [thread]{return !(thread->active);});
                        std::unique_lock<std::mutex> lock = WaitForThreadSafe(thread);

                        thread->active = true;

                        thread->functionreference = tf;
                        ThreadCalls.removeAt(j);

                        lock.unlock();
                        thread->cv.notify_all();

                        started = true;
                        CouldStart = true;
                        j--;
                        break;
                }
            }
        }


    }

    //ThreadListMutex.unlock();
    return started;
}



std::unique_lock<std::mutex> LisemThreadPool::WaitForThreadSafe(LisemThread *thread)
{

    bool done  = false;

    while(!done)
    {
        std::unique_lock<std::mutex> lock(thread->mutex_fr);

        thread->cv.wait(lock, [thread]{return !(thread->active);});

        if(thread->active == false)
        {
             done = true;
             return lock;
        }
        lock.unlock();
    }
    std::unique_lock<std::mutex> lock(thread->mutex_fr);
    return lock;
}

void LisemThreadPool::ThreadDone(LisemThread *thread)
{

}

void LisemThreadPool::RunCellCompute(std::function<void (int)> f)
{
    //ThreadListMutex.lock();
    for(int i = 0; i < TP_NumberOfCores; i++)
    {


        ThreadFunction *tf = new ThreadFunction();
        tf->type = ThreadFunction::THREAD_SIMPLE;
        tf->core = i;
        tf->f = f;
        tf->final = false;
        ThreadCalls.append(tf);


    }
    //ThreadListMutex.unlock();

    StartThread();
}

void LisemThreadPool::RunDynamicCompute(std::function<void (int)> f)
{

    time_since_call = std::chrono::high_resolution_clock::now();

    //ThreadListMutex.lock();
    for(int i = 0; i < TP_NumberOfCores; i++)
    {


        ThreadFunction *tf = new ThreadFunction();
        tf->type = ThreadFunction::THREAD_SIMPLE;
        tf->core = i;
        tf->f = f;
        tf->final = false;
        ThreadCalls.append(tf);


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
    //ThreadListMutex.unlock();

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
            if(ThreadList.at(i)->active)
            {
                index = i;
                tr = (ThreadList.at(i));

                wait = true;
                break;
            }
        }

        if(wait)
        {
            //std::unique_lock<std::mutex> lock(tr->mutex_fr);
            //tr->cv.wait(lock, [tr]{return !(tr->active);});
            std::unique_lock<std::mutex> lock = WaitForThreadSafe(tr);

            lock.unlock();

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
    ThreadCalls.clear();


    std::chrono::duration<double> timespan = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now()-time_since_call);

    //qDebug() << "TOTAL TIME:   " << timespan.count();
}


void LisemThreadPool::clear()
{

}
