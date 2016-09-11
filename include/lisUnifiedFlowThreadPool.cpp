
#include "lisUnifiedFlowThread.h"
#include "lisUnifiedFlowThreadPool.h"

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

    InverseDTMask = world->NewMap(0.0);
    CoreMask = world->NewMap(0.0);
    minusone = world->NewMap(0.0);
    /*Mask = NewMap(0.0);
    DTMask = NewMap(0.0);
    cellR = NewMap(0.0);
    cellC = NewMap(0.0);*/

    MaskList.clear();
    CellRList.clear();
    CellCList.clear();
    MaskListOrdered.clear();
    CellRListOrdered.clear();
    CellCListOrdered.clear();

    _nrRows = world->_nrRows;
    _nrCols = world->_nrCols;

    ThreadListMutex.lock();
    clear();

    Priv_Iterator = 0;

    for(int i = 0; i < TP_NumberOfCores; i++)
    {
        LisemThread *thread = new LisemThread();
        thread->CreateResources(world,this,i);
        thread->ThreadPool = this;
        ThreadList.append(thread);

        MaskList.append(world->NewMap(0.0));
        CellRList.append(world->NewMap(0.0));
        CellCList.append(world->NewMap(0.0));
        MaskListOrdered.append(world->NewMap(0.0));
        CellRListOrdered.append(world->NewMap(0.0));
        CellCListOrdered.append(world->NewMap(0.0));
    }
    MaskListOrdered.append(world->NewMap(0.0));
    CellRListOrdered.append(world->NewMap(0.0));
    CellCListOrdered.append(world->NewMap(0.0));

    ThreadListMutex.unlock();
}
void LisemThreadPool::SetMaskInitial(cTMap * _demmask)
{

    bool allinend = false;
    int minr = _nrRows;
    int maxr = 0;
    int minc = _nrCols;
    int maxc = 0;

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
    int dr = (maxr - minr);
    int dc = (maxc - minc);
    if(TP_NumberOfCores == 1 || (dr < TP_NumberOfCores * 4) || (dc < TP_NumberOfCores * 4))
    {
        allinend = true;
    }
    int c2 = TP_NumberOfCores/2;
    int rsize = std::max(0,(dr)/c2);

    int csplit = minc + std::floor((maxc-minc)/2.0);
    for(int rc = minr; rc < maxr+1; rc++)
    {
        for (int cc = minc; cc < maxc+1; cc++)
        {
            if(!pcr::isMV(_demmask->data[rc][cc]) )
            {
                bool inend = false;
                int rend = (rc - minr) % rsize;
                if((cc == csplit) || (cc == csplit + 1) || ((rc - minr > 2) && ((rc - minr < dr - 3)) ) &&((rend == 0) || (rend == 1)))
                {
                    inend = true;
                }

                if((!allinend) && (!inend))
                {
                    for(int i = 0; i < TP_NumberOfCores; i++)
                    {
                        MaskListOrdered.at(i)->data[rc][cc] = 0;
                    }

                    int core = std::floor((rc - minr)/rsize);
                    core = std::min( (TP_NumberOfCores/2) -1,core);

                    if(cc >(csplit + 1) )
                    {
                        core += c2;
                    }
                    core = std::max(0,std::min( TP_NumberOfCores - 1,core));
                    MaskListOrdered.at(core)->data[rc][cc] = 1;
                    CoreMask->data[rc][cc] = core;
                }else
                {
                    for(int i = 0; i < TP_NumberOfCores; i++)
                    {
                        MaskListOrdered.at(i)->data[rc][cc] = 0;
                    }
                    MaskListOrdered.at(TP_NumberOfCores)->data[rc][cc] = 1;
                    CoreMask->data[rc][cc] = TP_NumberOfCores;
                }
            }
        }
    }
    qDebug() << "setmaskinitial4" ;
}

void LisemThreadPool::SetMask(cTMap * _demmask, cTMap * _dtmask, cTMap * _cellR, cTMap * _cellC, bool create_subdivides)
{
    this->cellC = _cellC;
    this->cellR = _cellR;
    this->Mask = _demmask;
    this->DTMask = _dtmask;

    double count = 0;
    for(int r = 0; r < _nrRows; r++)
    {
        for (int c = 0; c < _nrCols; c++)
        {
            if(!pcr::isMV(_demmask->data[r][c]) )
            {
                 if(_dtmask > 0)
                 {
                     count ++;
                     this->InverseDTMask->Drc = 0;
                 }else
                 {
                     this->InverseDTMask->Drc = 1;
                 }
            }
        }
    }

    if(true)
    {
        int countpercore = std::floor(count/TP_NumberOfCores);

        int currentiterator = 0;
        int core = 0;
        int trc = 0;
        int tcc = 0;

        for(int rc = 0; rc < _nrRows; rc++)
        {
            for (int cc = 0; cc < _nrCols; cc++)
            {

                int r = (int) (_cellR->data[rc][cc]);
                int c = (int) (_cellC->data[rc][cc]);

                if(!INSIDE(r,c)){break;}

                CellRList.at(core)->data[trc][tcc] = r;
                CellCList.at(core)->data[trc][tcc] = c;

                tcc ++;
                if(cc == _nrCols)
                {
                    trc ++;
                    tcc = 0;
                }

                currentiterator ++;
                if(currentiterator == countpercore && !(core == TP_NumberOfCores -1))
                {
                    core ++;
                    currentiterator = 0;
                    trc = 0;
                    tcc = 0;
                }
                if(INSIDE(trc,tcc))
                {
                    CellRList.at(core)->data[trc][tcc] = -1;
                    CellCList.at(core)->data[trc][tcc] = -1;
                }

            }
        }
    }
    if(true)
    {

        QList<int> trcl;
        QList<int> tccl;

        for(int i = 0; i < TP_NumberOfCores + 1; i ++)
        {
            trcl.append(0);
            tccl.append(0);
        }

        for(int rc = 0; rc < _nrRows; rc++)
        {
            for (int cc = 0; cc < _nrCols; cc++)
            {

                int r = (int) (_cellR->data[rc][cc]);
                int c = (int) (_cellC->data[rc][cc]);

                if(!INSIDE(r,c)){break;}

                int core = CoreMask->Drc;

                int trc = trcl.at(core);
                int tcc = tccl.at(core);

                CellRListOrdered.at(core)->data[trc][tcc] = r;
                CellCListOrdered.at(core)->data[trc][tcc] = c;

                tcc ++;
                if(cc == _nrCols)
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
        }
    }

}

void LisemThreadPool::StartThread()
{
    ThreadListMutex.lock();

    bool CouldStart = true;
    while(CouldStart)
    {
        CouldStart = false;

        for(int i = 0; i <  ThreadList.length(); i++)
        {
            if(!(ThreadCalls.length() > 0))
            {
                break;
            }

            qDebug() << "search available thread" << ThreadList.at(i)->active << ThreadCalls.length();
            if(!(ThreadList.at(i)->active))
            {
                if(ThreadCalls.length() > 0)
                {

                    //We now know an unactive thread exist, and a threading call is available
                    //Start the thread!

                    LisemThread * thread = ThreadList.at(i);
                    int j = 0;
                    while( j < ThreadCalls.length())
                    {
                        thread->functionreference = ThreadCalls.at(j);
                        if(thread->functionreference->type == ThreadFunction::THREAD_SIMPLE)
                        {
                            thread->mask = this->Mask;
                            thread->cellC = this->cellC;
                            thread->cellR = this->cellR;
                            CouldStart = true;
                        }else if(thread->functionreference->type == ThreadFunction::THREAD_SPATIAL)
                        {
                            if(thread->functionreference->initial == true)
                            {
                                int index = thread->functionreference->maskindex;
                                thread->mask = this->Mask;
                                thread->cellC = minusone;
                                thread->cellR = minusone;
                                CouldStart = true;

                            }else if(thread->functionreference->initial == false && thread->functionreference->final == false)
                            {
                                if(Allow_Work(thread->functionreference->waitid))
                                {
                                    int index = thread->functionreference->maskindex;
                                    thread->mask = this->MaskList.at(index);
                                    thread->cellC = this->CellCList.at(index);
                                    thread->cellR = this->CellRList.at(index);
                                    CouldStart = true;
                                }
                            }

                        }else if(thread->functionreference->type == ThreadFunction::THREAD_SPATIAL_ORDERED)
                        {
                            if(thread->functionreference->initial == true)
                            {
                                int index = thread->functionreference->maskindex;
                                thread->mask = this->Mask;
                                thread->cellC = minusone;
                                thread->cellR = minusone;
                                CouldStart = true;

                            }else if(thread->functionreference->initial == false && thread->functionreference->final == false)
                            {
                                if(Allow_Work(thread->functionreference->waitid))
                                {
                                    int index = thread->functionreference->maskindex;
                                    thread->mask = this->MaskList.at(index);
                                    thread->cellC = this->CellCList.at(index);
                                    thread->cellR = this->CellRList.at(index);
                                    CouldStart = true;
                                }
                            }else if(thread->functionreference->final == true)
                            {
                                if(Allow_Final(thread->functionreference->waitid))
                                {
                                    int index = MaskList.length()-1;
                                    thread->mask = this->MaskList.at(index);
                                    thread->cellC = this->CellCList.at(index);
                                    thread->cellR = this->CellRList.at(index);
                                    CouldStart = true;
                                }
                            }
                        }
                        if(CouldStart)
                        {
                            ThreadCalls.removeAt(j);
                            thread->Start();
                        }else
                        {
                            j++;
                        }

                    }



                }else
                {
                    break;
                }

            }
        }
    }

    ThreadListMutex.unlock();
}

bool LisemThreadPool::Allow_Work(int id)
{

    for(int i = 0; i <  ThreadList.length(); i++)
    {
        LisemThread *tr = ThreadList.at(i);
        if(((tr->functionreference->type == ThreadFunction::THREAD_SPATIAL) || (tr->functionreference->type == ThreadFunction::THREAD_SPATIAL_ORDERED))
                && (tr->functionreference->initial == true) && (tr->functionreference->waitid == id) && (tr->active == true))
        {
            return false;
        }
    }

    for(int j = 0; j <  ThreadCalls.length(); j++)
    {
        if(((ThreadCalls.at(j)->type == ThreadFunction::THREAD_SPATIAL) || (ThreadCalls.at(j)->type == ThreadFunction::THREAD_SPATIAL_ORDERED))
                && (ThreadCalls.at(j)->initial == true) && (ThreadCalls.at(j)->waitid == id))
        {
            return false;
        }

    }

    return true;
}
bool LisemThreadPool::Allow_Final(int id)
{
    for(int i = 0; i <  ThreadList.length(); i++)
    {
        LisemThread *tr = ThreadList.at(i);
        if(((tr->functionreference->type == ThreadFunction::THREAD_SPATIAL) || (tr->functionreference->type == ThreadFunction::THREAD_SPATIAL_ORDERED))
                && (tr->functionreference->waitid == id) && (tr->active == true))
        {
            return false;
        }
    }

    for(int j = 0; j <  ThreadCalls.length(); j++)
    {
        if(((ThreadCalls.at(j)->type == ThreadFunction::THREAD_SPATIAL) || (ThreadCalls.at(j)->type == ThreadFunction::THREAD_SPATIAL_ORDERED))
                 && (ThreadCalls.at(j)->waitid == id))
        {
            return false;
        }

    }


}

void LisemThreadPool::ThreadDone(LisemThread *thread)
{
    ThreadListMutex.lock();
    thread->active = false;

    ThreadListMutex.unlock();

    StartThread();

}

void LisemThreadPool::RunFuctionOnThread(std::function<void (LisemThread*)> f)
{

    ThreadListMutex.lock();

    ThreadFunction *tf = new ThreadFunction();
    tf->type = ThreadFunction::THREAD_SIMPLE;
    tf->f = f;
    ThreadCalls.append(tf);

    ThreadListMutex.unlock();

    StartThread();
}

void LisemThreadPool::RunFunctionThreaded(std::function<void (LisemThread*)> f)
{
    ThreadListMutex.lock();

    ThreadFunction *tf = new ThreadFunction();
    tf->type = ThreadFunction::THREAD_SPATIAL;
    tf->f = f;
    tf->maskindex = 0;
    tf->initial = true;
    ThreadCalls.append(tf);

    for(int i = 0; i < ThreadList.length(); i++)
    {
        ThreadFunction *tf = new ThreadFunction();
        tf->type = ThreadFunction::THREAD_SPATIAL;
        tf->f = f;
        tf->maskindex = i;
        ThreadCalls.append(tf);
    }

    ThreadListMutex.unlock();

    StartThread();
}

void LisemThreadPool::RunFunctionThreadedOrdered(std::function<void (LisemThread*)> f)
{
    ThreadListMutex.lock();
    int waitid = Priv_Iterator++;
    ThreadFunction *tf = new ThreadFunction();
    tf->type = ThreadFunction::THREAD_SPATIAL_ORDERED;
    tf->f = f;
    tf->maskindex = 0;
    tf->waitid= waitid;
    tf->initial = true;
    tf->final = true;
    ThreadCalls.append(tf);

    for(int i = 0; i < ThreadList.length(); i++)
    {
        ThreadFunction *tf = new ThreadFunction();
        tf->type = ThreadFunction::THREAD_SPATIAL_ORDERED;
        tf->f = f;
        tf->maskindex = i;
        tf->waitid= waitid;
        tf->initial = true;
        tf->final = true;
        ThreadCalls.append(tf);
    }

    ThreadFunction *tf2 = new ThreadFunction();
    tf->type = ThreadFunction::THREAD_SPATIAL_ORDERED;
    tf->f = f;
    tf->maskindex = 0;
    tf->waitid= waitid;
    tf->initial = true;
    tf->final = true;
    ThreadCalls.append(tf);

    ThreadListMutex.unlock();

    StartThread();
}


void LisemThreadPool::WaitForAll()
{
        qDebug() << "wait threads";
    bool wait = true;
    while(wait == true)
    {
        wait = false;
        std::thread *threadobject;

        ThreadListMutex.lock();
        for(int i = 0; i < ThreadList.length(); i++)
        {
            if(ThreadList.at(i)->threadobject.joinable())
            {
                threadobject = &(ThreadList.at(i)->threadobject);

                wait = true;
                break;
            }
        }
        ThreadListMutex.unlock();
        if(wait)
        {
            threadobject->join();
        }
    }
    qDebug() << "threads done";
}

void LisemThreadPool::clear()
{
    for(int i = ThreadList.length()-1; i > -1; i--)
    {
        ThreadList.at(i)->CloseResources();
        delete ThreadList.at(i);
        ThreadList.removeAt(i);
    }
}
