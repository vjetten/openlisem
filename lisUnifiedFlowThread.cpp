
#include "lisUnifiedFlowThread.h"
#include "lisUnifiedFlowThreadPool.h"
#include "include/model.h"

void LisemThread::CreateResources(TWorld * world, LisemThreadPool * Pool, int id)
{
    ThreadPool = Pool;

    /*mask = NewMap(0.0);
    cellx = NewMap(0.0);
    celly = NewMap(0.0);*/
    //not needed, instead uses reference to Threadpool maps

    UF_t1 = world->NewMap(0.0);
    UF_t2 = world->NewMap(0.0);
    UF_t3 = world->NewMap(0.0);
    UF_t4 = world->NewMap(0.0);
    UF_t5 = world->NewMap(0.0);
    UF_t6 = world->NewMap(0.0);
    UF_t7 = world->NewMap(0.0);
    UF_t8 = world->NewMap(0.0);
    UF_t9 = world->NewMap(0.0);
    UF_t10 = world->NewMap(0.0);
    UF_t11 = world->NewMap(0.0);

    threadindex = id;

}

void LisemThread::Start()
{
    active = true;

    if(threadobject.joinable())
    {
        threadobject.detach();
    }
    this->threadobject = std::thread(&(LisemThread::Start_intern),this);
}

void LisemThread::Start_intern()
{
    if(functionreference->type == ThreadFunction::THREAD_SIMPLE)
    {
        functionreference->f(this);
    }else if(functionreference->type == ThreadFunction::THREAD_SPATIAL)
    {
        functionreference->f(this);
    }else if(functionreference->type == ThreadFunction::THREAD_SPATIAL_ORDERED)
    {
        functionreference->f(this);
    }

    this->ThreadPool->ThreadDone(this);

}
void LisemThread::CloseResources()
{


}

void LisemThread::Quit()
{


}

