
#include "lisUnifiedFlowThread.h"
#include "lisUnifiedFlowThreadPool.h"
#include "include/model.h"

void LisemThread::Start()
{

    this->active.store(false);
    this->quit.store(false);
    this->threadobject = std::thread(&(LisemThread::Start_intern),this);
}

void LisemThread::Start_intern()
{

    std::unique_lock<std::mutex> lock(mutex_fr);
    while(true)
    {
        cv.notify_all();
        this->cv.wait(lock, [this]{return (quit.load() || active.load());});


        if(quit.load() == true)
        {
            active.store(false);
            lock.unlock();
            cv.notify_all();
            break;

        }
        if(this->active.load() && functionreference != 0)
        {
            if(functionreference->type == ThreadFunction::THREAD_SIMPLE)
            {
                time_used_in_last_function = std::chrono::high_resolution_clock::now();

                functionreference->f(functionreference->core);

                std::chrono::duration<double> timespan = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now()-time_used_in_last_function);

                //qDebug() << "time used in thread" << functionreference->core << "  " << timespan.count();
            }
        }

        this->active.store(false);
    }
    return;
}


int LisemThread::is_active()
{
    return active.load();
}
