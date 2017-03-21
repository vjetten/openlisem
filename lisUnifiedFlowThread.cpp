
#include "lisUnifiedFlowThread.h"
#include "lisUnifiedFlowThreadPool.h"
#include "include/model.h"

void LisemThread::Start()
{

    this->active.store(false);
    this->done.store(true);
    this->quit.store(false);
    this->threadobject = std::thread(&(LisemThread::Start_intern),this);
}

void LisemThread::Start_intern()
{

    while(true)
    {
        std::unique_lock<std::mutex> lock(mutex_fr,std::defer_lock);
        this->cv.wait(lock, [this]{return (quit.load() || active.load());});

        if(quit.load() == true)
        {
            active.store(false);
            if(lock.owns_lock())
            {
                lock.unlock();
            }
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

                //qDebug() << "end thread work" <<timespan.count();

            }
        }

        this->active.store(false);
        this->done.store(true);
        if(lock.owns_lock())
        {
            lock.unlock();
        }
        cv.notify_all();

    }

    return;
}


int LisemThread::is_active()
{
    return active.load();
}
