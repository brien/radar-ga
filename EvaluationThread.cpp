
#include "EvaluationThread.h"
#include "CGA.h"

EvaluationThread::EvaluationThread() {}

int EvaluationThread::Start(void * arg)
{
   Arg(arg); // store user data
   int code = pthread_create( &ThreadId_, NULL, EvaluationThread::EntryPoint, (void*) this );
   return code;
}

void EvaluationThread::Join()
{
    pthread_join( ThreadId_, NULL );
}


int EvaluationThread::Run(void * arg)
{
   Setup();
   Execute( arg );
}

void * EvaluationThread::EntryPoint( void *pthis)
{
   EvaluationThread * pt = (EvaluationThread*)pthis;
   pt->Run( pt->Arg() );
}

void EvaluationThread::Setup()
{
    // Do any setup here
}

void EvaluationThread::Execute(void* arg)
{
    // Your code goes here
    pchild = (GAIndividual*)arg;
    pchild->world->Evaluate(*pchild);

    pthread_exit(0);
}
