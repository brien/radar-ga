//EvaluationThread.h
//Thread class interface. For use in a multithreaded evaluation GA.
//
//Brien Smith-Martinez
//Unversity of Kansas
//

#ifndef BSM_GA_EVAL_THREAD
#define BSM_GA_EVAL_THREAD

#include <pthread.h>
#include <iostream>
using namespace std;

#define THREADID pthread_t;

//#include "GAIndividual.h"

class GAIndividual;

class EvaluationThread
{
   public:
      EvaluationThread();
      int Start(void * arg);
      void Join();
   private:
      GAIndividual* pchild;
      int Run(void * arg);
      static void * EntryPoint(void*);
      void Setup();
      void Execute(void*);
      void * Arg() const {return Arg_;}
      void Arg(void* a){Arg_ = a;}
   private:
      pthread_t ThreadId_;
      void * Arg_;

};

#endif //BSM_GA_EVAL_THREAD

