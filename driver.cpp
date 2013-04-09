//driver.cpp
//Application entry point that uses the CGA class to run a genetic algorithm.
//
//
//
//Brien Smith-Martinez
//Unversity of Kansas
//
#define NDEBUG
#define VERSION_NUMBER 10
//#include <cassert>

#include <fstream>
#include <vector>
#include <complex>
#include <string>
#include <cstdlib>
#include <iostream>
#include <ctime>
// Required by fork routine
#include <sys/types.h>
#include <unistd.h>
using namespace std;



#include "FitnessTranslation.h"
#include "CGA.h"



//Application entry point. Opens files and sets the experiments in motion
int main(int argc, char* argv[] )
{
    //Variables used to configure the GA, initalized to defaults.
    string savefile("savestate_");
    string loadfile("");
    int randomseed = -1;
    unsigned int num_gens=100;
    unsigned int mu=100;
    unsigned int lambda=100;
    double mutation_rate=-1.0;
    bool useFitnessScaling=false;
    bool useUniqueness=false;
    int MAX_THREADS = 1;

    int PSK = 4;
    int txlength = 13;


    if( argc < 2)
    {
        cout << "Default Configuration" << endl;
    }

    int optind=1;
    // get arguments
    while ((optind < argc) && (argv[optind][0]=='-'))
    {
        string sw = argv[optind];
        if(sw=="-ss")
        {
            optind++;
            loadfile = argv[optind];
        }
        if(sw=="-sf")
        {
            optind++;
            savefile = argv[optind];
        }
        else if(sw=="-psk")
        {
            optind++;
            PSK = atoi(argv[optind]);
        }
        else if(sw=="-tx")
        {
            optind++;
            txlength = atoi(argv[optind]);
        }
        else if(sw=="-s")
        {
            optind++;
            randomseed = atoi(argv[optind]);
        }
        else if (sw=="-g")
        {
            optind++;
            num_gens = atoi(argv[optind]);
        }
        else if (sw=="-m")
        {
            optind++;
            mutation_rate = atof(argv[optind]);
        }
        else if (sw=="-p")
        {
            optind++;
            mu = atoi(argv[optind]);
        }
        else if (sw=="-o")
        {
            optind++;
            lambda = atoi(argv[optind]);
        }
        else if (sw=="-t")
        {
            optind++;
            MAX_THREADS = atoi(argv[optind]);
        }
        else if (sw=="-f")
        {
            useFitnessScaling=true;
        }
        else if (sw=="-u")
        {
            useUniqueness=true;
        }
        else
        {
            cout << "Unknown argument: " << argv[optind] << endl;
        }
        optind++;
    }
    if (MAX_THREADS == 0)
    {
        MAX_THREADS = 1;
    }
    else if (MAX_THREADS < 0 ) //|| MAX_THREADS > lambda)
    {
        MAX_THREADS = lambda;
    }

    if ( mutation_rate < 0)
    {
        mutation_rate = 1.0 / (100.0 * txlength);
    } 

    cout << "|| GA version: " << VERSION_NUMBER << endl;

    cout << "||Configuration:" << endl;

    cout << "|| tx code length = " << txlength << endl;
    cout << "|| PSK = " << PSK << endl;
    cout << "-------------" << endl;
    cout << "|| random seed = " << randomseed << endl;
    cout << "|| number of generations = " << num_gens << endl;
    cout << "|| population size = " << mu << endl;
    cout << "|| children per generation = " << lambda << endl;
    cout << "|| mutation rate = " << mutation_rate << endl;
    cout << "|| fitness scaling = " << useFitnessScaling << endl;
    cout << "|| use uniqueness = " << useUniqueness << endl;
    cout << "|| max number of eval threads = " << MAX_THREADS << endl;

    clock_t start_time = clock();
    clock_t gen_time;

    std::vector< std::complex<float> > seeder(51, complex<float>(0,0) );

    CGA myGA(txlength, PSK);

    GAIndividual elite(&myGA, txlength, PSK);
    myGA.MAX_THREADS = MAX_THREADS;
    myGA.Initialize(randomseed, mu, lambda, mutation_rate, useFitnessScaling, useUniqueness);

    cout << "Gen Evals  AvgF       BestF       Elite            Time" << endl;
    cout << "====================================" << endl;

    if( loadfile != "")
    {
        myGA.LoadState(loadfile);
        cout << "|| loaded from statefile = " << loadfile << endl;
    }

    for(unsigned int i = 0; i < num_gens; i++)
    {
        gen_time = clock();
        cout << i << " " << mu + lambda * i << "  ";
        myGA.RunGeneration();
        myGA.SaveState(savefile);

        gen_time = clock() - gen_time;
        cout << "  (" << gen_time / CLOCKS_PER_SEC  << ")"<< endl;
    }

    if( elite.fitness < myGA.elite.fitness)
    {
        elite = myGA.elite;
    }
    cout << "ELITE = " << elite.MSE << "  " << elite << endl;

    cout << "TOTAL TIME: " << (clock() - start_time) / CLOCKS_PER_SEC << endl;

    return(0);
}
