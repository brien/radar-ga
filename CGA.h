//CGA.h
//Defines the CGA class interface. 
//
//Brien Smith-Martinez
//Unversity of Kansas
//

#ifndef BSM_CGA_H_INCLUDED
#define BSM_CGA_H_INCLUDED

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

#include "GAIndividual.h"
#include "MersenneTwister.h"
#include "EvaluationThread.h"


class CGA
{
    public:

    std::vector<GAIndividual> population;
    GAIndividual elite;

    CGA(unsigned int cl, unsigned int al);

    MTRand mtRand;

    void Initialize(int randomseed, unsigned int m, unsigned int l, double mr, bool fs, bool uni);
    void RunGeneration();
    int SaveState(string filename);
    int LoadState(string filename);
    
    void Evaluate( GAIndividual& ind);

    unsigned int MAX_THREADS;

    TxCodeMSE& GetFitnessFunction()
    {
        return txcodemse;
    }
    
    private:
    
    TxCodeMSE txcodemse;


    bool useFitnessScaling;
    bool useUniqueness;

    //--PROBLEM SPECIFIC STUFF:
    //Chromosome length:
    const unsigned int CHROMO_LENGTH;
    //length of the genetic alphabet
    unsigned int ALPHA_LENGTH;
    //---

    //number of children to create every generation,
    //and thus the number of old individuals replaced:
    unsigned int lambda;

    //number of individuals in the population
    unsigned int mu;

    //Mutation rate (per gene)
    double mutation_rate;
    
    //Children storage:
    std::vector<GAIndividual> children;

    //Functions:
    GAIndividual GenerateRandomIndividual();

    void ParentSelection(unsigned int &p1, unsigned int &p2);
    void ParentSelection_FitnessProportional(unsigned int &p1, unsigned int &p2);

    void NPointCrossover(const unsigned int n, const unsigned int p1, const unsigned int p2, unsigned int c1, unsigned int c2);
    void Mutate(const unsigned int c);

    //Fitness scaling:
    void Prescale(double umax, double uavg, double umin);
    double Scale(double u);
    void Scalepop(double max, double avg, double min);
    //a and b are used with scaled fitness
    double a;
    double b;
    double totalFit;
};



#endif // BSM_CGA_H_INCLUDED
