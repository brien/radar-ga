//CGA.cpp
//Defines the CGA class implementation. 
//
//Brien Smith-Martinez
//Unversity of Kansas
//


#include "CGA.h"
#include <algorithm>
using namespace std;


CGA::CGA(unsigned int cl, unsigned int al) : CHROMO_LENGTH(cl), ALPHA_LENGTH(al), txcodemse(cl, al), elite(this, cl, al)
{

}

void CGA::Evaluate( GAIndividual& ind )
{
    if( !ind.isEvaled )
    {
        ind.MSE = txcodemse.CalculateMSE( ind.genotype );
        ind.fitness = 1.0/ind.MSE;
        ind.isEvaled = true;
    }

}

int CGA::SaveState(string filename)
{
    ofstream f;
    f.open( filename.c_str() );
    if(f.is_open() )
    {
        for(int i=0; i < mu; i++)
        {
            f << population[i] << endl;
        }
        f.close();
        return 1;
    }
    return -1;
}

int CGA::LoadState(string filename)
{
    //LoadState should open a file that looks like this:
    //  #Tx: (-1,0),(1,0),(0,-1),(0,1)
    //  0, 1, 2, 3 0.54323
    //  1, 3, 0, 2 0.42450
    //
    //  END
    ifstream f;
    f.open( filename.c_str() );

    if( !f.is_open() )
    {
        cerr << "Population Inputfile: '" << filename << "' could not be found." << endl;
        return -1;
    }

    string s;

    //TO DO: Fix this to translate between GA and TX alphabet spaces
    for(int i=0; i < mu && !f.eof(); i++)
    {
        getline(f,s);
        stringstream ss(s);
        GAIndividual& ni = population[i];
        ni.isEvaled = false;
        ni.genotype.clear();
        for(int k=0; k<CHROMO_LENGTH; k++)
        {
            unsigned int ctmp;
            ss >> ctmp;
            ni.genotype.push_back(ctmp);
        }
        //cout << "Input individiual: " << ni << endl;
        Evaluate( ni );
    }

    f.close();
    return 1;
}

void CGA::Initialize(int randomseed, unsigned int m, unsigned int l, double mr, bool fs, bool uni)
{
   if(randomseed < 0)
        mtRand.seed();
    else
        mtRand.seed(randomseed);
        
    lambda = l;
    mu = m;
    mutation_rate = mr;
    useFitnessScaling = fs;
    useUniqueness = uni;
    
    totalFit = 0.0;

    cout << "|| Initalization: Parameters OKAY" << endl;
    population.reserve(mu+lambda);
    GAIndividual newi(this, CHROMO_LENGTH, ALPHA_LENGTH);

    for( unsigned int i = 0; i < mu; i++)
    {
        newi.GenerateRandomIndividual();
        population.push_back(newi);
    }

    std::vector<EvaluationThread> evalthreads;
    evalthreads.reserve(mu);

    unsigned int ck = 0;
    while( ck < mu)
    {
        for( unsigned int tk = 0; tk < MAX_THREADS; tk++)
        {
            if( ck < mu)
            {
                EvaluationThread thread1;
                evalthreads.push_back(thread1);
                evalthreads[tk].Start ( (void*) &population[ck]);
                //cout << "Thread " << tk << " Started. Child = " << ck <<  endl;
                ck++;
            }
        }
        for( unsigned int tk = 0; tk < evalthreads.size(); tk++)
        {
            evalthreads[tk].Join();
            //cout << "Thread " << tk << " Joined. " << endl;
        }
        evalthreads.clear();
    }
    
    for( unsigned int i = 0; i < mu; i++)
    {
        totalFit += population[i].fitness;
    }

    sort(population.begin(), population.end());
    elite = population[0];
    children.resize(lambda, newi);

    cout << "|| Initalization: Initial Random Population OKAY" << endl;

}



void CGA::NPointCrossover(const unsigned int n, const unsigned int p1, const unsigned int p2, unsigned int c1, unsigned int c2)
{
    GAIndividual& child1 = children[c1];
    GAIndividual& child2 = children[c2];
    children[c1].genotype.clear();
    children[c2].genotype.clear();
    children[c1].isEvaled = false;
    children[c2].isEvaled = false;
    unsigned int crossoverPoint=0;
    unsigned int lastCrossoverPoint=0;
    for(unsigned int i = 0; i<=n; i++)
    {
        lastCrossoverPoint = crossoverPoint;
        crossoverPoint=mtRand.randInt( CHROMO_LENGTH - 1 - lastCrossoverPoint);
        crossoverPoint += lastCrossoverPoint;
        if(i==n)
            crossoverPoint = CHROMO_LENGTH;

        if(i%2)
        {
            child1.genotype.insert( child1.genotype.begin()+lastCrossoverPoint,
                        population[p1].genotype.begin()+lastCrossoverPoint,
                        population[p1].genotype.begin()+crossoverPoint);

            child2.genotype.insert( child2.genotype.begin()+lastCrossoverPoint,
                        population[p2].genotype.begin()+lastCrossoverPoint,
                        population[p2].genotype.begin()+crossoverPoint);
        }
        else
        {
            child1.genotype.insert( child1.genotype.begin()+lastCrossoverPoint,
                        population[p2].genotype.begin()+lastCrossoverPoint,
                        population[p2].genotype.begin()+crossoverPoint);

            child2.genotype.insert( child2.genotype.begin()+lastCrossoverPoint,
                        population[p1].genotype.begin()+lastCrossoverPoint,
                        population[p1].genotype.begin()+crossoverPoint);
        }
    }
    
    
    return;
}
void CGA::Mutate(const unsigned int c)
{
    for(unsigned int ic=0; ic<CHROMO_LENGTH; ic++)
    {
        if(mtRand.rand(1.0f) < mutation_rate )
        {
            int r = children[c].genotype[ic];
            while( r == children[c].genotype[ic] )
            {
                r = mtRand.randInt(ALPHA_LENGTH-1);
            }
            children[c].genotype[ic] = r;
        }
    }
    return;
}

void CGA::RunGeneration()
{
    //In a change from ga5 to ga6,
    //The stats are calculated and printed at the beginning of each generation, instead of at the end.

    double sumMSE = 0.0;
    totalFit = 0.0;
    for( unsigned int i = 0; i < mu; i++)
    {
        sumMSE += population[i].MSE;
        totalFit += population[i].fitness;
    }
    double averageMSE = sumMSE / mu;


    cout << averageMSE << "  " << population[0].MSE << "  " << elite.MSE << " " << elite ;


    std::vector<EvaluationThread> evalthreads;
    evalthreads.reserve(lambda);
    if( useFitnessScaling )
    {
        //Find max, min and average fitness
        double maxf=0;
        double avgf=0;
        double minf=0;
        double tfit=0;

        for(unsigned int i=0; i<mu; i++)
        {
            if( population[i].fitness > maxf)
                maxf = population[i].fitness;
            if( population[i].fitness < minf)
                minf = population[i].fitness;
            tfit += population[i].fitness;
        }
        avgf = tfit/mu;
        Scalepop(maxf, avgf, minf);
    }
    //Parent Selection
    //Find lambda # of pairs
    for( unsigned int i = 0; i < lambda; i+=2)
    {
        
        unsigned int p1, p2;
        unsigned int c1 = i, c2 = i+1;
        
        if(c2 > lambda-1)
        {
            cerr << "child index out of bounds" << endl;
        }
        ParentSelection_FitnessProportional(p1, p2);

        //Children Generation
        //Create lambda # of children from pairs
        

        NPointCrossover( 1, p1, p2, c1, c2);
    
        Mutate(c1);
        Mutate(c2);
    
    }
    if( useUniqueness )
    {
        //check to see if there are any duplicates. If so, mutate them until they aren't
        for(int u = 0; u < lambda; u++)
        {
            for(int k = 0; k < lambda; k++)
            {
                while( u != k && children[u] == children[k])
                {
                    children[u].genotype[mtRand.randInt(CHROMO_LENGTH-1)] = mtRand.randInt(ALPHA_LENGTH-1);
                }
            }
            for(int k = 0; k < mu; k++)
            {
                while( children[u] == population[k])
                {
                    children[u].genotype[mtRand.randInt(CHROMO_LENGTH-1)] = mtRand.randInt(ALPHA_LENGTH-1);
                }
            }

        }
    }
    
 
    
    unsigned int ck = 0;
    while( ck < lambda)
    {
        for( unsigned int tk = 0; tk < MAX_THREADS; tk++)
        {
            if( ck < lambda)
            {
                EvaluationThread thread1;
                evalthreads.push_back(thread1);
                evalthreads[tk].Start ( (void*) &children[ck]);
                //cout << "Thread " << tk << " Started. Child = " << ck <<  endl;
                ck++;
            }
        }
        for( unsigned int tk = 0; tk < evalthreads.size(); tk++)
        {
            evalthreads[tk].Join();
            //cout << "Thread " << tk << " Joined. " << endl;
        }
        evalthreads.clear();
    }
   
    for( unsigned int k = 0; k < lambda; k++)
    {
        double oldmse = children[k].MSE;
        if( oldmse != children[k].MSE)
        {
            cout << "ERROR in multithreaded fitness eval" << endl;
        } 
    }
    

    //Survival Selection

    //Sort by fitness
    sort(population.begin(), population.end());
    if( elite.fitness < population[0].fitness )
    {
        elite = population[0];
    }

    //REPLACE WORST:
    //Delete the offspringSize worst individuals so we're back to PopSize
    //population.erase( population.begin() + mu , population.end() );
    //population.insert(population.end(), children.begin(), children.begin() );
    //sort(population.begin(), population.end());

    //ELITIST
    //Add the children to the population, sort them, delete the lambda worst so we're back down to mu.
    population.insert(population.end(), children.begin(), children.end() );
    sort(population.begin(), population.end());
    population.erase( population.begin() + mu , population.end() );

    //GENERATIONAL:
    //population = children;
    //sort(population.begin(), population.end());
    
    //GENERATIONA GAP: (replace the worst portion of the population with the G best kids.)
    //unsigned int G = mu/2;
    //sort(children.begin(), children.end() );
    //population.erase( population.begin() + G , population.end() );
    //population.insert(population.end(), children.begin(), children.begin() + G );
    

}


void CGA::ParentSelection(unsigned int &p1, unsigned int &p2)
{
    p1 = mtRand.randInt(ALPHA_LENGTH - 1);
    p2 = mtRand.randInt(ALPHA_LENGTH - 1);
}

void CGA::ParentSelection_FitnessProportional(unsigned int &p1, unsigned int &p2)
{
    //Select a parent based on the proportion of fitness
    double r = mtRand.rand(totalFit);


    double runningTotal=0.0f;
    p1=0;
    runningTotal += population[p1].fitness;
    while ( runningTotal < r )
    {
        p1++;
        runningTotal += population[p1].fitness;
    }

    r = mtRand.rand(totalFit);
    runningTotal=0.0f;
    p2=0;
    runningTotal += population[p2].fitness;
    while ( runningTotal < r )
    {
        p2++;
        runningTotal += population[p2].fitness;
    }
    if(p1 == p2)
    {
        p2--;
    }
    if(p1 == 0)
    {
        p2 = 1;
    }
    
    //cerr << "Parent Selection: Selected " << "Parent 1: " << p1 << ", Parent 2: " << p2 << endl;

}

//Fitness Scaling code:
void CGA::Prescale(double umax, double uavg, double umin)
{
    const double fmultiple = 2.0;
    double delta;

    if( umin > (fmultiple * uavg - umax)/(fmultiple - 1.0))
    {
        delta = umax - uavg;
        a = (fmultiple - 1.0) * uavg / delta;
        b = uavg * (umax - fmultiple*uavg) / delta;
    }
    else
    {
        delta = uavg - umin;
        a = uavg / delta;
        b = -1.0*umin * uavg / delta;
    }

}

double CGA::Scale(double u)
{
    //cout << "Scaling: " << u << " -> " << a * u + b << endl;
    return ( a * u + b);
}

void CGA::Scalepop(double max, double avg, double min)
{
    Prescale(max, avg, min);
    for(unsigned int i=0; i < population.size(); i++)
    {
        population[i].fitness = Scale(1.0 / population[i].MSE);
    }

}
