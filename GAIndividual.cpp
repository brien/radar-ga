
#include "GAIndividual.h"
#include "CGA.h"

/*
   GAIndividual::GAIndividual()
   {
   isEvaled = false;
   world = NULL;
   }*/

GAIndividual::GAIndividual(CGA* w, unsigned int cl, unsigned int al) : CHROMO_LENGTH(cl), ALPHA_LENGTH(al)
{
    fitness = 0;
    isEvaled = false;
    world = w;
}

ostream & operator<< (ostream & lhs, const GAIndividual & rhs)
{
    //lhs << "{";
    for(unsigned i = 0; i <  rhs.genotype.size(); i++)
    {
        lhs << rhs.genotype[i] << " ";
    }
    //lhs << "}";
    return lhs;
}

void GAIndividual::GenerateRandomIndividual()
{
    genotype.clear();
    for(unsigned int i=0; i<CHROMO_LENGTH; i++)
    {
        int r;
        r = world->mtRand.randInt(ALPHA_LENGTH-1);
        genotype.push_back(r);
    }
    isEvaled = false;
}
