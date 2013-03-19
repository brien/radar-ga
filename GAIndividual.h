//GAIndividual.h
//Represents an individual candidate solution in the CGA class. 
//
//Brien Smith-Martinez
//Unversity of Kansas
//

#ifndef BSM_GA_INDIVIDUAL
#define BSM_GA_INDIVIDUAL


#include <string>
#include <complex>
#include <vector>
using namespace std;

#include "FitnessTranslation.h"


class CGA;

class GAIndividual
{

private:
	

public:
	//--PROBLEM SPECIFIC STUFF:
    //Chromosome length:
    unsigned int CHROMO_LENGTH;
    //length of the genetic alphabet
    unsigned int ALPHA_LENGTH;
	
	//GAIndividual();
	GAIndividual(CGA* w, unsigned int cl, unsigned int al );
	void GenerateRandomIndividual();
	CGA* world;
	
	double fitness;
	double MSE;
	bool isEvaled;

	//vector< complex<double> > genotype;
	std::vector< unsigned int > genotype;
	
	bool operator==(const GAIndividual &x ) const
	{
		for(int i=0; i<genotype.size(); i++)
		{
			if( genotype[i] != x.genotype[i] )
			{
				return false;
			}
		}
		
		return true;
	}
	
	//Used for sorting
	bool operator<( const GAIndividual &x ) const
	{
		return (x.fitness < fitness );
	}

    friend ostream & operator<< (ostream & lhs, const GAIndividual & rhs);

};




#endif //BSM_GA_INDIVIDUAL
