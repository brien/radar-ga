//FitnessTranslation.h
//
//Brien Smith-Martinez
//Unversity of Kansas
//

#ifndef BSM_MSE_TRANSLATION
#define BSM_MSE_TRANSLATION

#include <vector>
#include <complex>

#define BOOST_UBLAS_NDEBUG

#include <boost/math/special_functions/sinc.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace boost::math;
using namespace boost::numeric::ublas;

class TxCodeMSE
{
public:
    TxCodeMSE(unsigned int cl, unsigned int al);
    float CalculateMSE( const std::vector< unsigned int > &chrom);
    inline std::vector< std::complex<float> > GetAlpha()
    {
        return alpha;
    }
    //translates a GA genotype to a transmit code
    std::vector< std::complex<float> > GenToCode( const std::vector< unsigned int > &gen );

   //translates a transmit code to a GA genotype 
   std::vector< unsigned int > CodeToGen( std::vector< std::complex<float> > &code );

private:
    const unsigned int CODE_LENGTH;
    const unsigned int ALPHA_LENGTH;

    //the set of chips to build a code from:
    std::vector< std::complex<float> > alpha;

    //B MATRIX
    std::vector< matrix< std::complex<float> > > B;
    int K;
    unsigned int M, T;

};



#endif //BSM_MSE_TRANSLATION
