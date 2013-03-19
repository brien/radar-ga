// Brien Smith-Martinez
// University of Kansas
// 
// This is an attempt to translate matlab to C++
// What follows is nasty, nasty code.

#include "FitnessTranslation.h"

#include <iostream>
using namespace std;

#include <math.h>
#include <vector>
#include <complex>



#define PI 3.14159265358979323846

std::vector<float> mysinc( std::vector<float> a )
{
    for( unsigned int i=0; i < a.size(); i++)
    {
        a[i] = sinc_pi( a[i] * PI );
    }
    return a;
}


std::vector<float> mylinspace( float a, float b, int n)
{
    std::vector<float> ret(n, 0);

    ret[0] = a;
    for(int i=0; i<n; i++)
    {
        ret[i] = a + i * ( (b - a) / (n-1) );
        //ret[i] = (float(i)/float(n)) * b + a;
    }
    return ret;
}


TxCodeMSE::TxCodeMSE(unsigned int cl, unsigned int al) : CODE_LENGTH(cl), ALPHA_LENGTH(al)
{
    complex<float> newsymbol;
    float phase_change = 360.0/float(ALPHA_LENGTH);
    phase_change = phase_change * PI / 180.0;
    complex<float> eye(0,1);
    for( unsigned int i = 0; i < ALPHA_LENGTH; i++)
    {
        complex<float> iminus1 = complex<float>(i,0);
        
        newsymbol = exp( iminus1*eye*phase_change );
        if( abs( real(newsymbol) ) < 0.001 )
            newsymbol = complex<float>(0.0, imag(newsymbol));
        if( abs( imag(newsymbol) ) < 0.001 )
            newsymbol = complex<float>(real(newsymbol), 0.0);
            
        alpha.push_back(newsymbol);
    }
    
    unsigned int N = CODE_LENGTH;

    //%M is the number of range cells. We want this number larger than the number
    //%of chips in our Tx Signal. Changing this number will affect the value of
    //%the MSE, so we need to make sure we're both using the same number for each
    //%code length N. I need to double check with Dr. Stiles to make sure 2*N  is
    //%adeqaute.

    //New:
    M = 2*N+10;
    T = M;

    //cout << M << endl;

    //%K is the time bandwidth product. We can arbitrarily set this to 500 as
    //%long as we consistantly use the same number for comparison

    //K = 500;

    //this is the bandwidth that we want to truncate our sinc function at, and corresponds to the null i.e. 2nd, 3rd,..... null
    int BW = 3;

    //This accounts for signal having positive and neg. freqs.
    BW = 2*BW;

    K = T*BW;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //%The following section of code creates a data cube B. At the risk of
    //%oversimplifying, B is a frequency domain representation of our Rx signal.
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    //T = M;


    //The linspace function generates linearly spaced vectors. It is similar to the colon operator ":", but gives direct control over the number of points.
    //y = linspace(a,b) generates a row vector y of 100 points linearly spaced between and including a and b.
    //y = linspace(a,b,n) generates a row vector y of n points linearly spaced between and including a and b. For n < 2, linspace returns b.


    //m = linspace(0,M-1,M);
    //k = linspace(-ceil(K/2),floor(K/2),K+1);


    std::vector<float> m = mylinspace(0, M-1, M);

    std::vector<float> k = mylinspace( -1 * ceil(K/2), K/2, K+1);
    

    //freq_samples = k./T; %This is the k/T term. This is a vector of all the frequency terms in the B matrix
    //sinc_pulse = sinc(freq_samples); %This is our frequency domain representation of our square pulse.
    //sinc_pulse = repmat(sinc_pulse,length(m),1);

    std::vector<float> freq_samples;

    for(unsigned int i=0; i < k.size(); i++)
        freq_samples.push_back( k[i] / T );

    //vector<float> sinc_pulse = sinc_pi(freq_samples);
    //matrix<float> sinc_pulse(freq_samples.size() ) = sinc_pi(freq_samples);
    std::vector<float> sinc_pulse( freq_samples.size(), 0 );
    //sinc_pulse = mysinc(freq_samples);

    for(unsigned int i=0; i < sinc_pulse.size(); i++)
    {
        sinc_pulse[i] = sinc_pi( freq_samples[i] * PI );
    }

/*
    for( int i=0; i < sinc_pulse.size(); i++)
        cout << "f[" << i << "] = " << freq_samples[i] << "  /  "<< "s[" << i << "] = " << sinc_pulse[i] << endl;
*/



    //sinc_pulse = myrepmat(sinc_pulse, m.size(), 1);

    matrix< float > mx_sinc_pulse( m.size(), sinc_pulse.size() );

    //std::copy(sinc_pulse.begin(), sinc_pulse.end(), mx_sinc_pulse.begin2());

    for( unsigned int i=0; i<mx_sinc_pulse.size1(); i++)
    {
        for( unsigned int j=0; j<mx_sinc_pulse.size2(); j++)
        {
            mx_sinc_pulse(i, j) = sinc_pulse[j];
        }
    }

    //%The following loop creates our data cube.

    //for n = 0:N-1 %iterate over the number of chips in the Tx signal
    //    B(:,:,n+1) = (sinc_pulse'.*exp(-j*2*pi*freq_samples'*(n+m)));
    //end

    //std::vector< std::vector< std::vector<float> > > B;
    //cout << B.size() << endl;

    complex< float > cplx_j(0.0,1.0);

    matrix< float > mx_freq_samples( 1, freq_samples.size() );
    mx_freq_samples.clear();
    std::copy( freq_samples.begin(), freq_samples.end() , mx_freq_samples.begin2());

    matrix< float > mx_m(1, m.size());
    std::copy( m.begin(), m.end() , mx_m.begin2());

    mx_sinc_pulse = trans(mx_sinc_pulse);

    matrix< float > mx_pp(1, m.size());
    mx_pp.clear();
    matrix< float > trans_freq_samples( trans(mx_freq_samples) );
    complex< float > cplx_tmp  = complex< float >((2.0 * PI * -1.0),0) * cplx_j ;


    B.clear();
    
    for( unsigned int i=0; i<N; i++)
    {
        matrix< std::complex<float> > Bmat;
        Bmat.clear();
        mx_pp.clear();


        for (unsigned j = 0; j <mx_pp.size1 (); j++)
            for (unsigned jj = 0; jj <mx_pp.size2 (); jj++)
                mx_pp (j, jj) = i + m[jj];

        Bmat = prod(trans_freq_samples, mx_pp);
        Bmat *= cplx_tmp;


        for(unsigned int j=0; j < Bmat.size1(); j++)
        {
            for(unsigned int jj=0; jj < Bmat.size2(); jj++)
            {
                //cout << exp(Bmat(j, jj)) << endl;
                //cout << "[" << j << "][" << jj << "]" << mx_sinc_pulse(j, jj) << " * " << exp(Bmat(j, jj));
                Bmat(j, jj) = mx_sinc_pulse(j, jj) * exp(Bmat(j, jj));
                //cout << " = " << Bmat(i, jj) << endl;

            }
        }

        //cout << Bmat << endl;
        B.push_back( Bmat );
    }
    
}

//translates a GA genotype to a transmit code
std::vector< std::complex<float> > TxCodeMSE::GenToCode( const std::vector< unsigned int > &gen )
{
    std::vector< std::complex<float> > inputv(CODE_LENGTH);
    
    for(int i=0; i<CODE_LENGTH; i++)
    {
        inputv[i] = alpha[gen[i]];
    }
    
    return inputv;
}

//translates a transmit code to a GA genotype 
std::vector< unsigned int > TxCodeMSE::CodeToGen( std::vector< std::complex<float> > &code )
{
    std::vector< unsigned int > inputv(CODE_LENGTH);
    
    for(int i=0; i<CODE_LENGTH; i++)
    {
        if( abs( real(code[i]) ) < 0.001 )
            code[i] = complex<float>(0.0, imag(code[i]));
        if( abs( imag(code[i]) ) < 0.001 )
            code[i] = complex<float>(real(code[i]), 0.0);
        //complex <float> lowestdiff(100, 100);
        float rlowestdiff = 100;
        float ilowestdiff = 100;
        float alowestdiff = 100;
        for(int j=0; j < ALPHA_LENGTH; j++)
        {
            //if( code[i] == alpha[j] )
            //  inputv[i] = j;
            
            float rdiff( abs( real(code[i]) - real(alpha[j]) ) );
            float idiff( abs( imag(code[i])  - imag(alpha[j])  ) );
            float adiff( abs(arg(code[i]) - arg(alpha[j])) );

            //cout << j << " " << code[i]  << " , " << alpha[j] << endl;
            //cout << i << " " << arg(code[i])  << " , " << arg(alpha[j]) << endl;
            //cout << i << " " << rdiff  << " , " << idiff << endl;


            //if(rdiff < rlowestdiff && idiff < ilowestdiff  )
            if( adiff < alowestdiff)
            {
                alowestdiff = adiff;
                rlowestdiff = rdiff;
                ilowestdiff = idiff;
                inputv[i] = j;
            }
        }
    }
    
    return inputv;
}

//variables that depend on inputv:
//matrix< complex<float> > P(K+1,M);


//float TxCodeMSE::CalculateMSE( const std::vector< complex<float> > &inputv)
float TxCodeMSE::CalculateMSE( const std::vector< unsigned int > &chrom)
{
    //cout << "Calling CalculateMSE !" << endl;
    std::vector< std::complex<float> > inputv(CODE_LENGTH);
    
    //cout << "inputv length = " << inputv.size() << endl;
    
    for(int i=0; i<CODE_LENGTH; i++)
    {
        //inputv.push_back( alpha[chrom[i]] );
        inputv[i] = alpha[chrom[i]];
    }
    
    //cout << "inputv length = " << inputv.size() << endl;
    
    //cout << "Translated code !" << endl;

    //--------------------------------
    matrix< complex<float> > P(K+1,M);
    P.clear();
    
    //cout << "Declared P !" << endl;

    //cout << P.size1() << " , " << P.size2() << endl;
    
    for (unsigned int i=0; i < inputv.size(); i++)
    {
        //cout << P << endl;
        P = inputv[i] * B[i] + P;
        //cout << "[" << i << "]" << P.size1() << " , " << P.size2() << endl;
    }
    
    //cout << P.size1() << " , " << P.size2() << endl;
    
    //cerr << "Made P !" << endl;

    matrix< float > Kg(M, 1);
    Kg.clear();

    for(unsigned int i=0; i<Kg.size1(); i++)
        Kg(i, 0) = 1;

    //matrix< complex<float> > transP;

    matrix< complex<float> > tmpP = prod(herm(P) , P);

    matrix<float> mag_pi_sqrd(tmpP.size2(), 1);

    //tmpP = real(tmpP);
    for(unsigned int i=0; i<mag_pi_sqrd.size1(); i++)
    {
        //cout << "real(tmpP(i,i)) = " << real(tmpP(i,i)) << endl;
        mag_pi_sqrd(i, 0) = real(tmpP(i,i));
    }

    int unsigned i = int(round(M/2.0));

    matrix< complex<float> > p_i(P.size1(), 1);
    matrix< complex<float> > p_j(P.size1(), M-1);
    
    for(unsigned int j=0; j < P.size1(); j++)
    {
        p_i(j, 0) = P(j, i);
    }

    for(unsigned int j=0; j < M-1; j++)
    {
        int ix = j;
        if( j >= i )
            ix++;
        for(unsigned int jj=0; jj < P.size1(); jj++)
            p_j(jj, j) = P(jj, ix);
    }
    matrix< complex<float> > hp_i = herm(p_i);
    matrix< complex<float> > prodhpi_pj = prod( hp_i, p_j );

    matrix< complex<float> > hermprodhpi_pj = herm(prodhpi_pj);



    matrix< float > mag_rhoij = real( prod( prodhpi_pj , hermprodhpi_pj ) );

    //cout << mag_rhoij << endl;

    //MSE_tmp(ii) = Kg(ii)*mag_rhoij/(mag_pi_sqrd(ii)^2);
    float MSE = (Kg(i,0) * mag_rhoij(0,0)) / (mag_pi_sqrd(i,0)*mag_pi_sqrd(i,0));

    return MSE;


}


