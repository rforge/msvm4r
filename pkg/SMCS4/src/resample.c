#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP resample_systematic(SEXP logWeights)
{
    // uses logweights
    double dSumWeights, dLowerSum, dUpperSum, *us, *rLogWeights, u;
    SEXP out;
    int *iout;
    
    int N = LENGTH(logWeights);
    rLogWeights = REAL(logWeights);
    
    PROTECT(out = allocVector(INTSXP, N));
    iout = INTEGER(out);
        
    // compute cumulative sum of unnormalised weights
    dSumWeights = 0.0;
    for(int i = 0; i < N; i++)
        dSumWeights += exp(rLogWeights[i]);
    
    // Generate a random variable
    GetRNGstate();
    u = runif(0.0,dSumWeights/N);
    
    PutRNGstate();

    int j = 0, i = 0;
    dLowerSum = 0.0;
    dUpperSum = exp(rLogWeights[0]);
    double ut;
    while (i < N && j < N)
    {
        ut = u + ((double) i )/( (double) N) * dSumWeights;
        if (ut > dLowerSum && ut <= dUpperSum)
        {
            iout[i] = j + 1;
            i++;
        }
        else
        {
            j++;
            dLowerSum = dUpperSum;
            if (j < N)
                dUpperSum += exp(rLogWeights[j]);
        }
    }   
    UNPROTECT(1);
    return(out);
}

SEXP resample_stratified(SEXP logWeights)
{
    // function to perform stratified resampling
    // code adapted from Johansen's SMTC library
    // input: a vector with (unnormalised) log weights
    // returns: a vector with particle indices
    
    SEXP out;
    int *iout;
    double *rLogWeights;
    
    int N = LENGTH(logWeights);
    rLogWeights = REAL(logWeights);
    
    PROTECT(out = allocVector(INTSXP, N));
    iout = INTEGER(out);
    
    unsigned uRSCount[N];
    double dSumWeights = 0.0;
    
    // Procedure for stratified sampling
    double dCumulativeWeight = 0.0;
    // Calculate the normalising constant of the weight vector
    for(int i = 0; i < N; i++)
	    dSumWeights += exp(rLogWeights[i]);
    //Generate a random number between 0 and 1/N times the sum of the weights
    GetRNGstate();
    double u = runif(0.0, dSumWeights / ((double)N));

    int j = 0, k = 0;
    dCumulativeWeight = exp(rLogWeights[0]);
    while(j < N) {
	    while((dCumulativeWeight - u) > (dSumWeights*((double)j))/((double)N) && j < N) {
	        uRSCount[k]++;
	        iout[j] = k + 1;
	        j++;
	        u = runif(0,dSumWeights / ((double)N));
	    }
	    k++;
	    dCumulativeWeight += exp(rLogWeights[k]);
    }
    PutRNGstate();
    UNPROTECT(1);
    return out;
}

