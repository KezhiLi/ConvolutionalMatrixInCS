#include "mex.h"
#include <math.h>

/* The gateway routine. */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double *Uxr, *Uxi, *Uyr, *Uyi, *bxr, *bxi, *byr, *byi;
    double *Vr;
    double xr, yr, xi, yi;
    bool bComplex;
    mwSize rows, cols;
    mwSize i, j, pos;
    
    /* check for the proper number of arguments */
    if(nrhs != 4)
      mexErrMsgTxt("Exactly 4 input arguments are required.");
    if(nlhs != 1)
      mexErrMsgTxt("Need 1 output argument.");
    /*Check that both inputs are row vectors*/
    rows = mxGetM(prhs[0]); cols = mxGetN(prhs[0]);
    if ( rows != mxGetM(prhs[1]) || \
         rows != mxGetM(prhs[2]) || \
         rows != mxGetM(prhs[3]) || \
         cols != mxGetN(prhs[1]) || \
         cols != mxGetN(prhs[2]) || \
         cols != mxGetN(prhs[3]) )
        mexErrMsgTxt("All input must have the same size.");
         
    bComplex = mxIsComplex(prhs[0]);
    if ( bComplex != mxIsComplex(prhs[1]) || \
         bComplex != mxIsComplex(prhs[2]) || \
         bComplex != mxIsComplex(prhs[3]))
        mexErrMsgTxt("All input must be consistently real or complex.");
  
    /* get pointers to the real and imaginary parts of the inputs */
    Uxr = mxGetPr(prhs[0]);
    Uyr = mxGetPr(prhs[1]);
    bxr = mxGetPr(prhs[2]);
    byr = mxGetPr(prhs[3]);

    if (bComplex)
    {
        Uxi = mxGetPi(prhs[0]);    
        Uyi = mxGetPi(prhs[1]);
        bxi = mxGetPi(prhs[2]);
        byi = mxGetPi(prhs[3]);
    }
  
    /* create V */
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    Vr = mxGetPr(plhs[0]);

    /* MAIN COMPUTATION */
    if (bComplex)
    {
    /* Complex COMPUTATION */
        for (pos = 0; pos < rows*cols; pos++)
        {
            xr = Uxr[pos] + bxr[pos];
            xi = Uxi[pos] + bxi[pos];
            yr = Uyr[pos] + byr[pos];
            yi = Uyi[pos] + byi[pos];
            Vr[pos] = sqrt(xr*xr + xi*xi + yr*yr + yi*yi);
        }
    }
    else
    {
    /* REAL-only COMPUTATION */
        for (pos = 0; pos < rows*cols; pos++)
        {
            xr = Uxr[pos] + bxr[pos];
            yr = Uyr[pos] + byr[pos];
            Vr[pos] = sqrt(xr*xr + yr*yr);
        }
    }

    return;
}




