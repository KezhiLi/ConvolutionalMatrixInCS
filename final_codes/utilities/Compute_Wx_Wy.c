#include "mex.h"
#include <math.h>

/*
INPUT: Ux, Uy, bx, by, tau
OUTPUT: Wx, Wy

Equivalent to the following MATLAB code

    tic;
    UUx = Ux + bx; UUy = Uy + by;
    V = sqrt(UUx.*conj(UUx) + UUy.*conj(UUy));
    V = max(V - tau, 0) ./ max(V,eps);
    Wx = V.*UUx; Wy = V.*UUy;
    toc;

 */

/* The gateway routine. */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double *Uxr, *Uxi, *Uyr, *Uyi, *bxr, *bxi, *byr, *byi;
    double *Wxr, *Wxi, *Wyr, *Wyi;
    double tau, Vr;
    double xr, yr, xi, yi;
    bool bComplex;
    mwSize rows, cols;
    mwSize pos;
    
    /* check for the proper number of arguments */
    if(nrhs != 5)
      mexErrMsgTxt("Exactly 5 input arguments are required.");
    if(nlhs != 2)
      mexErrMsgTxt("Need 2 output argument.");
    /*Check that both inputs are row vectors*/
    rows = mxGetM(prhs[0]); cols = mxGetN(prhs[0]);
    if ( rows != mxGetM(prhs[1]) || \
         rows != mxGetM(prhs[2]) || \
         rows != mxGetM(prhs[3]) || \
         cols != mxGetN(prhs[1]) || \
         cols != mxGetN(prhs[2]) || \
         cols != mxGetN(prhs[3]) || \
         1    != mxGetM(prhs[4]) || \
         1    != mxGetN(prhs[4]) )
        mexErrMsgTxt("First 4 input must have the same size. 5th input must be a scalar.");
         
    bComplex = mxIsComplex(prhs[0]);
    if ( bComplex != mxIsComplex(prhs[1]) || \
         bComplex != mxIsComplex(prhs[2]) || \
         bComplex != mxIsComplex(prhs[3]) || \
         true     == mxIsComplex(prhs[4]) )
        mexErrMsgTxt("First 4 input must be consistently real or complex.");
  
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
    
    tau = mxGetScalar(prhs[4]);
    if (tau<=0) mexErrMsgTxt("5th input must be a scalar > 0.");
  
    /* create Wx Wy */
    if (bComplex)
    {
        plhs[0] = mxCreateDoubleMatrix(rows, cols, mxCOMPLEX);
        Wxr = mxGetPr(plhs[0]); Wxi = mxGetPi(plhs[0]);
        plhs[1] = mxCreateDoubleMatrix(rows, cols, mxCOMPLEX);
        Wyr = mxGetPr(plhs[1]); Wyi = mxGetPi(plhs[1]);
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
        Wxr = mxGetPr(plhs[0]);
        plhs[1] = mxCreateDoubleMatrix(rows, cols, mxREAL);
        Wyr = mxGetPr(plhs[1]);
    }

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
            Vr = sqrt(xr*xr + xi*xi + yr*yr + yi*yi);
            if (Vr <= tau)
            {
                Wxr[pos] = 0; Wxi[pos] = 0; Wyr[pos] = 0; Wyi[pos] = 0;
            }
            else
            {
                Vr = (Vr - tau) / Vr;
                Wxr[pos] = xr*Vr; Wxi[pos] = xi*Vr; Wyr[pos] = yr*Vr; Wyi[pos] = yi*Vr;
            }
        }
    }
    else
    {
    /* REAL-only COMPUTATION */
        for (pos = 0; pos < rows*cols; pos++)
        {
            xr = Uxr[pos] + bxr[pos];
            yr = Uyr[pos] + byr[pos];
            Vr = sqrt(xr*xr + yr*yr);
            if (Vr <= tau)
            {
                Wxr[pos] = 0; Wyr[pos] = 0;
            }
            else
            {
                Vr = (Vr - tau) / Vr;
                Wxr[pos] = xr*Vr; Wyr[pos] = yr*Vr;
            }
        }
    }

    return;
}




