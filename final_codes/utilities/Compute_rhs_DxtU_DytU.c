#include "mex.h"

/*
 INPUT: Wx, Wy, bx, by, tau
 OUTPUT: RHS
 
 Equivalent to the following code:
    
    RHS = tau*(DxtU(Wx-bx)+DytU(Wy-by));

    % compute D'_x(U)
    function dxtu = DxtU(U)
        dxtu = [U(:,end)-U(:, 1) U(:,1:end-1)-U(:,2:end)];
    end

    % compute D'_y(U)
    function dytu = DytU(U)
        dytu = [U(end,:)-U(1, :); U(1:end-1,:)-U(2:end,:)];
    end

*/

/* The gateway routine. */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double *Wxr, *Wxi, *Wyr, *Wyi, *bxr, *bxi, *byr, *byi;
    double *RHSr, *RHSi;
    double tau;
    bool bComplex;
    mwSize rows, cols;
    mwSize i, j, pos, colt, rowt;
    
    /* check for the proper number of arguments */
    if(nrhs != 5)
      mexErrMsgTxt("Exactly 5 input arguments are required.");
    if(nlhs != 1)
      mexErrMsgTxt("Need 1 output argument.");
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
    Wxr = mxGetPr(prhs[0]);
    Wyr = mxGetPr(prhs[1]);
    bxr = mxGetPr(prhs[2]);
    byr = mxGetPr(prhs[3]);

    if (bComplex)
    {
        Wxi = mxGetPi(prhs[0]);    
        Wyi = mxGetPi(prhs[1]);
        bxi = mxGetPi(prhs[2]);
        byi = mxGetPi(prhs[3]);
    }
    
    tau = mxGetScalar(prhs[4]);
    if (tau<=0) mexErrMsgTxt("5th input must be a scalar > 0.");

    /* create RHS */
    if (bComplex)
    {
        plhs[0] = mxCreateDoubleMatrix(rows, cols, mxCOMPLEX);
        RHSr = mxGetPr(plhs[0]); RHSi = mxGetPi(plhs[0]);
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
        RHSr = mxGetPr(plhs[0]);
    }
    
    /* MAIN COMPUTATION */
    if (bComplex)
    {
        /* Complex COMPUTATION */
        pos = 0;
        colt = rows*(cols-1);
        rowt = rows-1;
        /* for j = 1 */
            /* for i = 1*/
            {
                RHSr[pos] = tau*( Wxr[pos+colt]-bxr[pos+colt] - Wxr[pos]+bxr[pos] \
                                + Wyr[pos+rowt]-byr[pos+rowt] - Wyr[pos]+byr[pos] );
                RHSi[pos] = tau*( Wxi[pos+colt]-bxi[pos+colt] - Wxi[pos]+bxi[pos] \
                                + Wyi[pos+rowt]-byi[pos+rowt] - Wyi[pos]+byi[pos] );
                pos++;
            }
            for (i = 1; i < rows; i++)
            {
                RHSr[pos] = tau*( Wxr[pos+colt]-bxr[pos+colt] - Wxr[pos]+bxr[pos] \
                                + Wyr[pos -  1]-byr[pos -  1] - Wyr[pos]+byr[pos] );
                RHSi[pos] = tau*( Wxi[pos+colt]-bxi[pos+colt] - Wxi[pos]+bxi[pos] \
                                + Wyi[pos -  1]-byi[pos -  1] - Wyi[pos]+byi[pos] );
                pos++;
            }
        for (j = 1; j < cols; j++)
        {
            /* for i = 1*/
            {
                RHSr[pos] = tau*( Wxr[pos-rows]-bxr[pos-rows] - Wxr[pos]+bxr[pos] \
                                + Wyr[pos+rowt]-byr[pos+rowt] - Wyr[pos]+byr[pos] );
                RHSi[pos] = tau*( Wxi[pos-rows]-bxi[pos-rows] - Wxi[pos]+bxi[pos] \
                                + Wyi[pos+rowt]-byi[pos+rowt] - Wyi[pos]+byi[pos] );
                pos++;
            }
            for (i = 1; i < rows; i++)
            {
                RHSr[pos] = tau*( Wxr[pos-rows]-bxr[pos-rows] - Wxr[pos]+bxr[pos] \
                                + Wyr[pos -  1]-byr[pos -  1] - Wyr[pos]+byr[pos] );
                RHSi[pos] = tau*( Wxi[pos-rows]-bxi[pos-rows] - Wxi[pos]+bxi[pos] \
                                + Wyi[pos -  1]-byi[pos -  1] - Wyi[pos]+byi[pos] );
                pos++;
            }
        }
    }
    else
    {
    /* REAL-only COMPUTATION */
        pos = 0;
        colt = rows*(cols-1);
        rowt = rows-1;
        /* for j = 1 */
            /* for i = 1*/
            {
                RHSr[pos] = tau*( Wxr[pos+colt]-bxr[pos+colt] - Wxr[pos]+bxr[pos] \
                                + Wyr[pos+rowt]-byr[pos+rowt] - Wyr[pos]+byr[pos] );
                pos++;
            }
            for (i = 1; i < rows; i++)
            {
                RHSr[pos] = tau*( Wxr[pos+colt]-bxr[pos+colt] - Wxr[pos]+bxr[pos] \
                                + Wyr[pos -  1]-byr[pos -  1] - Wyr[pos]+byr[pos] );
                pos++;
            }
        for (j = 1; j < cols; j++)
        {
            /* for i = 1*/
            {
                RHSr[pos] = tau*( Wxr[pos-rows]-bxr[pos-rows] - Wxr[pos]+bxr[pos] \
                                + Wyr[pos+rowt]-byr[pos+rowt] - Wyr[pos]+byr[pos] );
                pos++;
            }
            for (i = 1; i < rows; i++)
            {
                RHSr[pos] = tau*( Wxr[pos-rows]-bxr[pos-rows] - Wxr[pos]+bxr[pos] \
                                + Wyr[pos -  1]-byr[pos -  1] - Wyr[pos]+byr[pos] );
                pos++;
            }
        }
    }

    return;
}




