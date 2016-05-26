// MILP with CBC library, mps
// Finds the solution by using CBC Library
// Code Authors: Akshay Miterani and Pranav Deshpande

#include <sci_iofunc.hpp>

// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"=
#include "CbcCutGenerator.hpp"
#include "CbcHeuristicLocal.hpp"
#include "OsiClpSolverInterface.hpp"
extern "C" {
#include <api_scilab.h>

int mps_cppintlinprog()
{
    OsiClpSolverInterface solver;  
    
    // Path to the MPS file
    char *mpsFilePath;

    // Options to set maximum iterations
    double *options;

    // Input - 1 or 2 arguments allowed.
    CheckInputArgument(pvApiCtx, 2, 2);

    // Get the MPS File Path from Scilab
    getStringFromScilab(1, &mpsFilePath);
    
    // Receive the options for setting the maximum number of iterations etc.
    if( getFixedSizeDoubleMatrixFromScilab(2, 1, 4, &options))
    {
        return 1;
    }
    
    // Read the MPS file
    solver.readMps(mpsFilePath);

    // Cbc Library used from here
    CbcModel model(solver);

    model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
    
    if((int)options[0]!=0)
            model.setIntegerTolerance(options[0]);
    if((int)options[1]!=0)
            model.setMaximumNodes((int)options[1]); 
    if((int)options[2]!=0)
            model.setMaximumSeconds(options[2]);
    if((int)options[3]!=0)
            model.setAllowableGap(options[3]);
    
    model.branchAndBound();

    int nVars = model.getNumCols();
    int nCons = model.getNumRows();
    
    const double *val = model.getColSolution();
    
    //Output the solution to Scilab
    
    //get solution for x
    double* xValue = model.getColSolution();

    //get objective value
    double objValue = model.getObjValue();

    //Output status
    double status_=-1;
    if(model.isProvenOptimal()){
        status_=0;
    }
    else if(model.isProvenInfeasible()){
        status_=1;
    }
    else if(model.isSolutionLimitReached()){
        status_=2;
    }
    else if(model. isNodeLimitReached()){
        status_=3;
    }
    else if(model.isAbandoned()){
        status_=4;
    }
    else if(model.isSecondsLimitReached()){
        status_=5;
    }
    else if(model.isContinuousUnbounded()){
        status_=6;
    }
    else if(model.isProvenDualInfeasible()){
        status_=7;
    }

    double nodeCount = model.getNodeCount();
    double nfps = model.numberIntegers();
    double U = model.getObjValue();
    double L = model.getBestPossibleObjValue();
    double iterCount = model.getIterationCount();

    returnDoubleMatrixToScilab(1 , nVars, 1 , xValue);
    returnDoubleMatrixToScilab(2 , 1 , 1 , &objValue);
    returnDoubleMatrixToScilab(3 , 1 , 1 , &status_);
    returnDoubleMatrixToScilab(4 , 1 , 1 , &nodeCount);
    returnDoubleMatrixToScilab(5 , 1 , 1 , &nfps);
    returnDoubleMatrixToScilab(6 , 1 , 1 , &L);
    returnDoubleMatrixToScilab(7 , 1 , 1 , &U);
    returnDoubleMatrixToScilab(8 , 1 , 1 , &iterCount);

    return 0;
}

}
