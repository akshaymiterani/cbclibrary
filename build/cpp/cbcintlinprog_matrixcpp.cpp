// MILP with CBC library, Matrix
// Code Authors: Akshay Miterani and Pranav Deshpande

#include <sci_iofunc.hpp>

// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristicLocal.hpp"
#include "OsiClpSolverInterface.hpp"
extern "C"{
#include <api_scilab.h>
#include "sciprint.h"

int matrix_cppintlinprog(){

    //Objective function
    double* obj;  
    //Constraint matrix coefficients
    double* conMatrix;
    //intcon Matrix
    double* intcon;  
    //Constraints upper bound
    double* conlb;
    //Constraints lower bound
    double* conub;
    //Lower bounds for variables
    double* lb;  
    //Upper bounds for variables
    double* ub;
    //options for maximum iterations and writing mps
    double* options;
    //Flag for Mps
    double flagMps;
    //mps file path
    char * mpsFile;
    //Error structure in Scilab  
    SciErr sciErr;
    //Number of rows and columns in objective function
    int nVars=0, nCons=0,temp1=0,temp2=0;
    int numintcons=0;
    double valobjsense;
    
    CheckInputArgument(pvApiCtx , 11 , 11);             //Checking the input arguments
    CheckOutputArgument(pvApiCtx , 8, 8);               //Checking the output arguments

    ////////// Manage the input argument //////////
    
    //Number of Variables
    if(getIntFromScilab(1,&nVars))
    {
        return 1;
    }

    //Number of Constraints
    if (getIntFromScilab(2,&nCons))
    {
        return 1;
    }

    //Objective function from Scilab
    temp1 = nVars;
    temp2 = nCons;
    if (getFixedSizeDoubleMatrixFromScilab(3,1,temp1,&obj))
    {
        return 1;
    }

    //intcon matrix
    if (getDoubleMatrixFromScilab(4,&numintcons,&temp2,&intcon))
    {
        return 1;
    }

    if (nCons!=0)
    {
        //conMatrix matrix from scilab
        temp1 = nCons;
        temp2 = nVars;

        if (getFixedSizeDoubleMatrixFromScilab(5,temp1,temp2,&conMatrix))
        {
            return 1;
        }

        //conLB matrix from scilab
        temp1 = nCons;
        temp2 = 1;
        if (getFixedSizeDoubleMatrixFromScilab(6,temp1,temp2,&conlb))
        {
            return 1;
        }

        //conUB matrix from scilab
        if (getFixedSizeDoubleMatrixFromScilab(7,temp1,temp2,&conub))
        {
            return 1;
        }

    }

    //lb matrix from scilab
    temp1 = 1;
    temp2 = nVars;
    if (getFixedSizeDoubleMatrixFromScilab(8,temp1,temp2,&lb))
    {
        return 1;
    }


    //ub matrix from scilab
    if (getFixedSizeDoubleMatrixFromScilab(9,temp1,temp2,&ub))
    {
        return 1;
    }

    //Object Sense
    if(getDoubleFromScilab(10,&valobjsense))
    {
        return 1;
    }

    //get options from scilab
    if(getFixedSizeDoubleMatrixFromScilab(11 , 1 , 4 , &options))
    {
        return 1;      
    }

    //------------Temporary Version to make coin packed matrix------
    OsiClpSolverInterface solver1;  
    
    CoinPackedMatrix *matrix =  new CoinPackedMatrix(false , 0 , 0);
    matrix->setDimensions(0 , nVars);
    for(int i=0 ; i<nCons ; i++)
    {
        CoinPackedVector row;
        for(int j=0 ; j<nVars ; j++)
        {
            row.insert(j, conMatrix[i+j*nCons]);
        }
        matrix->appendRow(row);
    }


    solver1.loadProblem(*matrix, lb, ub, obj, conlb, conub);
    
    for(int i=0;i<numintcons;i++)
        solver1.setInteger(intcon[i]-1);

    solver1.setObjSense(valobjsense);

    //-------------------------------------------------------------
    
    CbcModel model(solver1);

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

    double nodeCount=model.getNodeCount();
    double nfps=model.numberIntegers();
    double U=model.getObjValue();
    double L=model.getBestPossibleObjValue();
    double iterCount=model.getIterationCount();

    returnDoubleMatrixToScilab(1 , nVars, 1 , xValue);
    returnDoubleMatrixToScilab(2 , 1 , 1 , &objValue);
    returnDoubleMatrixToScilab(3 , 1 , 1 , &status_);
    returnDoubleMatrixToScilab(4 , 1 , 1 , &nodeCount);
    returnDoubleMatrixToScilab(5 , 1 , 1 , &nfps);
    returnDoubleMatrixToScilab(6 , 1 , 1 , &L);
    returnDoubleMatrixToScilab(7 , 1 , 1 , &U);
    returnDoubleMatrixToScilab(8 , 1 , 1 , &iterCount);

    //-------------------------------------------------------------
    
    return 0;
}
}
