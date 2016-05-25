
#include <sci_iofunc.hpp>
#include <cassert>
#include <iomanip>

// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristicLocal.hpp"
#include "OsiClpSolverInterface.hpp"
extern "C"{
#include <api_scilab.h>

int cppintlinprog(){

	//Objective function
	double* f;
	//Variables to be Integers
	double* intcon;
	//Constraint matrix coefficients
	double* A;  
	//Constraints upper bound
	double* b;
	//Equal Constraint matrix coefficients
	double* Aeq;  
	//Equal Constraints upper bound
	double* beq;
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

	/*
	 * Problem format: x = myintlinprog(number_of_variables,number_of_constraints,number_of_intcon,f,intcon,A,b,Aeq,beq,lb,ub,options)
	 */

	CheckInputArgument(pvApiCtx , 6 , 6);           //Checking the input arguments
	//CheckOutputArgument(pvApiCtx , 4, 4);             //Checking the output arguments

	// Number of Variables and Number of Constraints
	int nVars=0, nCons=0;
	int temp1=0, temp2=0;
	int numintcons=0;
	int numequals=0;

	//Objective function from Scilab
	temp1 = nVars;
	temp2 = nCons;
	if (getDoubleMatrixFromScilab(1,&nVars,&temp2,&f))
	{
		return 1;
	}


	//intcon matrix from scilab
	temp1 = 1;
	temp2 = nVars;
	if (getDoubleMatrixFromScilab(2,&temp1,&numintcons,&intcon))
	{
		return 1;
	}

	//A matrix from scilab
	temp1 = nCons;
	temp2 = nVars;
	if (getDoubleMatrixFromScilab(3,&nCons,&temp2,&A))
	{
		return 1;
	}

	//b matrix from scilab
	temp1 = nCons;
	temp2 = 1;
	if (getDoubleMatrixFromScilab(4,&temp1,&temp2,&b))
	{
		return 1;
	}

	//Aeq matrix from scilab
	temp1 = nCons;
	temp2 = nVars;
	if (getDoubleMatrixFromScilab(5,&numequals,&temp2,&Aeq))
	{
		return 1;
	}

	//b matrix from scilab
	temp1 = nCons;
	temp2 = 1;
	if (getDoubleMatrixFromScilab(6,&temp1,&temp2,&beq))
	{
		return 1;
	}

	//std::cout<<numintcons<<std::endl;
	//-----------------------Added Arrays--------------------------
	double rowLower[nCons+2*numequals];
	double columnLower[nVars];
	double columnUpper[nVars];
	lb=columnLower;
	ub=columnUpper;
	//Initializing it all
	for(int i=0; i<nCons+2*numequals ; i++){
		rowLower[i]= -COIN_DBL_MAX;
	}
	for(int i=0; i<nVars ; i++){
		columnLower[i]= -COIN_DBL_MAX;
	}
	for(int i=0; i<nVars ; i++){
		columnUpper[i]= COIN_DBL_MAX;
	}
	//------------Temporary Version to make coin packed matrix------
	OsiClpSolverInterface solver1;  
	OsiSolverInterface *solver = &solver1;

	CoinPackedMatrix *matrix =  new CoinPackedMatrix(false , 0 , 0);
	matrix->setDimensions(0 , nVars);
	for(int i=0 ; i<nCons ; i++)
	{
		CoinPackedVector row;
		for(int j=0 ; j<nVars ; j++)
		{
			row.insert(j, A[i+j*nCons]);
		}

		matrix->appendRow(row);
	}

	for(int i=0 ; i<numequals ; i++)
	{
		CoinPackedVector row;
		for(int j=0 ; j<nVars ; j++)
		{
			row.insert(j, Aeq[i+j*numequals]);
		}

		matrix->appendRow(row);
	}

	for(int i=0 ; i<numequals ; i++)
	{
		CoinPackedVector row;
		for(int j=0 ; j<nVars ; j++)
		{
			row.insert(j, -Aeq[i+j*numequals]);
		}

		matrix->appendRow(row);
	}

	double newb[nCons+2*numequals];

	for(int i=0;i<nCons;i++){
		newb[i]=b[i];
	}
	for(int i=0;i<numequals;i++){
		newb[i+nCons]=beq[i];
	}
	for(int i=0;i<numequals;i++){
		newb[i+nCons+numequals]=-beq[i];
	}

	solver1.loadProblem(*matrix, lb, ub, f, rowLower, newb);
	
	for(int i=0;i<numintcons;i++)
		solver1.setInteger(intcon[i]-1);
	
	solver1.setObjSense(1.0);

	//-------------------------------------------------------------
	
	CbcModel model(solver1);

	model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
	model.branchAndBound();
	
	const double *val = model.getColSolution();
	
	//Output the solution to Scilab
	//get solution for x
	double* xValue = model.getColSolution();

	//get objective value
	double objValue = model.getObjValue();


	returnDoubleMatrixToScilab(1 , 1 , nVars , xValue);
	returnDoubleMatrixToScilab(2 , 1 , 1 , &objValue);
	//printf("%f\n\nx0 = %f, x1 = %f, x2 = %f\n\n", model.getObjValue(), val[0], val[1],val[2]);
	
	//-------------------------------------------------------------
	
	return 0;
}
}


