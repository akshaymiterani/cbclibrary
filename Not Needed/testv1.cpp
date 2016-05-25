
#include <sci_iofunc.hpp>
extern "C"{
#include <api_scilab.h>

int myintlinprog(){

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
	 * Problem format: x = myintlinprog(number_of_variables,number_of_constraints,f,intcon,A,b,Aeq,beq,lb,ub,options)
	 */

	CheckInputArgument(pvApiCtx , 11 , 11);           //Checking the input arguments
	CheckOutputArgument(pvApiCtx , 4, 4);             //Checking the output arguments

	// Number of Variables and Number of Constraints
	int nVars=0, nCons=0;
	int temp1=0, temp2=0;
	int numintcons;
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
	if (getFixedSizeDoubleMatrixFromScilab(3,1,nVars,&f))
	{
		return 1;
	}


	//intcon matrix from scilab
	temp1 = 1;
	temp2 = nVars;
	if (getDoubleMatrixFromScilab(4,&temp1,&numintcons,&intcon))
	{
		return 1;
	}

	//A matrix from scilab
	temp1 = nCons;
	temp2 = nVars;
	if (getFixedSizeDoubleMatrixFromScilab(5,temp1,temp2,&A))
	{
		return 1;
	}

	//b matrix from scilab
	temp1 = nCons;
	temp2 = 1;
	if (getFixedSizeDoubleMatrixFromScilab(6,temp1,temp2,&b))
	{
		return 1;
	}

	//Aeq matrix from scilab
	temp1 = nCons;
	temp2 = nVars;
	if (getFixedSizeDoubleMatrixFromScilab(7,temp1,temp2,&Aeq))
	{
		return 1;
	}

	//beq matrix from scilab
	temp1 = nCons;
	temp2 = 1;
	if (getFixedSizeDoubleMatrixFromScilab(8,temp1,temp2,&beq))
	{
		return 1;
	}

	//lb matrix from scilab
	temp1 = 1;
	temp2 = nVars;
	if (getFixedSizeDoubleMatrixFromScilab(9,temp1,temp2,&lb))
	{
		return 1;
	}


	//ub matrix from scilab
	if (getFixedSizeDoubleMatrixFromScilab(10,temp1,temp2,&ub))
	{
		return 1;
	}

	//get options from scilab
	if(getFixedSizeDoubleMatrixInList(11 , 2 , 1 , 1 , &options))
	{
		return 1;      
	}

	//------------Temporary Version to make coin packed matrix------
	CoinPackedMatrix *matrix =  new CoinPackedMatrix(false , 0 , 0);
	matrix->setDimensions(0 , numVars_);
	for(int i=0 ; i<numCons_ ; i++)
	{
		CoinPackedVector row;
		for(int j=0 ; j<numVars_ ; j++)
		{
			row.insert(j, A[i+j*numCons_]);
		}
		matrix->appendRow(row);
	}
	//-------------------------------------------------------------
	
	//

	return 0;
}
}


