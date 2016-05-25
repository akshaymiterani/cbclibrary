#include <cassert>
#include <iomanip>

// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristicLocal.hpp"
#include "OsiClpSolverInterface.hpp"

extern "C"{

int main(){

	//Objective function
	double* f;
	//Variables to be Integers
	int* intcon;
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

	/*
	 * Problem format: x = myintlinprog(number_of_variables,number_of_constraints,f,intcon,A,b,Aeq,beq,lb,ub,options)
	 */

	OsiClpSolverInterface solver1;  
	OsiSolverInterface *solver = &solver1;

	int numVars_=0, numCons_=0;
	int temp1=0, temp2=0;
	int numintcons;

	//------------Temporary Initialization--------------------------

	double objValue[] = {-3,-2,-1};
	double columnLower[] = {0, 0, 0};
	double columnUpper[] = {COIN_DBL_MAX, COIN_DBL_MAX, 1};
	double rowLower[] = {-COIN_DBL_MAX,-COIN_DBL_MAX,-COIN_DBL_MAX};
	double rowUpper[] = {7,12,-12};
	double Aequality[] = {4,2,1};
	double bequality[] = {12};
	//int column[] = {0,1,0,1,0,1};
	//double element[] = {-1.0,-2.0,-4.0,-1.0,2.0,1.0};
	double conMatrix_[]= {1,4,-4,1,2,-2,1,1,-1};
	int whichInt[]={};
	

	f=objValue;
	A=conMatrix_;
	b=rowUpper;
	intcon=whichInt;
	Aeq=Aequality;
	beq=bequality;

	numVars_=(int)(sizeof(columnUpper)/sizeof(double));
	numCons_=(int)(sizeof(rowUpper)/sizeof(double));
	int numEquality_ = (int)(sizeof(bequality)/sizeof(double));

	//std::cout<<numVars_<<" "<<numCons_<<" "<<numEquality_<<std::endl;

	//--------------------------------------------------------------

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


	solver1.loadProblem(*matrix, columnLower, columnUpper, f, rowLower, b);

	for(int i=0;i<(int)(sizeof(whichInt)/sizeof(int));i++)
		solver1.setInteger(intcon[i]);
	
	solver1.setObjSense(1.0);

	//-------------------------------------------------------------
	
	CbcModel model(solver1);

	model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
	model.branchAndBound();
	
	const double *val = model.getColSolution();
	
	printf("%f\n\nx0 = %f, x1 = %f, x2 = %f\n\n", model.getObjValue(), val[0], val[1],val[2]);
	
	//-------------------------------------------------------------

	return 0;
}
}


