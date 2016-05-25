#include <cassert>
#include <iomanip>

// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristicLocal.hpp"
#include "OsiClpSolverInterface.hpp"

int main(char *fname)
{
  OsiClpSolverInterface solver1;  
  OsiSolverInterface *solver = &solver1;
  
  double objValue[] = {8.0, 1.0};
  double columnLower[] = {-COIN_DBL_MAX, -COIN_DBL_MAX};
  double columnUpper[] = {COIN_DBL_MAX, COIN_DBL_MAX};
  double rowLower[] = {-COIN_DBL_MAX, -COIN_DBL_MAX, -COIN_DBL_MAX};
  double rowUpper[] = {14, -33, 20};
  //int column[] = {0,1,0,1,0,1};
  //double element[] = {-1.0,-2.0,-4.0,-1.0,2.0,1.0};
  double conMatrix_[]= {-1.0, -4.0, 2.0, -2.0, -1.0, 1.0};
  int whichInt[]={1};
  int numberRows=(int) (sizeof(rowLower)/sizeof(double));
  int numberColumns=(int) (sizeof(columnLower)/sizeof(double));
  //int starts[]={0,2,4,6};
    
  //Adding the objective function and the bounds of the individual vars

  //------------------------
  int numVars_=(int)(sizeof(columnUpper)/sizeof(double));
  int numCons_=(int)(sizeof(rowUpper)/sizeof(double));
  CoinPackedMatrix *matrix =  new CoinPackedMatrix(false , 0 , 0);
  matrix->setDimensions(0 , numVars_);
   for(int i=0 ; i<numCons_ ; i++)
    {
      CoinPackedVector row;
      for(int j=0 ; j<numVars_ ; j++)
      {
          row.insert(j, conMatrix_[i+j*numCons_]);
      }
      
        matrix->appendRow(row);
    }
  //------------------------

  solver1.loadProblem(*matrix, columnLower, columnUpper, objValue, rowLower, rowUpper);
  for(int i=0;i<(int)(sizeof(whichInt)/sizeof(int));i++)
    solver1.setInteger(whichInt[i]);
  solver1.setObjSense(1.0);

  CbcModel model(solver1);
  
  model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  model.branchAndBound();
  const double *val = model.getColSolution();
  printf("%f\n\nx0 = %f, x1 = %f\n\n", model.getObjValue(), val[0], val[1]);
  return 0;
} 
