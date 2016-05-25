#ifdef __cplusplus
extern "C" {
#endif
#include <mex.h> 
#include <sci_gateway.h>
#include <api_scilab.h>
#include <MALLOC.h>
static int direct_gateway(char *fname,void F(void)) { F();return 0;};
extern Gatefunc matrix_cppintlinprog;
extern Gatefunc mps_cppintlinprog;
static GenericTable Tab[]={
  {(Myinterfun)sci_gateway,matrix_cppintlinprog,"sci_matrix_intlinprog"},
  {(Myinterfun)sci_gateway,mps_cppintlinprog,"sci_mps_intlinprog"},
};
 
int C2F(libFOSSEE_Scilab_Cbc)()
{
  Rhs = Max(0, Rhs);
  if (*(Tab[Fin-1].f) != NULL) 
  {
     if(pvApiCtx == NULL)
     {
       pvApiCtx = (StrCtx*)MALLOC(sizeof(StrCtx));
     }
     pvApiCtx->pstName = (char*)Tab[Fin-1].name;
    (*(Tab[Fin-1].f))(Tab[Fin-1].name,Tab[Fin-1].F);
  }
  return 0;
}
#ifdef __cplusplus
}
#endif
