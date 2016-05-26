// An example with variable bounds, inequality and equality constraints

c = [-1,-1/3]';
A = [1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1];
b = [2,1,2,1,-1,2];
Aeq = [1,1/4];
beq = [1/2];
lb = [-1,-0.5];
ub = [1.5,1.25];
intcon = [];

[xopt,fopt,exitflag,output] = cbcintlinprog(c, intcon, A, b, Aeq, beq, lb, ub)

// output  =
// 
//   relativegap: 0
//   absolutegap: 0
//   numnodes: 0
//   numfeaspoints: 0
//   numiterations: 1
//   constrviolation: 0
//   message: "Optimal Solution"
// exitflag  =
// 
//    0.  
// fopt  =
// 
// - 0.6041667  
// xopt  =
// 
//   0.1875  
//   1.25  

