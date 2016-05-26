// An example with inequality and equality constraints

c = [-1,-1/3]';
A = [1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1];
b = [2,1,2,1,-1,2];
Aeq = [1,1/4];
beq = [1/2];
intcon = [1, 2];

[xopt,fopt,exitflag,output] = cbcintlinprog(c, intcon, A, b, Aeq, beq)

//output  =
 
//   relativegap: 0
//   absolutegap: 0
//   numnodes: 0
//   numfeaspoints: 2
//   numiterations: 0
//   constrviolation: 0
//   message: "Optimal Solution"
// exitflag  =
// 
//    0.  
// fopt  =
// 
//  - 0.6666667  
// xopt  =
// 
//    0.  
//    2.  




