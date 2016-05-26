// Primal Infeasible Problem
c = [-1,-1,-1]';
A =[1,2,-1];
b = [-4];
Aeq = [1,5,3;1,1,0];
beq = [10,100];
lb = [0,0,0];
ub = [%inf,%inf,%inf];
intcon = []

[xopt,fopt,exitflag,output] = intlinprog(c, intcon, A, b, Aeq, beq, lb, ub)

//output  =
// 
//   relativegap: 0
//   absolutegap: 0
//   numnodes: 0
//   numfeaspoints: 0
//   numiterations: 0
//   constrviolation: 100
//   message: "Primal Infeasible"
// exitflag  =
// 
//    1.  
// fopt  =
// 
//    1.79D+308  
// xopt  =
// 
//    0.         
//    0.         
//    3.3333333  

// Note the message "Primal Infeasible" and the exitflag(which is 1)

