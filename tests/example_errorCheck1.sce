// Error Check
c = [-1,-1/3]';
A = [1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1];
b = [2,1,2,1,-1,2];
Aeq = [1,1/4];
beq = [1/2];

x = cbcintlinprog(c, [],A, b, Aeq, Aeq)

// !--error 6 
//Inconsistent row/column dimensions.
//at line     208 of function cbcmatrixintlinprog called by :  
//at line     186 of function cbcintlinprog called by :  
//x = cbcintlinprog(c, [],A, b, Aeq, Aeq)

