// Error Check Example
c = [-1,-1/3]';
A = [1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1];
intcon = [];

x = cbcintlinprog(c, intcon, A)

// !--error 10000 
//cbcintlinprog: Unexpected number of input arguments : 3 provided while should be in the set [4 6 8 9]
//at line     135 of function cbcintlinprog called by :  
//x = cbcintlinprog(c, intcon, A)
