// Solve an MILP problem given in an MPS file
// File p0033.mps taken from the CBC Library - http://www.coin-or.org/Cbc/

path = get_absolute_file_path('example_mps.sce');
[x, f] = cbcintlinprog(path + 'p0033.mps')
clear path;

// f  =
// 
//     3089.  
//  x  =
// 
//     1.  
//     0.  
//     0.  
//     0.  
//     0.  
//     0.  
//     1.  
//     1.  
//     0.  
//     1.  
//     0.  
//     0.  
//     0.  
//     1.  
//     0.  
//     0.  
//     0.  
//     1.  
//     1.  
//     0.  
//     0.  
//     1.  
//     0.  
//     1.  
//     1.  
//     1.  
//     1.  
//     1.  
//     1.  
//     1.  
//     0.  
//     0.  
//     0.  

