// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Guru Pradeep Reddy, Bhanu Priya Sayal
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in


function [xopt,fopt,exitflag,output] = cbcintlinprog (varargin)
	// Solves a linear programming problem.
	//
	//   Calling Sequence
	//   xopt = cbcintlinprog(c,A,b)
	//   xopt = cbcintlinprog(c,A,b,Aeq,beq)
	//   xopt = cbcintlinprog(c,A,b,Aeq,beq,lb,ub)
	//   xopt = cbcintlinprog(c,A,b,Aeq,beq,lb,ub,param)
	//   xopt = cbcintlinprog(file)
	//   xopt = cbcintlinprog(file,param)
	//   [xopt,fopt,exitflag,output,lambda] = cbcintlinprog( ... )
	//   
	//   Parameters
	//   c : a vector of double, contains coefficients of the variables in the objective 
	//   A : a matrix of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b. 
	//   b : a vector of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b.
	//   Aeq : a matrix of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   beq : a vector of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   lb : Lower bounds, specified as a vector or array of double. lb represents the lower bounds elementwise in lb ≤ x ≤ ub.
	//   ub : Upper bounds, specified as a vector or array of double. ub represents the upper bounds elementwise in lb ≤ x ≤ ub.
	//   options : a list containing the parameters to be set.
	//   file : a string describing the path to the mps file.
	//   xopt : a vector of double, the computed solution of the optimization problem.
	//   fopt : a double, the value of the function at x.
	//   status : status flag returned from symphony. See below for details.
	//   output : The output data structure contains detailed information about the optimization process. See below for details.
	//   lambda : The structure consist of the Lagrange multipliers at the solution of problem. See below for details.
	//   
	//   Description
	//   OSI-CLP is used for solving the linear programming problems, OSI-CLP is a library written in C++.
	//   Search the minimum of a constrained linear programming problem specified by :
	//
	//   <latex>
	//    \begin{eqnarray}
	//    &\mbox{min}_{x}
	//    & c^T⋅x  \\
	//    & \text{subject to} & A⋅x \leq b \\
	//    & & Aeq⋅x = beq \\
	//    & & lb \leq x \leq ub \\
	//    \end{eqnarray}
	//   </latex>
	//
	//   The routine calls Clp for solving the linear programming problem, Clp is a library written in C++.
	//
  	// The options allows the user to set various parameters of the Optimization problem. 
  	// It should be defined as type "list" and contains the following fields. In the current version it only contains maxiter.
	// <itemizedlist>
	//   <listitem>Syntax : options= list("MaxIter", [---]);</listitem>
	//   <listitem>MaxIter : a Scalar, containing the Maximum Number of Iteration that the solver should take.</listitem>
	//   <listitem>Default Values : options = list("MaxIter", [3000]);</listitem>
	// </itemizedlist>
	//
	// The exitflag allows to know the status of the optimization which is given back by CLP.
	// <itemizedlist>
	//   <listitem>exitflag=0 : Optimal Solution Found </listitem>
	//   <listitem>exitflag=1 : Primal Infeasible </listitem>
	//   <listitem>exitflag=2 : Dual Infeasible</listitem>
	//   <listitem>exitflag=3 : Maximum Number of Iterations Exceeded. Output may not be optimal.</listitem>
	//   <listitem>exitflag=4 : Solution Abandoned</listitem>
	//   <listitem>exitflag=5 : Primal objective limit reached.</listitem>
	//   <listitem>exitflag=6 : Dual objective limit reached.</listitem>
	// </itemizedlist>
	// 
	// The output data structure contains detailed informations about the optimization process. 
	// It has type "struct" and contains the following fields.
	// <itemizedlist>
	//   <listitem>output.iterations: The number of iterations performed during the search</listitem>
	//   <listitem>output.constrviolation: The max-norm of the constraint violation.</listitem>
	// </itemizedlist>
	//
	// The lambda data structure contains the Lagrange multipliers at the end 
	// of optimization. In the current version the values are returned only when the the solution is optimal. 
	// It has type "struct" and contains the following fields.
	// <itemizedlist>
	//   <listitem>lambda.lower: The Lagrange multipliers for variable lower bounds.</listitem>
	//   <listitem>lambda.eqlin: The Lagrange multipliers for the linear equality constraints.</listitem>
	//   <listitem>lambda.ineqlin: The Lagrange multipliers for the linear inequality constraints.</listitem>
	// </itemizedlist>
	//
	// Examples
	// //Optimal problems
	// //Linear program, linear inequality constraints
	// c=[-1,-1/3]'
	// A=[1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1]
	// b=[2,1,2,1,-1,2]
	// [xopt,fopt,exitflag,output,lambda]=cbcintlinprog(c, A, b)
	// // Press ENTER to continue
	//
	// Examples
	// //Linear program with Linear Inequalities and Equalities`
	// c=[-1,-1/3]'
	// A=[1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1]
	// b=[2,1,2,1,-1,2]  
	// Aeq=[1,1/4]
	// beq=[1/2]
	// [xopt,fopt,exitflag,output,lambda]=cbcintlinprog(c, A, b, Aeq, beq)
	// // Press ENTER to continue
	//
	// Examples
	// //Linear program with all constraint types 
	// c=[-1,-1/3]'
	// A=[1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1]
	// b=[2,1,2,1,-1,2]   
	// Aeq=[1,1/4]
	// beq=[1/2]
	// lb=[-1,-0.5]
	// ub=[1.5,1.25]
	// [xopt,fopt,exitflag,output,lambda]=cbcintlinprog(c, A, b, Aeq, beq, lb, ub)
	// // Press ENTER to continue
	//
	// Examples 
	// //Primal Infeasible Problem
	// c=[-1,-1,-1]'
	// A=[1,2,-1]
	// b=[-4]
	// Aeq=[1,5,3;1,1,0]
	// beq=[10,100]
	// lb=[0,0,0]
	// ub=[%inf,%inf,%inf]
	// [xopt,fopt,exitflag,output,lambda]= cbcintlinprog(c,A,b,Aeq,beq,lb,ub)
	// // Press ENTER to continue
	//
	// Examples
	// //Dual Infeasible Problem
	// c=[3,5,-7]'
	// A=[-1,-1,4;1,1,4]
	// b=[-8,5]
	// Aeq=[]
	// beq=[]
	// lb=[-%inf,-%inf,-%inf]
	// ub=[%inf,%inf,%inf]
	// [xopt,fopt,exitflag,output,lambda]= cbcintlinprog(c,A,b,Aeq,beq,lb,ub)
	// // Press ENTER to continue
	//
	// Examples
	// filepath = get_absolute_file_path('cbcintlinprog.dem.sce');
	// filepath = filepath + "exmip1.mps"
	// [xopt,fopt,exitflag,output,lambda] =cbcintlinprog(filepath)
	// Authors
	// Bhanu Priya Sayal, Guru Pradeep Reddy

	if(type(varargin(1))==1) then
      
      //To check the number of input and output argument
   	[lhs , rhs] = argn();
	
	//To check the number of argument given by user
	if ( rhs < 4 | rhs == 5 | rhs == 7 | rhs > 9 ) then
		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set [4 6 8 9]"), "cbcintlinprog", rhs);
		error(errmsg);
	end
   
	c = [];
	intcon = [];
	A = [];
	b = [];
	Aeq = [];
	beq = [];
	lb = [];
	ub = [];
	options = list();
	
	c = varargin(1)
	intcon = varargin(2)
	A = varargin(3)
	b = varargin(4)

	if(size(c,2) == 0) then
		errmsg = msprintf(gettext("%s: Cannot determine the number of variables because input objective coefficients is empty"),"cbcintlinprog");
		error(errmsg);
	end

   if (size(c,2)~=1) then
	errmsg = msprintf(gettext("%s: Objective Coefficients should be a column matrix"), "cbcintlinprog");
	error(errmsg);
   end

   nbVar = size(c,1);

   if ( rhs<5 ) then
      Aeq = []
      beq = []
   else
      Aeq = varargin(5);
      beq = varargin(6);
   end
   
   if ( rhs<7 ) then
      lb = repmat(-%inf,1,nbVar);
      ub = repmat(%inf,1,nbVar);
   else
      lb = varargin(7);
      ub = varargin(8);
   end
   
   if (rhs<9|size(varargin(9))==0) then
      options = list();
   else
      options = varargin(9);
   end
		[xopt,fopt,exitflag,output]=cbcmatrixintlinprog(c,intcon,A,b,Aeq,beq,lb,ub,options);	
	elseif(type(varargin(1))==10) then

		[lhs , rhs] = argn();

		//To check the number of argument given by user
	   	if ( rhs < 1 | rhs > 2) then
			errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set of [1 2]"),"cbcintlinprog",rhs);
			error(errmsg)
	   	end
	   	mpsFile = varargin(1);
		if ( rhs<2 | size(varargin(2)) ==0 ) then
		  param = list();
		else
		   param =varargin(2);
		end 
		[xopt,fopt,exitflag,output]=cbcmpsintlinprog(mpsFile,param);	
	end

endfunction
