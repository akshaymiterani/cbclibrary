exec doitagain.sce
f = [-3;-2;-1];
intcon = 3;
A = [1,1,1];
b = 7;
Aeq = [4,2,1];
beq = 12;
lb = zeros(3,1);
ub = [%inf;%inf;1];
options=list();
[xopt,fopt,exitflag,output]=cbcintlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options)
