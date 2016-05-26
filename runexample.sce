exec doitagain.sce
c=[-1,-1/3]'
opt=list('MaxNodes',200);
[xopt,fopt,exitflag,output]=cbcintlinprog('p0033.mps',opt)
