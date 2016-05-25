mode(-1)

toolboxTitle = 'FOSSEE_Scilab_Cbc';

path = get_absolute_file_path('builder.sce');
include_path = path + 'build/Cbc/include/coin';
lib_path = path + 'build/Cbc/lib';
sci_path = path + 'build/Scilab/';

C_Flags = [' -w -fpermissive -I ' + include_path + ' -I -Wl,-rpath=' + lib_path];
Linker_Flag = ['-L' + lib_path + 'libCoinUtils.so ' + lib_path + 'libClp.so ' + lib_path + 'libClpSolver.so ' + lib_path + 'libOsi.so ' + lib_path + 'libOsiClp.so ' + lib_path +  'libCgl.so ' + 'libCbc.so ' + lib_path + 'libCbcSolver.so ' + lib_path + 'libOsiCbc.so'];

files = ['cbcintlinprog_matrixcpp.cpp','cbcintlinprog_mpscpp.cpp', 'sci_iofunc.hpp', 'sci_iofunc.cpp'] ;
functionNames = ['sci_matrix_intlinprog', 'matrix_cppintlinprog',
		'sci_mps_intlinprog','mps_cppintlinprog'];

tbx_build_gateway(toolboxTitle, functionNames, files,path, [], Linker_Flag, C_Flags, [], 'g++');

clear WITHOUT_AUTO_PUTLHSVAR toolboxTitle functionNames files Linker_Flag C_Flags;

link(lib_path + '/libCoinUtils.so');
link(lib_path + '/libClp.so');
link(lib_path + '/libClpSolver.so');
link(lib_path + '/libOsi.so');
link(lib_path + '/libOsiClp.so');
link(lib_path + '/libCgl.so');
link(lib_path + '/libCbc.so');
link(lib_path + '/libCbcSolver.so');
link(lib_path + '/libOsiCbc.so');
    
exec loader.sce

disp('---------------------------------------------------------------------------------')
disp('Done')
