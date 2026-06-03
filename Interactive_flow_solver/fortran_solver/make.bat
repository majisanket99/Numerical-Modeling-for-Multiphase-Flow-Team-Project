cls
@echo off
g95 fftpack5.1.o module_globalvariables.F90 module_fluid.F90 module_immersedboundary.F90 initialize.F90 output.F90 solve.F90 deallocation.F90 IBM2D.F90 -O3 -w -o IBM2D.exe
echo\
echo  IBM2D.exe wurde erstellt...