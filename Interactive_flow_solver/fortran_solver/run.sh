## compile
#make fft
make ibm2d

## write cursor position to targetPos/targetPos.dat
## by system call using xdotool in getTargetPos() in module_immersedboundary
## OR by kinect python script
#xterm -e "cd targetPos && python writeTargetPosByKinect.py" &
#sleep 0.1

## execute solver and show result
xterm -e ./ibm2d.out &
sleep 0.1
xterm -e "cd postproc && gnuplot plot_solution.plt"
echo 'Done.'
