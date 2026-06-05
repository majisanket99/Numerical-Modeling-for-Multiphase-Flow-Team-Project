# Fluid-friendly visualization layout

unset key
unset grid

set size ratio -1
set xrange [0:N]
set yrange [0:N]

set xlabel "x"
set ylabel "y"

set colorbox
set cblabel "Fluid field"

# Color map
set palette rgbformulae 33,13,10

plot \
'../results/fluid.bin' binary array=(N,N) format='%float' with image, \
'../results/ib.dat' using 1:2 with lines lw 3 lc rgb "black"
#'../targetPos/targetPos.dat' using ($1/1920*N):((1.0-$2/1080)*N) with points pt 7 ps 1.5 lc rgb "white"
