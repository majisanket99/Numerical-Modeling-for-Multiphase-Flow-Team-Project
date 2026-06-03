#!/usr/bin/gnuplot -persist

stats '../results/info.dat'
N  = STATS_mean_x
Nk = STATS_mean_y

set xrange [1:N]
set yrange [-0.1:1]
set xtics 5
set ytics 0.1
set xzeroaxis

set title 'Line plot'
set xlabel 'grid index'
set ylabel '|u|'
show xlabel
show ylabel

plot '../results/slice.dat' with linespoints


c=1
while (c==1) {
pause 0.06
replot
}

# EOF
