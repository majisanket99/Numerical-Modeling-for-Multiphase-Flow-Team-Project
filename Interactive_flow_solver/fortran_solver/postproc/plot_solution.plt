#!/usr/bin/gnuplot -persist

set terminal x11 enhanced persist size 900,900

stats '../results/info.dat'
N  = STATS_mean_x
Nk = STATS_mean_y

#set terminal wxt size 1366,768

c=1
while (c==1) {

pause 0.06
#replot
load "layout.plt"

}

# EOF
