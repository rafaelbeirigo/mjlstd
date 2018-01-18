set terminal x11 enhanced

# set terminal postscript eps color enhanced font ",10"
# set output '../../fig/sam.eps

set datafile separator ','

set grid

# set yrange[-1:]

set xlabel "t*"
set ylabel "Max error F"
set style fill transparent solid 0.2 noborder
set title "lambda=0.1; Fs=Fopt*0.75; R=10; J=2; T=5; K=50; epsilon=1e-3;"

plot \
     "Fee_avg_std.csv" using 1:($2-$5):($2+$5) t "" w filledcurves lc 0,\
     "Fee_avg_std.csv" using 1:($2) smooth csplines t "" w lp pointinterval 12.5 lt 1 pt 2 lc 0,\
     "Fee_avg_std.csv" using 1:($3-$6):($3+$6) t "" w filledcurves lc 0,\
     "Fee_avg_std.csv" using 1:($3) smooth csplines t "" w lp pointinterval 12.5 lt 1 pt 2 lc 1,\
     "Fee_avg_std.csv" using 1:($4-$7):($4+$7) t "" w filledcurves lc 0,\
     "Fee_avg_std.csv" using 1:($4) smooth csplines t "" w lp pointinterval 12.5 lt 1 pt 2 lc 2
