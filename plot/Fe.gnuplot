set terminal x11 enhanced

# set terminal postscript eps color enhanced font ",10"
# set output '../../fig/sam.eps

set datafile separator ','

set grid

plot \
     "Fe.csv" using 1:2 t "" w l,\
     "Fe.csv" using 1:3 t "" w l,\
     "Fe.csv" using 1:4 t "" w l
     # "itTransf_avg_std.csv" using 1:($2-$3):($2+$3) t "" w filledcurves lc 0,\
     # "itTransf_avg_std.csv" smooth csplines t "With transfer" w lp pointinterval 12.5 lt 1 pt 2 lc 3,\
     # "itAlg_avg_std.csv" smooth csplines t "Without transfer" w lp pointinterval 12.5 lt 1 pt 1 lc 1
