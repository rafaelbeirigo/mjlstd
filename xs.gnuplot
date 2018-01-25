set terminal x11 enhanced

# set terminal postscript eps color enhanced font ',10'
# set output '../../fig/fig.eps

set datafile separator ','

set grid

set xlabel 't'
set ylabel 'x(2,:) (mean)'
set style fill transparent solid 0.25 noborder
set title ''

plot \
     "xs_opt.csv" using 1:($2-$3):($2+$3) t "" w filledcurves lc 0,\
     "xs.csv" using 1:($2-$3):($2+$3) t "" w filledcurves lc 1,\
     'xs_opt.csv' using 1:($2) smooth csplines t '' w lp pointinterval 1 lt 1 pt 1 lc 0,\
     'xs.csv' using 1:($2) smooth csplines t '' w lp pointinterval 1 lt 1 pt 2 lc 1
