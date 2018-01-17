set terminal x11 enhanced

# set terminal postscript eps color enhanced font ",10"
# set output '../../fig/sam.eps

set datafile separator ','

set grid

set yrange[-1:]

set xlabel "t*"
set ylabel "Max error F"

set title "lambda=0.1; Fs=Fopt*0.75; R=10; J=2; T=5; K=50; epsilon=1e-3;"

plot \
     "Fee.csv" using 1:2 t "1" w l,\
     "Fee.csv" using 1:3 t "2" w l,\
     "Fee.csv" using 1:4 t "3" w l
