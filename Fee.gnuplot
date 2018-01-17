set terminal x11 enhanced

# set terminal postscript eps color enhanced font ",10"
# set output '../../fig/sam.eps

set datafile separator ','

set grid

plot \
     "Fee.csv" using 1:2 t "1" w l,\
     "Fee.csv" using 1:3 t "2" w l,\
     "Fee.csv" using 1:4 t "3" w l
