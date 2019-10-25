# Use this as a template to graph results
set terminal pdf enhance font "palatino,10" size 4,4

set output 'pr_grains.pdf'
set title  "Sand grain count PDF"

set palette rgb 21,22,23

set autoscale xfix
set autoscale yfix
set autoscale cbfix

set datafile missing '0'

#xmin = 
#xmax = GPVAL_X_MAX
#set xrange [xmin:xmax] noreverse nowriteback

plot 'pr_grains.dat' using 1:($2 == 0 ? NaN : $2) with lines ls 1 notitle
