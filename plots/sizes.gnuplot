# Use this as a template to graph results
set terminal pdf enhance font "palatino,10" size 4,4

set output 'sizes.pdf'
set title  "Avalanche sizes at each iteration"

set palette rgb 21,22,23

set autoscale xfix
set autoscale yfix
set autoscale cbfix

plot 'timecount.dat' using 1:3 with lines ls 2 notitle
