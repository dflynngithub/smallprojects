# Use this as a template to graph results
set terminal pdf enhance font "palatino,10" size 4,4

set output 'grains.pdf'
set title  "Total sand grain counter over time"

set palette rgb 21,22,23

set autoscale xfix
set autoscale yfix
set autoscale cbfix

plot 'timecount.dat' using 1:2 with lines ls 1 notitle
