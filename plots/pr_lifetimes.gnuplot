# Use this as a template to graph results
set terminal pdf enhance font "palatino,10" size 4,4

set output 'pr_lifetimes.pdf'
set title  "Avalanche lifetime PDF (log y scale)"

set logscale y

set palette rgb 21,22,23

set autoscale xfix
set autoscale yfix
set autoscale cbfix

plot 'pr_lifetimes.dat' using 1:2 with lines ls 3 notitle
