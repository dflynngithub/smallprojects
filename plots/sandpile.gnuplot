# Use this as a template to graph results
set terminal pdf size 4,4
set output 'sandpile.pdf'
set title  "Total sand grain counter over time"

set palette rgb 21,22,23

set autoscale xfix
set autoscale yfix
set autoscale cbfix

plot 'sandpile.dat' matrix with image notitle
