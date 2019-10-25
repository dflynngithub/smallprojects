#regions.gnuplot
#
#  Usage:
#  gnuplot < regions.gnuplot
#
#
# Terminal output specs
set terminal pdf size 4,4
set output "plots/regions.pdf"
#
load "plots/pals/jet.pal"
#
# Axes and title
set title sprintf("regions")
set xrange [-0.5:  49.5]
set yrange [  49.5:-0.5] reverse
#
# Plot data to file
plot "plots/regions.dat" matrix with image notitle
