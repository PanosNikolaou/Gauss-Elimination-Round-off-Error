set xlabel "dimension"
set ylabel "error"
set logscale y
m="./data.dat"
set terminal qt size 1200, 700
set nokey
set grid
set title 'Gaussian Elimination Roundoff Error'
plot m using 2:1 with linespoints
replot
