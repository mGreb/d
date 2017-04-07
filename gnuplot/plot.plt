#!/gnuplot

set terminal png
set output 'plot.png'
set xlabel "x"
set ylabel "U"
plot 'out.txt' u 2:3 w l t "Solution development by time"

set terminal png
set output 'grunvald.png'
set xlabel "n"
set ylabel "coeff"
plot 'grunvald.txt' u 1:2 w l t "Grunvald-Letnikov coefficients", 'grunvald.txt' u 1:3 w l t "Grunvald-Letnikov coefficients absolute value"

set terminal png
set output 'hyper.png'
set xlabel "x"
set ylabel "func"
plot 'hyper.txt' u 1:2 w l t "2F1(alpha, alpha, alpha + 1, x)"
exit
