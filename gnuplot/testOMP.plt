#!/gnuplot

set terminal png
set output 'testOMP.png'
set xlabel 'T'
set ylabel 'x'
set zlabel 'U'
set view 60, 70
set xyplane 0
splot 'out.txt' u 1:2:3 w pm3d t "Solution development by time"
exit
