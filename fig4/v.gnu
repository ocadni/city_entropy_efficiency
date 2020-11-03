set term post color eps enhanced 'Helvetica' 25
set out 'v.eps'

set size 2,1
set origin 0,0

set multiplot


#######################

set size 1,1
set origin 0,0

#set title 'L=20' offset 11,-3
set xlabel '{/Symbol a}' font 'Helvetica,30'
set ylabel 'v_{OD}' font 'Helvetica,30'
set key top right samplen 1 spacing 1
plot[:][:] "avC-L20-m3-g05-o16-0.dat" u 1:8:9 with linespoints pt 4 ps 2 lw 3 t'{/Symbol l}=0: g=0.5',"avC-L20-m3-g1-o16-0.dat" u 1:8:9 with linespoints pt 6 ps 2 lw 3 t'g=1.0',"avC-L20-m3-g2-o16-0.dat" u 1:8:9 with linespoints pt 8 ps 2 lw 3 t'g=2.0',"avC-L20-m3-g4-o16-0.dat" u 1:8:9 with linespoints pt 10 ps 2 lw 3 t'g=4.0'

#######################

set size 1,1
set origin 1,0

#set title 'L=20' offset 11,-3
set xlabel '{/Symbol a}' font 'Helvetica,30'
set ylabel 'v_{OD}' font 'Helvetica,30'
set key top right samplen 1 spacing 1
plot[:0.2][:] "avC-L20-m3-g05-o16.dat" u 1:8:9 with linespoints pt 4 ps 2 lw 3 t'{/Symbol l}=0.5: g=0.5',"avC-L20-m3-g10-o16.dat" u 1:8:9 with linespoints pt 6 ps 2 lw 3 t'g=1.0',"avC-L20-m3-g20-o16.dat" u 1:8:9 with linespoints pt 8 ps 2 lw 3 t'g=2.0',"avC-L20-m3-g40-o16.dat" u 1:8:9 with linespoints pt 10 ps 2 lw 3 t'g=4.0'

