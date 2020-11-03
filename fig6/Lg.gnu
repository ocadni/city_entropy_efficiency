set term post color eps enhanced 'Helvetica' 25
set out 'Lg.eps'

set size 2,2
set origin 0,0

set multiplot


#######################

set size 1,1
set origin 0,1

#set title 'L=20' offset 11,-3
set xlabel 'g' font 'Helvetica,30'
set ylabel '{/Symbol h}_{OD}' font 'Helvetica,30'
set key top right samplen 1 spacing 1
plot[:][:] "avC-L20-m3-a0-o16-0.dat" u 1:2 with linespoints pt 4 ps 2 lw 3 t'{/Symbol l}={/Symbol a}=0: L=20',"avC-L30-m3-a0-o24-0.dat" u 1:2 with linespoints pt 4 ps 2 lw 3 t'L=30',"avC-L40-m3-a0-o32-0.dat" u 1:2 with linespoints pt 4 ps 2 lw 3 t'L=40'

#######################

set size 1,1
set origin 1,1

#set title 'L=20' offset 11,-3
set xlabel 'g' font 'Helvetica,30'
set ylabel '{/Symbol D}S_{OD}' font 'Helvetica,30'
set key bottom right samplen 1 spacing 1
plot[:][:] "avC-L20-m3-a0-o16-0.dat" u 1:14 with linespoints pt 4 ps 2 lw 3 t'{/Symbol l}={/Symbol a}=0: L=20',"avC-L30-m3-a0-o24-0.dat" u 1:14 with linespoints pt 4 ps 2 lw 3 t'L=30',"avC-L40-m3-a0-o32-0.dat" u 1:14 with linespoints pt 4 ps 2 lw 3 t'L=40'


#######################


set size 1,1
set origin 0,0

#set title 'L=20' offset 11,-3
set xlabel 'g' font 'Helvetica,30'
set ylabel 'v_{OD}' font 'Helvetica,30'
set key top right samplen 1 spacing 1
plot[:][:] "avC-L20-m3-a0-o16-0.dat" u 1:8 with linespoints pt 4 ps 2 lw 3 t'{/Symbol l}={/Symbol a}=0: L=20',"avC-L30-m3-a0-o24-0.dat" u 1:8 with linespoints pt 4 ps 2 lw 3 t'L=30',"avC-L40-m3-a0-o32-0.dat" u 1:8 with linespoints pt 4 ps 2 lw 3 t'L=40'


#######################


set size 1,1
set origin 1,0

#set title 'L=20' offset 11,-3
set xlabel 'g' font 'Helvetica,30'
set ylabel '{/Symbol P}_{t,t}' font 'Helvetica,30'
set key top right samplen 1 spacing 1
plot[:][:] "avC-L20-m3-a0-o16-0.dat" u 1:16 with linespoints pt 4 ps 2 lw 3 t'{/Symbol l}={/Symbol a}=0: L=20',"avC-L30-m3-a0-o24-0.dat" u 1:16 with linespoints pt 4 ps 2 lw 3 t'L=30',"avC-L40-m3-a0-o32-0.dat" u 1:16 with linespoints pt 4 ps 2 lw 3 t'L=40'


#######################

