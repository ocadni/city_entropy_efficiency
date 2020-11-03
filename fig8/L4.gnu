set term post color eps enhanced 'Helvetica' 25
set out 'L4.eps'

set size 2,2
set origin 0,0

set multiplot


#######################

set size 1,1
set origin 0,1

#set title 'L=20' offset 11,-3
set xlabel '{/Symbol a}' font 'Helvetica,30'
set ylabel '{/Symbol h}_{OD}' font 'Helvetica,30'
set key top right samplen 1 spacing 1
plot[:][:] "avC-L40-m3-g0-o32-0r.dat" u 1:2 with linespoints pt 4 ps 2 lw 3 t'{/Symbol l}=0:g=0.0',"avC-L40-m3-g05-o32-0r.dat" u 1:2 with linespoints pt 4 ps 2 lw 3 t'g=0.5',"avC-L40-m3-g10-o32-0r.dat" u 1:2 with linespoints pt 4 ps 2 lw 3 t'g=1.0',"avC-L40-m3-g40-o32-0r.dat" u 1:2 with linespoints pt 4 ps 2 lw 3 t'g=4.0'

#######################

set size 1,1
set origin 1,1

#set title 'L=20' offset 11,-3
set xlabel '{/Symbol a}' font 'Helvetica,30'
set ylabel '{/Symbol D}S_{OD}' font 'Helvetica,30'
set key top left samplen 1 spacing 1
plot[:][:] "avC-L40-m3-g0-o32-0r.dat" u 1:14 with linespoints pt 4 ps 2 lw 3 t'{/Symbol l}=0:g=0.0',"avC-L40-m3-g05-o32-0r.dat" u 1:14 with linespoints pt 4 ps 2 lw 3 t'g=0.5',"avC-L40-m3-g10-o32-0r.dat" u 1:14 with linespoints pt 4 ps 2 lw 3 t'g=1.0',"avC-L40-m3-g40-o32-0r.dat" u 1:14 with linespoints pt 4 ps 2 lw 3 t'g=4.0'


#######################


set size 1,1
set origin 0,0

#set title 'L=20' offset 11,-3
set xlabel '{/Symbol a}' font 'Helvetica,30'
set ylabel 'v_{OD}' font 'Helvetica,30'
set key top right samplen 1 spacing 1
plot[:][:] "avC-L40-m3-g0-o32-0r.dat" u 1:8 with linespoints pt 4 ps 2 lw 3 t'{/Symbol l}=0:g=0.0',"avC-L40-m3-g05-o32-0r.dat" u 1:8 with linespoints pt 4 ps 2 lw 3 t'g=0.5',"avC-L40-m3-g10-o32-0r.dat" u 1:8 with linespoints pt 4 ps 2 lw 3 t'g=1.0',"avC-L40-m3-g40-o32-0r.dat" u 1:8 with linespoints pt 4 ps 2 lw 3 t'g=4.0'


#######################


set size 1,1
set origin 1,0

#set title 'L=20' offset 11,-3
set xlabel '{/Symbol a}' font 'Helvetica,30'
set ylabel '{/Symbol P}_{t,t}' font 'Helvetica,30'
set key top right samplen 1 spacing 1
plot[:][:] "avC-L40-m3-g0-o32-0r.dat" u 1:16 with linespoints pt 4 ps 2 lw 3 t'{/Symbol l}=0:g=0.0',"avC-L40-m3-g05-o32-0r.dat" u 1:16 with linespoints pt 4 ps 2 lw 3 t'g=0.5',"avC-L40-m3-g10-o32-0r.dat" u 1:16 with linespoints pt 4 ps 2 lw 3 t'g=1.0',"avC-L40-m3-g40-o32-0r.dat" u 1:16 with linespoints pt 4 ps 2 lw 3 t'g=4.0'


#######################

