set term post color eps enhanced 'Helvetica' 25
set out 'Ti.eps'

set size 2,2
set origin 0,0

set multiplot

set logscale x
set logscale y


#a               = 27236.2          +/- 1752         (6.431%)
#b               = 2.43946          +/- 0.01167      (0.4783%)
f1(x)=27236.2/x**2.43946

#a               = 6.09583          +/- 0.1037       (1.702%)
#b               = 1.9645           +/- 0.003616     (0.1841%)
f2(x)=100*6.09583/x**1.9645

#a               = 13.037           +/- 0.2646       (2.03%)
#b               = 0.804117         +/- 0.003315     (0.4122%)
f3(x)=13.037/x**0.804117

#a               = 1213.5           +/- 64.48        (5.314%)
#b               = 1.74543          +/- 0.009365     (0.5365%)
f4(x)=1213.5/x**1.74543


#######################

set size 1,1
set origin 0,1

set title '(a) {/Symbol l}=0.0: T^{-2.43}' offset 11,-3
set xlabel 'T' font 'Helvetica,30'
set ylabel 'CPD(T)' font 'Helvetica,30'
set key bottom left samplen 1 spacing 1
plot[:2000][1e-3:1] "Ti-L40-a0-o16-0.dat" u 1:2 pt 4 ps 2 t'{/Symbol D}T_O=16: g=0.5',"Ti-L40-a0-o16-0.dat" u 1:3 pt 6 ps 2 t'g=2.0',f1(x) t'' lc 0 lw 6 

#######################

set size 1,1
set origin 1,1


set title '(b) {/Symbol l}=0.0: T^{-1.96}' offset 11,-3
set xlabel 'T' font 'Helvetica,30'
set ylabel 'CPD(T)' font 'Helvetica,30'
set key bottom left samplen 1 spacing 1
plot[:2000][1e-3:1] "Ti-L40-a0-o32-0.dat" u 1:2 pt 4 ps 2 t'{/Symbol D}T_O=32: g=0.5',"Ti-L40-a0-o32-0.dat" u 1:3 pt 6 ps 2 t'g=2.0',f2(x) t'' lc 0 lw 6 


#######################

set size 1,1
set origin 0,0


set title '(c) {/Symbol l}=0.5: T^{-0.80}' offset 11,-3
set xlabel 'T' font 'Helvetica,30'
set ylabel 'CPD(T)' font 'Helvetica,30'
set key bottom left samplen 1 spacing 1
plot[:20000][1e-3:1] "Ti-L40-a0-o16.dat" u 1:2 pt 4 ps 2 t'{/Symbol D}T_O=16: g=0.5',"Ti-L40-a0-o16.dat" u 1:3 pt 6 ps 2 t'g=2.0',f3(x) t'' lc 0 lw 6 



#######################

set size 1,1
set origin 1,0


set title '(d) {/Symbol l}=0.5: T^{-1.74}' offset 11,-3
set xlabel 'T' font 'Helvetica,30'
set ylabel 'CPD(T)' font 'Helvetica,30'
set key bottom left samplen 1 spacing 1
plot[:3000][1e-3:1] "Ti-L40-a0-o32.dat" u 1:2 pt 4 ps 2 t'{/Symbol D}T_O=32: g=0.5',"Ti-L40-a0-o32.dat" u 1:3 pt 6 ps 2 t'g=2.0',f4(x) t'' lc 0 lw 6 


