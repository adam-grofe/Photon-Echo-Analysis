set terminal pdf enhanced size 6in,6in #size 800,800
set output "2dir.pdf"

set xlabel "Frequency(cm-1)"

set palette color 
set loadpath '~/.gnuplot/palettes/'
load 'jet.pal'
# set view map
# set dgrid3d
set cbrange []
set xrange [1970:2020]
set yrange [1970:2020] 
#set arrow 1 from 1620,1620 to 1680,1680 nohead front lw 2
set mytics
set mxtics
set ytics scale 3,1
set xtics scale 3,1
unset key
do for [i=1:60] {set linetype i lc rgb "black" lw 1 }
set contour surface
set cntrparam levels incremental -1.0, 0.2, 1
set pm3d map  interpolate 5,5

set title "T_{w} = 0.00 ps"
splot  "2dir.prm.lsp" u 1:2:3 w pm3d t ""
