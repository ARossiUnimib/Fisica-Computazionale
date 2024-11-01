set terminal pngcairo size 800,400 enhanced font 'Verdana,15'
set output image

set size square 
set multiplot layout 1,2
set xrange [0:6*pi]
plot filename u 1:2 w l t "posizione", filename u 1:3 w l t "velocita'"

set yrange [-1.1:1.1]
set xrange [-1.1:1.1]

plot filename u 2:3 w l t "traiettoria"
