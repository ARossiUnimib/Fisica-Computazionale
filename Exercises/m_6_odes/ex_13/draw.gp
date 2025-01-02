set terminal pngcairo size 800,400 enhanced font 'Verdana,15'
set output "1.png"

filename = "data/3_0.5_0.5.dat"

set size square 
set multiplot layout 1,2
set xrange [0:3*pi]
plot filename u 1:2 w l t "posizione", filename u 1:3 w l t "velocita'"

set xrange [-5:5]

plot filename u 2:3 w l t "traiettoria"
