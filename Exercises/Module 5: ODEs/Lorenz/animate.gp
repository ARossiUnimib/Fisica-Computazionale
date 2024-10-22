set term gif animate delay 1 size 1280, 720
set output 'gravitation_animation.gif'  # Output file for the GIF

set xrange [-20:20]  # Adjust these ranges to better fit your data if necessary
set yrange [-30:30]
set zrange [0:50]

set grid
set xlabel "X"
set ylabel "Y"
set zlabel "Z"

# Number of data points (time steps)
n = system("wc -l < lorenz.dat")

# Loop through each time step and plot the positions of the bodies with bold lines
do for [i=0:int(n)-1] {
    splot "lorenz.dat" every ::0::i u 2:3:4 w l lw 2 lc rgb "blue" title 'Body 1', #\
#          "lorenz.dat" every ::0::i u 5:6:7 w l lw 2 lc rgb "green" title 'Body 2', \
#          "lorenz.dat" every ::0::i u 8:9:10 w l lw 2 lc rgb "red" title 'Body 3'
}
