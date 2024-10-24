set term gif animate delay 2.5 size 1920,1080 # Create an animated GIF with a 10ms delay
set output 'gravitation_animation.gif'  # Output file for the GIF

set xrange [-3:3]`` 
set yrange [-3:3]
set zrange [0:8]

set grid
set xlabel "X"
set ylabel "Y"
set zlabel "Z"

# Number of data points (time steps)
n = system("wc -l < gravitation_2.dat")

# Loop through each time step and plot the positions of the bodies with bold lines
do for [i=0:int(n)-1] {
    splot "gravitation_1.dat" every ::0::i u 2:3:4 w l lw 2 lc rgb "blue" title 'Body 1', \
          "gravitation_1.dat" every ::0::i u 5:6:7 w l lw 2 lc rgb "green" title 'Body 2', \
          "gravitation_1.dat" every ::0::i u 8:9:10 w l lw 2 lc rgb "red" title 'Body 3'
}
