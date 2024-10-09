set terminal pngcairo size 1280,960 enhanced font 'Verdana,12'
set output 'precision_error.png'

# Title and Labels
set title "Precision Error for Single (Float) and Double (Double) Precision"
set xlabel "Iteration"
set ylabel "Precision Error (log scale)"
set logscale y 10
set grid

# Skip lines with infinite values to avoid distorting the plot
set datafile missing "inf"

# Read the data file (with tab separator)
set datafile separator "\t"

# Plot the errors: f_div and d_div errors on the same graph
plot "output.txt" using 0:($6-1) title "Float Precision Error (1 + f_div - 1)" with linespoints lt rgb "red", \
     "output.txt" using 0:($4-1) title "Double Precision Error (1 + d_div - 1)" with linespoints lt rgb "blue"

