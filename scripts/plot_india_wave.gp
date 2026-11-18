set terminal pngcairo size 1200,700
set output 'results/india_second_wave.png'

set datafile separator ','
set xdata time
set timefmt '%Y-%m-%d'
set format x '%d %b'
set xtics rotate by -45
set title 'India COVID-19 Second Wave (March-June 2021)'
set xlabel 'Date'
set ylabel 'Smoothed daily reported cases'
set grid
set key off

plot 'data/india_second_wave.csv' using 1:3 every ::1 with lines lw 2 title 'Smoothed cases'

set output
