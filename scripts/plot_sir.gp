
set terminal pngcairo size 1200,700
set output 'results/sir_simulation.png'

set datafile separator ','
set title 'SIR Epidemic Simulation'
set xlabel 'Time (days)'
set ylabel 'Population'
set key left top
set grid

plot 'results/sir_simulation.csv' using 1:2 every ::1 with lines lw 2 title 'Susceptible', \
     '' using 1:3 every ::1 with lines lw 2 title 'Infected', \
     '' using 1:4 every ::1 with lines lw 2 title 'Recovered'

set output
