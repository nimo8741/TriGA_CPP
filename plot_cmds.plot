cd 'C:\\Users\\nickm\\Documents\\TriGA_CPP\\IO_files\\dat_files'
set size square
plot 'control.dat' with points linecolor rgb 'blue' pointtype 7 ps 0.2 notitle, \
     'edges.dat' with lines linecolor rgb 'red' lw 0.5 notitle
set size square
pause -1 "Hit return to continue"