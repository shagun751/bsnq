set terminal pngcairo
set output 'vol.png'
plot 'volumeData.dat' using 1:2 with lines
exit