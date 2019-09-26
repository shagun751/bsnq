set terminal pngcairo size 1280, 640 linewidth 3
#set yrange [-0.03:0.03]
set xrange [0:30]
#set ytics 0.96,0.005,1.04
#set ytics 0.02
set xtics 5
set grid ytics lt 0 lw 1
set grid xtics lt 0 lw 1
show grid
set palette model HSV
do for [ii=1:6] {
	fileN=sprintf("waveProbes.dat")
	fileN3=sprintf("waveProbes.dat")	
	fileN2=sprintf("WP%d.png",ii)
	set yrange [-0.03:0.03]
	set output fileN2
	iib=5*(ii-1)+1+2
	plot fileN3 using 1:iib with lines, fileN using 1:iib with lines
	#plot fileN using 1:iib with lines
}
