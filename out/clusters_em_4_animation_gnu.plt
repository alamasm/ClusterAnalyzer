set term gif animate optimize delay 100 size 600, 600 background "#ffeedf" crop
set output 'clusters_em_4_animation.gif'
set size square
do for [i=0:24] {
set palette model RGB defined (0 "red",1 "blue", 2 "green", 3 "yellow", 4 "orange", 5 "black", 6 "violet")
plot 'clusters_em_4.txt' using 1:2:3 notitle with points pt 2 palette, 'clusters_em_4_eighen.txt' using 1:2:3:4 with vectors filled head lw 3, 'clusters_em_4_animation_data.txt' every::i*4::(i+1)*4-1 using 1:2:3:4:5 with ellipses}
