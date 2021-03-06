set term gif \
    animate \
    optimize \
    delay 100 \
    size 300, 300 \
    background "#ffeedf" \
    crop \


do for [i=0:10] {
  set palette model RGB defined (0 "red",1 "blue", 2 "green", 3 "yellow", 4 "orange", 5 "black", 6 "violet")
    plot 'clusters_em_4.txt' using 1:2:3 notitle with points pt 2 palette, 'clusters_em_4_eighen.txt' using 1:2:3:4 with vectors filled head lw 3, 'clusters_em_4_em_animation_data.txt' every::i::(i+4) using 1:2:3:4:5 with ellipses
}

