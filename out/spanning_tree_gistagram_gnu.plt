width=0.01
bin(x, s) = s*int(x/s) + width/2
set boxwidth width
plot 'spanning_tree_gistagram.txt' u (bin($1,width)):(1.0)
s f w boxes fs solid 0.5 title 'spanning tree distances'