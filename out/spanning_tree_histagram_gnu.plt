binwidth=0.1
bin(x,width)=width*floor(x/width) + width/2.0
set boxwidth binwidth
plot 'spanning_tree_histagram.txt' using (bin($1,binwidth)):(1.0) smooth freq with boxes