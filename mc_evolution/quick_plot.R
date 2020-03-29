library(ggplot2)
data = read.csv('germ_alignment.csv')
ggplot(data, aes(update, mean_germ_alignment)) + 
    geom_line();
