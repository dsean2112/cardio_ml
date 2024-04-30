plot_func <- function(x = 1:5000/500,y,color) {

color[color == 1] <- 'p'
color[color == 2] <- 'N'
color[color == 3] <- 't'


frame <- data.frame(Time = x, Signal = y)
plot <-
  ggplot(frame, aes(Time, Signal, color = color)) + 
  geom_path(linewidth = 0.25, aes(group = 1)) + geom_point() + 
  scale_x_continuous(breaks = seq(0, 10, 1)) + theme(legend.position="none") + theme(legend.title = element_blank())

return(plot)
}