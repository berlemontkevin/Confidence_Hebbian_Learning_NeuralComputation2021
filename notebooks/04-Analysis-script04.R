library(ggplot2)
library(ggthemes)
library(scales)
library(viridis)
theme_Publication <- function(base_size=14, base_family="sans") {
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}


temptest <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\04-dataUniform-20Hz.csv')


p <- ggplot(data=subset(temptest,type_TC == "Uniform" & box_conv == 0.1 & Seuil_confidence ==20),aes(x=stim,y=Perf))
p <- p + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black")
p <- p + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N)))
p


p <- ggplot(data = subset(temptest,type_TC == "Uniform" & box_conv == 1 & Seuil_confidence ==15 ), aes(x=stim,y=N,fill=Perf))
p <- p + geom_tile() 
p <- p + theme_grey(base_size = 10) + 
  ggtitle("Heatmap (ggplot)") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p <- p +  scale_fill_distiller(palette = "Bu",limits=c(0.4,1.0))
p

temptestOpt <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\04-dataOptimized-20Hz.csv')
p <- ggplot(data = subset(temptestOpt,type_TC == "Optimized" & box_conv == 0.1 & Seuil_confidence ==15 ), aes(x=stim,y=N,fill=Perf))
p <- p + geom_tile() 
p <- p + theme_grey(base_size = 10) + 
  ggtitle("Heatmap (ggplot)") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p <- p +  scale_fill_distiller(palette = "Bu",limits=c(0.4,1.0))
p

