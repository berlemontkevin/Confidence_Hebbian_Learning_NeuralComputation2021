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


temptest <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\01-temptest.csv')



p <- ggplot(data=subset(temptest,Alpha == 0.25 & Type_TC == "Uniform" & Stim == 0.45),aes(x=Percentage,y=Delta_Perf))
p <- p + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black")
p <- p + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N)))
p

p <- ggplot(data=subset(temptest,Alpha == 0.25 & Type_TC == "Optimized" & Stim == 0.45),aes(x=Percentage,y=Delta_Perf))
p <- p + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
p <- p + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N)))
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black")
p


p <- ggplot(data=subset(temptest,Alpha == 0.25 & Type_TC == "Uniform" & N == 25),aes(x=Percentage,y=Delta_Perf))
p <- p + geom_point(aes(colour=as.factor(Stim)),size=4) + theme_Publication() + scale_color_viridis(option="A",discrete=TRUE)
p <- p + geom_line(alpha=0.5,size=2,aes(colour=as.factor(Stim)))
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black")
p


p <- ggplot(data=subset(temptest,Alpha == 0.25 & Type_TC == "Optimized" & N == 25),aes(x=Percentage,y=Delta_Perf))
p <- p + geom_point(aes(colour=as.factor(Stim)),size=4) + geom_line(alpha=0.5,size=2,aes(colour=as.factor(Stim)))+ theme_Publication() + scale_color_viridis(option="A",discrete=TRUE)
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black")
p


%%%%%%%%
  % % DM network with 15Hz as decision th.
%%%%%%%

temptest <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\01-temptest2.csv')



p <- ggplot(data=subset(temptest,Alpha == 0.25 & Type_TC == "Uniform" & Stim == 0.45),aes(x=Percentage,y=Delta_Perf))
p <- p + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black")
p <- p + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N)))
p

p <- ggplot(data=subset(temptest,Alpha == 0.25 & Type_TC == "Optimized" & Stim == 0.45),aes(x=Percentage,y=Delta_Perf))
p <- p + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
p <- p + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N)))
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black")
p


p <- ggplot(data=subset(temptest,Alpha == 0.25 & Type_TC == "Uniform" & N == 25),aes(x=Percentage,y=Delta_Perf))
p <- p + geom_point(aes(colour=as.factor(Stim)),size=4) + theme_Publication() + scale_color_viridis(option="A",discrete=TRUE)
p <- p + geom_line(alpha=0.5,size=2,aes(colour=as.factor(Stim)))
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black")
p


p <- ggplot(data=subset(temptest,Alpha == 0.25 & Type_TC == "Optimized" & N == 25),aes(x=Percentage,y=Delta_Perf))
p <- p + geom_point(aes(colour=as.factor(Stim)),size=4) + geom_line(alpha=0.5,size=2,aes(colour=as.factor(Stim)))+ theme_Publication() + scale_color_viridis(option="A",discrete=TRUE)
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black")
p

%% Creation heatmap
temptest <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\01-fulldata-20Hz.csv')
p <- ggplot(data = subset(temptest,Alpha==0.20 & Type_TC == "Optimized" & Theta ==20.0 ), aes(x=Stim,y=N,fill=Delta_Perf))
p <- p + geom_tile() 
p <- p + theme_grey(base_size = 10) + 
  ggtitle("Heatmap (ggplot)") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p <- p +  scale_fill_distiller(palette = "RdBu",limits=c(0.0,1.0))
p


t = subset(temptest,Alpha==0.2 & Type_TC == "Uniform" & Theta ==20.0 )
t2 = subset(temptest,Alpha==0.2 & Type_TC == "Optimized" & Theta ==20.0 )
q=t$Delta_Perf-t2$Delta_Perf
qdf = t
qdf$Delta_Perf = q

p <- ggplot(data = qdf, aes(x=Stim,y=N,fill=Delta_Perf))
p <- p + geom_tile() 
p <- p + theme_grey(base_size = 10) + theme(axis.ticks = element_blank(),panel.background = element_blank()) 
p <- p +  scale_fill_distiller(palette = "RdBu",limits=c(-0.1,0.1),oob=squish)
p


# Section 15Hz threshold  ---------------------------------
data15 <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\01-fulldata-15Hz.csv')


## Section 15Hz threshold  =================================
alpha = 0.22
theta = 11.0

tempUni = subset(data15,Alpha==alpha & Type_TC == "Uniform" & Theta ==15.0 )
tempOpti = subset(data15,Alpha==alpha & Type_TC == "Optimized" & Theta ==15.0 )

Delta_OptiUni = tempOpti
dummy=tempOpti$Delta_Perf-tempUni$Delta_Perf
Delta_OptiUni$Delta_Perf = dummy

p <- ggplot(data = Delta_OptiUni, aes(x=Stim,y=N,fill=Delta_Perf))
p <- p + geom_tile() 
p <- p + theme_grey(base_size = 10) + theme(axis.ticks = element_blank(),panel.background = element_blank()) 
p <- p +  scale_fill_distiller(palette = "RdBu",limits=c(-0.1,0.1),oob=squish)
p


tempUni = subset(data15,Alpha==alpha & Type_TC == "Uniform" & Theta ==theta )
tempOpti = subset(data15,Alpha==alpha & Type_TC == "Optimized" & Theta ==theta )

Delta_OptiUni = tempOpti
dummy=tempOpti$Delta_Perf-tempUni$Delta_Perf
Delta_OptiUni$Delta_Perf = dummy

p2 <- ggplot(data = Delta_OptiUni, aes(x=Stim,y=N,fill=Delta_Perf))
p2 <- p2 + geom_tile() 
p2 <- p2 + theme_grey(base_size = 10) + theme(axis.ticks = element_blank(),panel.background = element_blank()) 
p2 <- p2 +  scale_fill_distiller(palette = "RdBu",limits=c(-0.1,0.1),oob=squish)
p2

library(gridExtra)
grid.arrange(p,p2)




