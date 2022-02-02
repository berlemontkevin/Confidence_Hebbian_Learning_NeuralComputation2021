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


# Results from script 04

### Tuning cruves uniform
perf_Uni <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\06-Script04-TC-Uniform.csv')


p <- ggplot(data = subset(perf_Uni, box_conv == 0.1 & Seuil_confidence ==20 ), aes(x=Stim,y=N,fill=Perf))
p <- p + geom_tile() 
p <- p + theme_grey(base_size = 10) + 
  ggtitle("C = Uniform") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p <- p +  scale_fill_distiller(limits=c(0.4,1.0))
p

p4 <- ggplot(data = subset(perf_Uni, box_conv == 0.1 & Seuil_confidence ==14 ), aes(x=Stim,y=N,fill=Perf))
p4 <- p4 + geom_tile() 
p4 <- p4 + theme_grey(base_size = 10) + 
  ggtitle("C = Uniform") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p4 <- p4 +  scale_fill_distiller(limits=c(0.4,1.0))

# 100 ms ne change rien aux performances

perf_Opti <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\06-Script04-TC-Optimized.csv')


p2 <- ggplot(data = subset(perf_Opti, box_conv == 0.1 & Seuil_confidence ==20 ), aes(x=Stim,y=N,fill=Perf))
p2 <- p2 + geom_tile() 
p2 <- p2 + theme_grey(base_size = 10) + 
  ggtitle("C = Opti") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p2 <- p2 +  scale_fill_distiller(limits=c(0.4,1.0))

p5 <- ggplot(data = subset(perf_Opti, box_conv == 0.1 & Seuil_confidence ==14 ), aes(x=Stim,y=N,fill=Perf))
p5 <- p5 + geom_tile() 
p5 <- p5 + theme_grey(base_size = 10) + 
  ggtitle("C = Opti") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p5 <- p5 +  scale_fill_distiller(limits=c(0.4,1.0))


perf_H <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\06-Script04-TC-Hebbian.csv')
p3 <- ggplot(data = subset(perf_H, box_conv == 0.1 & Seuil_confidence ==20 ), aes(x=Stim,y=N,fill=Perf))
p3 <- p3 + geom_tile() 
p3 <- p3 + theme_grey(base_size = 10) + 
  ggtitle("C = Hebbian") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p3 <- p3 +  scale_fill_distiller(limits=c(0.4,1.0))

p6 <- ggplot(data = subset(perf_H, box_conv == 0.1 & Seuil_confidence ==14 ), aes(x=Stim,y=N,fill=Perf))
p6 <- p6 + geom_tile() 
p6 <- p6 + theme_grey(base_size = 10) + 
  ggtitle("C = Hebbian") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p6 <- p6 +  scale_fill_distiller(limits=c(0.4,1.0))

pf <- ggarrange(p,p2,p3,p4,p5,p6)
pf


perf_Delta_OU <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\06-Script04-TC-Delta-OU.csv')

p <- ggplot(data = subset(perf_Delta_OU, box_conv == 0.1 & Seuil_confidence ==20 ), aes(x=Stim,y=N,fill=Perf))
p <- p + geom_tile() 
p <- p + theme_grey(base_size = 10) + 
  ggtitle("C = Delta Opti - Uniform Sans confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p <- p +  scale_fill_distiller(palette="RdBu",limits=c(-0.2,0.2))
p

p <- ggplot(data = subset(perf_Delta_OU, box_conv == 0.1 & Seuil_confidence ==14 ), aes(x=Stim,y=N,fill=Perf))
p <- p + geom_tile() 
p <- p + theme_grey(base_size = 10) + 
  ggtitle("C = Delta Opti - Uniform") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p <- p +  scale_fill_distiller(palette="RdBu",limits=c(-0.2,0.2))
p



p <- ggplot(data=subset(perf_Uni, box_conv == 0.1 & Seuil_confidence ==20),aes(x=Stim,y=Perf))
p <- p + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black")
p <- p + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N))) + ggtitle("TC: Uniform") 
p

p <- ggplot(data=subset(perf_Opti, box_conv == 0.1 & Seuil_confidence ==14),aes(x=Stim,y=Perf))
p <- p + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black") + ggtitle("TC: Opti") 
p <- p + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N)))
p


p <- ggplot(data=subset(perf_Opti, box_conv == 0.1 & Stim == 0.45 & N %in% c(10,15,20,25) ),aes(x=percentageConfidence,y=Perf))
p <- p + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black") + ggtitle("TC: Opti Stim=0.45") 
p <- p + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N)))
p

p <- ggplot(data=subset(perf_Uni, box_conv == 0.1 & Stim == 0.45 & N %in% c(10,15,20,25)),aes(x=percentageConfidence,y=Perf))
p <- p + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black")
p <- p + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N))) + ggtitle("TC: Uniform, Stim=0.45") 
p

p <- ggplot(data=subset(perf_Uni, box_conv == 1 & Stim == 0.4),aes(x=percentageConfidence,y=Perf))
p <- p + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black") + ggtitle("TC: Uniform Stim = 0.4") 
p <- p + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N)))
p


p <- ggplot(data=subset(perf_Delta_OU, box_conv == 1 & Stim == 0.45),aes(x=percentageConfidence,y=Perf))
p <- p + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black")
p <- p + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N)))
p


#####
# Analyse donnees performances avec opti non lineaire
#####

perf_opti_U <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\06-Script03-NonLinear-Optimisation-TC-Uniform.csv')


p <- ggplot(data = subset(perf_opti_U, alpha == 0.2 ), aes(x=Stim,y=N,fill=Perf))
p <- p + geom_tile() 
p <- p + theme_grey(base_size = 10) + 
  ggtitle("C = Delta Opti - Uniform Sans confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p <- p +  scale_fill_distiller(palette="RdBu",limits=c(-0.4,1.0))
p



# Analyse script 03 ============================

### Tuning cruves uniform ###################
perf <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\06-Script03-Perf-20Hz.csv')


p <- ggplot(data = subset(perf, TC == "Uniform" & box_conv == 1.0 & Seuil_confidence ==15 & alpha == 0.25 ), aes(x=Stim,y=N,fill=Perf))
p <- p + geom_tile() 
p <- p + theme_grey(base_size = 10) + 
  ggtitle("C = Uniform") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p <- p +  scale_fill_distiller(limits=c(0.4,1.0))
p

### Optimized - Uniform ###################

perf_Delta_OU <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\06-Script03-Delta-Optimized-Uniform-20Hz.csv')

sc=14

#§ a ameliorer pour la fin
p <- ggplot(data = subset(perf_Delta_OU, box_conv == 1.0 & Seuil_confidence ==20 & alpha == 0.25 ), aes(x=Stim,y=N,fill=Perf))
p <- p + geom_tile() 
p <- p + theme_grey(base_size = 10) + 
  ggtitle("Delta OU : alpha = 0.25, Sans confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p <- p +  scale_fill_gradient2()
p


p2 <- ggplot(data = subset(perf_Delta_OU, box_conv == 1.0 & Seuil_confidence ==20 & alpha == 0.35 ), aes(x=Stim,y=N,fill=Perf))
p2 <- p2 + geom_tile() 
p2 <- p2 + theme_grey(base_size = 10) + 
  ggtitle("Delta OU : alpha = 0.35, Sans confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p2 <- p2 +  scale_fill_gradient2()
p2

p3 <- ggplot(data = subset(perf_Delta_OU, box_conv == 1.0 & Seuil_confidence ==sc & alpha == 0.25 ), aes(x=Stim,y=N,fill=Perf))
p3 <- p3 + geom_tile() 
p3 <- p3 + theme_grey(base_size = 10) + 
  ggtitle("Delta OU : alpha = 0.25, Avec confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p3 <- p3   +  scale_fill_gradient2()
p3


p4 <- ggplot(data = subset(perf_Delta_OU, box_conv == 1.0 & Seuil_confidence ==sc & alpha == 0.35 ), aes(x=Stim,y=N,fill=Perf))
p4 <- p4 + geom_tile() 
p4 <- p4 + theme_grey(base_size = 10) + 
  ggtitle("Delta OU : alpha = 0.35, Avec confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p4 <- p4  +  scale_fill_gradient2()
p4


p5 <- ggplot(data = subset(perf_Delta_OU, box_conv == 1.0 & Seuil_confidence ==20 & alpha == 0.15 ), aes(x=Stim,y=N,fill=Perf))
p5 <- p5 + geom_tile() 
p5 <- p5 + theme_grey(base_size = 10) + 
  ggtitle("Delta OU : alpha = 0.15, Sans confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p5 <- p5  +  scale_fill_gradient2()


p6 <- ggplot(data = subset(perf_Delta_OU, box_conv == 1.0 & Seuil_confidence ==sc & alpha == 0.15 ), aes(x=Stim,y=N,fill=Perf))
p6 <- p6 + geom_tile() 
p6 <- p6 + theme_grey(base_size = 10) + 
  ggtitle("Delta OU : alpha = 0.15, Avec confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p6 <- p6  +  scale_fill_gradient2()


library(ggpubr)
pf <- ggarrange(p,p2,p5,p3,p4,p6)
pf





### best - Uniform ###################

perf_Delta_BU <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\06-Script03-Delta-Best-Uniform-20Hz.csv')

#§ a ameliorer pour la fin
p <- ggplot(data = subset(perf_Delta_BU, box_conv == 1.0 & Seuil_confidence ==20 & alpha == 0.2 ), aes(x=Stim,y=N,fill=Perf))
p <- p + geom_tile() 
p <- p + theme_grey(base_size = 10) + 
  ggtitle("Delta Best - Uniform : alpha = 0.2, Sans confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p <- p +  scale_fill_gradient(low = "white", high = "red",limits=c(-0.03,0.4))
p


p2 <- ggplot(data = subset(perf_Delta_BU, box_conv == 1.0 & Seuil_confidence ==20 & alpha == 0.35 ), aes(x=Stim,y=N,fill=Perf))
p2 <- p2 + geom_tile() 
p2 <- p2 + theme_grey(base_size = 10) + 
  ggtitle("Delta Best - Uniform : alpha = 0.35, Sans confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p2 <- p2 +  scale_fill_gradient(low = "white", high = "red",limits=c(-0.03,0.4))
p2

p3 <- ggplot(data = subset(perf_Delta_BU, box_conv == 1.0 & Seuil_confidence ==15 & alpha == 0.2 ), aes(x=Stim,y=N,fill=Perf))
p3 <- p3 + geom_tile() 
p3 <- p3 + theme_grey(base_size = 10) + 
  ggtitle("Delta Best - Uniform : alpha = 0.2, Avec confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p3 <- p3 +  scale_fill_gradient(low = "white", high = "red",limits=c(-0.03,0.4))
p3


p4 <- ggplot(data = subset(perf_Delta_BU, box_conv == 1.0 & Seuil_confidence ==15 & alpha == 0.35 ), aes(x=Stim,y=N,fill=Perf))
p4 <- p4 + geom_tile() 
p4 <- p4 + theme_grey(base_size = 10) + 
  ggtitle("Delta Best - Uniform : alpha = 0.35, Avec confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p4 <- p4 +  scale_fill_gradient(low = "white", high = "red",limits=c(-0.03,0.4))
p4
  
library(ggpubr)
pf <- ggarrange(p,p2,p3,p4)
pf


### best - Optimized ###################

perf_Delta_BU <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\06-Script03-Delta-Best-Optimized-20Hz.csv')

#§ a ameliorer pour la fin
p <- ggplot(data = subset(perf_Delta_BU, box_conv == 1.0 & Seuil_confidence ==20 & alpha == 0.25 ), aes(x=Stim,y=N,fill=Perf))
p <- p + geom_tile() 
p <- p + theme_grey(base_size = 10) + 
  ggtitle("Delta Best - Optimized : alpha = 0.25, Sans confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p <- p +  scale_fill_gradient(low = "white", high = "red",limits=c(0.0,0.4))
p


p2 <- ggplot(data = subset(perf_Delta_BU, box_conv == 1.0 & Seuil_confidence ==20 & alpha == 0.35 ), aes(x=Stim,y=N,fill=Perf))
p2 <- p2 + geom_tile() 
p2 <- p2 + theme_grey(base_size = 10) + 
  ggtitle("Delta Best - Optimized : alpha = 0.35, Sans confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p2 <- p2 +  scale_fill_gradient(low = "white", high = "red",limits=c(-0.03,0.4))
p2

p3 <- ggplot(data = subset(perf_Delta_BU, box_conv == 1.0 & Seuil_confidence ==15 & alpha == 0.25 ), aes(x=Stim,y=N,fill=Perf))
p3 <- p3 + geom_tile() 
p3 <- p3 + theme_grey(base_size = 10) + 
  ggtitle("Delta Best - Optimized : alpha = 0.25, Avec confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p3 <- p3 +  scale_fill_gradient(low = "white", high = "red",limits=c(-0.03,0.4))
p3


p4 <- ggplot(data = subset(perf_Delta_BU, box_conv == 1.0 & Seuil_confidence ==15 & alpha == 0.35 ), aes(x=Stim,y=N,fill=Perf))
p4 <- p4 + geom_tile() 
p4 <- p4 + theme_grey(base_size = 10) + 
  ggtitle("Delta Best - Optimized : alpha = 0.35, Avec confiance") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 
p4 <- p4 +  scale_fill_gradient(low = "white", high = "red",limits=c(-0.03,0.4))
p4

library(ggpubr)
pf <- ggarrange(p,p2,p3,p4)
pf





### Analyse influence confiance ##########################
perf <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\06-Script03-Perf-20Hz.csv')

x = 0.46

p <- ggplot(data=subset(perf, alpha == 0.25 & TC == "Uniform" & box_conv == 1.0 & Stim == x & N %in% c(10,15,20,25,35)),aes(x=percentage_confidence,y=Perf))
p <- p + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
#p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black")
p <- p + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N))) + ggtitle("TC: Uniform, alpha = 0.25, Stim=0.45") 
p

p2 <- ggplot(data=subset(perf, alpha == 0.35 & TC == "Uniform" & box_conv == 1.0 & Stim == x & N %in% c(10,15,20,25,35)),aes(x=percentage_confidence,y=Perf))
p2 <- p2 + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
#p2 <- p2 + geom_hline(yintercept=00, linetype="dashed", color = "black")
p2 <- p2 + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N))) + ggtitle("TC: Uniform, alpha = 0.35, Stim=0.45") 
p

p3 <- ggplot(data=subset(perf, alpha == 0.25 & TC == "Optimized" & box_conv == 1.0 & Stim == x & N %in% c(10,15,20,25,35)),aes(x=percentage_confidence,y=Perf))
p3 <- p3 + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
#p3 <- p3 + geom_hline(yintercept=00, linetype="dashed", color = "black")
p3 <- p3 + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N))) + ggtitle("TC: Opt, alpha = 0.25, Stim=0.45") 
p

p4 <- ggplot(data=subset(perf, alpha == 0.35 & TC == "Optimized" & box_conv == 1.0 & Stim == x & N %in% c(10,15,20,25,35)),aes(x=percentage_confidence,y=Perf))
p4 <- p4 + geom_point(alpha=0.5,aes(colour=as.factor(N)),size=4) + theme_Publication() + scale_color_viridis(discrete=TRUE)
#p4 <- p4 + geom_hline(yintercept=00, linetype="dashed", color = "black")
p4 <- p4 + geom_line(alpha=0.5,size=2,aes(colour=as.factor(N))) + ggtitle("TC: Opt, alpha = 0.35,Stim=0.45") 
p

pf <- ggarrange(p,p2,p3,p4)
pf
