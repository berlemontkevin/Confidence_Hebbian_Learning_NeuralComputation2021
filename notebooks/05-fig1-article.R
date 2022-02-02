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

### Section : generation de la figure des tuning curves
N=20
TC_opti <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\05-TuningCurves.csv')

p <- ggplot(data=subset(TC_opti, N ==1 & type=="Optimized"),aes(x=stim,y=TC))
for (i in 1:N){
  p <- p + geom_line(data=subset(TC_opti,N==i & type=="Uniform"),size=1.5,alpha=1.0,
                      linetype="dashed",colour="firebrick",show.legend=FALSE)
}
p <- p  + theme_Publication()
p <- p + theme(axis.text.x = element_text(face="bold", 
                                          size=14),
               axis.text.y = element_text(face="bold", 
                                          size=14))
for (i in 1:N){
  p <- p + geom_line(data=subset(TC_opti,N==i & type=="Optimized"),size=1.5,colour="dodgerblue4",show.legend=FALSE)
  }

p
ggsave(file="C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\plots\\03-Tuning-curves-Information-theory\\fig1-TC.png", plot=p, width=10, height=8)


## Section génération figure Fisher information + Catégories
Delta_Fcode <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\05-FisherInformation.csv')

p <- ggplot(data=Delta_Fcode,aes(x=stim,y=Delta_FCode))
p <- p +  theme_Publication()
p <- p + geom_line(size=2,colour="dodgerblue4")
p <- p + geom_hline(yintercept=00, linetype="dashed", color = "black")
p <- p + theme(axis.text.x = element_text(face="bold", 
                                        size=14),
             axis.text.y = element_text(face="bold", 
                                        size=14))

Cat <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\05-Categories.csv')
for (i in 1:2){
p <- p + geom_line(data = subset(Cat, Categorie ==i), aes(x=stim,y=1000*prob_Cat),colour="black",size=1)
}
p
ggsave(file="C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\plots\\03-Tuning-curves-Information-theory\\fig1-Fcode.svg", plot=p, width=10, height=8)


### Figures Tuning curves individuelles

N=15
TC_opti <- read.csv(file = 'C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\data\\notebooks\\05-TuningCurves.csv')

p <- ggplot(data=subset(TC_opti, N ==1 & type=="Optimized"),aes(x=stim,y=TC))
p <- p  + theme_Publication()
p <- p + theme(axis.text.x = element_text(face="bold", 
                                          size=14),
               axis.text.y = element_text(face="bold", 
                                          size=14))
for (i in 1:N){
  p <- p + geom_line(data=subset(TC_opti,N==i & type=="Optimized"),size=1.5,colour="dodgerblue4",show.legend=FALSE)
}

p
ggsave(file="C:\\Documents\\1 - PhD Projects\\4-project_coding_categories\\plots\\03-Tuning-curves-Information-theory\\fig1-TC-O.svg", plot=p, width=10, height=8)

