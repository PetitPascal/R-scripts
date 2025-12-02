#------------------------------------------------------------------
#### Configurations  ####
#------------------------------------------------------------------

## Installing packages if needed
pack_needed<-c("tidyverse","ggridges")

for (i in 1:length(pack_needed)){
  if(pack_needed[i]%in%.packages(all.available=TRUE)){
    #ok
  }else{
    install.packages(pack_needed[i],repos = "http://cran.us.r-project.org")
  }
}

## Loading packages
library(tidyverse)
library(ggridges)

#------------------------------------------------------------------
#### Creating a simulated dataset  ####
#------------------------------------------------------------------

# Number of time periods or factors to consider (e.g., 5 in this example)
nb_periode<-5

# Creating a vector with the geometric mean of a continuous variable for each time period/factor
GeoMean<-c(0.34,0.23,0.17,0.10,0.06)

#  Creating a vector with the geometric standard deviation of a continuous variable for each time period/factor
GeoSD<-c(5.78,7.48,5.21,4.08,3.38)

# Creating a vector with the time period/factor names
Periode<-c("1987 - 1995","1996 - 2001", "2002 - 2007", "2008 - 2012", "2013 - 2019")

# Creating an empty table to fill with the upcoming for loops
tab_data<-c()

# Creating a table containing the continuous variable distributions for each time period/factor
for(i in 1:nb_periode){
  set.seed(123)                             # random seed for reproducibility
  distri<-rlnorm(as.numeric(100),           # number of values to generate
                 meanlog=log(GeoMean[i]),   # mean of the distribution on the log scale for each time period/factor
                 sdlog=log(GeoSD[i]))       # standard deviation of the distribution on the log scale for each time period/factor
  
  tab_data<-rbind(tab_data,data.frame(Valeur=as.numeric(distri),Periode=rep(Periode[i],100)))
  
}

# Transforming tab_data into a tibble
tab_data<-as_tibble(tab_data)

# Adding a new variable 
tab_data$Profession<-rep(c("Job1","Job2","Job3","Job4"),nrow(tab_data)/4)

#------------------------------------------------------------------
#### Graphical representation  ####
#------------------------------------------------------------------

# Graph 1 - density plots on the same window/figure
tab_data %>% ggplot()+ geom_density(aes(x=Valeur,group=factor(Periode),
                                        color=factor(Periode),fill=factor(Periode)),alpha=0.1) +
             theme_bw()+
             scale_x_log10("Concentration")+
             theme(axis.text = element_text(size=12,color="black"),
                   axis.title=element_text(size=14,face="bold",color="black"),
                   plot.title = element_text(size = 16, face = "bold"),
                   axis.line = element_line(color = "black",size = 0.1, linetype = "solid"))

# Graph 2 - density plots on the several windows
tab_data %>% ggplot()+ geom_density(aes(x=Valeur,group=factor(Periode),
                                        color=factor(Periode),fill=factor(Periode)),alpha=0.1,size=1) +
                      theme_bw()+
                      scale_x_log10("Concentration")+
                      theme(strip.text = element_text(size = 12, colour = "black", angle = 0, face = "bold"),
                            strip.background = element_rect(fill=alpha('#009E73', 0.2), colour="black", size=1),
                            axis.text = element_text(size=12,color="black",face="bold"),
                            axis.title=element_text(size=12,face="bold",color="black"),
                            legend.text = element_text(size = 12, face = "bold"),
                            legend.title=element_blank(),
                            legend.position = "none",
                            axis.line = element_line(color = "black",size = 0.1, linetype = "solid"))+
                      facet_wrap(.~Periode)

# Graph 3 - boxplot
tab_data %>% ggplot()+ geom_boxplot(aes(x=factor(Periode),y=Valeur)) + 
                       theme_bw()+
                       scale_x_discrete("Time period")+
                       scale_y_log10("Concentration")+
                        theme(strip.text = element_text(size = 12, colour = "black", angle = 0, face = "bold"),
                              strip.background = element_rect(fill=alpha('#009E73', 0.2), colour="black", size=1),
                              axis.text = element_text(size=12,color=grey(0.3),face="bold"),
                              axis.title=element_text(size=12,face="bold",color="black"),
                              legend.text = element_text(size = 12, face = "bold"),
                              legend.title=element_blank(),
                              legend.position = "none",
                              axis.line = element_line(color = "black",size = 0.1, linetype = "solid"))

# Graph 4 - violin plot
tab_data %>% ggplot(aes(x=factor(Periode),y=Valeur))+ 
                        geom_violin(alpha=0.2,size=1) + 
                        geom_boxplot(width=0.1,alpha=0.4)+
                        theme_light()+
                        scale_x_discrete("Period")+
                        scale_y_log10("Concentration")+
                        theme(strip.text = element_text(size = 12, colour = "black", angle = 0, face = "bold"),
                              strip.background = element_rect(fill=alpha('#009E73', 0.2), colour="black", size=1),
                              axis.text = element_text(size=12,color=grey(0.3),face="bold"),
                              axis.title=element_text(size=12,face="bold",color="black"),
                              legend.text = element_text(size = 12, face = "bold"),
                              legend.title=element_blank(),
                              legend.position = "none",
                              axis.line = element_line(color = "black",size = 0.1, linetype = "solid"))

# Graph 5 - simple ridgeline plot
tab_data %>% ggplot()+ 
  geom_density_ridges(data=tab_data,aes(y=Profession,x=Valeur)) + 
  scale_x_log10("Concentration")+
  theme_ridges()+
  theme(strip.text = element_text(size = 12, colour = "black", angle = 0, face = "bold"),
        strip.background = element_rect(fill=alpha('#009E73', 0.2), colour="black", size=1),
        axis.text = element_text(size=12,color=grey(0.3),face="bold"),
        axis.title=element_text(size=12,face="bold",color="black"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title=element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black",size = 0.1, linetype = "solid")) +
  facet_grid(Profession~Periode,drop = TRUE,scales='free')

# Graph 6 - Ridgeline plot with probability
ggplot(data=tab_data,aes(y=Profession,x=Valeur,fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.025, 0.975)
  ) +
  scale_fill_manual(
    name = "Probability", values = c("#00a064", grey(0.9), "#FF0000A0"),
    labels = c("]0 - 0.025]", "]0.025, 0.975]", "]0.975, 1]")
  )+
  scale_x_log10("Concentration")+
  theme_ridges()+
  theme(strip.text = element_text(size = 12, colour = "black", angle = 0, face = "bold"),
        strip.background = element_rect(fill=alpha('#009E73', 0.2), colour="black", size=1),
        axis.text = element_text(size=12,color=grey(0.3),face="bold"),
        axis.title=element_text(size=12,face="bold",color="black"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=12,face="bold",color="black"),
        axis.line = element_line(color = "black",size = 0.1, linetype = "solid")) +
  facet_grid(Profession~Periode,drop = TRUE,scales='free')
