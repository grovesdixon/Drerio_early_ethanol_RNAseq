#rose_diagrams.R

library(tidyverse)
library(readxl)
library(cowplot)
library(plotrix)
library(NISTunits)
theme_set(theme_cowplot())



# UPLOAD AND FORMAT -------------------------------------------------------

rm(list=ls())
#read in the data
dat =read_excel('figure_plotting/myvangl2_protrusionsV3.xlsx') %>% 
  dplyr::rename(embryo=`final embryo`) %>% 
  dplyr::select(genotype,
                treatment,
                angle,
                length,
                embryo) %>% 
  unite('id',
        genotype,
        treatment,
        embryo,
        remove=FALSE)



# SELECT BOUNDARIES FOR ROSE BARS -----------------------------------------

#set up bin coordinates
length = 15
mids = seq(0,(361-length),by=length)
lefts = mids - length/2
negLeft = lefts<0
lefts[negLeft]<-lefts[negLeft]+360
rights = mids + length/2
coords = data.frame(lefts,
                    rights,
                    mids,
                    span=length)
dat$bin = 'not assigned'
dat$binLeft = 'not assigned'
dat$binRight = 'not assigned'
for (i in 1:length(lefts)){
  l=lefts[i]
  r=rights[i]
  binDeg = mids[i]
  #deal with looping around
  if (i==1){
    above = dat$angle >= l & dat$angle <= 360
    below = dat$angle >= 0 & dat$angle < r
    inBin = above | below
  } else {
    inBin = dat$angle >= l & dat$angle < r
  }
  dat$bin[inBin] <- binDeg
  dat$binLeft[inBin] <- l
  dat$binRight[inBin] <- r
}
dat2 = dat %>% 
  mutate(id=factor(id),
         bin=factor(bin, levels=mids))



# PLOTTING FUNCTIONS -------------------------------------------------------

offset = NISTdegTOradian(90) + NISTdegTOradian(360/nrow(coords)*.5)


#function to add leading and trailing edge information
add_leading_trailing = function(df){
  df %>% 
    mutate(bin_num = as.numeric(as.character(bin)),
           leading = (bin_num >= 135 & bin_num <= 225),
           trailing = (bin_num >= 0 & bin_num <= 45) | (bin_num >= 315 & bin_num <= 360),
           edge = 'none',
           edge = if_else(leading,
                          'leading',
                          edge),
           edge = if_else(trailing,
                          'trailing',
                          edge),
           edge = factor(edge,
                         levels=c('leading', 'trailing', 'none')))
}

#make a dataframe to do stats on (avoid confusion with binning by using exact angles)
sdat = dat2 %>% 
  dplyr::select(-bin, -binLeft, -binRight) %>% 
  dplyr::rename(bin=angle) %>% 
  add_leading_trailing()

#function to build the rose diagram
plot_rose = function(df, Ycol='modMn', offset, title, border){
  df$border=border
  df %>% 
    ggplot() + 
    geom_bar(aes_string(x='bin', y=Ycol, fill='edge', color='edge'),
             width = 1, stat='identity') +
    scale_fill_manual(values = c('forestgreen', 'dodgerblue', 'firebrick')) +
    scale_color_manual(values = c('forestgreen', 'dodgerblue', 'firebrick')) +
    # geom_vline(xintercept = c(0.5,1,4)) +
    #second set of bars who's borders act as grid
    geom_bar(aes(x=bin,
                 y = border),
                 width = 1,
             position = "stack",
             stat = "identity",
             fill = NA,
             colour = "grey",
             size=0.1) +
    coord_polar(start=-offset ) + 
    labs(x = "", y = "",subtitle=title) + 
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks = element_blank())
}

# PLOT FOR EACH EMBRYO -----------------------------------------------------

#get means for embryos
mdat = dat2 %>% 
  group_by(id, genotype, treatment, bin, .drop=FALSE) %>% 
  summarize(N=n(),
            totLength = sum(length, na.rm=TRUE),
            med = median(length),
            mn = mean(length),
            se = std.error(length),
            mockSe = mn/5) %>% 
  arrange(as.numeric(bin)) %>% 
  mutate(modMn = if_else(is.na(mn),
                         0,
                         mn),
         modTotLength = if_else(is.na(totLength),
                                0,
                                totLength)) %>% 
  add_leading_trailing() %>% 
  ungroup() %>% 
  mutate(genotype = factor(genotype, levels=c('wthet', 'mutant')),
         treatment = factor(treatment))
  

#PLOT SUMMED LENGTHS FOR EACH EMBRYO
column_to_plot = 'totLength'
border = max(mdat[,column_to_plot])

embs = as.character(unique(mdat$id))
pltList=list()
for (e in embs){
  esub = mdat %>% 
    filter(id==e)
  plt=plot_rose(df=esub,
                Ycol=column_to_plot,
                offset=offset,
                title=e,
                border=border)
  pltList[[e]]=plt
}
label = ggdraw() + draw_label('Summed length')
plts = plot_grid(plotlist = pltList)
plot_grid(label, plts, nrow=2, rel_heights = c(1,15))


#PLOT TOTAL PROJECTIONS
column_to_plot = 'N'
border = max(mdat[,column_to_plot])

embs = as.character(unique(mdat$id))
pltList=list()
for (e in embs){
  esub = mdat %>% 
    filter(id==e)
  plt=plot_rose(df=esub,
                Ycol=column_to_plot,
                offset=offset,
                title=e,
                border=border)
  pltList[[e]]=plt
}
label = ggdraw() + draw_label('N projections')
# plts = plot_grid(plotlist = pltList)
# plot_grid(label, plts, nrow=2, rel_heights = c(1,15))


# RUN STATS FOR MEAN SUMMED LENGTH -----------------------------------------------------

#here compare leading, trailing, and none between the treatment groups
#boxplots
sdat %>% 
  ggplot(aes(x=treatment, y=length, fill=genotype)) +
  geom_boxplot() +
  labs(y='summed projection lengths') +
  facet_grid(~edge) 


#try a 2-way anova within each type of extension
edge_type = 'none'

esub = sdat %>% 
  filter(edge==edge_type) 
res_aov = aov(length ~
                genotype * treatment,
              data=esub)
summary(res_aov)

#barplots of means for significant treatment effect
sdat %>% 
  filter(edge==edge_type) %>% 
  group_by(treatment) %>% 
  summarize(mn = mean(length),
            se = std.error(length)) %>% 
  ggplot(aes(x=treatment, y=mn)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mn-se, ymax=mn+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  labs(y='summed non-leading/trailing projection length')


#run a 3-way anova for interactions
table(sdat$treatment, sdat$genotype)

res_aov <- aov(length ~ 
                  genotype + 
                  treatment + 
                  edge + 
                  genotype:edge + 
                  treatment:genotype + 
                  genotype:treatment + 
                  genotype:treatment:edge, 
                data = sdat)
summary(res_aov)
TukeyHSD(res_aov, which = "genotype:treatment:edge")


#barplots of means 
sdat %>% 
  group_by(genotype, treatment, edge) %>% 
  summarize(mn = mean(length),
            se = std.error(length)) %>% 
  ggplot(aes(x=treatment, y=mn, fill=genotype)) +
  geom_bar(position='dodge', stat="identity") +
  geom_errorbar(aes(ymin=mn-se, ymax=mn+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  labs(y='mean summed projection length') +
  facet_grid(~edge)


# RUN STATS FOR MEAN COUNT ----------------------------------------------------

sdat2 = sdat %>% 
  group_by(id, genotype, treatment, edge) %>% 
  summarize(N=n()) %>% 
  ungroup() %>% 
  mutate(genotype = factor(genotype, levels=c('wthet', 'mutant')))


#here compare leading, trailing, and none between the treatment groups
#boxplots
sdat2 %>% 
  ggplot(aes(x=treatment, y=N, fill=genotype)) +
  geom_boxplot() +
  labs(y='projection count') +
  facet_grid(~edge)


#try a 2-way anova within each type of extension
edge_type = 'none'
esub = sdat2 %>% 
  filter(edge==edge_type) 
res_aov = aov(N ~ genotype * treatment,
              data=esub)
summary(res_aov)
TukeyHSD(res_aov, which = "treatment")
TukeyHSD(res_aov, which = "genotype:treatment")


#build barplot for this (NOTE HARD-CODED TUKEY LETTERS!)
letter_add = 0.4
sdat2 %>% 
  filter(edge==edge_type) %>% 
  group_by(genotype, treatment) %>% 
  summarize(mn = mean(N),
            se = std.error(N)) %>% 
  ungroup() %>% 
  mutate(tukey = c('ab', 'a', 'ab', 'b')) %>% 
  ggplot(aes(x=treatment, y=mn, fill=genotype)) +
  geom_bar(position='dodge', stat="identity") +
  geom_errorbar(aes(ymin=mn-se, ymax=mn+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  labs(y = 'projection count',
       fill='') +
  scale_fill_manual(values = grey.colors(2),
                    labels = c(bquote('+/?'), bquote("-/-"))) +
  geom_text(aes(x=treatment, y=mn+se+letter_add, label=tukey),
            position=position_dodge(.9))


#run a 3-way anova for interactions
table(sdat2$treatment, sdat2$genotype)
res_aov3 <- aov(N ~ edge * genotype * treatment,
               data = sdat2)
summary(res_aov3)
TukeyHSD(res_aov3, which = "edge:treatment")
TukeyHSD(res_aov3, which = "genotype:treatment")
TukeyHSD(res_aov3, which = "edge:genotype:treatment")

#get Tukey letters
require(multcomp)
sdat_grp = sdat2 %>% 
  unite('group', genotype, treatment, edge) %>% 
  mutate(group = factor(group))
res_aovgrp <- aov(N ~ group,
                data = sdat_grp)
summary(res_aovgrp)
tm = glht(res_aovgrp, linfct = mcp(group = "Tukey"))
tlet = cld(tm)
my_let = tlet$mcletters$monospacedLetters
tukey_df = data.frame(group = names(my_let),
                      tukey = my_let)



#barplots of means 
letter_add = 0.4
sdat2 %>% 
  group_by(genotype, treatment, edge) %>% 
  summarize(mn = mean(N),
            se = std.error(N),
            N=n()) %>% 
  unite('group', genotype, treatment, edge, sep='_', remove=FALSE) %>% 
  left_join(tukey_df, by = 'group') %>% 
  ggplot(aes(x=treatment, y=mn, fill=genotype)) +
  geom_bar(position='dodge', stat="identity") +
  geom_errorbar(aes(ymin=mn-se, ymax=mn+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  labs(y = 'projection count',
       fill='') +
  scale_fill_manual(values = grey.colors(2),
                    labels = c(bquote(italic(vangl2)^'+/?'), bquote(italic(vangl2)^"m209/m209"))) +
  geom_text(aes(x=treatment, y=mn+se+letter_add, label=tukey),
            position=position_dodge(.9)) +
  facet_grid(~edge)

#doublecheck sample sizes
sdat_grp %>% 
  group_by(group) %>% 
  summarize(N=n())
  
  


# STATS FOR PROPORTION NONE -----------------------------------------------

#get proportions
pdat = sdat2 %>% 
  mutate(none_edge = if_else(edge=='none',
                             'anterior/posterior',
                             'leading/trailing')) %>% 
  group_by(id, genotype, treatment, none_edge) %>% 
  summarize(totN = sum(N)) %>% 
  pivot_wider(id_cols = c('id', 'genotype', 'treatment'),
              names_from = none_edge,
              values_from = totN) %>% 
  mutate(prop_none = `anterior/posterior` / (`leading/trailing` + `anterior/posterior`))


#stats
res_aov = aov(prop_none ~ genotype * treatment,
              data=pdat)
summary(res_aov)
TukeyHSD(res_aov, which = "treatment")
TukeyHSD(res_aov, which = "genotype:treatment")


#plot
pdat %>% 
  group_by(genotype, treatment) %>% 
  summarize(mn = mean(prop_none),
            se = std.error(prop_none)) %>% 
  ggplot(aes(x=treatment, y=mn, fill=genotype)) +
  geom_bar(position='dodge', stat="identity") +
  geom_errorbar(aes(ymin=mn-se, ymax=mn+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  labs(y = 'proportion anterior/posterior')
  


#------ BUILD THE PLOT
#get mean sums for groups
gdat = mdat %>% 
  group_by(genotype, treatment, bin) %>% 
  summarize(mean_sum = mean(totLength),
            mean_count = mean(N)) %>% 
  unite('group',
        genotype,
        treatment,
        sep='-') %>% 
  add_leading_trailing()


#build group level plots
#plot sums for each embryo
groups = as.character(unique(gdat$group))
g_pltList=list()
border = max(gdat$mean_sum)
for (g in groups){
  gsub = gdat %>% 
    filter(group==g)
  plt=plot_rose(df=gsub,
                Ycol='mean_sum',
                offset=offset,
                title=g,
                border=border)
  g_pltList[[g]]=plt
}
legend = cowplot::get_legend(g_pltList[[1]] + theme(legend.position = 'bottom',
                                                    legend.justification = 'center'))
label = ggdraw() + draw_label('Mean summed length')
plts = plot_grid(plotlist = g_pltList)
plot_grid(label, plts, legend, nrow=3, rel_heights = c(1,15, 2))



# PLOT MEAN COUNT PER GROUP -----------------------------------------------

border = max(gdat$mean_count)
for (g in groups){
  gsub = gdat %>% 
    filter(group==g)
  plt=plot_rose(df=gsub,
                Ycol='mean_count',
                offset=offset,
                title=g,
                border=border)
  g_pltList[[g]]=plt
}
legend = cowplot::get_legend(g_pltList[[1]] + theme(legend.position = 'bottom',
                                                    legend.justification = 'center'))
label = ggdraw() + draw_label('Mean projection count')
plts = plot_grid(plotlist = g_pltList)
plot_grid(label, plts, legend, nrow=3, rel_heights = c(1,15, 2))


# LEARN TO BUILD PLOT -----------------------------------------------------
#from here: https://learnr.wordpress.com/2010/08/16/consultants-chart-in-ggplot2/
#set up data
set.seed(9876)
DF <- data.frame(variable = 1:10, value = sample(10,replace = TRUE))
DF

#plot 
DF %>% 
ggplot(aes(factor(variable), value, fill = factor(variable))) + 
  geom_bar(width = 1, stat='identity') +
  scale_y_continuous(breaks = 0:10) +
  coord_polar() + 
  labs(x = "", y = "") + 
  theme(legend.position = "none",
       axis.text.x = element_blank(),
       axis.text.y = element_blank(),
       axis.ticks = element_blank())

#plot with grids
library(plyr)
library(dplyr)
DF <- ddply(DF, .(variable), transform, border = rep(1,value))
ggplot(DF, aes(factor(variable))) + 
  geom_bar(width = 1,
           aes(y = value, fill = factor(variable)),
           stat='identity') + 
  #second set of bars who's borders act as grid
  geom_bar(aes(y = border, width = 1),
           position = "stack",
           stat = "identity",
           fill = NA, colour = "white") +
  scale_y_continuous(breaks = 0:10) + 
  coord_polar() +
  labs(x = "", y = "") + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())




