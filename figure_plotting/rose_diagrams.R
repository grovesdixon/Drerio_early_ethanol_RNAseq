#rose_diagrams.R

rm(list=ls())
source('figure_plotting/rose_diagram_functions.R')

# UPLOAD AND FORMAT -------------------------------------------------------

#read in the data
dat =read_excel('figure_plotting/final_vangl2_protrusions.xlsx') %>% 
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


# SUMMARIZE DATASET -------------------------------------------------------

nrow(dat)  #177 total protrusions among the embryos
length(unique(dat$id)) # N = 20 total embroys
length(unique(dat$genotype)) # 2 genotypes
length(unique(dat$treatment)) # 2 treatments
dat %>% 
  group_by(genotype, treatment) %>% 
  summarize(N=n())            #total protrusions by group
dat %>% 
  group_by(genotype, treatment) %>% 
  summarize(N=length(unique(id)))    #total embroys by group



# SELECT BOUNDARIES FOR ROSE BARS -----------------------------------------

coords = make_bin_coords(length = 15)
dat2 = assign_bins(dat, coords)

#make a similarly formatted dataframe with bin as the exact angle for statistics
sdat = dat2 %>% 
  dplyr::select(-bin, -binLeft, -binRight) %>% 
  dplyr::rename(bin=angle) %>% 
  add_leading_trailing()

# PLOT FOR EACH EMBRYO -----------------------------------------------------

#set 90 degrees straight up as in the photos
offset = NISTdegTOradian(90) + NISTdegTOradian(360/nrow(coords)*.5)

#get binned counts per embryo
mdat = get_embryo_means(dat2)

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
all_emb_plts = plot_grid(plotlist = pltList)
# plot_grid(label, all_emb_plts, nrow=2, rel_heights = c(1,15))

# PLOT SUMMARY PROJECTIONS ------------------------------------------------

#get means for groups
gdat = mdat %>% 
  group_by(genotype, treatment, bin) %>% 
  summarize(mean_sum = mean(totLength),
            mean_count = mean(N)) %>% 
  unite('group',
        genotype,
        treatment,
        sep='-') %>% 
  add_leading_trailing()


#plot for small bins
border = max(gdat$mean_count)
groups = as.character(unique(gdat$group))
g_pltList=list()
for (g in groups){
  gsub = gdat %>% 
    filter(group==g)
  plt=plot_rose(df=gsub,
                Ycol='mean_count',
                offset=offset,
                title=g,
                border=border)
  g_pltList[[g]]=plt + labs(subtitle='')
}
legend = cowplot::get_legend(g_pltList[[1]] + theme(legend.position = 'bottom',
                                                    legend.justification = 'center'))
label = ggdraw() + draw_label('Mean projection count')
little_pies = plot_grid(plotlist = g_pltList)



# REPEAT FOR BIG BINS -----------------------------------------------------

coords = make_bin_coords(length = 90)
offset = NISTdegTOradian(90) + NISTdegTOradian(360/nrow(coords)*.5)
dat2 = assign_bins(dat, coords)
mdat = get_embryo_means(dat2)

#get means for groups
gdat = mdat %>% 
  group_by(genotype, treatment, bin) %>% 
  summarize(mean_sum = mean(totLength),
            mean_count = mean(N)) %>% 
  unite('group',
        genotype,
        treatment,
        sep='-') %>% 
  add_leading_trailing()

#plot 
border = max(gdat$mean_count)
groups = as.character(unique(gdat$group))
g_pltList=list()
for (g in groups){
  gsub = gdat %>% 
    filter(group==g)
  plt=plot_rose(df=gsub,
                Ycol='mean_count',
                offset=offset,
                title=g,
                border=border)
  g_pltList[[g]]=plt + labs(subtitle='')
}
legend = cowplot::get_legend(g_pltList[[1]] + theme(legend.position = 'bottom',
                                                    legend.justification = 'center'))
label = ggdraw() + draw_label('Mean projection count')
big_pies = plot_grid(plotlist = g_pltList)



# STATS FOR MEAN COUNT ----------------------------------------------------

#get the number of protrusions for each edge type for each embryo
sdat2 = sdat %>% 
  group_by(id, genotype, treatment, edge) %>% 
  summarize(N=n()) %>% 
  ungroup() %>% 
  mutate(genotype = factor(genotype, levels=c('wthet', 'mutant')))


#compare leading, trailing, and none between the treatment groups
#boxplots
sdat2 %>% 
  ggplot(aes(x=treatment, y=N, fill=genotype)) +
  geom_boxplot() +
  labs(y='projection count') +
  facet_grid(~edge)


#ANOVA FOR MEAN NON-LEADING/TRAILING PROTRUSIONS
#2-way anova for the mean number of protrusions that are not leading or trailing (edge_type=='none)

#subset 
edge_type = 'anterior-posterior'
esub = sdat2 %>% 
  filter(edge==edge_type) 

#anova
res_aov = aov(N ~ genotype * treatment,
              data=esub)
summary(res_aov)
TukeyHSD(res_aov, which = "treatment")
TukeyHSD(res_aov, which = "genotype:treatment") #double-check these against the Tukey letters from multcomp below

#view grouping summary
esub %>% 
  group_by(genotype, treatment) %>% 
  summarize(N=n())

#get Tukey letters
require(multcomp)
esub_grp = esub %>% 
  unite('group', genotype, treatment) %>% 
  mutate(group = factor(group))
res_aovgrp <- aov(N ~ group,
                  data = esub_grp)
summary(res_aovgrp)
tm = glht(res_aovgrp, linfct = mcp(group = "Tukey"))
tlet = cld(tm)
my_let = tlet$mcletters$monospacedLetters
tukey_df = data.frame(group = names(my_let),
                      tukey = my_let)


#build barplot for the mean number of non-leading/trailing protrusions by treatment group
letter_add = 0.4
bars_2way = sdat2 %>% 
  filter(edge==edge_type) %>% 
  group_by(genotype, treatment) %>% 
  summarize(mn = mean(N),
            se = std.error(N)) %>% 
  ungroup() %>% 
  unite('group', genotype, treatment, sep='_', remove=FALSE) %>% 
  left_join(tukey_df, by = 'group') %>% 
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
bars_2way


# ASSEMBLE FINAL PLOT FOR SIZE --------------------------------------------
little_pies
big_pies
bars_2way


# repeat above with 3-way anova -------------------------------------------

#run a 3-way anova for interactions
table(sdat2$treatment, sdat2$genotype)
res_aov3 <- aov(N ~ edge * genotype * treatment,
                data = sdat2)
summary(res_aov3)
TukeyHSD(res_aov3, which = "edge:treatment")
TukeyHSD(res_aov3, which = "genotype:treatment")
TukeyHSD(res_aov3, which = "edge:genotype:treatment")


#WRITE OUT P-VALUES FOR REVIEWERS WHO STRUGGLE WITH ABCs
pval_df = TukeyHSD(res_aov3, which = "edge:genotype:treatment")$`edge:genotype:treatment` %>% 
  data.frame()
pval_df %>% 
  rownames_to_column('group_pair') %>% 
  mutate(group_pair = sub('wthet', '+/?', group_pair),
         group_pair = sub('mutant', '-/-', group_pair) ) %>% 
  write_csv('figure_plotting/fig7_Tukey_pval_df.csv')

#get Tukey letters
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
                      tukey = trimws(my_let))


#barplots of means 
letter_add = 0.4
bplt = sdat2 %>% 
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
  labs(y = 'filopodia count',
       fill='') +
  scale_fill_manual(values = grey.colors(2),
                    labels = c(bquote('+/?'), bquote("-/-"))) +
  geom_text(aes(x=treatment, y=mn+se+letter_add, label=tukey),
            position=position_dodge(.9)) +
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(fill=c("red", 'grey')),
        strip.text = element_text(colour = 'white')) +
  facet_grid(~edge)

#change the rectangle stip colors to match other figures
g <- ggplot_gtable(ggplot_build(bplt))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c('forestgreen', 'dodgerblue', 'firebrick')
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
#final plot
ggdraw(g)



#match with letters?
tdf = TukeyHSD(res_aov3, which = "edge:genotype:treatment")$`edge:genotype:treatment`
ps = tdf[,4]
b = tdf[ps < 0.1,]
b
nrow(b) #matches the 11 'a' bars in figure

