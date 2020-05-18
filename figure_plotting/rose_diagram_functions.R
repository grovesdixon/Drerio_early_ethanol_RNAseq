#rose_diagram_functions.R

#LIBS
library(tidyverse)
library(readxl)
library(cowplot)
library(plotrix)
library(NISTunits)
theme_set(theme_cowplot())

#FUNCTIONS

make_bin_coords = function(length){
  #set up bin coordinates
  mids = seq(0,(361-length),by=length) #midpoints for bins in degrees
  lefts = mids - length/2              #left boundaries for bins
  negLeft = lefts<0                    #modify for going around cirlce
  lefts[negLeft]<-lefts[negLeft]+360
  rights = mids + length/2             #right boundaries for bins
  coords = data.frame(lefts,
                      rights,
                      mids,
                      span=length)
  return(coords)
}

#function to assign projections to bins based on their angles
assign_bins = function(input_df, coords){
  #loop through bins and assign each protrusion to one of them
  dat = input_df
  dat$bin = 'not assigned'
  dat$binLeft = 'not assigned'
  dat$binRight = 'not assigned'
  for (i in 1:nrow(coords)){
    l=coords$lefts[i]
    r=coords$rights[i]
    binDeg = coords$mids[i]
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
  #switch to factors
  dat2 = dat %>% 
    mutate(id=factor(id),
           bin=factor(bin, levels=coords$mids))
  return(dat2)
}

#function to get embryo projection counts and means for the assigned bins
get_embryo_means = function(dat2){
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
  return(mdat)
}


#function to add leading and trailing edge information
add_leading_trailing = function(df){
  df %>% 
    mutate(bin_num = as.numeric(as.character(bin)),
           leading = (bin_num >= 135 & bin_num <= 225),
           trailing = (bin_num >= 0 & bin_num <= 45) | (bin_num >= 315 & bin_num <= 360),
           edge = 'anterior-posterior',
           edge = if_else(leading,
                          'leading',
                          edge),
           edge = if_else(trailing,
                          'trailing',
                          edge),
           edge = factor(edge,
                         levels=c('leading', 'trailing', 'anterior-posterior')))
}


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