schematic_plot_b <- function() {
  library(ggplot2)
  library(ggpubr)
  library(data.table)
  # viral load
  vl=splinefunH(c(0,5,20), c(0,1,0), c(0,0,0))
  vl=vl(0:(20/0.01)*0.01)
  vl=data.table(x=0:(20/0.01)*0.01, y=vl)
  vl.ci.lwr=splinefunH(c(1,6,17), c(10,35,10), c(0,0,0))
  vl.ci.lwr=vl.ci.lwr(1:(17/0.01)*0.01)
  vl.ci.lwr=data.table(x=1:(17/0.01)*0.01, y=vl.ci.lwr)
  vl.ci.upp=splinefunH(c(0,3,6,23), c(10,40,40,10), c(0,0,0,0))
  vl.ci.upp=vl.ci.upp(0:(20/0.01)*0.01)
  vl.ci.upp=data.table(x=0:(20/0.01)*0.01, y=vl.ci.upp)
  vl.ci = data.table(x=c(vl.ci.lwr$x, 20, rev(vl.ci.upp$x)), y=c(vl.ci.lwr$y, 10, rev(vl.ci.upp$y)))
  # ct trajectory
  ct = data.table(x=c(0,4.3,18.5), y=c(40,17.5,40))
  ct.ci.lwr = splinefunH(c(1,4.6,17), c(40,24,40), c(0,0,0))
  ct.ci.lwr=ct.ci.lwr(2:(17/0.01)*0.01)
  ct.ci.lwr=data.table(x=2:(17/0.01)*0.01, y=ct.ci.lwr)
  ct.ci.upp = splinefunH(c(-2,4.5,24), c(40,15,40), c(0,0,0))
  ct.ci.upp = ct.ci.upp(0:(20/0.01)*0.01)
  ct.ci.upp = data.table(x=0:(20/0.01)*0.01, y=ct.ci.upp)
  ct.ci = data.table(x=c(ct.ci.lwr$x, 20, rev(ct.ci.upp$x)), y=c(ct.ci.lwr$y, 40, rev(ct.ci.upp$y)))
  ct.obs = data.table(x=c(2,3,6,7,9,10,15), y=c(33,24,23.5,20,27,27,33))
  # reverse ct
  ct$y=50-ct$y
  ct.ci$y=50-ct.ci$y
  ct.obs$y=50-ct.obs$y
  # symptoms onset
  symp.obs=data.table(x=5, y=36)
  # colours
  CONFIG = list() 
  CONFIG$cols = c('#E64B35FF', '#4DBBD5FF', '#00A087FF', '#3C5488FF', '#F39B7FFF', 
                  '#8491B4FF', '#91D1C2FF', '#DC0000FF', '#7E6148FF', '#B09C85FF') # Nature
  lightup = function(c, alpha)
  {
    z=col2rgb(c)/255
    return(rgb(z[1],z[2],z[3],alpha))  # 0.125 for col3
  }
  CONFIG$colsLight1 = c();for(i in seq(CONFIG$cols))CONFIG$colsLight1[i] = lightup(CONFIG$cols[i], alpha = 0.25)
  CONFIG$colsLight2 = c();for(i in seq(CONFIG$cols))CONFIG$colsLight2[i] = lightup(CONFIG$cols[i], alpha = 0.1)
  # sample plot
  test.a=ggplot() +
    # geom_polygon(aes(x=c(2,4,12,15,4,-1), y=c(1,0.5,1,1,0.3,1)), fill=c(CONFIG$colsLight1[6])) +
    
    # geom_polygon(data=vl.ci,aes(x=x, y=y), fill=c(CONFIG$colsLight2[5])) +
    # geom_line(data=vl, aes(x=x, y=y*30+10, colour='Viral load'), size=0.25) +
    
    geom_polygon(data=ct.ci,aes(x=x, y=y), fill=c(CONFIG$colsLight2[6])) +
    # geom_line(data=ct, aes(x=x, y=y, colour='Ct trajectory'), size=0.25) +
    geom_line(data=ct, aes(x=x, y=y), colour=CONFIG$cols[4], size=0.25) +
    
    # scale_color_manual(values=c(CONFIG$cols[4], CONFIG$cols[1]),
    #                    guide = guide_legend(reverse = TRUE),
    #                    name = 'Lines') +
    geom_hline(yintercept=10,linetype=2, colour='gray60') +
    geom_segment(aes(x=0,xend=4.3,y=50-17.5,yend=50-17.5),linetype=2, colour='gray60') +
    geom_segment(aes(x=4.3,xend=4.3,y=10,yend=50-17.5),linetype=2, colour='gray60') +
    
    # geom_point(aes(x=0, y=0.2*30+10, shape='01infection'), colour=CONFIG$cols[9]) +
    # 
    # geom_line(aes(x=c(2,2,4,4,4,6,6),y=c(0.2,0.175,0.175,0.15,0.175,0.175,0.2)*30+10), size=0.1) +
    # annotate(geom="text", x=4, y=0.14*30+10, label='symptoms window', size=1.5) +
    # geom_line(aes(x=c(2,6), y=c(0.2,0.2)*30+10), colour=CONFIG$cols[5],size=0.25) +
    # geom_point(aes(x=c(2,4,6), y=0.2*30+10, shape='02symptoms'), colour=c(CONFIG$cols[5],CONFIG$cols[1],CONFIG$cols[5])) +
    # 
    # geom_line(aes(x=c(0,0,2,2,2,4,4),y=c(0.125,0.1,0.1,0.075,0.1,0.1,0.125)*30+10), size=0.1) +
    # annotate(geom="text", x=2, y=0.065*30+10, label='incubation period', size=1.5) +
    # 
  # 
  geom_point(data=ct.obs, aes(x=x, y=y, shape='03ct'), colour=CONFIG$cols[4]) +
    
    geom_point(data=symp.obs, aes(x=x, y=y, shape='04onset'), colour=CONFIG$cols[1]) +
    
    annotate(geom="text", x=2.5, y=39, label='incubation period', size=2) +
    geom_segment(aes(x = 0,
                     y = 38,
                     xend = 5,
                     yend = 38), size=0.25,
                 arrow = arrow(length = unit(0.1, "cm"),
                               type = "closed",
                               ends='last')) +
    # 
    # annotate(geom="text", x=2.5, y=24, label='time to peak Ct', size=1.5) +
    # geom_segment(aes(x = 0,
    #                  y = 23.5,
    #                  xend = 5,
    #                  yend = 23.5), size=0.25,
    #              arrow = arrow(length = unit(0.1, "cm"), 
    #                            type = "closed")) +
    annotate(geom="text", x=5, y=37, label='symptoms onset window', size=2) +
    geom_segment(aes(x = 3,
                     y = 36,
                     xend = 7,
                     yend = 36), size=0.25,
                 arrow = arrow(length = unit(0.1, "cm"),
                               type = "closed",
                               ends='both')) +
    # scale_shape_manual(values=c(17,16,4),
    #                    guide = guide_legend(override.aes = list(colour = c(CONFIG$cols[9], CONFIG$cols[1], CONFIG$cols[4]))),
    #                    name = 'Event',
    #                    labels =c('infection', 'symptoms', 'lab ct')) +
    scale_shape_manual(values=c(16, 17),
                       guide = guide_legend(override.aes = list(colour = c(CONFIG$cols[4], CONFIG$cols[1]))),
                       name = 'Observation',
                       labels =c('lab ct', 'symptoms onset')) +
    
    scale_x_continuous(expand = c(0, 0),
                       limits = c(0,20),
                       breaks = sort(c(seq(0,20,1), 4.3, 18.5)),
                       labels = c(expression(t[e]==0),'','','','',expression(t[p]),5,'','','','',10,'','','','',15,'','','',expression(t[LOD]),'',20)) + 
    
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(7,40),
      breaks = sort(c(seq(10,40,5), 32.5)),
      # labels = c('LOD=40','',30,'',20,expression(c[p]),'',10),
      "Ct"#,
      # sec.axis = sec_axis(~ . /30 - (1/3), 
      #                     name = 'Relative probability of transmission (%)', 
      #                     breaks = seq(0, 1, by = 0.125),
      #                     labels = c(0,'',25,'',50,'',75,'',100))
    ) +
    
    labs(x='Days since exposure')+
    # annotate(geom="text", x=0.5, y=41, label='Limit of detection',
    #          color="gray60", size=2) +
    
    # annotate(geom="text", x=13, y=38, label='Model parameters', size=2) +
    # annotate(geom="text", x=14, y=36.5, label=expression(t[e]=='unobserved time of exposure'), size=2) +
    # annotate(geom="text", x=13.9, y=35, label=expression(t[p]=='unobserved time to peak Ct'), size=2) +
    # annotate(geom="text", x=14.75, y=33.5, label=expression(t[LOD]=='unobserved time until LOD is reached'), size=2) +
    # annotate(geom="text", x=12.65, y=32, label=expression(c[p]=='peak Ct'), size=2) +
    
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_rect(colour = "black", fill=NA),
          axis.text.y   = element_text(size=8),
          axis.text.x   = element_text(size=8),
          axis.title.y  = element_text(size=8),
          axis.title.x  = element_text(size=8),
          # legend.position = c(0.65, 0.63),
          legend.position = "none",
          legend.background=element_blank(),
          legend.title = element_text(size=5),
          legend.text = element_text(size=5),
          legend.key = element_rect(color = NA, fill = NA),
          legend.key.size = unit(0.25, "cm"),
          legend.spacing.y = unit(0.01, "cm")) +
    
    labs(tag = 'A')
  
}
