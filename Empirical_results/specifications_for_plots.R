##------------------------------------------------------------------------------
## Plot specifications
## Citations: 
## - taken from Yun and Barber[2023]'s codes from 
##   https://github.com/yjyun97/cluster_inf_unknown_var/blob/main/functions.R

TITLE_SIZE = 45
AXIS_TITLE_SIZE = 35
AXIS_TEXT_SIZE = 25
LEGEND_TITLE_SIZE = 40
LEGENT_TEXT_SIZE = 40

plot_details = theme(plot.title = element_text(hjust = 0.5, 
                                               size = TITLE_SIZE),
                     axis.title = element_text(size = AXIS_TITLE_SIZE, 
                                               color = "black"),
                     axis.text = element_text(size = AXIS_TEXT_SIZE, 
                                              color = "black"),
                     axis.line = element_line(color = "black"), 
                     panel.background = element_rect(fill = 'white', 
                                                     color = 'white'),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     legend.position = "right", 
                     legend.justification = "left",
                     legend.direction = "vertical",
                     legend.background = element_rect(fill = "transparent"), 
                     legend.title = element_text(size = LEGEND_TITLE_SIZE),
                     legend.text = element_text(size = LEGENT_TEXT_SIZE), 
                     legend.key = element_rect(color = NA, fill = NA))

##------------------------------------------------------------------------------