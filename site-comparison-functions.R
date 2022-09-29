#Site comparison functions

ks_test_recycled_vs_not = function(data, positions = c(175, 175, 175, 5000)) {
  data = data %>%
    mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
           Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
           Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
    filter(Length <= 175 & Width <= 175 & Thickness <= 175) #size of calipers
  
  pd = data %>% gather(key = "measurement", value = "value", Weight, Length, Width, Thickness)
  pd$measurement = factor(pd$measurement, levels = c("Length", "Width", "Thickness", "Weight"))
  
  stat.test1 = pd %>% group_by(measurement) %>% t_test(value ~ recycled, ref.group = "FALSE") %>%
    add_xy_position(x = "recycled")
  stat.test1$ks_p = signif(ks.test((data %>% filter(recycled == T))$Weight, (data %>% filter(recycled != T))$Weight)$p)
  stat.test1[2,]$ks_p = signif(ks.test((data %>% filter(recycled == T))$Width, (data %>% filter(recycled != T))$Width)$p)
  stat.test1[3,]$ks_p = signif(ks.test((data %>% filter(recycled == T))$Thickness, (data %>% filter(recycled != T))$Thickness)$p)
  stat.test1[1,]$ks_p = signif(ks.test((data %>% filter(recycled == T))$Length, (data %>% filter(recycled != T))$Length)$p)
  
  stat.test1 = stat.test1 %>% add_significance(
    p.col = "ks_p",
    output.col = "ks_p_sigf",
    cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
    symbols = c("****", "***", "**", "*", "ns")
  )
  
  
  plot = ggplot(pd, aes(x = as.factor(recycled), y = value, color = as.factor(recycled))) +
    geom_boxplot() +
    facet_wrap(. ~ measurement, scales = "free_y") +
    labs(x = "Recycled?", y = "Value", color = "Recycled?") +
    theme(legend.title = element_blank()) +
    stat_pvalue_manual(stat.test1, label = "ks_p_sigf", tip.length = 0.01,
                       label.size = 3, y.position = positions) +
    theme_tufte() +
    scale_color_colorblind()
  return(plot)
}

dist_comp_recycled_vs_all = function(data) {
  data = data %>%
    mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
           Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
           Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
    filter(Length <= 200 & Width <= 200 & Thickness <= 200) #size of calipers
  
  pd = data %>% gather(key = "measurement", value = "value", Weight, Length, Width, Thickness)
  pd$measurement = factor(pd$measurement, levels = c("Length", "Width", "Thickness", "Weight"))
  
  # plot = ggplot(data = pd) +
  #   geom_density_ridges(mapping = aes(x = log(value + 1), 
  #                                     y = recycled, 
  #                                     fill = measurement), 
  #                       color = "white") + 
  #   geom_density(aes(x = log(value + 1), fill = measurement), color = "white") +
  #   scale_fill_manual(values = alpha(c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"), 0.7)) +
  #   facet_wrap(~measurement)
  
  p1 = ggplot(data = pd %>% filter(recycled == T)) + 
    geom_density(mapping = aes(x = log(value), 
                               color = recycled, fill = recycled), 
                 size = 0.5) + 
    geom_density(data = pd, mapping = aes(x = log(value)), 
                 color = "gray20", size = 0.5) + 
    facet_wrap(~measurement) +
    guides(color = "none", fill = "none") +
    scale_fill_manual(values = alpha("#56B4E9", 0.7)) + 
    scale_color_manual(values = alpha("#56B4E9", 1)) +
    labs(title = "Recycled") +
    ylim(0, 1) +
    xlim(-1, 8) +
    theme(axis.text=element_text(size=5),
          axis.title=element_text(size=6,face="bold"), 
          title = element_text(size=8), 
          strip.text.x = element_text(size = 6))
  
  p2 = ggplot(data = pd %>% filter(recycled == F)) + 
    geom_density(mapping = aes(x = log(value), 
                               color = recycled, fill = recycled), 
                 size = 0.5) + 
    geom_density(data = pd, mapping = aes(x = log(value)), 
                 color = "gray20", size = 0.5) + 
    facet_wrap(~measurement) +
    guides(color = "none", fill = "none") +
    scale_fill_manual(values = alpha("#E69F00", 0.7)) + 
    scale_color_manual(values = alpha("#E69F00", 1)) +
    labs(title = "Not recycled") +
    ylim(0, 1) +
    xlim(-1, 8) +
    theme(axis.text=element_text(size=5),
          axis.title=element_text(size=7,face="bold"), 
          title = element_text(size=8), 
          strip.text.x = element_text(size = 6))
  
  return(cowplot::plot_grid(p1, p2))
}

dist_comp_only_recycled = function(data) {
  data = data %>%
    mutate(Thickness = ifelse(is.na(Flake_thickness), Maximum_core_thickness, Flake_thickness),
           Length = ifelse(is.na(Flake_length), Maximum_core_length, Flake_length),
           Width = ifelse(is.na(Flake_width), Maximum_core_width, Flake_width)) %>%
    filter(Length <= 200 & Width <= 200 & Thickness <= 200) #size of calipers
  
  pd = data %>% gather(key = "measurement", value = "value", Weight, Length, Width, Thickness)
  pd$measurement = factor(pd$measurement, levels = c("Length", "Width", "Thickness", "Weight"))
  
  # plot = ggplot(data = pd) +
  #   geom_density_ridges(mapping = aes(x = log(value + 1), 
  #                                     y = recycled, 
  #                                     fill = measurement), 
  #                       color = "white") + 
  #   geom_density(aes(x = log(value + 1), fill = measurement), color = "white") +
  #   scale_fill_manual(values = alpha(c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"), 0.7)) +
  #   facet_wrap(~measurement)
  
  p1 = ggplot(data = pd %>% filter(recycled == T)) + 
    geom_density(data = pd, mapping = aes(x = log(value)), 
                 color = "gray20", size = 0.5) + 
    geom_density(mapping = aes(x = log(value), 
                               color = recycled, fill = recycled), 
                 size = 0.5) + 
    facet_wrap(~measurement) +
    guides(color = "none", fill = "none") +
    scale_fill_manual(values = alpha("#56B4E9", 0.7)) + 
    scale_color_manual(values = alpha("#56B4E9", 1)) +
    ylim(0, 1) +
    xlim(-1, 8) +
    theme(axis.text=element_text(size=5),
          axis.title=element_text(size=6,face="bold"), 
          title = element_text(size=8), 
          strip.text.x = element_text(size = 6))
  
  return(p1)
}


