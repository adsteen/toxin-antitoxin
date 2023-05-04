# Makes boxplot of raw data
raw_data_boxplot <- function(df) {
  
  
  df_summ <- df %>%
    group_by(category) %>%
    summarise(mean = mean(gene.count, na.rm=TRUE),
              sd = sd(gene.count, na.rm=TRUE)) #%>%
    #mutate(category = fct_reorder(category, gene.count, mean, .desc=TRUE))
  
  p <- ggplot(df_summ, aes(x=category)) + 
    geom_boxplot(aes(lower = mean-sd, middle=mean, upper=mean+sd,
                     ymin = mean-1.5*sd, ymax = mean+1.5*sd), stat="identity") + 
    #geom_point(position=position_jitter(height = 0.3), alpha = 0.5) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle=-45, hjust=0))
  #p
 p
}


raw_line_plot <- function(df) {
  p <- ggplot(df, aes(x=gene.count, y=freq, color=category)) + 
    geom_point() + 
    geom_line()  + 
    scale_color_brewer(type="qual")
  p
}

raw_bar_plot <- function(df) {
  p <- ggplot(df, aes(x=gene.count, y=freq, fill=category)) + 
    geom_bar(stat="identity", position="dodge") + 
    scale_fill_brewer(type="qual")
  p
}
