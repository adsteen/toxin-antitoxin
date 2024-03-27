# function to recreate raw data from the frequency data
recreate_raw <- function(d) {
  #browser()
  # Initialize a zero-row 
  recreated_raw_data <- data.frame("category"=character(0), 
                                   "count"=double(0), 
                                   "gene.count"=double(0)) # create empty data frame to bind to
  #browser()
  # This is a really inelegant way to recreate the raw count data from the summary counts 
  # "count" is the number of genomes that has a given number of relevant genes
  # gene.count is the number of relevant genes that genome has
  # So this loop runs through each row of the original data frame
  # and appends a number of rows equal to the "count" to the new df
  # Thus giving us a "recreated" data frame with a row for every genome observed
  
  # Note how I'm growing the df in the loop, very cool stuff
  for(i in 1:nrow(d)) {
    this_row <- d[i, ]
    
    for(j in 1:this_row$count) {
      if(this_row$count > 0) {
        recreated_raw_data <- rbind(recreated_raw_data, this_row)
      }
    }
  }
  recreated_raw_data 
}