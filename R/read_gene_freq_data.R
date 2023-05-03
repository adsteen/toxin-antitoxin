# Read in raw data, removes "total" row
# Makes the df long

read_gene_freq_data <- function(path, rename.vec = NULL) {
  d <- read_csv(path=path, sheet=sheet, range=range) 
  
  if(!is.null(rename.vec)) {
    d <- rename(d, rename.vec)
  }
  browser()
  mutate(n.genes = as.numeric(n.genes)) %>%
    pivot_longer(cols = -1, names_to = "category", values_to = "count") %>%
    group_by(category) %>%
    mutate(freq = count / sum(count, na.rm = TRUE)) %>%
    arrange(category, count)
}
