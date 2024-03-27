# Monte Carlo analysis of test statistic distributions

# Introduction

The previous analysis ([normality_analysis.md](normality_analysis.md))
showed that parametric ANOVA is inappropriate for these datasets,
because they’re pretty badly unbalanced and because they’re pretty badly
heteroskedastic.

Here I’m taking a resampling/Monte Carlo approach, in which I
repeatedly, randomly shuffle the labels on our data set, thus creating a
data set for which the null hypothesis (no significant differences among
treatments) is known to be true. Comparing the distribution of test
statistics (the t value) under the null hypothesis to the observed test
statistic. If the observed t statistic is more extreme than (1-$\alpha$;
typically 95%) of the test statistics under the null hypothesis, we
would call the result statistically significant.

# zor-orz

We’ll start with the full analysis for zor-orz.

``` r
library(tidyverse)
library(MASS)
library(furrr) # purrr for parallel processes
library(tictoc)
library(rlang) # this is maybe needed by calculate_mean_diff
plan(multisession) # seems like multicore is way faster but I'm working in RStudio which doesn't support multicore
source("helper_funs.R")

theme_set(theme_classic()) # set the graphical theme
```

Load the data, which are grouped counts but need to be recast as
original observations.

``` r
# Load data - 
# Read in raw zor-orz data
region_zor <- readxl::read_xlsx("data/FINAL data for Steen.xlsx", 
                            sheet = "zor-orz gene number",
                            range = "A11:G18") %>%
  rename(gene.count = Continent) %>%
  pivot_longer(-1, names_to = "category", values_to = "count") %>%
  group_by(category) %>%
  mutate(freq = count / sum(count, na.rm = TRUE))
  
path_zor <- readxl::read_xlsx("data/FINAL data for Steen.xlsx", 
                            sheet = "zor-orz gene number",
                            range = "A1:E8") %>%
  rename(gene.count = Pathotype) %>%
  pivot_longer(-1, names_to = "category", values_to = "count") %>%
  group_by(category) %>%
  mutate(freq = count / sum(count, na.rm = TRUE))

# recreate raw data based on frequencies of gene abundance
raw_region_zor <- recreate_raw(region_zor) %>%
  arrange(category) # this appears to have worked
raw_path_zor <- recreate_raw(path_zor)
```

# Monte Carlo simulations

Based on the results of the `normality_analysis.qmd` workbook, it is
inappropriate to assume that f-values are distributed according to the
f-distribution.

So first we’ll determine the distribution of f values under the null
hypothesis. We do 10,000 simulations, which (as we’ll see) is plenty to
get a smooth distribution and to get a reasonably precise idea of where
the 5% cutoff lies.

``` r
# Set the number of monte carlo replicates
n <- 10000 # can change back to 10000

#tic()
reg.f.vec <- 1:n |> 
  future_map_dbl(.f = shuf_calc_f,
                 df = raw_region_zor,
                 .options = furrr_options(seed = 2112)) # great, if you're *really* into prog rock
path.f.vec <- 1:n |> 
  future_map_dbl(.f = shuf_calc_f,
                 df = raw_path_zor,
                 .options = furrr_options(seed = 444)) # aesthetically pleasing, dontcha think?
#toc() # Runs in about 12 seconds on 8 core macbook pro
#... or about 6 seconds with multicore, but that's not an option

# put the simulated f values in a data frame
f_vals <- data.frame(reg.sim.f = reg.f.vec,
                     path.sim.f = path.f.vec)
```

Now we measure the two actual f values:

``` r
# Pull out actual f values
region_model <- lm(gene.count ~ category, data = raw_region_zor)
reg.f.real <- summary(aov(region_model))[[1]][1,4] # 9.335
path_model <- lm(gene.count ~ category, data = raw_path_zor)
path.f.real <- summary(aov(path_model))[[1]][1,4] # 46.035
```

How do the real f values compare to the simulated, null-hypothesis
values?

``` r
p_reg_hist <- ggplot(f_vals, aes(x=reg.sim.f)) + 
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = reg.f.real, color="red")  + 
  ggtitle("data by region")
print(p_reg_hist)
```

![](monte-carlo_files/figure-commonmark/plot_by_regions-1.png)

``` r
p_path_hist <- ggplot(f_vals, aes(x=path.sim.f)) + 
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = path.f.real, color="red") + 
  ggtitle("data by pathology")
print(p_path_hist)
```

![](monte-carlo_files/figure-commonmark/plot_by_pathology-1.png)

So: I have simulated 10,000 and found that, for each case, the actual
measured *f* values are much, much larger than they would be likely to
be if the null hypothesis were true - so much larger that we can’t
calculate a p value, because none of our 10,000 simulations captured a
*f* value that big. We can say conservatively say that, in each case, p
\< 10^{-4}.

In summary:

| f-value      | Simulated maximum | Observed |
|--------------|-------------------|----------|
| by pathology | 6.7               | 46       |
| by region    | 6.1               | 9.3      |

# Monte Carlo Tukey Test

A Tukey test works by comparing means of all possible combinations of
populations (in this case, regions or pathotypes) and then comparing to
a studentized range distribution. I’m going to do exactly this, except
that the studentized range distribution is replaced with the observed
distribution of mean differences in shuffled data.

``` r
# Let's make a function to calculate actual means and then simulated means
tic()
n.tukey <- 1e4
path_diffs <- monte_carlo_tukey(raw_path_zor, n.tukey)
region_diffs <- monte_carlo_tukey(raw_region_zor, n.tukey)
toc()
```

    35.471 sec elapsed

# Results

## zor-orz results

``` r
knitr::kable(path_diffs)
```

| diff.id                            | cutoff.diff | mean.diff | sig.diff |
|:-----------------------------------|------------:|----------:|:---------|
| bacteremia-healthy                 |   0.5558217 | 0.2684063 | FALSE    |
| bacteremia-intestinal disease      |   0.4466198 | 1.3412415 | TRUE     |
| bacteremia-urinary disease         |   0.4985302 | 0.0912472 | FALSE    |
| healthy-intestinal disease         |   0.3769255 | 1.0728352 | TRUE     |
| healthy-urinary disease            |   0.4319432 | 0.1771591 | FALSE    |
| intestinal disease-urinary disease |   0.2842389 | 1.2499943 | TRUE     |

### How to interpret this table

This is a table comparing differences in the mean gene number between
each pair of groups. For instance, the top row is `Africa-Asia`,
`mean.diff` indicates that the absolute value of the difference in the
mean number of zor/orz genes between Africa and Asia is 0.091.
`cutoff.diff`, the “cutoff” above which a difference would be
statistically significant (p \< 0.05), is 0.286. 0.091 is not greater
than 0.286, so there is no significant difference. Thus, the `sig.diff`
entry is FALSE.

`Asia-North America`, does have a significant difference: `cutoff.diff`
is 0.087, `mean.diff` is 0.183, so `sig.diff` is TRUE.

``` r
knitr::kable(region_diffs)
```

| diff.id                     | cutoff.diff | mean.diff | sig.diff |
|:----------------------------|------------:|----------:|:---------|
| Africa-Asia                 |   0.5852703 | 0.0325988 | FALSE    |
| Africa-Europe               |   0.5825097 | 0.3538880 | FALSE    |
| Africa-North America        |   0.5878287 | 0.3771936 | FALSE    |
| Africa-Oceania              |   0.7363064 | 0.8914286 | TRUE     |
| Africa-South America        |   0.8505411 | 0.7085714 | FALSE    |
| Asia-Europe                 |   0.1637049 | 0.3864868 | TRUE     |
| Asia-North America          |   0.1708636 | 0.4097924 | TRUE     |
| Asia-Oceania                |   0.4922427 | 0.9240273 | TRUE     |
| Asia-South America          |   0.6270674 | 0.6759727 | TRUE     |
| Europe-North America        |   0.1564865 | 0.0233056 | FALSE    |
| Europe-Oceania              |   0.4854046 | 0.5375405 | TRUE     |
| Europe-South America        |   0.6309310 | 1.0624595 | TRUE     |
| North America-Oceania       |   0.4927248 | 0.5142350 | TRUE     |
| North America-South America |   0.6277024 | 1.0857650 | TRUE     |
| Oceania-South America       |   0.7886147 | 1.6000000 | TRUE     |

# Human vs non-human animal

Based on a reviewer request, we need to do the same test with humans vs
non-human animals, for zor-orz and tis-istR.

Again, the first step is to recreate the raw data:

``` r
human_tis <- readxl::read_xlsx("data/FINAL data for Steen after reviewer comments.xlsx", 
                            sheet = "to analyze based on reviewers c",
                            range = "A2:C4") %>%
  rename(gene.count = `gene copy number`) %>%
  pivot_longer(-1, names_to = "category", values_to = "count") %>%
  group_by(category) %>%
  mutate(freq = count / sum(count, na.rm = TRUE))

human_zor <- readxl::read_xlsx("data/FINAL data for Steen after reviewer comments.xlsx", 
                            sheet = "to analyze based on reviewers c",
                            range = "A10:C17") %>%
  rename(gene.count = `gene copy number`) %>%
  pivot_longer(-1, names_to = "category", values_to = "count") %>%
  group_by(category) %>%
  mutate(freq = count / sum(count, na.rm = TRUE))

raw_human_tis <- recreate_raw(human_tis) #%>%
raw_human_orz <- recreate_raw(path_zor)
```

Now we have to compare the actual f-values for the ANOVA to the
distributions of simulated f-values.

``` r
# Make up artificial f values for the human/nonhuman tis data
# Not even sure whether this is different, but it works after fixing some silly errors
shuf_calc_f_alt <- function(throwaway, df) { # This is some hacky shit, but
  # map iterates over an iterator, which I don't actually use. 
  # So I've got a throwaway parameter that doesn't get used. I'm sure 
  # there's a better way to do this, but so it goes.
  
  df.rows <- nrow(df)
  
  df <- df %>% 
    ungroup() %>%
    mutate(shuf.data = sample(gene.count, size = df.rows, replace=FALSE)) 
  
  m <- aov(shuf.data ~ category, data = df)
  f <- summary(m)[[1]][1,4]
  f
}

# Iterate the shuffle-and-count-f function n times; place data in a vector
human.tis.f.vec <- future_map_dbl(seq_along(1:n), 
                            shuf_calc_f_alt, 
                            df=raw_human_tis, #nrow=nrow(raw_human_tis), 
                             .options = furrr_options(seed = 55)) # a speed at which I can't drive
human.orz.f.vec <- future_map_dbl(seq_along(1:n), 
                             shuf_calc_f_alt, 
                             df=path_zor, #nrow=nrow(raw_human_orz), 
                             .options = furrr_options(seed = 99991)) # allegedly the largest 5-digit prime number 
#toc() # Runs in about 14 seconds on 6 core macbook pro; pretty sweet

# put the simulated f values in a data frame
human_f_vals <- data.frame(human.tis.f = human.tis.f.vec,
                     human.orz.f = human.orz.f.vec)
```

Now calculate the actual f values:

``` r
human_tis_model <- lm(gene.count ~ category, data = raw_human_tis)
human.tis.real <-  summary(aov(human_tis_model))[[1]][1,4] # 171.1295
human_orz_model <- lm(gene.count ~ category, data = raw_human_orz)
human.orz.real <-  summary(aov(human_orz_model))[[1]][1,4] # 46.03
```

``` r
p_human_tis <- ggplot(human_f_vals, aes(x=human.tis.f)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = human.tis.real, color="red") +
  ggtitle("human-vs-nonhuman tis")
print(p_human_tis)
```

![](monte-carlo_files/figure-commonmark/unnamed-chunk-4-1.png)

``` r
p_human_orz <- ggplot(human_f_vals, aes(x=human.orz.f)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = human.orz.real, color="red") +
  ggtitle("human-vs-nonhuman orz")
print(p_human_orz)
```

![](monte-carlo_files/figure-commonmark/unnamed-chunk-5-1.png)

Looks like these are both highly significant, again! As before, we can’t
assign a p value because the observed f value is (way) more extreme than
any observed f value in the simulations, but we can say that p \<
10^{-4}.
