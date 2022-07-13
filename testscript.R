
# Packages load/install -------------------------------------------------
if (!require("pacman")) install.packages("pacman") #check pacman is loaded and install if not
pacman::p_load(readr, #bring in data
               tidyverse, # for data manipulation 
               naniar, #replace na's
               gt,#nice table
               gtsummary,
               #changepoint,
               CPAT) #finding peaks



# functions ---------------------------------------------------------------

cleanup <- function(x){
  x %>% 
    select(c("(Volts)","(Amps)")) %>% 
    rename(Volts = "(Volts)", Amps = "(Amps)") %>% 
    mutate(lg_Amps = log(abs(Amps))) 
}

bp_calc_I <- function(x){
  
  bp_index <- function(x){
    y<- x %>% pull(Amps)
    bp <-  CUSUM.test(y)
    bp <- bp[["estimate"]][["t*"]]
     return(bp)
  }
  x %>% 
    summarise(bp_idx = bp_index(.),
              bp_v = nth(Volts,bp_idx[[1]]),
              bp_i = nth(Amps,bp_idx[[1]]),
              pbp_v = nth(Volts,(bp_idx[[1]]+1)),
              pbp_i = nth(Amps,(bp_idx[[1]]+1))) 
 }

bp_log_I <- function(x){
  
 # bp_index_lg <- function(x){
    y<- x %>% pull(lg_Amps)
    bp <-  CUSUM.test(y)
    bp <- bp[["estimate"]][["t*"]]
    # return(bp)
 # }
  x %>% 
    summarise(bp_lg_idx = bp_index(.),
              bp_lg_v = nth(Volts,bp_lg_idx[[1]]),
              bp_lg_i = nth(lg_Amps,bp_lg_idx[[1]])) 
}

plots <- function(files,dat,flag,bp_v,bp_i){
  dat %>% 
    ggplot(.,aes(x = Volts, y = Amps, color=flag)) + 
    geom_point() +
    geom_line() +
    geom_vline(xintercept = bp_v, linetype="dotted", color = "orange", size=0.75) +
    geom_hline(yintercept = bp_i, linetype="dotted", color = "orange", size=0.75) +
    scale_y_continuous(trans='log10')+
    scale_color_manual(values = c("ok" = "green",
                                  "broken" ="red")) +
    labs(title = files,
         y = "Log(Amps)")
}

# Data import -------------------------------------------------------------

data_path <- "./raw_data/"   
files <- dir(data_path, pattern = "*.csv", recursive = T) 

df <- tibble(files = files) %>%
  mutate( dat = map(files, ~ read_csv(file.path(data_path, .),skip = 46,show_col_types = F) ) ) %>% 
  mutate( die = sub("\\/.*", "", files),
          die = str_remove_all(die, "Die "),
          field = sub(".*\\/(.)-.*", "\\1", files),
          device = sub(".*-", "", files),
          device = str_remove_all(device, ".csv")) %>% 
  mutate( dat = map(dat, cleanup)) %>% 
  mutate(files = str_remove_all(files,".csv"))

df_bp <- df %>% 
  mutate( #bp_lg = map(dat,bp_log_I),
          bp = map(dat,bp_calc_I))

df_flag <- df_bp %>% 
  unnest(bp) %>% 
  mutate(flag = case_when(
  #  bp_v <= 5 | bp_v > 55 ~ "broken",
    (bp_v > 5 & bp_v < 50) & (bp_i <= 0.0001 & pbp_i > 0.000001) ~ "ok",
    TRUE ~ "broken"
  )) %>% 
  mutate(grp = case_when(
    startsWith(die,"F") ~ "Flat",
    TRUE ~ "Round"))


  
# plotting ----------------------------------------------------------------

df_plots <- df_flag %>% 
  mutate(plots = pmap(list(files,dat,flag,bp_v,bp_i),plots)) 

#df_plots$plots[[1]]
#print(df_plots$plot)

plots_save <- df_flag %>% 
  mutate(plot = pmap(list(files,dat,flag,bp_v,bp_i),plots),
         filename = paste0(files, ".png"),
         filename = str_replace_all(filename,"/","")) %>% 
  select(filename, plot)

pwalk(plots_save, ggsave, path = getwd())



#plot histogram of the breakpoint voltage per group
df_flag  %>% 
  filter(bp_v >= 20)%>% 
  filter(., flag == 'ok') %>% 
  select(bp_v,grp) %>% 
  ggplot( aes(x=bp_v, fill=grp)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=15) +
  labs(fill="")

  
df_flag %>% 
  filter(bp_v >= 20)%>% 
  filter(., flag == 'ok') %>% 
  select(bp_v,bp_i,grp) %>% 
  group_by(grp) %>% 
  summarise(min = min(bp_v)
            ,max = max(bp_v)
            ,mean = mean(bp_v)
            ,sd = sd(bp_v)
            ,n = n()
            ,median = median(bp_v)
            ,q25 = quantile(bp_v, .25)
            ,q75 = quantile(bp_v, .75)) 
  

  
df_flag %>% 
 # filter(bp_v >= 20)%>% 
  filter(., flag == 'ok') %>% 
  select(bp_v,grp) %>% 
  tbl_summary(by= grp,
              statistic = list(all_continuous() ~ "{mean} ({sd})",
                               all_categorical() ~ "{n} / {N} ({p}%)"),
              digits = all_continuous() ~ 2) %>% add_p()

library(tidyverse)
library(broom)
library(knitr)
library(ggfortify)

"{mean} ({sd})"
"{median} ({p25}, {p75})"
pc_fit <- df_flag %>% 
  select(files,bp_v,bp_i) %>% 
  nest() %>% 
  mutate(pca = map(data, ~ prcomp(.x %>% select(-files), 
                                  center = TRUE, scale = TRUE)),
         pca_aug = map2(pca, data, ~augment(.x, data = .y)))

pc_fit %>% 
  unnest(pca_aug) 

%>% 
  summarize(across(contains("PC")), ~ list(var)) %>% 
  gather(key = pc, value = variance) %>% 
  mutate(var_exp = variance/sum(variance),
         cum_var_exp = cumsum(var_exp),
         pc = str_replace(pc, ".fitted", ""))



df_flag %>% 
  select(files,dat) %>% 
  unnest(dat) %>% 
  select(Volts,Amps) %>% 
  mutate(Amps = as.numeric(Amps),
         Volts = as.numeric(Volts))  %>% 
  prcomp(scale=TRUE) # do P
  
  group_by(files) %>% 
  select(where(is.numeric)) %>% 
  prcomp(scale=TRUE) # do PCA


# PCA testing -------------------------------------------------------------
df_flag %>% 
  select(bp_v,bp_i)


pca_fit <- biopsy %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  scale() %>% # scale data
  prcomp() # d


# Junk --------------------------------------------------------------------

# df_plots <- df_flag %>% 
#   mutate(plot = map2(dat,files, function(.x,.y) { 
#     .x %>% 
#       ggplot(.,aes(x = Volts, y = log(Amps))) + geom_point() +
#       labs(title = .y)  }))

# df_plots <- df_flag %>% 
#   mutate(plot = map2(dat,files, function(.x,.y) { 
#     .x %>% 
#       ggplot(.,aes(x = Volts, y = Amps)) + geom_point() +
#       labs(title = .y)  })) 




# 
# i_max <- max(x$Amps)
# i_min <- min(x$Amps)
# Imid <-  (i_max - i_min)*0.5 + i_min
# br_pos <-  which.max(x$Amps > Imid) 
# return(br_pos)

# library(CPAT)
# 
# y <- c(rnorm(10, mean = 0), rnorm(90, mean = 2))
# 
# plot(z)  # If you want to visualize the data
# CUSUM.test(z)  # The CUSUM test
# DE.test(z)     # The Darling-Erdös test
# HR.test(z)     # The Rényi-type test introduced in a forthcoming paper
