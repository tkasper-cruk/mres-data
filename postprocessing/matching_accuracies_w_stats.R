library(tidyverse)

get_p_vs_random <- function(data){
  tests <- tibble(ident= unique(data$type))
  scores <- list()
  p_values <- list()
  for (src in tests$ident){
    t <- wilcox.test(
      x = filter(data,type==src)$accuracy,
      y = filter(data,type=="other_random")$accuracy,
      paired = TRUE,
      alternative = "greater",
      exact = FALSE
    )
    p_values[src] <- t$p.value
    scores[src] <- t$statistic
  }
  tests$statistic <- scores
  tests$p <- p_values
  tests$p_corr <- p.adjust(p_values, method = "bonferroni")
  return(tests)
}

get_p_normality <- function(data){
  tests <- tibble(ident = unique(data$type))
  p_values <- list()
  scores <- list()
  for (src in tests$ident){
    t <-  shapiro.test(filter(data,type==src)$accuracy)
    p_values[src] <-t$p.value
    scores[src] <- t$statistic
  }
  tests$statistic <- scores
  tests$p <- p_values
  tests$p_corr <- p.adjust(p_values, method = "bonferroni")
  return(tests)
}


data_dir <- "accuracies/"
# Normality tests and pairwise wilcoxon test for preliminary tests
pretest_data <- read_csv(paste0(data_dir,"line2line_acc_pretest.csv"))
pretest_data <- pretest_data %>% 
  mutate(
    type = str_replace_all(paste(pair,method,sep="-"),c("random_random-random"="other_random"))
  )

##Normality
prelim_test_norm <- get_p_normality(pretest_data)


## get median and sd 
pretest_stats <- pretest_data %>%
  group_by(type) %>%
  summarise(
    median = median(accuracy),
    sd = sd(accuracy)
  )

pretest_stats %>%
  write_csv("statistics_and_summaries/line2line_prelim_summary.csv")
##comparisons to random (>)
pretest_v_random <- get_p_vs_random(pretest_data)

##comparisons between modalities (!=)
pretest_pairwise_p = list()
pretest_pairwise_stat = list()
for (src1 in unique(pretest_data$type)){
  if (src1 == "other_random"){ #ignore random for power
    next
  }
  for (src2 in unique(pretest_data$type)){
    if (src2 == "other_random"){ #ignore random for power
      next
    }
    if (str_split_i(src1,"-",2)!=str_split_i(src2,"-",2)){
      next #ignore pearson vs cosine comparisons to keep power
    }
    if (src1==src2){
      next #ignore self comparisons for power
    }
    pair <- paste(src1,src2,sep="|")
    t <- wilcox.test(
      x = filter(pretest_data,type==src1)$accuracy,
      y = filter(pretest_data,type==src2)$accuracy,
      paired = TRUE,
      alternative = "two.sided",
      exact = FALSE
    )
    pretest_pairwise_p[pair] <- t$p.value
    pretest_pairwise_stat[pair] <- t$statistic
  }
}
pretest_pairwise <- tibble(
  ident = names(pretest_pairwise_p),
  statistic=pretest_pairwise_stat,
  p=pretest_pairwise_p,
  p_corr = p.adjust(pretest_pairwise_p, method = "bonferroni")
  )


# plot results

pretest_stats %>%
  filter(
    type != "other_random"
  ) %>%
  mutate(
    median = round(median,2),
    metric = str_replace_all(
      str_split_i(type,"-",2),
      c("cosine"="Cosine similarity","pearson"="PCC")
    ),
    dna = str_replace_all(
      str_split_i(str_split_i(type,"-",1),"_",1),
      c ("cor"="Batch\ncorrected","norm"="Normalised","pre"="Raw")
    ),
    rna = str_replace_all(
      str_split_i(str_split_i(type,"-",1),"_",2),
      c ("cor"="Batch\ncorrected","norm"="Normalised","pre"="Raw")
    )
  ) %>%
  select(
    c(-type,-sd)
  ) %>%
  ggplot(aes(x=dna,y=rna,fill = median,label=median)) +
  facet_wrap(vars(metric))+
  geom_tile(colour="white",lwd = 0.75, linetype = 1) +
  geom_text(colour="white",size = 4) +
  xlab("DNA Copy Number Profile type") +
  ylab("Gene expression profile type") + 
  labs(
    fill = "Median\naccuracy"
  ) + 
  theme_bw()
ggsave("figures/prelim_test_accuracies.png",width=6,height = 3,dpi=1200)


# Normality tests and pairwise wilcoxon test for line2line
all_data <- read_csv(paste0(data_dir,"line2line_acc_all.csv"))
all_data$type <- paste(all_data$tool,all_data$method,sep="_")


all_norm <- get_p_normality(filter(all_data,type!="other_clonealign"))

## get median and sd 
all_stats <- all_data %>%
  group_by(type) %>%
  summarise(
    median = median(accuracy),
    sd = sd(accuracy)
  )
all_stats %>%
  write_csv("statistics_and_summaries/line2line_all_summary.csv")
all_v_random <- get_p_vs_random(filter(all_data,type!="other_clonealign"))

###plotting
all_data %>%
  mutate(
    tool = factor(
      str_replace_all(
        tool,
        c(
          "copykat"="CopyKAT","copyvae"="copyVAE","rna"="Gene Expression",
          "numbat_cat"="numbat\nCategorical","numbat_cn"="numbat\nInteger",
          "other"="Other"
        )
      ),
      levels = c(
        "Gene Expression","CopyKAT","copyVAE","numbat\nCategorical","numbat\nInteger","Other"
      )
    ),
    method = factor(
      str_replace_all(
        method,
        c(
          "cosine"="Cosine similarity",
          "pearson"="Pearson's r",
          "spearman"=paste("Spearman's",expression(rho)),
          "kendall"=paste("Kendall's",expression(tau)),
          "clonealign"="clonealign",
          "macrodna"="MaCroDNA",
          "random" = "Random\npermutation"
          
        )
      ),
      levels = c(
        "Cosine similarity","Pearson's r",paste("Spearman's",expression(rho)),
        paste("Kendall's",expression(tau)),"clonealign","MaCroDNA","Random\npermutation"
      )
    )
  ) %>%
  ggplot(aes(x=tool,y=accuracy,colour = method)) +
  geom_boxplot()+
  labs(
    colour = "Similarity metric"
  ) + 
  xlab("") + 
  ylim(0,1) + 
  ylab("Cellline assignment accuracy") + 
  #ggtitle("Performance of the tested cell-matching methods in cellline matching") +
  theme_bw() +
  theme(legend.position="bottom")

ggsave("figures/line2line_all_accuracies.png",width=6.5,height = 4.5,dpi=1200)

#subsets of line2line
peo_data <- read_csv(paste0(data_dir,"line2line_acc_PEO_only.csv"))
peo_data$type <- paste(peo_data$tool,peo_data$method,sep="_")

## get median and sd 
peo_stats <- peo_data %>%
  group_by(type) %>%
  summarise(
    median = median(accuracy),
    sd = sd(accuracy)
  )
peo_stats %>%
  write_csv("statistics_and_summaries/line2line_peo_summary.csv")
##Statistical tests
peo_norm <- get_p_normality(peo_data)
peo_v_random <- get_p_vs_random(peo_data)

##Plotting
peo_data %>%
  mutate(
    tool = factor(
      str_replace_all(
        tool,
        c(
          "copykat"="CopyKAT","copyvae"="copyVAE","rna"="Gene Expression",
          "numbat_cat"="numbat\nCategorical","numbat_cn"="numbat\nInteger",
          "other"="Other"
        )
      ),
      levels = c(
        "Gene Expression","CopyKAT","copyVAE","numbat\nCategorical","numbat\nInteger","Other"
      )
    ),
    method = factor(
      str_replace_all(
        method,
        c(
          "cosine"="Cosine similarity",
          "pearson"="PCC",
          "spearman"= "SCC",
          "kendall"= "KCC",
          "random" = "Random\npermutation"
          
        )
      ),
      levels = c(
        "Cosine similarity","PCC","SCC","KCC","clonealign","MaCroDNA","Random\npermutation"
      )
    )
  ) %>%
  ggplot(aes(x=tool,y=accuracy,colour = method)) +
  geom_boxplot()+
  labs(
    colour = "Similarity metric"
  ) + 
  xlab("") + 
  ylim(0,1) + 
  ylab("Cellline assignment accuracy") + 
  #ggtitle("Performance of the tested cell-matching methods for cell line\nmatching on PEO1, PEO1-Missense and PEO1-Stop") +
  theme_bw() +
  theme(legend.position="bottom")

ggsave("figures/line2line_peo_accuracies.png",width=6.5,height = 4.5,dpi=900)


#PEO1 vs NA12878
easy_data <- read_csv(paste0(data_dir,"line2line_acc_PEO_NA12878.csv"))
easy_data$type <- paste(easy_data$tool,easy_data$method,sep="_")
##median and SD
easy_stats <- easy_data %>%
  group_by(type) %>%
  summarise(
    median = median(accuracy),
    sd = sd(accuracy)
  )
easy_stats %>%
  write_csv("statistics_and_summaries/line2line_peo_na12878_summary.csv")
##Statistical tests
easy_norm <- get_p_normality(easy_data)
easy_v_random <- get_p_vs_random(easy_data)

easy_data %>%
  mutate(
    tool = factor(
      str_replace_all(
        tool,
        c(
          "copykat"="CopyKAT","copyvae"="copyVAE","rna"="Gene Expression",
          "numbat_cat"="numbat\nCategorical","numbat_cn"="numbat\nInteger",
          "other"="Other"
        )
      ),
      levels = c(
        "Gene Expression","CopyKAT","copyVAE","numbat\nCategorical","numbat\nInteger","Other"
      )
    ),
    method = factor(
      str_replace_all(
        method,
        c(
          "cosine"="Cosine\nsimilarity",
          "pearson"="PCC",
          "spearman"="SCC",
          "kendall"="KCC",
          "random" = "Random\npermutation"
          
        )
      ),
      levels = c(
        "Cosine\nsimilarity","PCC","SCC","KCC","clonealign","MaCroDNA","Random\npermutation"
      )
    )
  ) %>%
  ggplot(aes(x=tool,y=accuracy,colour = method)) +
  geom_boxplot()+
  labs(
    colour = "Similarity metric"
  ) + 
  xlab("") + 
  ylim(0,1) + 
  ylab("Cellline assignment accuracy") + 
  #ggtitle("Performance of the tested cell-matching methods in cellline matching on PEO1 and NA12878") +
  theme_bw() +
  theme(legend.position="bottom")

ggsave("figures/line2line_easy_accuracies.png",width=6.5,height = 4.5,dpi=1200)


## cell2cell
cell_data <- read_csv(paste0(data_dir,"cell2cell_acc.csv"))
cell_data$type <- paste(cell_data$tool,cell_data$method,sep="_")

## get median and sd 
cell_stats <- cell_data %>%
  group_by(type) %>%
  summarise(
    median = median(accuracy),
    sd = sd(accuracy)
  )
cell_stats %>%
  write_csv("statistics_and_summaries/cell2cell_acc_summary.csv")
##Statistical tests
cell_norm <- get_p_normality(filter(cell_data,type!="other_clonealign"))
cell_v_random <- get_p_vs_random(filter(cell_data,type!="other_clonealign"))


###plotting
cell_data %>%
  mutate(
    tool = factor(
      str_replace_all(
        tool,
        c(
          "copykat"="CopyKAT","copyvae"="copyVAE","rna"="Gene Expression",
          "other"="Other"
        )
      ),
      levels = c(
        "Gene Expression","CopyKAT","copyVAE","Other"
      )
    ),
    method = factor(
      str_replace_all(
        method,
        c(
          "cosine"="Cosine\nsimilarity",
          "pearson"="PCC",
          "spearman"= "SCC",
          "kendall"= "KCC",
          "clonealign"="clonealign",
          "macrodna"="MaCroDNA",
          "random" = "Random\npermutation"
          
        )
      ),
      levels = c(
        "Cosine\nsimilarity","PCC","SCC","KCC","clonealign","MaCroDNA","Random\npermutation"
      )
    )
  ) %>%
  ggplot(aes(x=tool,y=accuracy,colour = method)) +
  geom_boxplot()+
  labs(
    colour = "Similarity metric"
  ) + 
  xlab("") + 
  #ylim(0,1) + 
  ylab("Cell assignment accuracy") + 
  #ggtitle("Performance of the tested cell-matching methods in cell matching") +
  theme_bw() +
  theme(legend.position="bottom")
ggsave("figures/cell2cell_accuracies.png",width=6.5,height = 4.5,dpi=1200)

#recall
recall_data <- read_csv(paste0(data_dir,"cell2cell_rec.csv")) %>%
  select(1:4) %>%
  rename(
    "accuracy" = "recall"
  )
recall_data$type <- paste(recall_data$tool,recall_data$metric,sep="_")
##median and sd
recall_stats <- recall_data %>%
  group_by(type,top_n) %>%
  summarise(
    median = median(accuracy),
    sd = sd(accuracy)
  )
recall_stats %>%
  mutate(
    profile = str_split_i(type,"_",1),
    metric = str_split_i(type,"_",2),
    median = round(median,4),
    sd = round(sd,4)
  ) %>%
  write_csv("statistics_and_summaries/cell2cell_rec_summary.csv")

##statistical test -> grouped by top_n
norm_tibbles <- recall_data %>%
  filter(type!="other_random") %>%
  group_split(top_n) %>%
  map(get_p_normality)
names(norm_tibbles) <- unique(recall_data$top_n)
recall_norm <- bind_rows(norm_tibbles,.id="top_n")
recall_norm$p_corr <- p.adjust(recall_norm$p,method="bonferroni")

random_tibbles <- recall_data %>%
  group_split(top_n) %>%
  map(get_p_vs_random)
names(random_tibbles) <- unique(recall_data$top_n)
recall_random <- bind_rows(random_tibbles,.id="top_n")
recall_random$p_corr <- p.adjust(recall_random$p,method="bonferroni")


recall_data %>%
  mutate(
    tool = factor(
      str_replace_all(
        tool,
        c(
          "copykat"="CopyKAT","copyvae"="copyVAE","rna"="Gene\nExpression",
          "other"="Other"
        )
      ),
      levels = c(
        "Gene\nExpression","CopyKAT","copyVAE","Other"
      )
    ),
    method = factor(
      str_replace_all(
        metric,
        c(
          "cosine"="Cosine\nsimilarity",
          "pearson"="PCC",
          "spearman"="SCC",
          "kendall"="KCC",
          "random" = "Random\npermutation"
          
        )
      ),
      levels = c(
        "Cosine\nsimilarity","PCC","SCC","KCC","Random\npermutation"
      )
    )
  ) %>%
  ggplot(aes(x=tool,y=accuracy,colour = method)) +
  geom_boxplot()+
  labs(
    colour = "Similarity metric"
  ) + 
  facet_wrap(vars(top_n),scales = "free")+
  xlab("") + 
  #ylim(0,1) + 
  ylab("Matched Cell recall") + 
  #ggtitle("Performance of the tested cell-matching methods in cell matching") +
  theme_bw() +
  theme(legend.position="bottom")
ggsave("figures/cell2cell_recall.png",width=8,height = 5,dpi=1200)


#collate stats for export
norm_tests <- bind_rows(
  list(
    cell2cell_acc = cell_norm,
    cell2cell_recall = select(recall_norm,-top_n),
    line2line_prelim = prelim_test_norm,
    line2line_all = all_norm,
    line2line_peo = peo_norm,
    line2line_easy = easy_norm
  ),
  .id = "Experiment"
)  %>% mutate(across(where(is.list), unlist)) 

random_tests <- bind_rows(
  list(
    cell2cell_acc = cell_v_random,
    cell2cell_recall = select(recall_random,-top_n),
    line2line_prelim = pretest_v_random,
    line2line_prelim_pairwise = pretest_pairwise,
    line2line_all = all_v_random,
    line2line_peo = peo_v_random,
    line2line_easy = easy_v_random
  ),
  .id = "Experiment"
) %>% mutate(across(where(is.list), unlist)) 

norm_tests %>%
  write_csv("statistics_and_summaries/shapiro.csv")
random_tests %>%
  write_csv("statistics_and_summaries/wilcoxon.csv")