library(tidyverse)

args <- commandArgs(T)
work_dir = args[1]
prefix = args[2]

df = read_table(paste0(work_dir, prefix, ".clumped")) %>% 
	 mutate(NUM=S0001+S001) %>% 
	 arrange(CHR, BP) %>% 
	 select(CHR, SNP, BP, P, NUM, SP2)

df_pos_ID = read_delim(paste0(work_dir, prefix, ".rsID.bp.txt"), " ",col_names = c("CHR","BP","SNP"))

get_interval = function(x) {
  if (x=="NONE") {return(tibble(interval="-",start=NA,end=NA))}
  list_SNP = tibble(SNP=str_split(x, "\\,") %>% unlist() %>% str_replace_all("\\(1\\)", "")) %>%
    left_join(df_pos_ID, by="SNP")
  interval_start = min(list_SNP$BP)
  interval_end = max(list_SNP$BP)
  tibble(interval=paste0(interval_start,"-",interval_end),
    start=sprintf("%.2f",interval_start/1e6),
    end=sprintf("%.2f",interval_end/1e6))
}

df_unusual_regions = df %>%
  rowwise %>%
  mutate(get_interval(SP2)) %>%
  ungroup() %>%
  select(-SP2) %>%
  arrange(CHR,BP)

write_tsv(df_unusual_regions, paste0(work_dir, prefix, ".clump.cleaned.txt"), na="")

