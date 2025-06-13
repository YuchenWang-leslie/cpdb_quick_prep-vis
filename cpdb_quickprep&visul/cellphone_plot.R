library(purrr)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(scales)
INpvalue <- read.delim('Infant_non_AMD/statistical_analysis_pvalues_Infant_non_AMD_results.txt',check.names = F)
INmean <- read.delim('Infant_non_AMD/statistical_analysis_means_Infant_non_AMD_results.txt',check.names = F)

ONpvalue <- read.delim('Old_non_AMD/statistical_analysis_pvalues_Old_non_AMD_results.txt',check.names = F)
ONmean <- read.delim('Old_non_AMD/statistical_analysis_means_Old_non_AMD_results.txt',check.names = F)

OApvalue <- read.delim('Old_AMD/statistical_analysis_pvalues_Old_AMD_results.txt',check.names = F)
OAmean <- read.delim('Old_AMD/statistical_analysis_means_Old_AMD_results.txt',check.names = F)

unprocessed_p_list <- list(INpvalue,ONpvalue,OApvalue)
unprocessed_means_list <- list(INmean,ONmean,OAmean)

filter_columns <- function(df, sender = NULL, receptor = NULL) {
  col_names <- colnames(df)
  if (!all(grepl("\\|", col_names))) {
    stop("部分列名不符合'sender|receptor'的格式要求")
  }
  split_names <- strsplit(col_names, "\\|")
  sender_part <- sapply(split_names, function(x) x[1])
  receptor_part <- sapply(split_names, function(x) x[2])
  keep <- rep(TRUE, length(col_names))
  if (!is.null(sender)) {
    keep <- keep & (sender_part %in% sender)
  }
  if (!is.null(receptor)) {
    keep <- keep & (receptor_part %in% receptor)
  }
  df[, keep, drop = FALSE]
}

process_cpdb_pvals <- function(unprocessed_list, 
                                       senders = NULL, 
                                       receptors = NULL, 
                                       pval_threshold = 0.05) {
  # 初始化结果列表
  processed_list <- vector("list", length(unprocessed_list))
  
  for (i in seq_along(unprocessed_list)) {
    df <- unprocessed_list[[i]]
    rownames(df) <- df$interacting_pair
    df <- df[, 14:ncol(df), drop = FALSE]
    df <- filter_columns(df, sender = senders, receptor = receptors)
    df <- df[, sapply(df, function(col) {
      if(is.numeric(col)) {
        min(col, na.rm = TRUE) <= pval_threshold
      } else {
        TRUE
      }
    }), drop = FALSE]
    df <- df %>% 
      filter(if_any(where(is.numeric), ~ .x < pval_threshold))
    processed_list[[i]] <- df
  }
  if (!is.null(names(unprocessed_list))) {
    names(processed_list) <- names(unprocessed_list)
  }
  
  return(processed_list)
}

process_cpdb_means <- function(unprocessed_list, #  对于means，只需要筛选细胞对
                               senders = NULL, 
                               receptors = NULL) {
  # 初始化结果列表
  processed_list <- vector("list", length(unprocessed_list))
  
  for (i in seq_along(unprocessed_list)) {
    df <- unprocessed_list[[i]]
    rownames(df) <- df$interacting_pair
    df <- df[, 14:ncol(df), drop = FALSE]
    df <- filter_columns(df, sender = senders, receptor = receptors)
    processed_list[[i]] <- df
  }
  if (!is.null(names(unprocessed_list))) {
    names(processed_list) <- names(unprocessed_list)
  }
  
  return(processed_list)
}

get_common_rownames <- function(df_list) {
  df_list %>% 
    map(rownames) %>% 
    reduce(intersect)
}

processed_p_list <- process_cpdb_pvals(
  unprocessed_list = unprocessed_p_list,
  senders = c('B_cell','RPE','Neuron','Fibroblast','Pyramidal'),
  receptors = c('Mac_LipidMet','Mac_ProInf_1','Mac_ProInf_2','Mac_AntiInf_1','Mac_AntiInf_2','CD8_T_GZMK','CTL'),
  pval_threshold = 0.05
)


processed_means_list <- process_cpdb_means(
  unprocessed_list = unprocessed_means_list,
  senders = c('B_cell','RPE','Neuron','Fibroblast','Pyramidal'),
  receptors = c('Mac_LipidMet','Mac_ProInf_1','Mac_ProInf_2','Mac_AntiInf_1','Mac_AntiInf_2','CD8_T_GZMK','CTL')
)

common_interaction <- get_common_rownames(processed_p_list)

for (p in 1:length(processed_p_list)) {
  df <- processed_p_list[[p]]
  df <- df[rownames(df)%in%common_interaction,]
  processed_p_list[[p]] <- df
}

for (p in 1:length(processed_means_list)) {
  df <- processed_means_list[[p]]
  df <- df[rownames(df)%in%common_interaction,]
  processed_means_list[[p]] <- df
}

stopifnot(
  length(processed_p_list) == length(processed_means_list),
  all(sapply(1:length(processed_p_list), function(i) 
    identical(rownames(processed_p_list[[i]]), rownames(processed_means_list[[i]])) &&
      identical(colnames(processed_p_list[[i]]), colnames(processed_means_list[[i]])))))


######################################################-
##############   气泡折线图   ########################
######################################################-


combined_data <- map_dfr(1:length(processed_p_list), function(i) {
  # 转换p值数据框
  p_df <- processed_p_list[[i]] %>%
    rownames_to_column("Ligand_Receptor") %>%
    pivot_longer(-Ligand_Receptor, 
                 names_to = "Cell_Pair", 
                 values_to = "p_value")
  
  # 转换mean值数据框
  mean_df <- processed_means_list[[i]] %>%
    rownames_to_column("Ligand_Receptor") %>%
    pivot_longer(-Ligand_Receptor, 
                 names_to = "Cell_Pair", 
                 values_to = "mean_value")
  
  # 合并数据并添加组标识
  inner_join(p_df, mean_df, by = c("Ligand_Receptor", "Cell_Pair")) %>%
    mutate(
      Group = paste0("Group", i),
      Group_Factor = factor(i, levels = 1:length(processed_p_list)),
      log_p = -log10(p_value + 1e-10)  # 避免log(0)
    )
})
combined_data$Cell_Pair <- gsub('\\|',' -> ',combined_data$Cell_Pair)

#####   按照重要程度挑通讯   #######
# top_pathways <- combined_data %>%
#   group_by(Ligand_Receptor, Cell_Pair) %>%
#   summarize(
#     max_mean = max(mean_value),
#     min_p = min(p_value)
#   ) %>%
#   ungroup() %>%
#   filter(min_p < 0.05 & max_mean > quantile(max_mean, 0.75)) %>%
#   arrange(desc(max_mean)) %>%
#   slice_head(n = 50)  # 选择前20个最重要的通路
# 
# filtered_data <- combined_data %>%
#   semi_join(top_pathways, by = c("Ligand_Receptor", "Cell_Pair"))

########   按照关键词挑通讯   ######
keyword <- c('CCL','CXCL','TNF','TGF')
keypath <- paste0(keyword,collapse = '|')
filtered_data <- combined_data[grep(keypath,combined_data$Ligand_Receptor),]


ggplot(filtered_data, aes(x = Group_Factor, y = mean_value)) +
  geom_line(aes(group = interaction(Ligand_Receptor, Cell_Pair)), 
            color = "grey70", alpha = 0.6) +
  geom_point(aes(size = mean_value, color = log_p), alpha = 0.8) +
  geom_text(aes(label = ifelse(p_value < 0.05, "*", "")), 
            size = 5, vjust = 0.7, color = "black") +
  # 修改分面设置：去掉y轴刻度标签
  facet_grid(Ligand_Receptor ~ Cell_Pair, scales = "free_y", switch = "y") +
  
  # 调整视觉元素
  scale_size_continuous(
    name = "Strength",
    range = c(2, 8),
    breaks = pretty_breaks(4)
  ) +
  scale_color_gradientn(
    name = "-log10(p-value)",
    colours = c("navy",'white',"firebrick"),
    values = rescale(c(0,5,10))
  ) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  
  # 修改主题设置
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_blank(),  # 移除Y轴刻度标签
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_blank(),  # 移除Y轴标题
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9, face = "bold"),
    strip.text.x = element_text(size = 9, face = "bold"),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),  # 移除水平网格线
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "grey40"),
    strip.placement = "outside"
  ) +
  
  labs(
    x = "Experimental Group",
    title = "Cell-Cell Communication Analysis",
    subtitle = "Lines connect the same ligand-receptor pair across groups\nBubble size = Interaction strength | Color = Significance level"
  )


###################################################-
#############   分组气泡图   ######################
###################################################-
# 
# combined_data <- map_dfr(1:length(processed_p_list), function(i) {
#   # 转换p值数据框
#   p_df <- processed_p_list[[i]] %>%
#     rownames_to_column("Ligand_Receptor") %>%
#     pivot_longer(-Ligand_Receptor, 
#                  names_to = "Cell_Pair", 
#                  values_to = "p_value")
#   
#   # 转换mean值数据框
#   mean_df <- processed_means_list[[i]] %>%
#     rownames_to_column("Ligand_Receptor") %>%
#     pivot_longer(-Ligand_Receptor, 
#                  names_to = "Cell_Pair", 
#                  values_to = "mean_value")
#   
#   # 合并数据并添加组标识
#   inner_join(p_df, mean_df, by = c("Ligand_Receptor", "Cell_Pair")) %>%
#     mutate(
#       Group = paste0("Group", i),
#       log_p = -log10(p_value + 1e-10)  # 避免log(0)
#     )
# })
# 
# combined_data$Cell_Pair <- gsub('\\|',' -> ',combined_data$Cell_Pair)
#   
# ########   筛选需要的通路  #########
# ###   按means筛选
# filtered_data <- combined_data %>%
#   group_by(Ligand_Receptor, Cell_Pair) %>%
#   filter(mean(mean_value) > quantile(combined_data$mean_value, 0.9)) %>%
#   ungroup()
# ###   筛选特定关键词
# keyword <- c('CCL','CXCL','TNF','TGF')
# keypath <- paste0(keyword,collapse = '|')
# filtered_data <- combined_data[grep(keypath,combined_data$Ligand_Receptor),]
# 
# ggplot(filtered_data, aes(x = Cell_Pair, y = Ligand_Receptor)) +
#   geom_point(aes(
#     size = mean_value,
#     color = log_p,
#     alpha = ifelse(p_value < 0.05, 0.9, 0.4)  # 显著交互更明显
#   )) +
#   scale_size_continuous(
#     name = "Communication\nStrength (mean)",
#     range = c(1, 8),  # 调整气泡大小范围
#     breaks = seq(0, max(filtered_data$mean_value), length.out = 5)
#   ) +
#   scale_color_gradientn(
#     name = "-log10(p-value)",
#     colours = c("navy",'white',"firebrick"),
#     values = scales::rescale(c(0,5,10))
#   ) +
#   scale_alpha_identity() +  # 使用alpha表示显著性
#   facet_wrap(~ Group, ncol = 2) +  # 按组分面展示
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
#     axis.text.y = element_text(size = 8),
#     axis.title = element_blank(),
#     panel.grid.major = element_line(color = "grey90"),
#     legend.position = "right",
#     strip.text = element_text(face = "bold", size = 10)
#   ) +
#   labs(
#     title = "Cell-Cell Communication Analysis",
#     subtitle = "Bubble size: Interaction strength | Color: Significance level"
#   ) +
#   guides(
#     size = guide_legend(override.aes = list(color = "grey50")),
#     color = guide_colorbar(barwidth = 1, barheight = 10)
#   )


