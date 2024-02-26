# 0. Load packages --------------------------------------------------------
library(tidyverse)
library(readxl)
library(writexl)
library(RColorBrewer)
library(ggpubr)
library(ggpmisc)
library(ggsci)
library(gg.gap)
library(ggbreak)
library(ComplexHeatmap)
library(patchwork)

# 1. Set working directory ------------------------------------------------------------
setwd("D:/OneDrive - hust.edu.cn/【00-doing】项目/【09】菌菌互作关系/【R】15 - Paper")
dir.create("Figure2", showWarnings = F)
dir.create("Figure2/input", showWarnings = F)
setwd("D:/OneDrive - hust.edu.cn/【00-doing】项目/【09】菌菌互作关系/【R】15 - Paper/Figure2")

# 2. Load data ------------------------------------------------------------
taxa <- as.data.frame(read_xlsx("../Figure1/input/20221113_Total_clean_taxa.xlsx"))
taxa$Species <- gsub("\\[", "", taxa$Species)
taxa$Species <- gsub("\\]", "", taxa$Species)

interactions <- as.data.frame(read_xlsx("../Figure1/input/20221113_Total_clean_interactions.xlsx"))

# 3. Set colors & sign.comp -----------------------------------------------------------
colors.6 <- c("#22B7B2", "#ACE6E9", "#ECECED", "#FFE4C3", "#FFBB7F", "#FF9148")
names(colors.6) <- c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")
colors.3 <- c("#22B7B2", "#E5E5E5", "#FF9148")
names(colors.3) <- c("1", "0", "-1")

phylum.color <- brewer.pal(5, "Set2")
names(phylum.color) <- c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Verrucomicrobia")

phylum.color.1 <- data.frame(Phylum = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Verrucomicrobia"), color = brewer.pal(5, "Set2"))

sign.comp <- function(list) {
  pheno <- unique(list)
  result <- list()
  if (length(pheno) != 2) {
    for (i in 1:(length(pheno) - 1)) {
      for (j in (i + 1):length(pheno)) {
        comp <- list(c(pheno[i], pheno[j]))
        result <- append(result, comp)
      }
    }
  } else {
    result <- append(result, list(c(pheno[1], pheno[2])))
  }

  return(result)
}


# 4. Species from same or different taxa ---------------------------------------------
interactions.and.taxa <- interactions
interactions.and.taxa <- merge(interactions.and.taxa, taxa[, c(1:8)], by.x = "A", by.y = "Name", all.x = T)
interactions.and.taxa <- merge(interactions.and.taxa, taxa[, c(1:8)], by.x = "B", by.y = "Name", all.x = T)
colnames(interactions.and.taxa)[7:13] <- c("Kingdom.A", "Phylum.A", "Class.A", "Order.A", "Family.A", "Genus.A", "Species.A")
colnames(interactions.and.taxa)[14:20] <- c("Kingdom.B", "Phylum.B", "Class.B", "Order.B", "Family.B", "Genus.B", "Species.B")
interactions.and.taxa <- interactions.and.taxa[, c(2, 1, 3:ncol(interactions.and.taxa))]

# 分组
tax.list <- c("Phylum", "Class", "Order", "Family", "Genus")

for (level in tax.list) {
  col.a.name <- paste0(level, ".A")
  col.b.name <- paste0(level, ".B")
  interactions.and.taxa$flag <- ifelse(interactions.and.taxa[, col.a.name] == interactions.and.taxa[, col.b.name], "Same", "Different")
  colnames(interactions.and.taxa)[ncol(interactions.and.taxa)] <- level
}

rm(tax.list, col.a.name, col.b.name, level)

# 5. Figure 2A/B ---------------------------------------------------------
{ ## Figure2A ----

  dist <- read.table("input/phylogenetic_distance.txt", sep = "\t", header = F)
  # dist <- rbind(dist,dist[,c(2,1,3)])
  # dist <- dist[dist$V2>=dist$V1,] %>% unique()

  dist$Species.A <- apply(dist, 1, function(x) {
    paste(str_split(x[1], "_|\\.")[[1]][1], str_split(x[1], "_|\\.")[[1]][2])
  })
  dist$Species.B <- apply(dist, 1, function(x) {
    paste(str_split(x[2], "_|\\.")[[1]][1], str_split(x[2], "_|\\.")[[1]][2])
  })

  dist <- dist[, c(4, 5, 3)]
  colnames(dist)[3] <- "distance"

  interactions.1 <- interactions.and.taxa[, c("Species.A", "Species.B", "Interaction")]
  interactions.1 <- merge(interactions.1, dist, by = c("Species.A", "Species.B"))
  interactions.1$Interaction <- factor(interactions.1$Interaction, levels = c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition"))
  interactions.1$Interaction1[interactions.1$Interaction == "Mutualism" | interactions.1$Interaction == "Commensalism"] <- "Positive"
  interactions.1$Interaction1[interactions.1$Interaction == "Exploitation" | interactions.1$Interaction == "Amensalism" | interactions.1$Interaction == "Competition"] <- "Negative"
  interactions.1$Interaction1[interactions.1$Interaction == "Neutralism"] <- "Neutralism"

  write_xlsx(interactions.1, "Table S8.xlsx")

  plot.data <- interactions.1[, c(3, 4)]
  plot.data$type <- "six"
  plot.data.1 <- interactions.1[, c(5, 4)]
  plot.data.1$type <- "three"
  colnames(plot.data.1) <- colnames(plot.data)
  plot.data <- rbind(plot.data, plot.data.1)

  colors.3 <- c("#42B540B2", "#ED0000B2", "#ECECED")
  names(colors.3) <- c("Positive", "Negative", "Neutralism")

  colors.9 <- c(colors.3, colors.6)
  plot.data$Interaction <- factor(plot.data$Interaction, levels = c("Mutualism", "Commensalism", "Positive", "Neutralism", "Negative", "Amensalism", "Exploitation", "Competition"))

  figure2a <- plot.data %>% ggplot(aes(x = Interaction, y = distance, fill = Interaction)) +
    geom_boxplot(width = 0.3) +
    facet_grid(~type, space = "free", scales = "free") +
    scale_fill_manual(values = colors.9, guide = guide_legend(reverse = T)) +
    theme_bw() +
    guides(fill = "none") +
    geom_signif(
      # comparisons = list(c("Mutualism","Commensalism"),c("Commensalism","Neutralism"),c("Neutralism","Amensalism"),c("Amensalism","Exploitation"),c("Exploitation","Competition")),
      # comparisons = sign.comp(c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")),
      comparisons = list(
        c("Mutualism", "Competition"), c("Commensalism", "Neutralism"),
        c("Commensalism", "Amensalism"), c("Commensalism", "Exploitation"), c("Commensalism", "Competition"),
        c("Amensalism", "Exploitation"), c("Amensalism", "Competition")
      ),
      step_increase = 0.1,
      tip_length = 0.015,
      map_signif_level = T,
      vjust = -0.01
    ) +
    labs(x = "", y = "Genetic distance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  figure2a

  ggsave(figure2a, file = "Figure2A.pdf", width = 4, height = 5)
}

{ ## Figure 2B ----

  { ## 在不同的分类学水平分为same 和 different ----

    prepare.plot.data.1 <- function(data) {
      tax <- c("Phylum", "Class", "Order", "Family", "Genus")

      final.plot.data <- data.frame()
      final.plot.data.1 <- data.frame()
      for (level in tax) {
        col.a.name <- paste0(level, ".A")
        col.b.name <- paste0(level, ".B")
        data$flag <- ifelse(data[, col.a.name] == data[, col.b.name], "Same", "Different")

        plot.data <- as.data.frame(table(data[, c("Interaction", "flag")]))
        plot.data$percent[plot.data$flag == "Different"] <- plot.data$Freq[plot.data$flag == "Different"] / sum(plot.data$Freq[plot.data$flag == "Different"]) * 100
        plot.data$percent[plot.data$flag == "Same"] <- plot.data$Freq[plot.data$flag == "Same"] / sum(plot.data$Freq[plot.data$flag == "Same"]) * 100
        plot.data$level <- level
        final.plot.data <- rbind(final.plot.data, plot.data)

        temp <- data[, c(col.a.name, col.b.name, "Pheno.A", "Pheno.B", "Interaction", "flag")]
        colnames(temp)[1:2] <- c("A", "B")
        temp$level <- level
        final.plot.data.1 <- rbind(final.plot.data.1, temp)
      }
      final.plot.data$interaction <- factor(final.plot.data$Interaction, levels = rev(c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")))
      final.plot.data$level <- factor(final.plot.data$level, levels = tax)
      final.plot.data.1$interaction <- factor(final.plot.data.1$Interaction, levels = rev(c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")))
      final.plot.data.1$group[final.plot.data.1$Interaction == "Mutualism" | final.plot.data.1$Interaction == "Commensalism"] <- "Positive"
      final.plot.data.1$group[final.plot.data.1$Interaction == "Exploitation" | final.plot.data.1$Interaction == "Amensalism" | final.plot.data.1$Interaction == "Competition"] <- "Negative"
      final.plot.data.1$group[final.plot.data.1$Interaction == "Neutralism"] <- "Neutralism"
      final.plot.data$group[final.plot.data$Interaction == "Mutualism" | final.plot.data$Interaction == "Commensalism"] <- "Positive"
      final.plot.data$group[final.plot.data$Interaction == "Exploitation" | final.plot.data$Interaction == "Amensalism" | final.plot.data$Interaction == "Competition"] <- "Negative"
      final.plot.data$group[final.plot.data$Interaction == "Neutralism"] <- "Neutralism"

      return(list(stats = final.plot.data, total = final.plot.data.1))
    }

    taxa.same.diff <- prepare.plot.data.1(data = interactions.and.taxa)

    plot.data <- taxa.same.diff[["stats"]]
    plot.data$Interaction <- factor(plot.data$Interaction, levels = rev(c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")))
    plot.data$flag <- factor(plot.data$flag, levels = c("Same", "Different"))
    plot.data <- plot.data %>%
      group_by(level, flag) %>%
      mutate(sum = sum(Freq))
    ggplot(plot.data, aes(x = flag, y = percent, fill = Interaction)) +
      geom_bar(stat = "identity", position = "stack", alpha = 1, width = 0.4) +
      guides(fill = guide_legend(ncol = 1, reverse = TRUE)) +
      scale_fill_manual(name = "Interaction", values = colors.6) +
      scale_y_continuous(labels = scales::percent_format(scale = 1)) +
      labs(x = "", y = "", title = "") +
      facet_grid(~level) +
      theme_bw() +
      theme(
        text = element_text(size = 15, family = "sans"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, vjust = 1, size = 20),
        axis.text = element_text(size = 15),
        strip.background.x = element_rect(fill = "#FFFFFF")
      ) +
      geom_text(aes(label = sum, y = 100), size = 3, vjust = -0.5)
  }

  { ## 卡方 ----

    { ### 计算差异 ----
      tax <- c("Phylum", "Class", "Order", "Family", "Genus")
      chi.result <- data.frame()
      comp.result <- data.frame()
      ratio.result <- data.frame()
      for (level in tax) {
        data <- taxa.same.diff[["total"]][taxa.same.diff[["total"]]$level == level, ]
        # 总体
        mytable <- xtabs(~ flag + group, data = data)
        mytable
        p.0.1 <- chisq.test(mytable)[["p.value"]]

        # 总体
        mytable <- xtabs(~ flag + interaction, data = data)
        mytable
        p.0.2 <- chisq.test(mytable)[["p.value"]]

        # 合作
        data$group1[data$group == "Positive"] <- "Positive"
        data$group1[data$group != "Positive"] <- "None"
        mytable <- xtabs(~ flag + group1, data = data)
        p.1 <- chisq.test(mytable)[["p.value"]]

        # 竞争
        data$group2[data$group == "Negative"] <- "Negative"
        data$group2[data$group != "Negative"] <- "None"
        mytable <- xtabs(~ flag + group2, data = data)
        p.2 <- chisq.test(mytable)[["p.value"]]

        # 竞争
        data$group3[data$Interaction == "Competition"] <- "Competition"
        data$group3[data$Interaction != "Competition"] <- "None"
        mytable <- xtabs(~ flag + group3, data = data)
        p.3 <- chisq.test(mytable)[["p.value"]]
        # exploitation
        data$group4[data$Interaction == "Exploitation"] <- "Exploitation"
        data$group4[data$Interaction != "Exploitation"] <- "None"
        mytable <- xtabs(~ flag + group4, data = data)
        p.4 <- chisq.test(mytable)[["p.value"]]
        # Amensalism
        data$group5[data$Interaction == "Amensalism"] <- "Amensalism"
        data$group5[data$Interaction != "Amensalism"] <- "None"
        mytable <- xtabs(~ flag + group5, data = data)
        p.5 <- chisq.test(mytable)[["p.value"]]
        # Neutralism
        data$group6[data$Interaction == "Neutralism"] <- "Neutralism"
        data$group6[data$Interaction != "Neutralism"] <- "None"
        mytable <- xtabs(~ flag + group6, data = data)
        p.6 <- chisq.test(mytable)[["p.value"]]
        # Commensalism
        data$group7[data$Interaction == "Commensalism"] <- "Commensalism"
        data$group7[data$Interaction != "Commensalism"] <- "None"
        mytable <- xtabs(~ flag + group7, data = data)
        p.7 <- chisq.test(mytable)[["p.value"]]

        # Mutualism
        data$group8[data$Interaction == "Mutualism"] <- "Mutualism"
        data$group8[data$Interaction != "Mutualism"] <- "None"
        mytable <- xtabs(~ flag + group8, data = data)
        p.8 <- chisq.test(mytable)[["p.value"]]

        # comp signal
        six.comp <- as.data.frame(table(data[, c("flag", "interaction")]))
        six.comp <- spread(six.comp, flag, Freq)
        six.comp$Different <- six.comp$Different / sum(six.comp$Different)
        six.comp$Same <- six.comp$Same / sum(six.comp$Same)
        six.comp$flag <- ifelse(six.comp$Different > six.comp$Same, 1, -1)

        three.comp <- as.data.frame(table(data[, c("flag", "group")]))
        three.comp <- spread(three.comp, flag, Freq)
        three.comp$Different <- three.comp$Different / sum(three.comp$Different)
        three.comp$Same <- three.comp$Same / sum(three.comp$Same)
        three.comp$flag <- ifelse(three.comp$Different > three.comp$Same, 1, -1)
        colnames(three.comp) <- colnames(six.comp)

        ratio.result <- rbind(ratio.result, six.comp %>% mutate(Level = level))
        ratio.result <- rbind(ratio.result, three.comp %>% mutate(Level = level))

        chi.result <- rbind(chi.result, c(level, p.0.1, p.0.2, p.1, p.2, p.3, p.4, p.5, p.6, p.7, p.8))
        comp.result <- rbind(comp.result, c(level, three.comp$flag[c(2, 1)], six.comp$flag))
      }
      colnames(chi.result) <- c("Level", "Three types", "Six types", "Positive", "Negative", "Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism")
      colnames(comp.result) <- c("Level", "Positive", "Negative", "Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism")
      rownames(comp.result) <- comp.result$Level
      comp.result <- comp.result[, -1]
      for (i in 1:ncol(comp.result)) {
        comp.result[, i] <- as.numeric(comp.result[, i])
      }

      for (i in 2:ncol(chi.result)) {
        chi.result[, i] <- as.numeric(chi.result[, i])
      }
      chi.result.signif <- chi.result
      for (j in 2:ncol(chi.result.signif)) {
        for (i in 1:nrow(chi.result.signif)) {
          num <- as.numeric(chi.result.signif[i, j])
          chi.result.signif[i, j] <- ifelse(num < 0.05 & num > 0.01, "*",
            ifelse(num > 0.05 & num < 0.1, ".",
              ifelse(num < 0.01, "**", "")
            )
          )
        }
      }

      rownames(chi.result) <- chi.result$Level
      chi.result <- chi.result[, -1]
      rownames(chi.result.signif) <- chi.result.signif$Level
      chi.result.signif <- chi.result.signif[, -1]
    }

    { ### 折线图 ----

      plot.data <- ratio.result
      plot.data$interaction <- as.character(plot.data$interaction)
      plot.data <- gather(plot.data[, -4], Group, Ratio, -c(interaction, Level))
      plot.data$Group <- factor(plot.data$Group, levels = c("Same", "Different"))
      plot.data$interaction <- factor(plot.data$interaction, levels = c("Negative", "Positive", "Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism"))
      plot.data$Level <- factor(plot.data$Level, levels = c("Phylum", "Class", "Order", "Family", "Genus"))
      plot.data$Ratio <- plot.data$Ratio * 100

      plot.data.1 <- plot.data[plot.data$interaction %in% c("Negative", "Positive", "Neutralism"), ]

      df1 <- chi.result.signif[, c("Negative", "Positive", "Neutralism")]
      df1$Level <- rownames(df1)
      df1 <- gather(df1, interaction, sig, -Level)
      df1 <- df1[df1$sig != "", ]
      df1$x <- 1.5
      df1$y <- apply(df1, 1, function(x) {
        num <- mean(plot.data.1$Ratio[plot.data.1$interaction == x[2] & plot.data.1$Level == x[1]])
      })

      df1$Level <- factor(df1$Level, levels = c("Phylum", "Class", "Order", "Family", "Genus"))
      colors.6.1 <- colors.6
      colors.6.1[["Neutralism"]] <- "grey"
      colors.6.1[["Negative"]] <- "#ED0000B2"
      colors.6.1[["Positive"]] <- "#42B540B2"
      colors.6.1 <- colors.6.1[c("Negative", "Positive", "Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism")]

      plot.data.1 %>% ggplot(aes(x = Group, y = Ratio)) +
        geom_point(aes(color = interaction)) +
        geom_line(aes(group = interaction, color = interaction)) +
        facet_grid(~Level) +
        geom_text(data = df1, aes(x = x, y = y, label = sig)) +
        scale_color_manual(name = "Intearction", values = colors.6.1, labels = c("Negative", "Positive", "Neutralism")) +
        theme_bw() +
        # guides(color = guide_legend(ncol = 1, reverse = TRUE))
        theme(
          text = element_text(size = 15, family = "sans"),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, vjust = 1, size = 20),
          axis.text = element_text(size = 15),
          strip.background.x = element_rect(fill = "#FFFFFF"),
        ) +
        labs(x = "", y = "") +
        scale_y_break(c(25, 60), # 截断位置及范围
          space = 0.3, # 间距大小
          scales = 1, # 上下显示比例，大于1上面比例大，小于1下面比例大
          expand = c(0, 0)
        ) +
        scale_y_continuous(labels = scales::percent_format(scale = 1, accuracy = 1)) +
        theme(
          axis.text.y.right = element_blank(),
          axis.title.y.right = element_blank(),
          axis.ticks.y.right = element_blank()
        )
      # scale_y_continuous()
    }
  }

  { ## Figure 2B ----

    { ### 折线图 ----

      plot.data <- ratio.result
      plot.data$interaction <- as.character(plot.data$interaction)
      plot.data <- gather(plot.data[, -4], Group, Ratio, -c(interaction, Level))
      plot.data$Group <- factor(plot.data$Group, levels = c("Same", "Different"))
      plot.data$interaction <- factor(plot.data$interaction, levels = c("Negative", "Positive", "Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism"))
      plot.data$Level <- factor(plot.data$Level, levels = c("Phylum", "Class", "Order", "Family", "Genus"))
      plot.data$Ratio <- plot.data$Ratio * 100

      plot.data.1 <- plot.data[plot.data$interaction %in% c("Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism"), ]

      df1 <- chi.result.signif[, c("Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism")]
      df1$Level <- rownames(df1)
      df1 <- gather(df1, interaction, sig, -Level)
      df1 <- df1[df1$sig != "", ]
      df1$x <- 1.5
      df1$y <- apply(df1, 1, function(x) {
        num <- mean(plot.data.1$Ratio[plot.data.1$interaction == x[2] & plot.data.1$Level == x[1]])
      })

      df1$Level <- factor(df1$Level, levels = c("Phylum", "Class", "Order", "Family", "Genus"))
      colors.6.1 <- colors.6
      colors.6.1[["Neutralism"]] <- "grey"
      colors.6.1[["Negative"]] <- "#ED0000B2"
      colors.6.1[["Positive"]] <- "#42B540B2"
      colors.6.1 <- colors.6.1[c("Negative", "Positive", "Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism")]

      figure2b <- plot.data.1 %>% ggplot(aes(x = Group, y = Ratio)) +
        geom_point(aes(color = interaction), size = 3.5) +
        geom_line(aes(group = interaction, color = interaction)) +
        facet_grid(~Level) +
        geom_text(data = df1, aes(x = x, y = y, label = sig)) +
        scale_color_manual(name = "Intearction", values = colors.6.1, labels = c("Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism")) +
        theme_bw() +
        # guides(color = guide_legend(ncol = 1, reverse = TRUE))
        theme(
          text = element_text(size = 15, family = "sans"),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, vjust = 1, size = 20),
          axis.text = element_text(size = 15),
          strip.background.x = element_rect(fill = "#FFFFFF"),
        ) +
        labs(x = "", y = "") +
        # scale_y_break(c(25,60),#截断位置及范围
        #               space = 0.3,#间距大小
        #               scales = 1,#上下显示比例，大于1上面比例大，小于1下面比例大
        #               expand = c(0,0))+
        scale_y_continuous(labels = scales::percent_format(scale = 1, accuracy = 1)) +
        theme(
          axis.text.y.right = element_blank(),
          axis.title.y.right = element_blank(),
          axis.ticks.y.right = element_blank()
        )
      # scale_y_continuous()
      figure2b
      ggsave(figure2b, file = paste0("Figure2B.pdf"), width = 11, height = 4)
    }
  }
}

{ ## Figure S5 ----

  inter.6.compute <- function(data, level) {
    { # 统计每个类内部的关系 ----
      same.data <- data[data[, level] == "Same", c(paste0(level, ".A"), "Interaction")]
      colnames(same.data)[1] <- "Taxa"
      same.table <- same.data %>%
        group_by(Taxa, Interaction) %>%
        summarise(Freq = n())
      same.table <- spread(same.table, Interaction, Freq)
      same.table[is.na(same.table)] <- 0
      same.table <- gather(same.table, Interaction, Freq, -Taxa)
      same.table <- same.table %>%
        group_by(Taxa) %>%
        mutate(Sum = sum(Freq))
      same.table$Ratio <- same.table$Freq / same.table$Sum * 100
      same.table$Interaction <- factor(same.table$Interaction, levels = rev(c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")))
    }

    { # 统计每个类跟其他类的关系----
      diff.data <- data[data[, level] == "Different", c(paste0(level, ".A"), paste0(level, ".B"), "Interaction")]
      colnames(diff.data)[1:2] <- c("Taxa.A", "Taxa.B")
      diff.table <- diff.data %>%
        group_by(Taxa.A, Taxa.B, Interaction) %>%
        summarise(Freq = n())
      diff.table <- spread(diff.table, Interaction, Freq)
      diff.table[is.na(diff.table)] <- 0
      diff.table <- gather(diff.table, Interaction, Freq, -c(Taxa.A, Taxa.B))

      diff.table.1 <- diff.table[, c(2, 1, 3, 4)]
      colnames(diff.table.1) <- colnames(diff.table)

      diff.table <- rbind(diff.table, diff.table.1)
      diff.table <- aggregate(diff.table$Freq, by = list(diff.table$Taxa.A, diff.table$Taxa.B, diff.table$Interaction), FUN = sum)
      diff.table <- rbind(diff.table, diff.table[, c(2, 1, 3, 4)]) %>%
        subset(Group.2 >= Group.1) %>%
        unique()
      colnames(diff.table) <- colnames(diff.table.1)
      diff.table <- diff.table %>%
        group_by(Taxa.A, Taxa.B) %>%
        mutate(Sum = sum(Freq))
      diff.table$Ratio <- diff.table$Freq / diff.table$Sum * 100
      diff.table$Interaction <- factor(diff.table$Interaction, levels = rev(c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")))
    }
    { # 不分内外统计总体 ----
      total.data <- data[, c(paste0(level, ".A"), paste0(level, ".B"), "Interaction")]
      colnames(total.data)[1:2] <- c("Taxa.A", "Taxa.B")
      total.table <- total.data %>%
        group_by(Taxa.A, Taxa.B, Interaction) %>%
        summarise(Freq = n())
      total.table <- spread(total.table, Interaction, Freq)
      total.table[is.na(total.table)] <- 0
      total.table <- gather(total.table, Interaction, Freq, -c(Taxa.A, Taxa.B))

      total.table.1 <- total.table[, c(2, 1, 3, 4)]
      colnames(total.table.1) <- colnames(total.table)

      total.table <- rbind(total.table, total.table.1)
      total.table <- aggregate(total.table$Freq, by = list(total.table$Taxa.A, total.table$Taxa.B, total.table$Interaction), FUN = sum)
      total.table <- rbind(total.table, total.table[, c(2, 1, 3, 4)]) %>%
        subset(Group.2 >= Group.1) %>%
        unique()
      colnames(total.table) <- colnames(total.table.1)
      total.table <- total.table %>%
        group_by(Taxa.A, Taxa.B) %>%
        mutate(Sum = sum(Freq))
      total.table$Ratio <- total.table$Freq / total.table$Sum * 100
      total.table$Interaction <- factor(total.table$Interaction, levels = rev(c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")))
    }
    { # 合并内部外部数据 ----

      bind.table <- data.frame()
      for (i in unique(c(data[, paste0(level, ".A")], data[, paste0(level, ".B")]))) {
        # for(i in unique(c(data$Class.A,data$Class.B))){

        same.temp <- same.table[same.table$Taxa == i, ]
        same.temp$Group <- "Inner"

        diff.temp1 <- diff.table[diff.table$Taxa.A == i, ]
        colnames(diff.temp1)[1:2] <- c("Taxa", "Group")
        diff.temp2 <- diff.table[diff.table$Taxa.B == i, ]
        colnames(diff.temp2)[1:2] <- c("Group", "Taxa")
        diff.temp <- rbind(diff.temp1, diff.temp2)
        diff.temp3 <- aggregate(diff.temp$Freq, list(diff.temp$Taxa, diff.temp$Interaction), sum)
        colnames(diff.temp3) <- c("Taxa", "Interaction", "Freq")
        diff.temp3$Sum <- sum(diff.temp3$Freq)
        diff.temp3$Ratio <- diff.temp3$Freq / diff.temp3$Sum * 100
        diff.temp3$Group <- "Outer"

        total.temp1 <- total.table[total.table$Taxa.A == i, ]
        colnames(total.temp1)[1:2] <- c("Taxa", "Group")
        total.temp2 <- total.table[total.table$Taxa.B == i, ]
        colnames(total.temp2)[1:2] <- c("Group", "Taxa")
        total.temp <- rbind(total.temp1, total.temp2)
        total.temp <- unique(total.temp)
        total.temp3 <- aggregate(total.temp$Freq, list(total.temp$Taxa, total.temp$Interaction), sum)
        colnames(total.temp3) <- c("Taxa", "Interaction", "Freq")
        total.temp3$Sum <- sum(total.temp3$Freq)
        total.temp3$Ratio <- total.temp3$Freq / total.temp3$Sum * 100
        total.temp3$Group <- "Total"

        bind.table <- rbind(bind.table, same.temp, diff.temp, diff.temp3, total.temp, total.temp3)
      }
    }

    { # 按照某种关系进行排序 ----

      order <- bind.table[bind.table$Interaction == "Competition" | bind.table$Interaction == "Exploitation" | bind.table$Interaction == "Amensalism", ]
      order <- order[order$Group == "Total", ]
      order <- order %>%
        group_by(Taxa) %>%
        summarise(sum = sum(Ratio)) %>%
        arrange(desc(sum))
      bind.table$Taxa <- factor(bind.table$Taxa, levels = order$Taxa)
    }

    { # 合并门水平信息 ----

      bind.table <- merge(bind.table, taxa[, c("Phylum", level)], by.x = "Taxa", by.y = level)
      bind.table <- unique(bind.table)
      colnames(bind.table)[ncol(bind.table)] <- "Phylum"
    }
    return(bind.table)
  }

  phylum.inter.6 <- inter.6.compute(data = interactions.and.taxa, level = "Phylum")

  { ## 内外部一起画图(根据内外部分面) ----
    # bind.table$Taxa <- factor(bind.table$Taxa,levels = c( "Verrucomicrobia","Actinobacteria","Bacteroidetes","Firmicutes",   "Proteobacteria"))
    plot.data <- phylum.inter.6
    text.col <- unique(plot.data[, c("Taxa", "Phylum")])
    text.col <- merge(text.col, phylum.color.1, by = "Phylum", all.x = T)
    rownames(text.col) <- text.col$Taxa
    text.col <- text.col[levels(plot.data$Taxa), ]

    plot.data$Taxa <- factor(plot.data$Taxa, levels = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Verrucomicrobia"))
    plot.data.1 <- plot.data %>% subset(Group == "Outer" | Group == "Inner")
    plot.data.1 <- plot.data.1[, c(1, 2, 3, 6)]
    plot.data.1 <- spread(plot.data.1, Group, Freq, fill = 0)
    plot.data.1 <- gather(plot.data.1, Group, Freq, -c(Taxa, Interaction))
    plot.data.1 <- plot.data.1 %>%
      group_by(Taxa, Group) %>%
      mutate(Sum = sum(Freq))
    plot.data.1$Ratio <- plot.data.1$Freq / plot.data.1$Sum * 100
    plot.data.1 %>%
      # subset(Taxa != "Verrucomicrobia") %>%
      ggplot(aes(
        x = Group, # 将第一列转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
        y = Ratio, # 判断分组情况，将两个柱子画在0的两侧
        fill = Interaction
      )) +
      facet_grid(. ~ Taxa, scales = "free_x", space = "free_x") +
      guides(fill = guide_legend(ncol = 1, reverse = F)) +
      geom_bar(stat = "identity", width = 0.4) + # 画柱形图
      scale_fill_manual(name = "Interaction", values = colors.6) +
      scale_y_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0, 0.05))) + # 在y轴的两侧，留下一部分的空白位置，防止加标签的时候，显示不全
      theme_bw() +
      geom_hline(yintercept = 0, colour = "black") +
      labs(x = "", y = "") +
      geom_text(aes(label = Sum, y = 100), size = 3, vjust = -0.5) +
      theme(strip.background.x = element_rect(fill = "#FFFFFF")) +
      theme(
        text = element_text(size = 15, family = "serif"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, vjust = 1, size = 20),
        axis.text = element_text(size = 15)
      )
  }

  { # 卡方 ----

    # tax <- c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria")
    tax <- unique(plot.data$Taxa[plot.data$Group == "Inner"])
    chi.result <- data.frame()
    for (level in tax) {
      data <- plot.data[plot.data$Taxa == level, ]
      mytable <- spread(data[, c(2, 3, 6)], Interaction, Freq)
      rownames(mytable) <- mytable$Group
      mytable <- mytable[-3, -1]

      # 总体
      mytable.1 <- mytable
      mytable.1$Positive <- mytable.1$Mutualism + mytable.1$Commensalism
      mytable.1$Negative <- mytable.1$Competition + mytable.1$Exploitation + mytable.1$Amensalism
      mytable.1 <- mytable.1[, c("Positive", "Negative", "Neutralism")]
      p.0.1 <- chisq.test(mytable.1)[["p.value"]]

      # 总体
      p.0.2 <- chisq.test(mytable)[["p.value"]]

      # 合作
      mytable.1 <- mytable
      mytable.1$Positive <- mytable.1$Mutualism + mytable.1$Commensalism
      mytable.1$No <- mytable.1$Competition + mytable.1$Exploitation + mytable.1$Amensalism + mytable.1$Neutralism
      mytable.1 <- mytable.1[, c("Positive", "No")]
      p.1 <- chisq.test(mytable.1)[["p.value"]]

      # 竞争
      mytable.1 <- mytable
      mytable.1$Negative <- mytable.1$Competition + mytable.1$Exploitation + mytable.1$Amensalism
      mytable.1$No <- mytable.1$Mutualism + mytable.1$Commensalism + mytable.1$Neutralism
      mytable.1 <- mytable.1[, c("Negative", "No")]
      p.2 <- chisq.test(mytable.1)[["p.value"]]

      # 竞争
      mytable.1 <- mytable
      mytable.1$No <- rowSums(mytable.1[, -1])
      mytable.1 <- mytable.1[, c("Competition", "No")]
      p.3 <- chisq.test(mytable.1)[["p.value"]]

      # Exploitation
      mytable.1 <- mytable
      mytable.1$No <- rowSums(mytable.1[, -2])
      mytable.1 <- mytable.1[, c("Exploitation", "No")]
      p.4 <- chisq.test(mytable.1)[["p.value"]]
      # Amensalism
      mytable.1 <- mytable
      mytable.1$No <- rowSums(mytable.1[, -3])
      mytable.1 <- mytable.1[, c("Amensalism", "No")]
      p.5 <- chisq.test(mytable.1)[["p.value"]]
      # Neutralism
      mytable.1 <- mytable
      mytable.1$No <- rowSums(mytable.1[, -4])
      mytable.1 <- mytable.1[, c("Neutralism", "No")]
      p.6 <- chisq.test(mytable.1)[["p.value"]]
      # Commensalism
      mytable.1 <- mytable
      mytable.1$No <- rowSums(mytable.1[, -5])
      mytable.1 <- mytable.1[, c("Commensalism", "No")]
      p.7 <- chisq.test(mytable.1)[["p.value"]]

      # Mutualism
      mytable.1 <- mytable
      mytable.1$No <- rowSums(mytable.1[, -6])
      mytable.1 <- mytable.1[, c("Mutualism", "No")]
      p.8 <- chisq.test(mytable.1)[["p.value"]]

      chi.result <- rbind(chi.result, c(level, p.0.1, p.0.2, p.1, p.2, p.3, p.4, p.5, p.6, p.7, p.8))
    }
    colnames(chi.result) <- c("Level", "Three types", "Six types", "Positive", "Negative", "Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism")

    for (i in 2:ncol(chi.result)) {
      chi.result[, i] <- as.numeric(chi.result[, i])
    }

    chi.result.signif <- chi.result
    for (j in 2:ncol(chi.result.signif)) {
      for (i in 1:nrow(chi.result.signif)) {
        num <- as.numeric(chi.result.signif[i, j])
        chi.result.signif[i, j] <- ifelse(num < 0.05 & num > 0.01, "*",
          ifelse(num > 0.05 & num < 0.1, ".",
            ifelse(num < 0.01, "**", "")
          )
        )
      }
    }

    rownames(chi.result) <- chi.result$Level
    chi.result <- chi.result[, -1]
    rownames(chi.result.signif) <- chi.result.signif$Level
    chi.result.signif <- chi.result.signif[, -1]


    chi.result.2 <- chi.result[, -c(1, 2)]

    comp.result <- chi.result.2
    for (i in 1:nrow(comp.result)) {
      for (j in unique(plot.data.1$Interaction)) {
        temp <- plot.data.1[plot.data.1$Taxa == rownames(comp.result)[i] & plot.data.1$Interaction == j, ]
        if (temp$Ratio[temp$Group == "Inner"] >= temp$Ratio[temp$Group == "Outer"]) {
          comp.result[rownames(comp.result)[i], j] <- -1
        } else {
          comp.result[rownames(comp.result)[i], j] <- 1
        }
      }
    }
    plot.data <- chi.result.2 * comp.result
    plot.data$Level <- rownames(plot.data)
    plot.data <- gather(plot.data, type, p.value, -Level)
    plot.data$p.abs <- abs(plot.data$p.value)
    plot.data$signal <- ifelse(plot.data$p.value > 0, "Down", "Up")
    plot.data$signif <- ifelse(plot.data$p.abs < 0.05, "Signif.", "Not")

    plot.data$Level <- factor(plot.data$Level, levels = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Verrucomicrobia"))
    plot.data$type <- factor(plot.data$type, levels = rev(c("Positive", "Negative", "Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism")))

    # 折线图
    plot.data.2 <- plot.data.1[, c(1, 2, 3, 4)]
    plot.data.2$Interaction <- as.character(plot.data.2$Interaction)
    plot.data.2$Interaction[plot.data.2$Interaction == "Competition" | plot.data.2$Interaction == "Exploitation" | plot.data.2$Interaction == "Amensalism"] <- "Negative"
    plot.data.2$Interaction[plot.data.2$Interaction == "Commensalism" | plot.data.2$Interaction == "Mutualism"] <- "Positive"
    plot.data.2 <- aggregate(plot.data.2$Freq, by = list(plot.data.2$Taxa, plot.data.2$Interaction, plot.data.2$Group), sum)
    colnames(plot.data.2) <- c("Taxa", "Interaction", "Group", "Freq")
    plot.data.2 <- plot.data.2 %>%
      group_by(Taxa, Group) %>%
      mutate(Sum = sum(Freq))
    plot.data.2$Ratio <- plot.data.2$Freq / plot.data.2$Sum * 100

    plot.data.3 <- rbind(plot.data.1, plot.data.2[plot.data.2$Interaction != "Neutralism", ])
    plot.data.3$Interaction <- factor(plot.data.3$Interaction, levels = c("Negative", "Positive", "Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism"))
    colors.6.1 <- colors.6
    colors.6.1[["Neutralism"]] <- "grey"
    colors.6.1[["Negative"]] <- "#ED0000B2"
    colors.6.1[["Positive"]] <- "#42B540B2"

    colors.6.1 <- colors.6.1[c("Negative", "Positive", "Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism")]

    plot.data.4 <- plot.data.3[plot.data.3$Interaction %in% c("Negative", "Positive", "Neutralism"), ]

    df1 <- chi.result.signif[, c("Negative", "Positive", "Neutralism")]
    df1$Taxa <- rownames(df1)
    df1 <- gather(df1, interaction, sig, -Taxa)
    df1 <- df1[df1$sig != "", ]
    df1$x <- 1.5
    df1$y <- apply(df1, 1, function(x) {
      num <- mean(plot.data.4$Ratio[plot.data.4$Interaction == x[2] & plot.data.4$Taxa == x[1]])
    })

    df1$Taxa <- factor(df1$Taxa, levels = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Verrucomicrobia"))

    plot.data.4 %>% ggplot(aes(x = Group, y = Ratio)) +
      geom_point(aes(color = Interaction)) +
      geom_line(aes(group = Interaction, color = Interaction)) +
      facet_grid(~Taxa) +
      geom_text(data = df1, aes(x = x, y = y, label = sig)) +
      scale_color_manual(name = "Intearction", values = colors.6.1, labels = c("Negative", "Positive", "Neutralism")) +
      theme_bw() +
      theme(
        text = element_text(size = 15, family = "serif"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, vjust = 1, size = 20),
        axis.text = element_text(size = 15),
        strip.background.x = element_rect(fill = "#FFFFFF"),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
      ) +
      labs(x = "", y = "") +
      scale_y_break(c(30, 60), # 截断位置及范围
        space = 0.3, # 间距大小
        scales = 1, # 上下显示比例，大于1上面比例大，小于1下面比例大
        expand = c(0, 0)
      ) +
      scale_y_continuous(labels = scales::percent_format(scale = 1, accuracy = 1))
  }

  { ## Figure S5 ----

    plot.data.5 <- plot.data.3[!(plot.data.3$Interaction %in% c("Negative", "Positive")), ]

    df1 <- chi.result.signif[, c("Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism")]
    df1$Taxa <- rownames(df1)
    df1 <- gather(df1, interaction, sig, -Taxa)
    df1 <- df1[df1$sig != "", ]
    df1$x <- 1.5
    df1$y <- apply(df1, 1, function(x) {
      num <- mean(plot.data.4$Ratio[plot.data.4$Interaction == x[2] & plot.data.4$Taxa == x[1]])
    })
    df1$Taxa <- factor(df1$Taxa, levels = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Verrucomicrobia"))

    figures5 <- ggplot(plot.data.5, aes(x = Group, y = Ratio)) +
      geom_point(aes(color = Interaction), size = 3.5) +
      geom_line(aes(group = Interaction, color = Interaction)) +
      facet_grid(. ~ Taxa) +
      geom_text(data = df1, aes(x = x, y = y, label = sig)) +
      scale_color_manual(name = "Intearction", values = colors.6.1, labels = c("Competition", "Exploitation", "Amensalism", "Neutralism", "Commensalism", "Mutualism")) +
      theme_bw() +
      scale_y_continuous(labels = scales::percent_format(scale = 1)) +
      # guides(color = guide_legend(ncol = 1, reverse = TRUE))
      theme(
        text = element_text(size = 15, family = "serif"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, vjust = 1, size = 20),
        axis.text = element_text(size = 15),
        strip.background.x = element_rect(fill = "#FFFFFF"),
        # strip.text.x = element_text(color = "#FFFFFF")
      ) +
      labs(x = "", y = "")
    figures5
    ggsave(figures5, file = paste0("Supp5.pdf"), width = 11, height = 4)
  }
}

# 6. Figure 2C ------------------------------------------------------------
{ ## Load data ----

  load("../Figure1/input/pair.combined.RData")
  load("../Figure1/input/single.cont.RData")
  load("../Figure1/input/single.dis.RData")
  load("../Figure1/inputimp.index.RData")
  single.cont$tmRNA[is.na(single.cont$tmRNA)] <- 0
  single.cont$SM.total[is.na(single.cont$SM.total)] <- 0

  load("../Figure1/input/species.inter.3.RData")
  species.inter.3[["outgoing"]]$Species <- gsub("\\[", "", species.inter.3[["outgoing"]]$Species)
  species.inter.3[["outgoing"]]$Species <- gsub("\\]", "", species.inter.3[["outgoing"]]$Species)
  species.inter.3[["incoming"]]$Species <- gsub("\\[", "", species.inter.3[["incoming"]]$Species)
  species.inter.3[["incoming"]]$Species <- gsub("\\]", "", species.inter.3[["incoming"]]$Species)
  incoming <- species.inter.3[["incoming"]]
  outgoing <- species.inter.3[["outgoing"]]
  rm(species.inter.3)

  incoming <- incoming[, c(12, 2, 3)]
  incoming <- incoming %>%
    group_by(Species, Group) %>%
    summarise(Freq = sum(Freq))
  incoming <- incoming %>%
    group_by(Species) %>%
    mutate(Sum = sum(Freq))
  incoming$Ratio <- incoming$Freq / incoming$Sum * 100
  incoming <- spread(incoming[, -c(3, 4)], Group, Ratio) %>% as.data.frame()
  rownames(incoming) <- incoming$Species
  incoming <- incoming[, -1]
  colnames(incoming) <- c("Incoming_Harmful", "Incoming_Neutral", "Incoming_Beneficial")
  incoming <- incoming[rownames(incoming) != "Klebsiella pneumoniae(sp)", ]

  outgoing <- outgoing[, c(12, 2, 3)]
  outgoing <- outgoing %>%
    group_by(Species, Group) %>%
    summarise(Freq = sum(Freq))
  outgoing <- outgoing %>%
    group_by(Species) %>%
    mutate(Sum = sum(Freq))
  outgoing$Ratio <- outgoing$Freq / outgoing$Sum * 100
  outgoing <- spread(outgoing[, -c(3, 4)], Group, Ratio) %>% as.data.frame()
  rownames(outgoing) <- outgoing$Species
  outgoing <- outgoing[, -1]
  colnames(outgoing) <- c("Outgoing_Harmful", "Outgoing_Neutral", "Outgoing_Beneficial")
  outgoing <- outgoing[rownames(outgoing) != "Klebsiella pneumoniae(sp)", ]

  total.3 <- merge(incoming, outgoing, by = "row.names")
  rownames(total.3) <- total.3$Row.names
  total.3 <- total.3[, -1]
  total.3 <- total.3[rownames(total.3) != "Klebsiella pneumoniae(sp)", ]

  inter.6to3 <- data.frame(
    Interaction = c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition"),
    Interaction.1 = c("Positive", "Positive", "Neutral", "Negative", "Negative", "Negative")
  )

  interactions <- as.data.frame(read_xlsx("../Figure1/input/20221113_Total_clean_interactions.xlsx"))

  interactions.and.taxa <- interactions
  interactions.and.taxa <- merge(interactions.and.taxa, taxa[, c(1:8)], by.x = "A", by.y = "Name", all.x = T)
  interactions.and.taxa <- merge(interactions.and.taxa, taxa[, c(1:8)], by.x = "B", by.y = "Name", all.x = T)
  colnames(interactions.and.taxa)[7:13] <- c("Kingdom.A", "Phylum.A", "Class.A", "Order.A", "Family.A", "Genus.A", "Species.A")
  colnames(interactions.and.taxa)[14:20] <- c("Kingdom.B", "Phylum.B", "Class.B", "Order.B", "Family.B", "Genus.B", "Species.B")
  interactions.and.taxa <- interactions.and.taxa[, c(2, 1, 3:ncol(interactions.and.taxa))]

  # 分组
  tax.list <- c("Phylum", "Class", "Order", "Family", "Genus")

  for (level in tax.list) {
    col.a.name <- paste0(level, ".A")
    col.b.name <- paste0(level, ".B")
    interactions.and.taxa$flag <- ifelse(interactions.and.taxa[, col.a.name] == interactions.and.taxa[, col.b.name], "Same", "Different")
    colnames(interactions.and.taxa)[ncol(interactions.and.taxa)] <- level
  }

  rm(tax.list, col.a.name, col.b.name, level)
}

{ ## dbRDA/envfit ----

  pheno <- merge(single.cont, single.dis, by = "Species")
  pheno <- pheno[, c("Species", imp.index)]
  colnames(pheno)[str_detect(colnames(pheno), "VFC")] <- gsub(" \\(", "_", colnames(pheno)[str_detect(colnames(pheno), "VFC")])
  colnames(pheno)[str_detect(colnames(pheno), "VFC")] <- gsub("\\)", "", colnames(pheno)[str_detect(colnames(pheno), "VFC")])
  colnames(pheno)[which(colnames(pheno) == "Antimicrobial activity/Competitive advantage_VFC0325")] <- "Antimicrobial_activity_Competitive_advantage_VFC0325"
  colnames(pheno)[str_detect(colnames(pheno), "VFC")] <- gsub(" ", "_", colnames(pheno)[str_detect(colnames(pheno), "VFC")])
  colnames(pheno)[which(colnames(pheno) == "Nutritional/Metabolic_factor_VFC0272")] <- "Nutritional_Metabolic_factor_VFC0272"
  colnames(pheno)[str_detect(colnames(pheno), "VFC")] <- gsub("-", "_", colnames(pheno)[str_detect(colnames(pheno), "VFC")])
  pheno[, -1] <- apply(pheno[, -1], 2, as.numeric)
  rownames(pheno) <- pheno$Species
  pheno <- pheno[, -1]
  pheno <- pheno[, -c(which(colnames(pheno) == "Gram.negative"))]


  { ## incoming vif 删除共现性指标 ----

    data <- decostand(t(incoming), "hellinger")
    env <- pheno

    sel <- decorana(data)
    sel
    temp1 <- rda(t(data) ~ ., env, scale = FALSE)
    vif.cca(temp1) %>%
      as.data.frame() %>%
      View()
    temp1 <- rda(t(data) ~ ., env[, -c(4, 2, 3, 13, 19, 38, 39, 14, 24, 8, 10, 11)])
    vif.cca(temp1) %>%
      as.data.frame() %>%
      View()
  }

  { ## outgoing vif 删除共现性指标 ----

    data <- decostand(t(outgoing), "hellinger")
    env <- pheno

    sel <- decorana(data)
    sel
    temp1 <- rda(t(data) ~ ., env)
    vif.cca(temp1) %>%
      as.data.frame() %>%
      View()
    temp1 <- rda(t(data) ~ ., env[, -c(4, 2, 3, 13, 19, 38, 39, 14, 24, 8, 10, 11)])
    vif.cca(temp1) %>%
      as.data.frame() %>%
      View()
  }

  { ## total vif 删除共现性指标 ----

    data <- decostand(t(total.3), "hellinger")
    env <- pheno

    sel <- decorana(data)
    sel

    temp1 <- rda(t(data) ~ ., env)
    vif.cca(temp1) %>%
      as.data.frame() %>%
      View()
    temp1 <- rda(t(data) ~ ., env[, -c(4, 2, 3, 13, 19, 38, 39, 14, 24, 8, 10, 11)])
    vif.cca(temp1) %>%
      as.data.frame() %>%
      View()
  }
}

{ ## Figure 2C ----

  { ## incoming ----
    env.i <- env[, -c(4, 2, 3, 13, 19, 38, 39, 14, 24, 8, 10, 11)]
    rda <- rda(incoming ~ ., env.i)
    sum <- summary(rda)
    # total:45.3%

    fit <- envfit(rda, env.i, permutations = 999)

    variance <- data.frame(
      "Explained_Variance" = c(fit$factors$r, fit$vectors$r),
      "Pval" = c(fit$factors$pvals, fit$vectors$pvals)
    )

    variance$"Value" <- c(names(fit$factors$pvals), names(fit$vectors$pvals))
    variance$label <- "ns."
    variance$label[which(variance$Pval < 0.05)] <- "*"
    variance$label[which(variance$Pval < 0.01)] <- "**"
    variance$label[which(variance$Pval < 0.001)] <- "***"
    variance$label[which(variance$Pval < 0.1 & variance$Pval > 0.05)] <- "."

    plot.data <- variance %>% mutate(Value = fct_reorder(Value, (Explained_Variance)))
    order <- levels(plot.data$Value) %>%
      rev() %>%
      head(20)
    plot.data <- plot.data[plot.data$Value %in% order, ]

    p.i <- plot.data %>%
      ggplot() +
      geom_col(aes(x = Value, y = Explained_Variance, fill = Value)) +
      geom_text(aes(x = Value, y = Explained_Variance + 0.01, label = label)) +
      scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Reds"))(20)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_y_continuous(labels = scales::percent) +
      guides(fill = "none") +
      labs(y = "Explained Variance", x = "", title = "Incoming") +
      coord_flip()

    p.i

    plot.data$Type <- "Incoming"
    table.s12 <- plot.data
  }

  { ## outgoing ----
    env.o <- env[, -c(4, 2, 3, 13, 19, 38, 39, 14, 24, 8, 10, 11)]
    rda <- rda(outgoing ~ ., env.o)
    sum <- summary(rda)
    sum
    # total:46.6%

    fit <- envfit(rda, env.o, permutations = 999)

    variance <- data.frame(
      "Explained_Variance" = c(fit$factors$r, fit$vectors$r),
      "Pval" = c(fit$factors$pvals, fit$vectors$pvals)
    )

    variance$"Value" <- c(names(fit$factors$pvals), names(fit$vectors$pvals))
    variance$label <- "ns."
    variance$label[which(variance$Pval < 0.05)] <- "*"
    variance$label[which(variance$Pval < 0.01)] <- "**"
    variance$label[which(variance$Pval < 0.001)] <- "***"
    variance$label[which(variance$Pval < 0.1 & variance$Pval > 0.05)] <- "."

    plot.data <- variance %>% mutate(Value = fct_reorder(Value, (Explained_Variance)))
    order <- levels(plot.data$Value) %>%
      rev() %>%
      head(20)
    plot.data <- plot.data[plot.data$Value %in% order, ]

    p.o <- plot.data %>%
      ggplot() +
      geom_col(aes(x = Value, y = Explained_Variance, fill = Value)) +
      geom_text(aes(x = Value, y = Explained_Variance + 0.01, label = label)) +
      scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Reds"))(20)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_y_continuous(labels = scales::percent) +
      guides(fill = "none") +
      labs(y = "Explained Variance", x = "", title = "Outgoing") +
      coord_flip()

    p.o

    plot.data$Type <- "Outgoing"
    table.s12 <- rbind(table.s12, plot.data)
  }

  { ## total ----

    env.t <- env[, -c(4, 2, 3, 13, 19, 38, 39, 14, 24, 8, 10, 11)]
    rda <- rda(total.3 ~ ., env.t)
    sum <- summary(rda)
    sum
    # total:46.0%

    fit <- envfit(rda, env.t, permutations = 999)

    variance <- data.frame(
      "Explained_Variance" = c(fit$factors$r, fit$vectors$r),
      "Pval" = c(fit$factors$pvals, fit$vectors$pvals)
    )

    variance$"Value" <- c(names(fit$factors$pvals), names(fit$vectors$pvals))
    variance$label <- "ns."
    variance$label[which(variance$Pval < 0.05)] <- "*"
    variance$label[which(variance$Pval < 0.01)] <- "**"
    variance$label[which(variance$Pval < 0.001)] <- "***"
    variance$label[which(variance$Pval < 0.1 & variance$Pval > 0.05)] <- "."

    plot.data <- variance %>% mutate(Value = fct_reorder(Value, (Explained_Variance)))
    order <- levels(plot.data$Value) %>%
      rev() %>%
      head(20)
    plot.data <- plot.data[plot.data$Value %in% order, ]

    p.t <- plot.data %>%
      ggplot() +
      geom_col(aes(x = Value, y = Explained_Variance, fill = Value)) +
      geom_text(aes(x = Value, y = Explained_Variance + 0.01, label = label)) +
      scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Reds"))(20)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_y_continuous(labels = scales::percent) +
      guides(fill = "none") +
      labs(y = "Explained Variance", x = "", title = "Total") +
      coord_flip()

    p.t

    plot.data$Type <- "Total"
    table.s12 <- rbind(table.s12, plot.data)
  }

  p <- p.i + p.o + p.t
  p
  ggsave(p, filename = "Figure2C.pdf", width = 18, height = 8)

  table.s12 <- table.s12[, c("Type", "Value", "Explained_Variance", "Pval", "label")]
  colnames(table.s12) <- c("Type", "Index", "Explained Variance", "Pvalue", "Label")
  write_xlsx(table.s12, "Table S12.xlsx")
}


# 7. Figure 2D ------------------------------------------------------------
{ ## Figure 2D ----

  data <- pair.combined[, c(1, 2, 4, 5, 3, 16, 17, 21, 22, 26, 27)]
  colnames(data)[1:4] <- c("A.Species", "B.Species", "A.Pheno", "B.Pheno")

  a.single.cont <- single.cont
  colnames(a.single.cont) <- paste0("A.", colnames(a.single.cont))
  b.single.cont <- single.cont
  colnames(b.single.cont) <- paste0("B.", colnames(b.single.cont))
  data <- merge(data, a.single.cont, by = "A.Species")
  data <- merge(data, b.single.cont, by = "B.Species")
  data <- data[, c(2, 1, 3:ncol(data))]
  data <- merge(data, inter.6to3, by = "Interaction")
  data <- data[, c(1, ncol(data), 2:(ncol(data) - 1))]

  data$Interaction <- factor(data$Interaction, levels = c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition"))
  data$Interaction.1 <- factor(data$Interaction.1, levels = c("Positive", "Neutral", "Negative"))
  data[, -c(1:4)] <- apply(data[, -c(1:4)], 2, as.numeric)

  colnames(data)
  index <- "PM.total"
  # index <- "Size"
  temp <- data[, c("Interaction", paste0("A.", index), paste0("B.", index))]

  temp$mean <- rowMeans(temp[, -1])
  p4 <- temp %>% ggplot(aes(x = Interaction, y = mean, fill = Interaction)) +
    geom_boxplot(width = 0.3) +
    scale_fill_manual(values = colors.6) +
    theme_bw() +
    geom_signif(
      comparisons = list(c("Mutualism", "Commensalism"), c("Commensalism", "Neutralism"), c("Neutralism", "Amensalism"), c("Amensalism", "Exploitation"), c("Exploitation", "Competition")),
      step_increase = 0.1,
      tip_length = 0.015,
      map_signif_level = T,
      vjust = -0.01
    ) +
    labs(x = "", y = paste0("Mean ", index)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  colors.3 <- c("#42B540B2", "#ED0000B2", "#ECECED")
  names(colors.3) <- c("Positive", "Negative", "Neutral")

  temp <- data[, c("Interaction.1", paste0("A.", index), paste0("B.", index))]
  temp$mean <- rowMeans(temp[, -1])
  figure2d <- temp %>% ggplot(aes(x = Interaction.1, y = mean, fill = Interaction.1)) +
    geom_boxplot(width = 0.3) +
    scale_fill_manual(values = colors.3) +
    theme_bw() +
    geom_signif(
      comparisons = sign.comp(c("Positive", "Neutral", "Negative")),
      step_increase = 0.1,
      tip_length = 0.015,
      map_signif_level = T,
      vjust = -0.01
    ) +
    labs(x = "", y = paste0("Mean ", index))
  # theme(axis.text.x = element_text(angle=45,hjust=1))

  ggsave(figure2d, filename = paste0("Figure2D.pdf"), width = 4, height = 6)
}

# 8. Figure 2E ------------------------------------------------------------
{ ## flagellum ----

  flagellum <- single.dis[, c("Species", "Flagellum")]
  colnames(flagellum) <- c("Species", "Flagellum")

  inter <- interactions.and.taxa[, c("Species.A", "Species.B", "Pheno.A", "Pheno.B", "Interaction")]

  inter$flagellum.A <- ""
  inter$flagellum.B <- ""

  for (i in 1:nrow(inter)) {
    inter$flagellum.A[i] <- ifelse(flagellum$Flagellum[flagellum$Species == inter$Species.A[i]] == 0, "No", "Have")

    inter$flagellum.B[i] <- ifelse(flagellum$Flagellum[flagellum$Species == inter$Species.B[i]] == 0, "No", "Have")
  }

  inter$group <- ""
  for (i in 1:nrow(inter)) {
    if (inter$flagellum.A[i] == "No" & inter$flagellum.B[i] == "No") {
      inter$group[i] <- "Zero"
    } else if (inter$flagellum.A[i] == "Have" & inter$flagellum.B[i] == "Have") {
      inter$group[i] <- "Two"
    } else {
      inter$group[i] <- "One"
    }
  }

  chi.data <- as.data.frame(table(inter[, c("Interaction", "group")]))

  chi.data.1 <- inter[, c(5, 8)]
  chi.data.1$Interaction <- factor(chi.data.1$Interaction, levels = c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition"))
  chi.data.1$group <- factor(chi.data.1$group, levels = c("Zero", "One", "Two"))

  chi.data <- spread(chi.data, group, Freq)
  rownames(chi.data) <- chi.data$Interaction
  chi.data <- chi.data[, -1]
  chi <- chisq.test(chi.data)
  chi$p.value

  colors.6 <- c("#22B7B2", "#ACE6E9", "#ECECED", "#FFE4C3", "#FFBB7F", "#FF9148")
  names(colors.6) <- c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")

  p1 <- ggstatsplot::ggbarstats(chi.data.1, Interaction,
    group,
    results.subtitle = F,
    xlab = "Group",
    title = "p=2.69e-36"
  ) +
    scale_fill_manual(values = colors.6) +
    theme(plot.title = element_text(face = "italic", hjust = 0.5))
  p1

  colors.3 <- c("#42B540B2", "#ED0000B2", "#ECECED")
  names(colors.3) <- c("Positive", "Negative", "Neutralism")

  chi.data.1$inter.3 <- ""
  chi.data.1$inter.3[chi.data.1$Interaction %in% c("Amensalism", "Exploitation", "Competition")] <- "Negative"
  chi.data.1$inter.3[chi.data.1$Interaction %in% c("Neutralism")] <- "Neutralism"
  chi.data.1$inter.3[chi.data.1$Interaction %in% c("Mutualism", "Commensalism")] <- "Positive"
  chi.data.1$inter.3 <- factor(chi.data.1$inter.3, levels = c("Positive", "Neutralism", "Negative"))

  figure2e <- ggstatsplot::ggbarstats(chi.data.1, inter.3, group,
    results.subtitle = F,
    xlab = "Group",
    title = "p=1.39e-28"
  ) +
    scale_fill_manual(values = colors.3) +
    theme(plot.title = element_text(face = "italic", hjust = 0.5))

  figure2e
  ggsave(figure2e, filename = "Figure2E.pdf", width = 5, height = 6)
}
