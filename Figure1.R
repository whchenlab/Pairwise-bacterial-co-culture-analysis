# 0. Load packages --------------------------------------------------------
library(tidyverse)
library(readxl)
library(writexl)
library(ggsci)
library(ggpubr)
library(RColorBrewer)
library(drc)
library(ggpmisc)
library(patchwork)

# 1. Set working directory ------------------------------------------------------------
setwd("D:/OneDrive - hust.edu.cn/【00-doing】项目/【09】菌菌互作关系/【R】15 - Paper")
dir.create("Figure1", showWarnings = F)
dir.create("Figure1/input", showWarnings = F)
setwd("D:/OneDrive - hust.edu.cn/【00-doing】项目/【09】菌菌互作关系/【R】15 - Paper/Figure1")

# 2. Load data ------------------------------------------------------------
total.taxa <- as.data.frame(read_xlsx("input/20221113_Total_clean_taxa.xlsx"))
total.taxa$Species <- gsub("\\[", "", total.taxa$Species)
total.taxa$Species <- gsub("\\]", "", total.taxa$Species)
total.interactions <- as.data.frame(read_xlsx("input/20221113_Total_clean_interactions.xlsx"))

interactions.and.taxa <- total.interactions
interactions.and.taxa <- merge(interactions.and.taxa, total.taxa[, c(1:8)], by.x = "A", by.y = "Name", all.x = T)
interactions.and.taxa <- merge(interactions.and.taxa, total.taxa[, c(1:8)], by.x = "B", by.y = "Name", all.x = T)
colnames(interactions.and.taxa)[7:13] <- c("Kingdom.A", "Phylum.A", "Class.A", "Order.A", "Family.A", "Genus.A", "Species.A")
colnames(interactions.and.taxa)[14:20] <- c("Kingdom.B", "Phylum.B", "Class.B", "Order.B", "Family.B", "Genus.B", "Species.B")
interactions.and.taxa <- interactions.and.taxa[, c(2, 1, 3:ncol(interactions.and.taxa))]
interactions.and.taxa$Species.A <- gsub("\\[", "", interactions.and.taxa$Species.A)
interactions.and.taxa$Species.A <- gsub("\\]", "", interactions.and.taxa$Species.A)
interactions.and.taxa$Species.B <- gsub("\\[", "", interactions.and.taxa$Species.B)
interactions.and.taxa$Species.B <- gsub("\\]", "", interactions.and.taxa$Species.B)


colors.6 <- c("#22B7B2", "#ACE6E9", "#ECECED", "#FFE4C3", "#FFBB7F", "#FF9148")
names(colors.6) <- c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")
colors.3 <- c("#42B540B2", "#ED0000B2", "#ECECED")
names(colors.3) <- c("Positive", "Negative", "Neutral")

inter.6to3 <- data.frame(
  Interaction = c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition"),
  Interaction.1 = c("Positive", "Positive", "Neutral", "Negative", "Negative", "Negative")
)

interactions.and.taxa <- merge(interactions.and.taxa, inter.6to3, by = "Interaction")
interactions.and.taxa <- interactions.and.taxa[, c(2, 3, 4, 5, 1, 6:ncol(interactions.and.taxa))]

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

load("species.inter.3.RData")

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

total.3 <- merge(incoming, outgoing, by = "row.names")
rownames(total.3) <- total.3$Row.names
total.3 <- total.3[, -1]
total.3 <- total.3[rownames(total.3) != "Klebsiella pneumoniae(sp)", ]



# 3. Figure 1A & Figure S1 --------------------------------------------------
{ # Load mapping rate: PRJNA835720 ----

  project.data <- list()
  taxa.temp <- as.data.frame(read_xlsx(paste0("input/", "PRJNA835720", "_taxa.xlsx")))
  meta.temp <- as.data.frame(read_xlsx(paste0("input/", "PRJNA835720", "_meta.xlsx")))
  map.temp <- as.data.frame(read.table(paste0("input/", "PRJNA835720", "_mapping.txt"), header = F, quote = "", sep = "\t"))

  # 对mapping数据进行处理
  colnames(map.temp) <- c("mapping.rate", "bac", "sample")
  map.temp$bac <- apply(map.temp, 1, function(x) {
    if (x[2] != "*") {
      paste(str_split(x[2], "_|\\.")[[1]][1], str_split(x[2], "_|\\.")[[1]][2])
    } else {
      x[2]
    }
  })
  map.temp <- spread(map.temp, sample, mapping.rate)
  rownames(map.temp) <- map.temp$bac
  map.temp <- map.temp[, -1]

  # 对mapping数据进行rate计算
  map.temp <- apply(map.temp, 2, function(x) {
    x / sum(x) * 100
  })
  map.temp <- as.data.frame(map.temp)
  map.temp <- map.temp[-1, ]

  # abun和meta数据保存至项目对应列表中
  project.data[["PRJNA835720"]][["metadata"]] <- meta.temp
  project.data[["PRJNA835720"]][["taxa"]] <- taxa.temp
  project.data[["PRJNA835720"]][["mapping"]] <- map.temp

  rm(taxa.temp, meta.temp, map.temp)
  save(project.data, file = "project.data.RData")
}

{ # Functions ----

  { ## Set colors ----

    phylum.color <- brewer.pal(5, "Set2")
    names(phylum.color) <- c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria", "Verrucomicrobia")

    country.color <- c("#FB8072", "#80B1D3", "#FDB462", "#BEBADA", "#8DD3C7") # Set3
    names(country.color) <- c("China", "USA", "Denmark", "Sweden", "Netherlands")
  }

  { ## Enzyme data ----

    uhgg.ec <- as.data.frame(read_xlsx("input/Total.EC.xlsx"))
    ec.file <- list.files(path = "input/BBI_enzyme/", pattern = "_")
    ec.data <- data.frame()
    for (i in 1:length(ec.file)) {
      ec.data.temp <- read.csv(paste0("input/BBI_enzyme/", ec.file[i]), sep = "\t", header = T)
      sum(ec.data.temp$EC_number != "")
      name <- paste(str_split(ec.file[i], "_|\\.")[[1]][1], str_split(ec.file[i], "_|\\.")[[1]][2])
      ec.data <- rbind(ec.data, data.frame(Taxa = rep(name, length(ec.data.temp$EC_number[ec.data.temp$EC_number != ""])), EC = ec.data.temp$EC_number[ec.data.temp$EC_number != ""]))
      colnames(ec.data) <- c("Taxa", "EC")
    }
    rm(i, name, ec.file, ec.data.temp)

    ec.data$value <- 1
    ec.data <- ec.data %>%
      group_by(Taxa, EC) %>%
      summarise(sum = sum(value))
    total.enzyme <- spread(ec.data, Taxa, sum, fill = 0)
    total.enzyme <- merge(uhgg.ec, total.enzyme, by.x = "final.ec", by.y = "EC", all.x = T, all.y = T) # 合并uhgg的酶
    total.enzyme[is.na(total.enzyme)] <- 0
    rownames(total.enzyme) <- total.enzyme$final.ec
    total.enzyme <- total.enzyme[, -1]
    total.enzyme <- as.data.frame(t(total.enzyme))
    rownames(total.enzyme)[rownames(total.enzyme) == "Clostridium leptum"] <- "[Clostridium] leptum"
    rownames(total.enzyme)[rownames(total.enzyme) == "Ruminococcus gnavus"] <- "[Ruminococcus] gnavus"
    total.enzyme <- total.enzyme[rownames(total.enzyme) != "Klebsiella pneumoniae(sp)", ]

    rm(ec.data, uhgg.ec)
  }

  abundance.boxplot <- function(data, sample.list, level, type) {
    data <- data[, sample.list]
    data <- merge(unique(total.taxa[, c(1:7)]), data, by.x = "Species", by.y = "row.names")
    data <- aggregate(data[-c(1:7)], by = list(data[, level]), sum)
    colnames(data)[1] <- "Taxa"
    data <- gather(data, sample, abun, -Taxa)
    data <- merge(unique(total.taxa[, c("Phylum", level)]), data, by.x = level, by.y = "Taxa")
    colnames(data)[1:2] <- c("Taxa", "Phylum")

    median <- data %>%
      group_by(Taxa) %>%
      summarise(median = median(abun)) %>%
      arrange(median)
    x.order <<- median$Taxa
    data$Taxa <- factor(data$Taxa, levels = median$Taxa)


    # 画图
    p <- data %>%
      ggplot(aes(x = Taxa, y = abun, fill = Phylum)) +
      geom_boxplot(outlier.alpha = 0, alpha = 0.9, size = 0.3) +
      scale_fill_manual(name = "Phylum", values = phylum.color) +
      coord_flip() +
      theme_bw() +
      # scale_x_discrete(limits = levels(data$taxa), labels = paste0(levels(data$taxa), " (", round(mean$mean, 4), "%)")) +
      scale_y_log10(position = "left", labels = scales::label_percent(scale = 1), limits = c(0.0001, 15)) +
      labs(x = level, y = "") +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(face = "italic")
      )

    if (level == "Genus") {
      num <- as.data.frame(table(total.taxa$Genus))
      num$label <- paste0(num$Var1, " [", num$Freq, "]")
      data <- merge(data, num, by.x = "Taxa", by.y = "Var1")
      x.label <- num$label
      names(x.label) <- num$Var1
      p <- p + scale_x_discrete(label = x.label) + theme(axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20))
    }

    if (type == "mapping") {
      p <- p + labs(y = "Mapping rate (Log-scaled)")
    } else {
      p <- p + labs(y = "Relative abundance (Log-scaled)")
    }

    p

    return(p)
  }
  prevalence.barplot <- function(data, sample.list, level, prev.cut) {
    data <- data[, sample.list]
    data <- merge(unique(total.taxa[, c(1:7)]), data, by.x = "Species", by.y = "row.names")
    data <- aggregate(data[-c(1:7)], by = list(data[, level]), sum)
    colnames(data)[1] <- "Taxa"
    data <- merge(unique(total.taxa[, c("Phylum", level)]), data, by.x = level, by.y = "Taxa")
    colnames(data)[1:2] <- c("Taxa", "Phylum")

    # prevalence计算
    data$prevalence <- apply(data[, 3:ncol(data)], 1, function(x) {
      sum(x > prev.cut) / length(x) * 100
    })

    # table s2
    # supp <- data[,c(1,2,ncol(data),3:(ncol(data)-1))]
    # supp <- supp %>% arrange(desc(prevalence))
    # # 还有肠型数据也一起输出
    # write_xlsx(supp,"01-代表性计算/02-output/Table S2.xlsx")

    data <- data[, c(1, 2, ncol(data))]

    # 按照abundance图中的x轴进行排序
    data$Taxa <- factor(data$Taxa, levels = x.order)

    # 画图
    p <- data %>%
      ggplot() +
      geom_bar(aes(x = Taxa, y = 100), alpha = 0.1, fill = "#FFFFFF", color = "grey", stat = "identity", size = 0.3, width = 0.7) +
      geom_bar(aes(x = Taxa, y = prevalence, fill = Phylum), stat = "identity", width = 0.65, show.legend = F) +
      scale_fill_brewer(name = "Phylum", palette = "Set2") +
      coord_flip() +
      theme_bw() +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(size = rel(0.5))) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      scale_y_continuous(position = "left", labels = scales::label_percent(scale = 1)) +
      labs(x = "", y = "Prevalence")
    p

    return(p)
  }
  cumulated.abun.coverage <- function(data, sample.list, level, type, meta) {
    data <- data[, sample.list]
    data <- merge(unique(total.taxa[, c(1:7)]), data, by.x = "Species", by.y = "row.names")
    data <- aggregate(data[-c(1:7)], by = list(data[, level]), sum)
    colnames(data)[1] <- "Taxa"

    data <- gather(data, sample, abun, -Taxa)
    data <- merge(data, meta[, c("sample", "country")], by = "sample")

    # 计算丰度平均值
    data <- data %>%
      group_by(country, Taxa) %>%
      summarise(mean = mean(abun))
    data <- as.data.frame(data)

    # 按照第一张图的y轴顺序计算累积值
    rownames(data) <- data$Taxa
    data <- data[x.order, ]
    data$Cumsum <- rev(cumsum(rev(data$mean)))
    data$Taxa <- factor(data$Taxa, levels = x.order)

    # 画图
    p <- data %>%
      ggplot(aes(x = Taxa, y = Cumsum, group = 1, color = country)) +
      geom_step(size = 0.8) +
      geom_point() +
      coord_flip() +
      scale_color_manual(name = "Country", values = country.color) +
      theme_bw() +
      labs(y = "", x = "") +
      # scale_y_continuous(limits = c(0,100),position = "left", labels = scales::label_percent(scale = 1))+
      scale_y_continuous(position = "left", labels = scales::label_percent(scale = 1)) +
      theme(legend.position = "bottom") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    # theme(text = element_text(size = 20),axis.title=element_text(size=15))

    if (type == "mapping") {
      p <- p + labs(y = "Cumulative mapping rate")
    } else {
      p <- p + labs(y = "Cumulative relative abundance")
    }

    p

    return(p)
  }
  cumulated.enzyme.coverage <- function(data, level) {
    data <- merge(unique(total.taxa[, c(1:7)]), data, by.x = "Species", by.y = "row.names")
    data <- aggregate(data[-c(1:7)], by = list(data[, level]), sum)
    colnames(data)[1] <- "Taxa"
    rownames(data) <- data$Taxa
    data <- data[, -1]

    # 根据前面生成的x.order的顺序依次计算累计的酶
    cumulated.data <- data.frame()
    x.order.1 <- rev(x.order)
    for (i in 1:length(x.order.1)) {
      plot.data <- data[x.order.1[1:i], ]
      cumulated.data <- rbind(cumulated.data, c(x.order.1[i], sum(colSums(plot.data) != 0)))
    }
    colnames(cumulated.data) <- c("Taxa", "Freq")
    cumulated.data$Freq <- as.numeric(cumulated.data$Freq)
    cumulated.data$Ratio <- cumulated.data$Freq / ncol(data) * 100
    cumulated.data$Taxa <- factor(cumulated.data$Taxa, levels = x.order)

    # 画图
    p <- cumulated.data %>%
      ggplot(aes(x = Taxa, y = Ratio, group = 1)) +
      geom_step(size = 0.8) +
      geom_point(color = "navyblue") +
      coord_flip() +
      theme_bw() +
      labs(y = "Cumulative enzyme coverage", x = "") +
      scale_y_continuous(position = "left", labels = scales::label_percent(scale = 1)) +
      theme(legend.position = "bottom") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    # theme(text = element_text(size = 20),axis.title=element_text(size=15))

    p

    return(p)
  }
  draw.combined.plot <- function(abun, meta, enzyme, sample.list, level, type, prev.cut) {
    p1 <- abundance.boxplot(
      data = abun,
      sample.list = sample.list,
      level = level,
      type = type
    )

    p2 <- prevalence.barplot(
      data = abun,
      sample.list = sample.list,
      level = level,
      prev.cut = prev.cut
    )

    p3 <- cumulated.abun.coverage(
      data = abun,
      sample.list = sample.list,
      level = level,
      type = type,
      meta = meta
    ) + guides(color = "none")

    p4 <- cumulated.enzyme.coverage(
      data = enzyme,
      level = level
    )

    p.all <- ggarrange(p1, p2, p4, p3, nrow = 1, widths = c(0.6, 0.2, 0.2, 0.2), common.legend = T, legend = "right")

    return(p.all)
  }
}

{ # Figure 1A: PRJNA835720 & genus & enterotypes ----

  p.genus <- draw.combined.plot(
    abun = project.data[["PRJNA835720"]][["mapping"]],
    meta = project.data[["PRJNA835720"]][["metadata"]],
    enzyme = total.enzyme,
    sample.list = project.data[["PRJNA835720"]][["metadata"]]$sample,
    level = "Genus",
    type = "mapping",
    prev.cut = 0.001
  )
  p.genus
  ggsave(p.genus, file = "Figure1A.pdf", width = 15, height = 18)
}

{ # Figure S1: all projects & species ----

  p.species <- draw.combined.plot(
    abun = project.data[["PRJNA835720"]][["mapping"]],
    meta = project.data[["PRJNA835720"]][["metadata"]],
    enzyme = total.enzyme,
    sample.list = project.data[["PRJNA835720"]][["metadata"]]$sample,
    level = "Species",
    type = "mapping",
    prev.cut = 0.001
  )
  p.species
  ggsave(p.species, file = "Supp1.pdf", width = 15, height = 18)
}


# 4. Figure 1D -------------------------------------------------------------
{ ## Load data ----

  load("input/project.data.RData")
  project.data[["PRJNA835720"]][["mapping"]] <- project.data[["PRJNA835720"]][["mapping"]][rownames(project.data[["PRJNA835720"]][["mapping"]]) != "Klebsiella pneumoniae(sp)", ]

  load("input/increase.ratio.rda")

  compute.interaction <- function(cutoff) {
    type <- c("west_anaero", "fiber_anaero", "west_aero", "fiber_aero")

    result <- list()
    for (i in type) {
      data <- get(i)
      data$strain1.ratio <- as.numeric(data$strain1.ratio)
      data$strain2.ratio <- as.numeric(data$strain2.ratio)
      data$strain1.pheno <- ifelse(data$strain1.ratio >= cutoff, 1, ifelse(data$strain1.ratio <= cutoff * (-1), -1, 0))
      data$strain2.pheno <- ifelse(data$strain2.ratio >= cutoff, 1, ifelse(data$strain2.ratio <= cutoff * (-1), -1, 0))

      data$sum <- data$strain1.pheno + data$strain2.pheno
      data$interaction <- ifelse(data$sum == 2, "Mutualism",
        ifelse(data$sum == -2, "Competition",
          ifelse(data$sum == 1, "Commensalism",
            ifelse(data$sum == -1, "Amensalism",
              ifelse(data$strain1.pheno == 0 & data$strain2.pheno == 0, "Neutralism", "Exploitation")
            )
          )
        )
      )

      data <- data[, -c(ncol(data) - 1)]
      colnames(data)[ncol(data)] <- "Interaction"

      data <- data[, c("Strain1", "Strain2", "strain1.pheno", "strain2.pheno", "Interaction")]
      colnames(data) <- c("A", "B", "Pheno.A", "Pheno.B", "Interaction")

      result[[i]] <- data
    }
    return(result)
  }
  cut10 <- compute.interaction(cutoff = 10)

  rm(west_anaero, fiber_anaero, west_aero, fiber_aero)
}
{ ## Set colors ----
  # 设置各组颜色
  colors.6 <- c("#22B7B2", "#ACE6E9", "#ECECED", "#FFE4C3", "#FFBB7F", "#FF9148")
  names(colors.6) <- c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")
}
{ ## Figure1D: 6种互作关系 柱状图 各项研究比较 ----

  # NBT
  plot.data <- data.frame()
  for (i in names(cut10)) {
    temp <- cut10[[i]]
    temp <- as.data.frame(table(temp$Interaction))
    temp$group <- i
    plot.data <- rbind(plot.data, temp)
  }
  plot.data$group[plot.data$group == "west_anaero"] <- "Western, anaerobic"
  plot.data$group[plot.data$group == "fiber_anaero"] <- "High fiber, anaerobic"
  plot.data$group[plot.data$group == "west_aero"] <- "Western, aerobic"
  plot.data$group[plot.data$group == "fiber_aero"] <- "High fiber, aerobic"

  # our
  # 如果要去掉反向点样的，就只剩1600对
  # temp <- interactions.and.taxa[,c("Species.A","Species.B","Interaction")]
  # temp1 <- rbind(temp,temp[,c(2,1,3)]) %>% subset(Species.A>Species.B) %>% unique()
  # our <- as.data.frame(table(temp1$Interaction))
  our <- as.data.frame(table(total.interactions$Interaction))

  our$group <- "PairInteraX"
  plot.data <- rbind(plot.data, our)
  plot.data <- plot.data %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(sum = sum(Freq))

  # human 66 pairs
  human66 <- data.frame(
    Var1 = c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition"),
    Freq = c(1, 2, 1, 24, 17, 21),
    group = rep("Venturelli_MSB_2018", 6),
    sum = rep(66, 6)
  )
  plot.data <- rbind(plot.data, human66)
  plot.data$ratio <- plot.data$Freq / plot.data$sum * 100

  # 设置画图顺序
  plot.data$Var1 <- factor(plot.data$Var1, levels = rev(c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")))
  plot.data <- plot.data[rev(order(plot.data$Var1)), ]

  plot.data$group <- factor(plot.data$group, levels = c("PairInteraX", "Venturelli_MSB_2018", "Western, anaerobic", "High fiber, anaerobic", "Western, aerobic", "High fiber, aerobic"))

  plot.data$facet <- ""
  plot.data$facet[plot.data$group == "PairInteraX" | plot.data$group == "Venturelli_MSB_2018"] <- "Experimentally validated"
  plot.data$facet[!(plot.data$group == "PairInteraX" | plot.data$group == "Venturelli_MSB_2018")] <- "In Silico simulated (Magnúsdóttir_NBT_2017)"
  label <- plyr::ddply(plot.data, "group", transform, label_y = cumsum(ratio) - 0.5 * ratio)

  figure1d <- ggplot(label, aes(x = group, y = ratio, fill = Var1)) +
    geom_bar(stat = "identity", position = "stack", width = 0.3) +
    scale_fill_manual(name = "Interaction", values = colors.6) +
    theme_bw() +
    facet_grid(~facet, scales = "free", space = "free") +
    labs(x = "Group", y = "Ratio") +
    geom_text(aes(label = sum, y = 100), size = 3, vjust = -0.5) +
    theme(strip.background.x = element_rect(fill = "#FFFFFF", colour = "#000000")) +
    scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100))
  figure1d
  ggsave(figure1d, file = "Figure1D.pdf", width = 10, height = 6)

  rm(temp, our, human66)
}
{ ## Table S4 ----

  table.s4 <- label[, 1:5]
  table.s4$Num <- paste0(table.s4$Freq, " (", round(table.s4$ratio, 1), "%)")
  table.s4 <- table.s4[, c(1, 3, 6)]
  table.s4 <- spread(table.s4, group, Num)
  rownames(table.s4) <- table.s4$Var1
  table.s4 <- table.s4[rev(c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")), c("PairInteraX", "Venturelli_MSB_2018", "Western, anaerobic", "High fiber, anaerobic", "Western, aerobic", "High fiber, aerobic")]
  table.s4 <- rbind(table.s4, c(3233, 66, 298378, 298378, 298378, 298378))
  rownames(table.s4)[nrow(table.s4)] <- "Total"
  table.s4$Interaction <- rownames(table.s4)
  table.s4 <- table.s4[, c(ncol(table.s4), 1:(ncol(table.s4) - 1))]
  write_xlsx(table.s4, "Table S4.xlsx")
}
# 5. Figure S3 ------------------------------------------------------------
{ ## Figure S3 A/B ----

  dry <- cut10[["west_anaero"]]
  wet <- interactions.and.taxa

  match.result <- function(dry, wet) {
    wet <- wet[, c("Species.A", "Species.B", "Pheno.A", "Pheno.B", "Interaction")]
    wet$Species.A <- gsub("\\[", "", wet$Species.A)
    wet$Species.A <- gsub("\\]", "", wet$Species.A)
    wet$Species.B <- gsub("\\[", "", wet$Species.B)
    wet$Species.B <- gsub("\\]", "", wet$Species.B)

    species.list <- unique(total.taxa$Species)
    species.list <- gsub("\\[", "", species.list)
    species.list <- gsub("\\]", "", species.list)

    colnames(wet) <- colnames(dry)
    dry$A <- apply(dry, 1, function(x) {
      paste(str_split(x[1], "_")[[1]][1], str_split(x[1], "_")[[1]][2])
    })
    dry$B <- apply(dry, 1, function(x) {
      paste(str_split(x[2], "_")[[1]][1], str_split(x[2], "_")[[1]][2])
    })
    dry <- dry[dry$A %in% total.taxa$Species | dry$B %in% total.taxa$Species, ]

    match.result <- data.frame()
    for (i in 1:nrow(wet)) {
      temp <- dry[(dry$A == wet[i, 1]) & (dry$B == wet[i, 2]), ]
      temp1 <- wet[i, ]
      if (nrow(temp) != 0) {
        temp2 <- merge(temp, temp1, by = c("A", "B"))
        match.result <- rbind(match.result, temp2)
      }
    }
    colnames(match.result) <- c("A", "B", "Pheno.A.dry", "Pheno.B.dry", "Interaction.dry", "Pheno.A.wet", "Pheno.B.wet", "Interaction.wet")

    matrix <- as.data.frame(t(table(match.result[, c("Interaction.dry", "Interaction.wet")])))
    matrix <- spread(matrix, Interaction.wet, Freq)
    rownames(matrix) <- matrix$Interaction.dry
    matrix <- matrix[, -1]
    # return(matrix)
    return(match.result)
  }
  west_anaero <- match.result(dry = cut10[["west_anaero"]], wet = interactions.and.taxa)
  west_aero <- match.result(dry = cut10[["west_aero"]], wet = interactions.and.taxa)
  fiber_anaero <- match.result(dry = cut10[["fiber_anaero"]], wet = interactions.and.taxa)
  fiber_aero <- match.result(dry = cut10[["fiber_aero"]], wet = interactions.and.taxa)

  matched.cut10 <- list()
  matched.cut10[["west_anaero"]] <- west_anaero
  matched.cut10[["west_aero"]] <- west_aero
  matched.cut10[["fiber_anaero"]] <- fiber_anaero
  matched.cut10[["fiber_aero"]] <- fiber_aero

  { ## Supp3A: 6种互作关系 匹配的pair之间进行比较 ----

    # NBT
    plot.data <- data.frame()
    for (i in names(matched.cut10)) {
      temp <- matched.cut10[[i]]
      temp <- as.data.frame(table(temp$Interaction.dry))
      temp$group <- i
      plot.data <- rbind(plot.data, temp)
    }
    plot.data$group[plot.data$group == "west_anaero"] <- "Western, anaerobic"
    plot.data$group[plot.data$group == "fiber_anaero"] <- "High fiber, anaerobic"
    plot.data$group[plot.data$group == "west_aero"] <- "Western, aerobic"
    plot.data$group[plot.data$group == "fiber_aero"] <- "High fiber, aerobic"

    # our
    our <- as.data.frame(table(west_anaero$Interaction.wet))
    our$group <- "PairInteraX"
    plot.data <- rbind(plot.data, our)
    plot.data <- plot.data %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(sum = sum(Freq))
    plot.data$ratio <- plot.data$Freq / plot.data$sum * 100

    # 设置画图顺序
    plot.data$Var1 <- factor(plot.data$Var1, levels = rev(c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")))
    plot.data <- plot.data[rev(order(plot.data$Var1)), ]

    plot.data$group <- factor(plot.data$group, levels = c("PairInteraX", "Western, anaerobic", "High fiber, anaerobic", "Western, aerobic", "High fiber, aerobic"))

    plot.data$facet <- ""
    plot.data$facet[plot.data$group == "PairInteraX"] <- "Experimentally validated"
    plot.data$facet[!(plot.data$group == "PairInteraX")] <- "In Silico simulated (Magnúsdóttir_NBT_2017)"
    label <- plyr::ddply(plot.data, "group", transform, label_y = cumsum(ratio) - 0.5 * ratio)

    figure.s3a <- ggplot(label, aes(x = group, y = ratio, fill = Var1)) +
      geom_bar(stat = "identity", position = "stack", width = 0.3) +
      scale_fill_manual(name = "Interaction", values = colors.6) +
      theme_bw() +
      facet_grid(~facet, scales = "free", space = "free") +
      labs(x = "Group", y = "Ratio") +
      geom_text(aes(label = sum, y = 100), size = 3, vjust = -0.5) +
      theme(strip.background.x = element_rect(fill = "#FFFFFF", colour = "#000000")) +
      scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100))
    figure.s3a
    ggsave(figure.s2a, file = "Supp2A.pdf", width = 10, height = 6)

    rm(temp, our, human66)
  }

  { ## Supp3B: 6种互作关系 dry vs. wet heatmap ----

    # 画热图展示dry和wet的区别
    west.anaero.data <- as.data.frame(table(west_anaero[, c(5, 8)]))
    west.aero.data <- as.data.frame(table(west_aero[, c(5, 8)]))
    fiber.anaero.data <- as.data.frame(table(fiber_anaero[, c(5, 8)]))
    fiber.aero.data <- as.data.frame(table(fiber_aero[, c(5, 8)]))

    west.anaero.data <- spread(west.anaero.data, Interaction.dry, Freq)
    west.aero.data <- spread(west.aero.data, Interaction.dry, Freq)
    fiber.anaero.data <- spread(fiber.anaero.data, Interaction.dry, Freq)
    fiber.aero.data <- spread(fiber.aero.data, Interaction.dry, Freq)

    rownames(west.anaero.data) <- west.anaero.data$Interaction.wet
    west.anaero.data <- west.anaero.data[, -1]
    rownames(west.aero.data) <- west.aero.data$Interaction.wet
    west.aero.data <- west.aero.data[, -1]
    rownames(fiber.anaero.data) <- fiber.anaero.data$Interaction.wet
    fiber.anaero.data <- fiber.anaero.data[, -1]
    rownames(fiber.aero.data) <- fiber.aero.data$Interaction.wet
    fiber.aero.data <- fiber.aero.data[, -1]

    order <- c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition")
    west.anaero.data <- west.anaero.data[order, order]
    west.aero.data <- west.aero.data[order, order]
    fiber.anaero.data <- fiber.anaero.data[order, order]
    fiber.aero.data <- fiber.aero.data[order, order]

    pdf("Supp3B.west.anaero.pdf", width = 8, height = 8)
    pheatmap::pheatmap(west.anaero.data,
      display_numbers = T,
      number_format = "%.0f",
      cluster_cols = F,
      cluster_rows = F,
      cellheight = 20,
      cellwidth = 20,
      color = colorRampPalette(c("white", "firebrick3"))(50)
    )
    dev.off()

    pdf("Supp3B.west.aero.pdf", width = 8, height = 8)
    pheatmap::pheatmap(west.aero.data,
      display_numbers = T,
      number_format = "%.0f",
      cluster_cols = F,
      cluster_rows = F,
      cellheight = 20,
      cellwidth = 20,
      color = colorRampPalette(c("white", "firebrick3"))(50)
    )
    dev.off()

    pdf("Supp3B.fiber.anaero.pdf", width = 8, height = 8)
    pheatmap::pheatmap(fiber.anaero.data,
      display_numbers = T,
      number_format = "%.0f",
      cluster_cols = F,
      cluster_rows = F,
      cellheight = 20,
      cellwidth = 20,
      color = colorRampPalette(c("white", "firebrick3"))(50)
    )
    dev.off()

    pdf("Supp3B.fiber.aero.pdf", width = 8, height = 8)
    pheatmap::pheatmap(fiber.aero.data,
      display_numbers = T,
      number_format = "%.0f",
      cluster_cols = F,
      cluster_rows = F,
      cellheight = 20,
      cellwidth = 20,
      color = colorRampPalette(c("white", "firebrick3"))(50)
    )
    dev.off()
  }
}
{ ## Table S5 ----

  table.s5 <- label[, 1:5]
  table.s5$Num <- paste0(table.s5$Freq, " (", round(table.s5$ratio, 1), "%)")
  table.s5 <- table.s5[, c(1, 3, 6)]
  table.s5 <- spread(table.s5, group, Num)
  rownames(table.s5) <- table.s5$Var1
  table.s5 <- table.s5[, -1]
  table.s5$Interaction <- rownames(table.s5)
  table.s5 <- table.s5[, c(ncol(table.s5), 1:(ncol(table.s5) - 1))]
  write_xlsx(table.s5, "Table S5.xlsx")
}
{ ## Figure S3 C/D ----

  { ## Set colors ----

    # 设置各组颜色
    colors.3 <- c("#42B540B2", "#ED0000B2", "#ECECED")
    names(colors.3) <- c("Positive", "Negative", "Neutralism")
  }

  { ## spearman/pearson 相关性计算 ----

    PRJNA835720.mapping <- as.data.frame(t(project.data[["PRJNA835720"]][["mapping"]]))

    spearman <- psych::corr.test(PRJNA835720.mapping, method = "spearman")
    pearson <- psych::corr.test(PRJNA835720.mapping, method = "pearson")

    process.corr.data <- function(data) {
      cor <- data[["r"]] %>% as.data.frame()
      p <- data[["p"]] %>% as.data.frame()

      cor$Species.x <- rownames(cor)
      cor <- gather(cor, Species.y, correlation, -Species.x)
      cor <- cor[cor$Species.x != cor$Species.y, ]

      p$Species.x <- rownames(p)
      p <- gather(p, Species.y, p.value, -Species.x)
      p <- p[p$Species.x != p$Species.y, ]

      cor <- merge(cor, p, by = c("Species.x", "Species.y"))

      return(cor)
    }

    spearman <- process.corr.data(spearman)
    pearson <- process.corr.data(pearson)

    match.result.2 <- function(dry, wet) {
      wet <- wet[, c("Species.A", "Species.B", "Pheno.A", "Pheno.B", "Interaction")]
      wet$Species.A <- gsub("\\[", "", wet$Species.A)
      wet$Species.A <- gsub("\\]", "", wet$Species.A)
      wet$Species.B <- gsub("\\[", "", wet$Species.B)
      wet$Species.B <- gsub("\\]", "", wet$Species.B)
      colnames(wet)[1:2] <- c("A", "B")
      colnames(dry)[1:2] <- c("A", "B")


      match.result <- data.frame()
      for (i in 1:nrow(wet)) {
        dry.subset <- dry[(dry$A %in% c(wet[i, 1:2])) & (dry$B %in% c(wet[i, 1:2])), ]
        if (nrow(dry.subset) != 0) {
          if (dry.subset$A[1] == wet[i, 1]) {
            match.result <- rbind(match.result, c(c(dry.subset[1, ]), c(wet[i, c(3, 4, 5)])))
          } else {
            match.result <- rbind(match.result, c(c(dry.subset[1, ]), c(wet[i, c(4, 3, 5)])))
          }
        }
      }

      match.result$Interaction.dry[match.result$p.value > 0.05] <- "Neutralism"
      # match.result$Interaction.dry[match.result$p.value<0.05&match.result$correlation>0.4] <- "Positive"
      # match.result$Interaction.dry[match.result$p.value<0.05&match.result$correlation<(-0.4)] <- "Negative"
      # match.result$Interaction.dry[match.result$p.value<0.05&match.result$correlation<0.4&match.result$correlation>-0.4] <- "Neutralism"
      match.result$Interaction.dry[match.result$p.value < 0.05 & match.result$correlation > 0] <- "Positive"
      match.result$Interaction.dry[match.result$p.value < 0.05 & match.result$correlation < 0] <- "Negative"

      match.result$Interaction.wet[match.result$Interaction %in% c("Amensalism", "Exploitation", "Competition")] <- "Negative"
      match.result$Interaction.wet[match.result$Interaction %in% c("Neutralism")] <- "Neutralism"
      match.result$Interaction.wet[match.result$Interaction %in% c("Mutualism", "Commensalism")] <- "Positive"

      return(match.result)
    }
    spearman <- match.result.2(dry = spearman, wet = interactions.and.taxa)
    pearson <- match.result.2(dry = pearson, wet = interactions.and.taxa)
  }

  { ## 画图 ----

    { ## Supp3C: 3种互作关系 匹配的pair之间进行比较 ----

      # our
      plot.data <- data.frame()
      our <- as.data.frame(table(spearman$Interaction.wet))
      our$group <- "PairInteraX"
      plot.data <- rbind(plot.data, our)

      # pearson
      pearson.stats <- as.data.frame(table(pearson$Interaction.dry))
      pearson.stats$group <- "Pearson"
      plot.data <- rbind(plot.data, pearson.stats)

      # spearman
      spearman.stats <- as.data.frame(table(spearman$Interaction.dry))
      spearman.stats$group <- "Spearman"
      plot.data <- rbind(plot.data, spearman.stats)

      plot.data <- plot.data %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(sum = sum(Freq))
      plot.data$ratio <- plot.data$Freq / plot.data$sum * 100

      # 设置画图顺序
      plot.data$Var1 <- factor(plot.data$Var1, levels = rev(c("Positive", "Neutralism", "Negative")))
      plot.data <- plot.data[rev(order(plot.data$Var1)), ]

      plot.data$group <- factor(plot.data$group, levels = c("PairInteraX", "Spearman", "Pearson"))

      plot.data$facet <- ""
      plot.data$facet[plot.data$group == "PairInteraX"] <- "Experimentally validated"
      plot.data$facet[!(plot.data$group == "PairInteraX")] <- "In Silico simulated"
      label <- plyr::ddply(plot.data, "group", transform, label_y = cumsum(ratio) - 0.5 * ratio)

      figure.s3c <- ggplot(label, aes(x = group, y = ratio, fill = Var1)) +
        geom_bar(stat = "identity", position = "stack", width = 0.3) +
        scale_fill_manual(name = "Interaction", values = colors.3) +
        theme_bw() +
        facet_grid(~facet, scales = "free", space = "free") +
        labs(x = "Group", y = "Ratio") +
        geom_text(aes(label = sum, y = 100), size = 3, vjust = -0.5) +
        theme(strip.background.x = element_rect(fill = "#FFFFFF", colour = "#000000")) +
        scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100))
      figure.s3c
      ggsave(figure.s3c, file = "Supp2C.pdf", width = 6, height = 6)
    }

    { ## Supp3D: 3种互作关系 dry vs. wet heatmap ----

      # 画热图展示dry和wet的区别
      spearman.data <- as.data.frame(table(spearman[, c(8, 9)]))
      pearson.data <- as.data.frame(table(pearson[, c(8, 9)]))

      spearman.data <- spread(spearman.data, Interaction.dry, Freq)
      pearson.data <- spread(pearson.data, Interaction.dry, Freq)

      rownames(spearman.data) <- spearman.data$Interaction.wet
      spearman.data <- spearman.data[, -1]
      rownames(pearson.data) <- pearson.data$Interaction.wet
      pearson.data <- pearson.data[, -1]

      pdf("Supp3D.spearman.pdf", width = 8, height = 8)
      pheatmap::pheatmap(spearman.data,
        display_numbers = T,
        number_format = "%.0f",
        cluster_cols = F,
        cluster_rows = F,
        cellheight = 20,
        cellwidth = 20,
        color = colorRampPalette(c("white", "firebrick3"))(50)
      )
      dev.off()

      pdf("Supp3D.pearson.pdf", width = 8, height = 8)
      pheatmap::pheatmap(pearson.data,
        display_numbers = T,
        number_format = "%.0f",
        cluster_cols = F,
        cluster_rows = F,
        cellheight = 20,
        cellwidth = 20,
        color = colorRampPalette(c("white", "firebrick3"))(50)
      )
      dev.off()
    }

    { ## Table s6 ----

      table.s6 <- merge(spearman[, c(1:4, 8)], pearson[, c(1:4, 8)], by = c("A", "B")) %>% unique()
      colnames(table.s6) <- c("A", "B", "Spearman correlation", "Spearman P value", "Spearman result", "Pearson correlation", "Pearson P value", "Pearson result")
      write_xlsx(table.s6, "Table S6.xlsx")
    }

    { ## Table s7 ----

      table.s7 <- label[, 1:5]
      table.s7$Num <- paste0(table.s7$Freq, " (", round(table.s7$ratio, 1), "%)")
      table.s7 <- table.s7[, c(1, 3, 6)]
      table.s7 <- spread(table.s7, group, Num)
      rownames(table.s7) <- table.s7$Var1
      table.s7 <- table.s7[, -1]
      table.s7$Interaction <- rownames(table.s7)
      table.s7 <- table.s7[, c(ncol(table.s7), 1:(ncol(table.s7) - 1))]
      write_xlsx(table.s7, "Table S7.xlsx")
    }
  }
}


# 6. Figure 1E/f/G ----------------------------------------------------------
# load("input/project.data.RData")
# mapping.rate <- project.data[["PRJNA835720"]][["mapping"]][rownames(project.data[["PRJNA835720"]][["mapping"]])!="Klebsiella pneumoniae(sp)",]
load("input/PRJNA422434_China_mapping_rate.RData")
mapping.rate <- PRJNA422434.data
load("input/PRJCA016454_China_mapping_rate.RData")
mapping.rate <- PRJCA016454.data


{ ## Figure 1E ----

  compute.meanmapping <- function(a, b) {
    mapping.subset <- mapping.rate[rownames(mapping.rate) %in% c(a, b), ]
    mapping.subset[mapping.subset < 0.01] <- 0
    mapping.subset <- as.data.frame(t(mapping.subset))
    mapping.subset <- mapping.subset[, c(a, b)]
    result <- mean(colMeans(mapping.subset))

    return(result)
  }

  pairs <- interactions.and.taxa[, c("Species.A", "Species.B")] %>% unique()
  pairs.result <- data.frame()
  for (i in 1:nrow(pairs)) {
    print(i)
    a <- pairs[i, 1]
    b <- pairs[i, 2]
    result <- compute.meanmapping(a = a, b = b)
    pairs.result <- rbind(pairs.result, c(a, b, result))
  }
  colnames(pairs.result) <- c("Species.A", "Species.B", "Mean_mapping_rate")
  pairs.result <- merge(interactions.and.taxa[, c("Species.A", "Species.B", "Pheno.A", "Pheno.B", "Interaction")],
    pairs.result,
    by = c("Species.A", "Species.B")
  )
  pairs.result$Mean_mapping_rate <- as.numeric(pairs.result$Mean_mapping_rate)
  pairs.result$Interaction <- factor(pairs.result$Interaction, levels = c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition"))

  p1 <- pairs.result %>% ggplot(aes(x = Interaction, y = Mean_mapping_rate, fill = Interaction)) +
    geom_boxplot(width = 0.3) +
    geom_signif(
      comparisons = list(
        c("Mutualism", "Commensalism"),
        c("Commensalism", "Neutralism"),
        c("Neutralism", "Amensalism"),
        c("Amensalism", "Exploitation"),
        c("Exploitation", "Competition")
      ),
      step_increase = 0.05,
      tip_length = 0.005,
      map_signif_level = T
    ) +
    scale_fill_manual(values = colors.6) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_log10(labels = scales::percent_format(scale = 1)) +
    labs(y = "Mean mapping rate (Log-scaled)")
  p1

  pairs.result <- merge(pairs.result, inter.6to3, by = "Interaction")
  pairs.result$Interaction.1 <- factor(pairs.result$Interaction.1, levels = c("Positive", "Neutral", "Negative"))

  p2 <- pairs.result %>% ggplot(aes(x = Interaction.1, y = Mean_mapping_rate, fill = Interaction.1)) +
    geom_boxplot(width = 0.3) +
    geom_signif(
      comparisons = sign.comp(c("Positive", "Negative", "Neutral")),
      step_increase = 0.2,
      map_signif_level = T
    ) +
    scale_fill_manual(values = colors.3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_log10(labels = scales::percent_format(scale = 1)) +
    labs(y = "Mean mapping rate (Log-scaled)")
  p2

  figure1e <- p1 + p2 + plot_layout(widths = c(2, 1))
  figure1e
  ggsave(figure1e, filename = "Figure1E.pdf", width = 13, height = 8)
}

{ ## Figure 1F ----

  compute.cooccur <- function(a, b, cutoff) {
    mapping.subset <- mapping.rate[rownames(mapping.rate) %in% c(a, b), ]
    mapping.subset[mapping.subset > cutoff] <- 1
    mapping.subset[mapping.subset < cutoff] <- 0

    mapping.subset <- as.data.frame(t(mapping.subset))
    mapping.subset <- mapping.subset[, c(a, b)]
    colnames(mapping.subset) <- c("a", "b")

    a.freq <- sum(mapping.subset$a == 1) / nrow(mapping.subset) * 100
    b.freq <- sum(mapping.subset$b == 1) / nrow(mapping.subset) * 100
    ab.freq <- sum(mapping.subset$a == 1 & mapping.subset$b == 1) / nrow(mapping.subset) * 100
    a0.freq <- sum(mapping.subset$a == 1 & mapping.subset$b == 0) / nrow(mapping.subset) * 100
    b0.freq <- sum(mapping.subset$a == 0 & mapping.subset$b == 1) / nrow(mapping.subset) * 100

    result <- c(a.freq, b.freq, ab.freq, a0.freq, b0.freq)
    names(result) <- c("a.freq", "b.freq", "ab.freq", "a0.freq", "b0.freq")

    return(result)
  }

  pairs <- interactions.and.taxa[, c("Species.A", "Species.B")] %>% unique()
  pairs.result <- data.frame()
  for (i in 1:nrow(pairs)) {
    print(i)
    a <- pairs[i, 1]
    b <- pairs[i, 2]
    result <- compute.cooccur(a = a, b = b, cutoff = 0.001)
    pairs.result <- rbind(pairs.result, c(a, b, result))
  }
  colnames(pairs.result) <- c("Species.A", "Species.B", "A.freq", "B.freq", "AB.freq", "A0.freq", "B0.freq")
  pairs.result <- merge(interactions.and.taxa[, c("Species.A", "Species.B", "Pheno.A", "Pheno.B", "Interaction")],
    pairs.result,
    by = c("Species.A", "Species.B")
  )
  pairs.result[, 6:10] <- apply(pairs.result[, 6:10], 2, as.numeric)
  pairs.result$Interaction <- factor(pairs.result$Interaction, levels = c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition"))

  # plot

  p1 <- pairs.result %>% ggplot(aes(x = Interaction, y = AB.freq, fill = Interaction)) +
    geom_boxplot(width = 0.3) +
    geom_signif(
      comparisons = list(
        c("Mutualism", "Commensalism"),
        c("Commensalism", "Neutralism"),
        c("Neutralism", "Amensalism"),
        c("Amensalism", "Exploitation"),
        c("Exploitation", "Competition"),
        c("Amensalism", "Commensalism"),
        c("Neutralism", "Exploitation"),
        c("Neutralism", "Competition")
      ),
      step_increase = 0.05,
      map_signif_level = T
    ) +
    scale_fill_manual(values = colors.6) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Co-occurrence (%)") +
    ylim(0, 150)
  p1

  pairs.result <- merge(pairs.result, inter.6to3, by = "Interaction")
  pairs.result$Interaction.1 <- factor(pairs.result$Interaction.1, levels = c("Positive", "Neutral", "Negative"))

  p2 <- pairs.result %>% ggplot(aes(x = Interaction.1, y = AB.freq, fill = Interaction.1)) +
    geom_boxplot(width = 0.3) +
    geom_signif(
      comparisons = sign.comp(c("Positive", "Negative", "Neutral")),
      step_increase = 0.05,
      map_signif_level = T
    ) +
    scale_fill_manual(values = colors.3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Co-occurrence (%)") +
    ylim(0, 150)
  p2

  figure1f <- p1 + p2 + plot_layout(widths = c(2, 1))
  figure1f
  ggsave(figure1f, filename = "Figure1F.pdf", width = 13, height = 8)
}

{ ## Figure 1G ----

  { # Figure 1G: 分段取丰度Rank ----

    # temp <- interactions.and.taxa[,c("Species.A","Species.B","Interaction")]
    # temp1 <- rbind(temp,temp[,c(2,1,3)]) %>% subset(Species.A>Species.B) %>% unique()
    # interactions.and.taxa <- temp1
    observed.species <- colSums(mapping.rate > 0.001) %>% as.data.frame()
    split.group <- data.frame(
      lower = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90),
      upper = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 97)
    )
    result <- data.frame()
    for (j in 1:nrow(split.group)) {
      lower <- split.group[j, 1]
      upper <- split.group[j, 2]

      for (i in 1:ncol(mapping.rate)) {
        data <- mapping.rate[, i, drop = F]
        colnames(data) <- "sample"
        sample <- colnames(mapping.rate)[i]

        data.subset <- data[data$sample >= 0.001, , drop = F] %>% arrange(sample)
        # data.subset <- data %>% arrange(sample)

        data.subset <- head(data.subset, upper)
        if (lower != 0) {
          data.subset <- data.subset[-c(1:lower), , drop = F]
        }

        species <- rownames(data.subset)

        inter <- interactions.and.taxa[(interactions.and.taxa$Species.A %in% species) & (interactions.and.taxa$Species.B %in% species), ]
        inter.num <- nrow(inter)

        if (inter.num != 0) {
          prop <- as.data.frame(table(inter$Interaction) / inter.num * 100)
          # prop <- as.data.frame(table(inter$Interaction.1)/inter.num*100)
          prop$sample <- sample
          prop$lower <- lower
          prop$upper <- upper
          result <- rbind(result, prop)
        }
      }
    }

    result$group <- paste0("[", result$lower, ",", result$upper, "]")
    result <- result[, -c(4, 5)]

    result <- spread(result, Var1, Freq, fill = 0)
    result <- gather(result, Interaction, Prop, -c(sample, group))
    result$group <- factor(result$group, levels = c("[0,10]", "[10,20]", "[20,30]", "[30,40]", "[40,50]", "[50,60]", "[60,70]", "[70,80]", "[80,90]", "[90,97]"))
    result$Interaction <- factor(result$Interaction, levels = c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition"))

    p <- result %>% ggplot(aes(x = group, y = Prop, fill = Interaction)) +
      geom_boxplot(width = 0.5, linewidth = 0.3) +
      theme_bw() +
      scale_fill_manual(values = colors.6) +
      stat_summary(fun = median, geom = "line", aes(group = Interaction), linewidth = 0.5) +
      # stat_summary(fun=median, geom="point",shape=5)+
      facet_grid(. ~ Interaction, space = "free", scales = "free") +
      theme(strip.background = element_rect(fill = "white"), axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Abundance Rank (Top N; Low to High)", y = "Proportion") +
      # labs(x="Abundance Rank (Top N; High to Low)",y="Proportion")+
      scale_y_continuous(labels = scales::percent_format(scale = 1))

    p

    ggsave(p, filename = "Figure1G_丰度分段_cutoff0.001.pdf", width = 15, height = 6)
    # ggsave(p,filename="Figure1G_丰度分段_nocutoff.pdf",width=15,height = 6)
  }

  { ## Low to high ----

    result <- data.frame()
    for (topn in c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 97)) {
      for (i in 1:ncol(mapping.rate)) {
        data <- mapping.rate[, i, drop = F]
        colnames(data) <- "sample"
        # data <- data %>% arrange(desc(sample))
        # data <- data %>% arrange(sample)

        sample <- colnames(mapping.rate)[i]

        data.subset <- data[data$sample >= 0.001, , drop = F] %>% arrange(sample)
        data.subset <- head(data.subset, topn)
        species <- rownames(data.subset)

        inter <- interactions.and.taxa[(interactions.and.taxa$Species.A %in% species) & (interactions.and.taxa$Species.B %in% species), ]
        inter.num <- nrow(inter)

        if (inter.num != 0) {
          prop <- as.data.frame(table(inter$Interaction) / inter.num * 100)
          # prop <- as.data.frame(table(inter$Interaction.1)/inter.num*100)
          prop$sample <- sample
          prop$topn <- topn
          result <- rbind(result, prop)
        }
      }
    }
    result <- spread(result, Var1, Freq, fill = 0)
    result <- gather(result, Interaction, Prop, -c(sample, topn))
    result$topn <- as.character(result$topn)
    result$topn <- factor(result$topn, levels = c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 97))
    result$Interaction <- factor(result$Interaction, levels = c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition"))
    # result$Interaction <- factor(result$Interaction,levels=c("Positive", "Negative", "Neutral"))

    p <- result %>% ggplot(aes(x = topn, y = Prop, fill = Interaction)) +
      geom_boxplot(width = 0.5, linewidth = 0.3) +
      theme_bw() +
      scale_fill_manual(values = colors.6) +
      stat_summary(fun = median, geom = "line", aes(group = Interaction), linewidth = 0.5) +
      # stat_summary(fun=median, geom="point",shape=5)+
      facet_grid(. ~ Interaction, space = "free", scales = "free") +
      theme(strip.background = element_rect(fill = "white")) +
      labs(x = "Abundance Rank (Top N; Low to High)", y = "Proportion") +
      # labs(x="Abundance Rank (Top N; High to Low)",y="Proportion")+
      scale_y_continuous(labels = scales::percent_format(scale = 1))

    p
    ggsave(p, filename = "Figure1G.pdf", width = 15, height = 6)
  }

  { ## High to low ----

    result <- data.frame()
    for (topn in c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 97)) {
      for (i in 1:ncol(mapping.rate)) {
        data <- mapping.rate[, i, drop = F]
        colnames(data) <- "sample"
        # data <- data %>% arrange(desc(sample))
        # data <- data %>% arrange(sample)

        sample <- colnames(mapping.rate)[i]

        data.subset <- data[data$sample >= 0.001, , drop = F] %>% arrange(desc(sample))
        data.subset <- head(data.subset, topn)
        species <- rownames(data.subset)

        inter <- interactions.and.taxa[(interactions.and.taxa$Species.A %in% species) & (interactions.and.taxa$Species.B %in% species), ]
        inter.num <- nrow(inter)

        if (inter.num != 0) {
          prop <- as.data.frame(table(inter$Interaction) / inter.num * 100)
          # prop <- as.data.frame(table(inter$Interaction.1)/inter.num*100)
          prop$sample <- sample
          prop$topn <- topn
          result <- rbind(result, prop)
        }
      }
    }
    result <- spread(result, Var1, Freq, fill = 0)
    result <- gather(result, Interaction, Prop, -c(sample, topn))
    result$topn <- as.character(result$topn)
    result$topn <- factor(result$topn, levels = c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 97))
    result$Interaction <- factor(result$Interaction, levels = c("Mutualism", "Commensalism", "Neutralism", "Amensalism", "Exploitation", "Competition"))
    # result$Interaction <- factor(result$Interaction,levels=c("Positive", "Negative", "Neutral"))

    p <- result %>% ggplot(aes(x = topn, y = Prop, fill = Interaction)) +
      geom_boxplot(width = 0.5, linewidth = 0.3) +
      theme_bw() +
      scale_fill_manual(values = colors.6) +
      stat_summary(fun = median, geom = "line", aes(group = Interaction), linewidth = 0.5) +
      # stat_summary(fun=median, geom="point",shape=5)+
      facet_grid(. ~ Interaction, space = "free", scales = "free") +
      theme(strip.background = element_rect(fill = "white")) +
      # labs(x="Abundance Rank (Top N; Low to High)",y="Proportion")+
      labs(x = "Abundance Rank (Top N; High to Low)", y = "Proportion") +
      scale_y_continuous(labels = scales::percent_format(scale = 1))

    p
    ggsave(p, filename = "Figure1G_reverse.pdf", width = 15, height = 6)
  }
}

# 7. Figure 1H ---------------------------------------------------------------
{ # Figure 1H：各个关系丰度 和 shannon ----
  shannon <- as.data.frame(vegan::diversity(t(mapping.rate), index = "shannon", base = exp(1)))
  colnames(shannon) <- "Index"
  simpson <- as.data.frame(vegan::diversity(t(mapping.rate), index = "simpson"))
  colnames(simpson) <- "Index"

  index <- shannon %>% arrange(Index)

  result <- data.frame()
  for (inter in c("Positive", "Neutral", "Negative")) {
    for (i in 1:ncol(mapping.rate)) {
      data <- mapping.rate[, i, drop = F]
      colnames(data)[1] <- "sample"
      data <- data[data$sample >= 0.001, , drop = F] %>% arrange(sample)
      # data <- data%>% arrange(sample)

      inter.subset <- interactions.and.taxa[, c("Species.A", "Species.B", "Interaction")]
      inter.subset <- merge(inter.subset, inter.6to3, by = "Interaction")

      species <- inter.subset[inter.subset$Interaction.1 == inter, ]
      species <- unique(c(species$Species.A, species$Species.B))

      data.subset <- data[rownames(data) %in% species, , drop = F]
      sum <- sum(data.subset$sample)

      data.subset.1 <- data[rownames(data) %in% unique(total.taxa$Species), , drop = F]
      sum.1 <- sum(data.subset.1$sample)

      # result <- rbind(result,c(colnames(mapping.rate)[i],sum))

      result <- rbind(result, c(colnames(mapping.rate)[i], sum / sum.1 * 100, inter))
    }
  }

  colnames(result) <- c("Sample", "Accumulated.abundance", "Interaction")
  plot.data <- merge(result, index, by.x = "Sample", by.y = "row.names")
  plot.data$Accumulated.abundance <- as.numeric(plot.data$Accumulated.abundance)

  colors.3 <- c("#42B540B2", "#ED0000B2", "grey")
  names(colors.3) <- c("Positive", "Negative", "Neutral")

  p1 <- plot.data %>% ggplot(aes(x = Index, y = Accumulated.abundance, color = Interaction, fill = Interaction)) +
    geom_point() +
    geom_smooth(method = "lm") +
    stat_poly_eq(
      formula = y ~ x,
      label.x.npc = "left", label.y.npc = "bottom",
      aes(
        label = paste(..eq.label..,
          ..rr.label..,
          ..p.value.label..,
          sep = "~~~"
        ),
        color = Interaction
      ),
      parse = TRUE
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    scale_fill_manual(values = colors.3) +
    scale_color_manual(values = colors.3) +
    labs(x = "Shannon", y = "Abundance of the species in interactions/total abundance of 97 species") +
    theme_bw()
  p1
  ggsave(p1, filename = "Supp4.pdf", width = 6, height = 5)
}
