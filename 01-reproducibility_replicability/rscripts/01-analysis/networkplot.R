library(network)
library(sna)
library(ggplot2)
library(ggnetwork)


# test_wal <- pivot_longer(Walsh_dada2_networkplot, cols = `(Intercept)`:ethnicityRecodedwhite)
# test_wal <- test_wal %>% filter(value == TRUE)
# test_wal <- test_wal %>% dplyr::select(-c(value))
# test_wal$study <- "Walsh"
# 
# test_ant <- pivot_longer(Antonio_dada2_networkplot, cols = `(Intercept)`:age)
# test_ant <- test_ant %>% filter(value == TRUE)
# test_ant <- test_ant %>% dplyr::select(-c(value))
# 
# test_ant$study <- "Antonio"
# 
# test_tse <- pivot_longer(Tsementzi_dada2_networkplot, cols = `(Intercept)`:age)
# test_tse <- test_tse %>% filter(value == TRUE)
# test_tse <- test_tse %>% dplyr::select(-c(value))
# test_tse$study <- "Tsementzi"
# 
# test_cha <- pivot_longer(Chao_dada2_networkplot, cols = `(Intercept)`:age)
# test_cha <- test_cha %>% filter(value == TRUE)
# test_cha <- test_cha %>% dplyr::select(-c(value))
# test_cha$study <- "Chao"
# 
# test_gre <- pivot_longer(Gressel_dada2_networkplot, cols = `(Intercept)`:histologyEC)
# test_gre <- test_gre %>% filter(value == TRUE)
# test_gre <- test_gre %>% dplyr::select(-c(value))
# test_gre$study <- "Gressel"


test_all <- rbind(Walsh_dada2_networkplot, Antonio_dada2_networkplot,
                  Tsementzi_dada2_networkplot, Chao_dada2_networkplot, Gressel_dada2_networkplot)
test_all <- test_all %>% dplyr::filter(name!='(Intercept)')
#test_all <- test_all %>% dplyr::filter(str_detect(name, "ethnicity|menopa")==FALSE)
test_all$taxon <- tolower(test_all$taxon)
#test_all <- test_all %>% filter(str_detect(taxon, 'peptoniphilus|ruminococcus|porphyromonas|prevotella|prophyromonas|varibaculum|fastidiosipila|anaeroglobus|criibacterium|lactobacillaceae|lactobacillus|finegoldia|campylobacter|anaerococcus|dialister|pseudomonas|brevundimonas|streptococcus|dnf00809|
#                                            firmicutes|spirochaetes|actinobacteria|atopobium|proteobacteria|bacteroides|blautia'))
test_all$taxon <- gsub(".*:","", test_all$taxon)
test_all$taxon <- gsub(" .*","", test_all$taxon)
test_all$taxon <- gsub("_.*","", test_all$taxon)

#test_all <- test_all %>% filter(expression=="positive")
test_all <- select(test_all, -c(value_lfc, value_diff_abun, expression))
test_all <- unique(test_all)

test_all1 <- test_all %>% dplyr::count(taxon, name , sort = TRUE)
freq <- test_all$n
test_all$name <- gsub('histologyEC', 'histology',as.character(test_all$name))
#test_all$name <- gsub('pHRecoded>4.5', 'pH',as.character(test_all$name))
test_all$name <- gsub('age', 'Age',as.character(test_all$name))
test_all$name <- gsub('ethnicityRecodeOther', 'Ethnicity',as.character(test_all$name))
test_all$name <- gsub('ethnicityRecodewhite', 'Ethnicity',as.character(test_all$name))
test_all$name <- gsub('ethnicityRecodeasian', 'Ethnicity',as.character(test_all$name))

test_all <- dplyr::filter(test_all, !(name == "histology" & n==1))

## Remove the taxa only associated with age, BMI, or ethnicity 
hist_only <- dplyr::filter(test_all, (name == "histology"))
age_only <- dplyr::filter(test_all, (name == "Age" & n==1))
age_only <- age_only[age_only$taxon %in% hist_only$taxon ,]
bmi_only <- dplyr::filter(test_all, (name == "BMI" & n==1))
bmi_only <- bmi_only[bmi_only$taxon %in% hist_only$taxon ,]
et_only <- dplyr::filter(test_all, (name == "Ethnicity" & n==1))
et_only <- et_only[et_only$taxon %in% hist_only$taxon ,]
test_all <- rbind(hist_only, age_only, bmi_only, et_only)

#expression <- test_all$expression
test_all <- test_all %>% dplyr::select(-c(n))


n_all <- network(as.matrix(test_all), matrix.type= "edgelist")
network::set.edge.attribute(n_all, "Number_of_studies", freq)
#network::set.edge.attribute(n_all, "expression", expression)


df_test <- data.frame(names=network.vertex.names(n_all))
df_test$type <- with(df_test, ifelse((names %in% c("(Intercept)", "Age", "BMI", 
                                                   "Ethnicity", "histology")), 
                                     "clinical covariate", "taxa"))
n_all%v%'vertex.names'<- df_test$names
n_all%v%'type'<- df_test$type

pdf("~/Desktop/species_4.pdf",  width=30, height=15.5)
det <- ggplot(ggnetwork(n_all, layout = "fruchtermanreingold", cell.jitter = 0), aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(linewidth = as.factor(Number_of_studies)), alpha = 0.2, curvature = 0.2) +
  geom_nodes(aes(x, y, size = type, color = type)) + 
  geom_nodelabel_repel(aes(label = ifelse(type == "clinical covariate", "", vertex.names)), size = 12) +
  geom_nodetext(aes(label = ifelse(type == "clinical covariate", vertex.names, "" )), size = 12, fontface = "bold") + 
  scale_color_manual("type", values = c("whitesmoke", "black"),  guide = "none") +
  scale_size_manual("type", values = c(70, 1), guide = "none") +
  scale_discrete_manual("linewidth", values = c(1, 2, 5, 7, 11)) +
 theme_blank() + theme(legend.text=element_text(size=18), plot.title = element_text(size = 18, hjust = 0.5)) + ggtitle("Differentially expressed taxa")
print(det)
dev.off()

ggplot(ggnetwork(n_all, layout = "fruchtermanreingold", cell.jitter = 1), aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(linewidth = Number_of_studies), alpha = 0.5, curvature = 0.2) +
  geom_nodelabel_repel(aes(label = ifelse(type == "clinical covariate", "", vertex.names)), size = 9) +
  geom_nodes(aes(fill =type, size = ifelse(type == "clinical covariate", 10, 0))) +
  geom_nodetext(aes(label = ifelse(type == "clinical covariate", vertex.names, "" ), size = type),
                fontface = "bold") + scale_size_manual("type", values = c(10, 0), guide = "none") + 
  scale_shape_manual("type", values = c(0, 0), guide = "none") + 
  scale_fill_manual("type", values = c("white", "black"), guide = "none") +theme_blank()
