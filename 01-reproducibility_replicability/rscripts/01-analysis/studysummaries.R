obj  <- Tsementzi_SOTAphyloseq_tree_raw
sample_data(obj) %>% dplyr::count(histology)

data.frame(sample_data(obj)) %>%
  mutate_at(c('age'), as.numeric)%>%
  group_by(histology) %>%
  get_summary_stats(age, type = "median_mad")

test <- data.frame(sample_data(obj)) %>%
  mutate_at(c('age'), as.numeric) %>% na.omit()
mad(test$age)

sample_data(obj) %>%
  group_by(histology) %>%
  mutate_at(c('age'), as.numeric) %>% na.omit() %>%
  identify_outliers(age)

sample_data(obj) %>%
  group_by(histology) %>%
  shapiro_test(age)

data.frame(sample_data(obj)) %>% levene_test(age ~ histology)

stat.test <-data.frame(sample_data(obj)) %>% 
  mutate_at(c('age'), as.numeric) %>% na.omit() %>%
  t_test(age ~ histology) %>%
  add_significance()
stat.test

sample_data(obj) %>% dplyr::count(histology, pHRecoded)

dat <- data.frame(
  "<=4.5" = c(4, 1),
  ">4.5" = c(21, 3),
  "NA" = c(3, 4),
  row.names = c("Benign", "EC"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("<=4.5", ">4.5", "NA")

dat
chisq.test(dat)$expected
test <- chisq.test(dat)
test

sample_data(obj) %>% dplyr::count(histology, ethnicityRecode)

dat <- data.frame(
  "White" = c(17, 2),
  "Other" = c(11, 6),
  row.names = c("Benign", "EC"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("White", "Other")

dat
chisq.test(dat)$expected
test <- chisq.test(dat)
test
