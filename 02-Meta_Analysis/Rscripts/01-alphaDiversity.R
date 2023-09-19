library(meta)

my_Shannon <- function(x){
  x <- x[x > 0]
  S <- length(x)
  p <- x/sum(x)
  (-sum(p * log(p)))
}

all_cohorts <- c("Antonio", "Walsh", "Tsementzi", "Gressel", "Chao")

for(cohort in all_cohorts){
  print(paste0("cohort_pipeline", cohort, '_', pipeline))
  phylo <- eval(parse(text = paste0(cohort, '_', pipeline, 'phyloseq_tree_raw')))
  otu <- abundances(phylo)
  to_keep <- which((colSums(otu))!=0)
  otu <- otu[, to_keep]
  ev <- apply(otu, 2, function(x) {
    my_Shannon(x)
  })
  ps1.meta <- meta(phylo)
  ps1.meta <- ps1.meta[to_keep,]
  ps1.meta$shannon <- as.data.frame(ev)$ev
  assign(paste0(cohort, "_SOTA_alphaDiversity"),ps1.meta,.GlobalEnv)
  assign(paste0(cohort, "_", pipeline, "alphaDiversity_sig"),test.sig,.GlobalEnv)
}

all_cohorts <- plyr::rbind.fill(eval(parse(text = paste0('Antonio_SOTA_alphaDiversity'))),
                                eval(parse(text = paste0('Chao_SOTA_alphaDiversity'))),
                                eval(parse(text = paste0('Gressel_SOTA_alphaDiversity'))),
                                eval(parse(text = paste0('Tsementzi_SOTA_alphaDiversity'))),
                                eval(parse(text = paste0('Walsh_SOTA_alphaDiversity'))))
all_cohorts_long <- gather(all_cohorts, metric, value, shannon, factor_key=TRUE)

for_metaAnalysis <- all_cohorts %>%
                      group_by(cohort, histology) %>%
                      summarise_at(vars(shannon), list(n=length,mean=mean, sd=sd))

for_metaAnalysis_wide <- for_metaAnalysis %>%
                              pivot_wider(names_from = histology,
                                          names_glue = "{histology}_{.value}",
                                          values_from = c(n, mean, sd))
for_metaAnalysis_wide$Year <- c("2016", "2022", "2021", "2020", "2019")
rownames(for_metaAnalysis_wide) <- for_metaAnalysis_wide$cohort
m1 <- metacont(EC_n, EC_mean, EC_sd, Benign_n, Benign_mean, Benign_sd,
               data = for_metaAnalysis_wide, common = FALSE, 
               random = TRUE, outclab = "shannon index", label.e = "EC",
               label.c = "Benign",
               sm = "SMD", studlab = paste0(cohort, " et al., ", Year))
m1
forest(m1)


all_cohorts_ft <- all_cohorts %>%
  dplyr::select(-c(X, shannon))



