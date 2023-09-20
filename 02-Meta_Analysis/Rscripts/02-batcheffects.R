phylo_use <- eval(parse(text = paste0(cohort, "_SOTAphyloseq_tree_raw")))
phylo_use <- prune_samples(sample_sums(phylo_use)> 0, phylo_use)
p1 <- phylo_use %>%
  tax_fix() %>%
  tax_transform("clr") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(axes = c(1, 2), color = "cohort", size = 2.5)

dist_obj <- phylo_use %>%
  tax_fix() %>%
  tax_transform("identity",) %>%
  dist_calc("aitchison")

mul_var_res <- mul_var_analysis(cohort = cohort, dist_mat = dist_obj@dist, phylo_to_use = phylo_use)
assign(paste0(cohort, "_SOTA_beta"), p1,.GlobalEnv)
print(mul_var_res)

p1 <- phylo_use %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("wunifrac") %>% 
  ord_calc("NMDS")%>%
  ord_plot(axes = c(1, 2), color = "cohort", size = 2.5) 
