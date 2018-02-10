library(grid)

create.configurations <- function(base.dir=".") {
  sizes <- c(500, 1000, 5000)
  graph.settings <- expand.grid(graph.type=c("small-world"), degree=5, p=c(0.03, 0.05, 0.1), power=NA, mu=NA, fw=NA, bw=NA, size=sizes)
  graph.settings <- rbind(graph.settings, expand.grid(graph.type=c("barabasi-albert"), power=c(0.1, 0.3, 0.5), degree=NA, p=NA, mu=NA, fw=NA, bw=NA, size=sizes))
  graph.settings <- rbind(graph.settings, expand.grid(graph.type=c("sbm"), mu=c(0.1, .2, .3), degree=NA, p=NA, power=NA, fw=NA, bw=NA, size=sizes))
  
  ff.settings <- expand.grid(graph.type=c("forest-fire"), mu=NA, degree=NA, p=NA, power=NA, size=sizes)
  ff.settings$fw <- c(.32, 0.37, 0.37)
  ff.settings$bw <- c(.33, 0.33, 0.35)
  graph.settings <- rbind(graph.settings, ff.settings)
  graph.settings$graph.no <- 1:nrow(graph.settings)
  exp.settings.red <- expand.grid(lambda_0=-1.5, lambda_1=.75, lambda_2=.5)
  exp1 <- merge(graph.settings, exp.settings.red)
  write.csv(exp1, file.path(base.dir, "binary_networks/all_graph_configurations.csv"))
  
  exp.settings <- expand.grid(lambda_0=c(-1.5), lambda_1=c(0.25, 0.5, 0.75, 1), lambda_2=c(0, 0.1, 0.5, 1.0))
  
  graph.settings.red <- expand.grid(graph.type="small-world", degree=5, p=0.05, power=NA, mu=NA, fw=NA, bw=NA, size=1000)
  graph.settings.red <- rbind(graph.settings.red, expand.grid(graph.type="barabasi-albert", power=0.3, degree=NA, p=NA, mu=NA, fw=NA, bw=NA, size=1000))
  graph.settings.red <- rbind(graph.settings.red, expand.grid(graph.type="sbm", mu=.2, degree=NA, p=NA, power=NA, fw=NA, bw=NA, size=1000))
  graph.settings.red <- rbind(graph.settings.red, subset(ff.settings, size==1000))
  
  graph.settings.red <- merge(graph.settings, graph.settings.red, by=colnames(graph.settings.red))
  exp2 <- merge(graph.settings.red, exp.settings)
  write.csv(exp2, file.path(base.dir, "all_adv_configurations.csv"))
  #exp2 <- subset(exp2, !(graph.no %in% exp1$graph.no & lambda_1 == 0.75 & lambda_2 == 0.5))
  
  #all.settings <- rbind(exp1, exp2)
  all.settings <- merge(graph.settings, exp.settings)
  write.csv(all.settings, file.path(base.dir, "all_adv_configurations_comb.csv"))
}

create.graph.inf.configurations <- function(base.dir=".") {
  sizes <- c(500)
  graph.settings <- expand.grid(graph.type=c("small-world"), degree=5, p=c(0.03, 0.04, 0.05, 0.06, 0.07), power=NA, mu=NA, size=c(500))
  graph.settings <- rbind(graph.settings, expand.grid(graph.type=c("barabasi-albert"), power=c(0.1, 0.2, 0.3, 0.4, 0.5), degree=NA, p=NA, mu=NA, size=sizes))
  graph.settings <- rbind(graph.settings, expand.grid(graph.type=c("sbm"), mu=c(0.1, .2, .3), degree=NA, p=NA, power=NA, size=c(500)))
  graph.settings$graph.no <- 1:nrow(graph.settings)
  exp.settings.red <- expand.grid(lambda_0=-1.5, lambda_1=.75, lambda_2=.5)
  all.settings <- merge(graph.settings, exp.settings.red)
  
  write.csv(all.settings, file.path(base.dir, "all_inf_graph_configurations2.csv"))
}

create.fb.configurations <- function(base.dir=".") {
  graph.settings <- expand.grid(graph.type=c("facebook"), degree=NA, p=NA, power=NA, mu=NA, size=NA)
  graph.settings$graph.no <- 1:nrow(graph.settings)
  exp.settings <- expand.grid(lambda_0=c(0), lambda_1=c(1.845e-4, 1.929e-04, 1.995e-4), lambda_2=c(1.086e-3, 1.111e-3, 1.136e-3))
  all.settings <- merge(graph.settings, exp.settings)
  
  write.csv(all.settings, file.path(base.dir, "all_adv_configurations_rw.csv"))
}