source("adversary-experiment.R")

test.single.config <- function(idx, configs, trials, all=FALSE) { 
  cat("Running", idx, "\n")
  print(configs[idx,])
  
  results <- data.frame(index=numeric(), size.of.dom=logical(), method=character(), 
                        pt.uncovered=numeric(), adversary.influence=numeric(), ATE.true=numeric(), 
                        variable=numeric(), value=numeric(), pt.covered=numeric(), n=numeric(), 
                        graph.type=character(), power=numeric(), degree=numeric(), p=numeric(), 
                        mu=numeric(), ncoms=numeric(), maxc=numeric(), minc=numeric(), 
                        lambda_0=numeric(), lambda_1=numeric(), lambda_2=numeric(), stringsAsFactors=FALSE)
  
  graph.params <- build.graph.params(configs, idx)
  adversary.params <- list()
  adversary.params$model <- reduction.adv.model
  adversary.params$all <- all
  outcome.params <- build.outcome.params(configs[idx,"lambda_0"], configs[idx,"lambda_1"], configs[idx,"lambda_2"], configs[idx,"sd.noise"])
  clustering <- "infomap"
  
  for(i in 1:trials) {
    graph.params$ind <- i
    
    cat("trial", i, "\n")
    bias.behavior.ATE <- adversary.experiment(graph.params, clustering, adversary.params, outcome.params)
    bias.behavior.ATE$adversary.influence <- as.numeric(bias.behavior.ATE$adversary.influence)
    bias.behavior.ATE$gui.beta <- as.numeric(bias.behavior.ATE$gui.beta)
    bias.behavior.ATE$gui.gamma <- as.numeric(bias.behavior.ATE$gui.gamma)
    
    bias.behavior.ATE <- add.graph.params(bias.behavior.ATE, graph.params)
    bias.behavior.ATE <- add.outcome.params(bias.behavior.ATE, outcome.params)
    bias.behavior.ATE$graph.id <- configs[idx,"graph.no"]
    bias.behavior.ATE$adv.bias <- bias.behavior.ATE$nonadv.ATE - bias.behavior.ATE$ATE.adv.gui
    
    results <- rbind(results, bias.behavior.ATE)
    write.csv(results, paste0("results/adversary-results-", graph.params$graph.type, "-", idx, ".csv"))
  }
}

test <- function() { 
  test.all(100)  
}

test.all <- function(trials, all=FALSE) { 
  configs <- read.csv("all_adv_configurations.csv")
  
  for(idx in 1:length(configs[[1]])) { 
    test.single.config(idx, configs, trials, all)
  }
}

test.inf.distr <- function(trials=500) { 
  configs <- read.csv("binary_networks/all_graph_configurations.csv")
  graph.configs <- configs
  graph.configs[,c("X", "lambda_0", "lambda_1", "lambda_2")] <- NULL
  graph.configs <- unique(graph.configs)
  graph.ty <- ""
  
  all.infs <- data.frame(graph.type=character(), degree=numeric(), p=numeric(), power=numeric(), 
                         mu=numeric(), size=numeric(), graph.no=numeric(), infs=numeric())    
  
  for(idx in 1:length(graph.configs[[1]])) { 
    #if(configs[idx, "graph.type"] == "sbm") { 
    cat("Running", idx, "\n")
    graph.params <- build.graph.params(configs, idx)
    
    infs <- c()
    ex.trials <- trials
    if(trials > 100 & graph.params$graph.type == "sbm") ex.trials <- 100
    
    for(i in 1:ex.trials) {
      graph.params$ind <- i
      g <- generate.graph(graph.params)
      n <- length(V(g))
      A <- get.adjacency(g)
      transition <- (1/degree(g)) * A
      
      inf <- colSums(transition)
      inf <- inf / graph.params$n
      inf <- as.data.frame(inf)
      infs <- c(infs, inf)
    }
    infs <- unlist(infs)
    infs <- as.data.frame(infs)
    infs <- merge(infs, graph.configs[idx,])
    all.infs <- rbind(all.infs, infs)
  }
  
  ggplot(subset(all.infs, graph.type=="barabasi-albert"), aes(infs, color=as.factor(power))) + 
    stat_density(geom="line") + guides(color=guide_legend(title="power")) + xlab("Influence") + ylab("Density") + theme_bw()+ theme(text = element_text(size = 15)) + theme(legend.position="bottom") 
  ggplot(subset(all.infs, graph.type=="small-world"), aes(infs, color=as.factor(p))) +
    stat_density(geom="line")  + guides(color=guide_legend(title="p")) + xlab("Influence") + ylab("Density") + theme_bw()+ theme(legend.position="bottom") + theme(text = element_text(size = 15)) 
  ggplot(subset(all.infs, graph.type=="sbm"), aes(infs, color=as.factor(mu))) + 
    stat_density(geom="line") + guides(color=guide_legend(title=expression(mu))) + xlab("Influence") + ylab("Density") + theme_bw() + theme(legend.position="bottom") + theme(text = element_text(size = 15))
  ggplot(subset(all.infs, graph.type=="forest-fire"), aes(infs)) + 
    stat_density(geom="line") + xlab("Influence") + ylab("Density") + theme_bw() + theme(legend.position="bottom") + theme(text = element_text(size = 15))
  
  
  ggplot(subset(all.infs, (graph.type=="barabasi-albert" & power==.3) | (graph.type=="sbm" & mu==.2) | (graph.type == "small-world" & p == 0.05) | (graph.type == "forest-fire")), aes(infs, color=graph.type)) + 
    stat_density(geom="line") + guides(color=guide_legend(title="graph type")) + xlab("Influence") + ylab("Density") + theme_bw()+ theme(text = element_text(size = 15)) + theme(legend.position="bottom") + facet_wrap(~graph.type)
}

test.small.world <- function(trials) { 
  graph.params <- list()
  graph.params$graph.type <- "small-world"
  graph.params$degree <- 5
  graph.params$n <- 225
  graph.params$p <- .05
  
  adversary.params <- list()
  adversary.params$num.adv <- 5
  adversary.params$max <- FALSE
  adversary.params$model <- reduction.adv.model
  
  outcome.params <- list()
  outcome.params$lambda_0 <- -1.5
  outcome.params$lambda_1 <- 0.75
  outcome.params$lambda_2 <- 0.5
  outcome.params$sd.noise <- 1
  
  #set.seed(1337)
  bias.behavior.ATE <- adversary.experiment(graph.params, "infomap", adversary.params, outcome.params)
  for(i in 1:(trials-1)) {
    bias.behavior.ATE <- rbind(bias.behavior.ATE, adversary.experiment(graph.params, "infomap", adversary.params, outcome.params))
  }
  
  bias.behavior.ATE$value <- as.numeric(bias.behavior.ATE$value)
  bias.behavior.ATE$adversary.influence <- as.numeric(bias.behavior.ATE$adversary.influence)
  
  library(ggplot2)
  plot1 <- ggplot(subset(bias.behavior.ATE, variable=="ATE.adv.gui"), aes(adversary.influence, value, shape=variable, color=method)) + geom_point() + geom_abline(intercept = bias.behavior.ATE$ATE.true, slope=0)
  plot(plot1)
}