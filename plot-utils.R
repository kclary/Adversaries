library(ggplot2)

prepare.for.plots <- function(g, adversaries, adversary.exposure, treatment.assignments="", labels=FALSE) { 
  V(g)$color <- "lightblue"
  V(g)$color[which(adversary.exposure[1,] > 0)] <- "orange"
  V(g)$color[which(adversaries > 0)] <- "red"
  
  bdr <- rep("black", length(V(g)))
  if(length(treatment.assignments) > 2) bdr <- ifelse(treatment.assignments, "green", "black")
  
  if(labels == FALSE) plot(g, vertex.frame.color=bdr, vertex.label=NA, layout=layout_with_fr) 
  else plot(g, vertex.frame.color=bdr) 
}

plot.realworld.ATE.bias <- function() { 
  res <- read.csv("adversary-results-polyblog.csv")  
  res$X <- NULL
  res.1 <- subset(res, variable=="ATE.adv.gui")
  res.1 <- rename(res.1, c("value"="ATE.adv.gui"))
  res.1$variable <- NULL
  
  res.2 <- subset(res, variable=="ATE.adv.est")
  res.2 <- rename(res.2, c("value"="ATE.adv.est"))
  res.2$variable <- NULL
  
  res2 <- merge(res.1, res.2, by=c("index", "method", "size.of.dom", "pt.uncovered", 
                                   "adversary.influence",  "ATE.true", "pt.covered", "n", "graph.type", "power", "p", "mu", "ncoms", "maxc", "minc", "lambda_0", "lambda_1", "lambda_2"))
  res2$bias <- res2$ATE.true - res2$ATE.adv.gui
  res.plot <- subset(res2, method != "influence" & size.of.dom==FALSE)
  ggplot(res.plot, aes(index, bias, color=method)) + geom_point()
  
}

plot.increase.ATE.bias <- function(res, g.type) {
  #res <- read.csv("adversary-results-revised.csv")
  #res2 <- read.csv("adversary-results-revised-sbm.csv")
  #res <- rbind(res, res2)
  res <- read.csv("results/all-results.csv")
  
  res$bias <- res$ATE.true - res$ATE.adv.gui
  res$est.diff <- res$nonadv.ATE - res$ATE.adv.gui
  res$bias.norm <- res$bias / res$ATE.true
  res$diff.norm <- res$est.diff / res$nonadv.ATE
  
  #res <- subset(res, method != "degree")
  res$method <- ifelse(res$method == "random", "random", "dominating")
  res$graph.type <- ifelse(res$graph.type == "barabasi-albert", "scale-free", as.character(res$graph.type))
  res$graph.type <- ifelse(res$graph.type == "sbm", "SBM", as.character(res$graph.type))
  res$lambda_1_lab <- paste0("\u03BB_1 = ", as.character(res$lambda_1))
  res$lambda_2_lab <- paste0("\u03BB_2 = ", as.character(res$lambda_2))
  
  #remove late indices
  res <- subset(res, !(graph.type == "forest-fire" & index > 226))
  #res <- subset(res, !(graph.type == "small-world" & index > 64))
  #res <- subset(res, !(graph.type == "scale-free" & index > 186))
  res <- subset(res, !(graph.type == "SBM" & index > 70))
  res$pt.adversaries <- res$index / res$n
  
  # remove either degree or influence 
  
  df <- subset(res, graph.type == "scale-free")
  #stats.barabasi <- df[,c("method", "lambda_1", "lambda_2", "index", "power", "diff.norm")] %>% group_by(method, lambda_1, lambda_2, index, power) %>% summarize_each_(funs(mean, sd, median, n), vars="diff.norm") 
  plot1 <- ggplot(subset(df, lambda_1 == 0.75 & lambda_2 == 0.5), aes(pt.adversaries, diff.norm, color=method)) + 
    geom_smooth() + facet_wrap(~power) + xlab("Adversarial fraction of network") + 
    ylab("Bias in estimated ATE / Estimated nonadversarial ATE") + geom_abline(slope=0) + 
    theme_bw()+ theme(text = element_text(size = 15)) + theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1))
  plot(plot1)
  
  plot4 <- ggplot(df, aes(pt.adversaries, diff.norm, color=method)) + geom_smooth() + facet_grid(lambda_1_lab ~ lambda_2_lab) + 
    xlab("Adversarial fraction of network") + ylab("Bias in Estimated ATE / Estimated nonadversarial ATE") + 
    geom_abline(slope=0) + theme_bw()+ theme(text = element_text(size = 15)) + theme(legend.position="bottom") +
    guides(color=guide_legend(override.aes=list(fill=NA))) + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1))
  plot(plot4) 
  
  df <- subset(res, size.of.dom==FALSE & graph.type == "small-world")
  #stats.sw <- df[,c("method", "lambda_1", "lambda_2", "index", "p", "diff.norm")] %>% group_by(method, lambda_1, lambda_2, index, p) %>% summarize_each_(funs(mean, sd, median), vars="diff.norm") 
  #stats.sw <- subset(stats.sw, index < 31)
  plot2 <- ggplot(subset(df, lambda_1 == 0.75 & lambda_2 == 0.5), aes(pt.adversaries, diff.norm, color=method)) + 
    geom_smooth() + facet_wrap(~p) + xlab("Adversarial fraction of network") + ylim(c(0,1)) +
    ylab("Bias in Estimated ATE / Estimated nonadversarial ATE") + geom_abline(slope=0) + 
    theme_bw() + theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1))
  plot(plot2)
  
  plot3 <- ggplot(df, aes(pt.adversaries, diff.norm, color=method)) + geom_smooth() + facet_grid(lambda_1_lab ~ lambda_2_lab) + 
    xlab("Adversarial fraction of network") + ylab("Bias in Estimated ATE / Estimated nonadversarial ATE") + 
    geom_abline(slope=0) + theme_bw() + theme(text = element_text(size = 15)) + 
    theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1))
  plot(plot3) 
  
  df <- subset(res, size.of.dom==FALSE & graph.type == "SBM")
  #stats.sw <- df[,c("method", "lambda_1", "lambda_2", "index", "p", "diff.norm")] %>% group_by(method, lambda_1, lambda_2, index, p) %>% summarize_each_(funs(mean, sd, median), vars="diff.norm") 
  #stats.sw <- subset(stats.sw, index < 31)
  plot7 <- ggplot(subset(df, lambda_1 == 0.75 & lambda_2 == 0.5), aes(pt.adversaries, diff.norm, color=method)) + 
    geom_smooth() + facet_wrap(~mu) + xlab("Adversarial fraction of network") + ylim(c(0,1)) +
    ylab("Bias in Estimated ATE / Estimated nonadversarial ATE") + geom_abline(slope=0) + 
    theme_bw() + theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1))
  plot(plot7)
  
  plot8 <- ggplot(df, aes(pt.adversaries, diff.norm, color=method)) + geom_smooth() + facet_grid(lambda_1_lab ~ lambda_2_lab) + 
    xlab("Adversarial fraction of network (SBM)") + ylab("Bias in Estimated ATE / Estimated nonadversarial ATE") + 
    geom_abline(slope=0) + theme_bw() + theme(text = element_text(size = 15)) + 
    theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1))
  plot(plot8) 
  
  df <- subset(res, size.of.dom==FALSE & graph.type == "forest-fire")
  plot8 <- ggplot(df, aes(pt.adversaries, diff.norm, color=method)) + geom_smooth() + facet_grid(lambda_1_lab ~ lambda_2_lab) + 
    xlab("Adversarial fraction of network (forest fire)") + ylab("Bias in Estimated ATE / Estimated nonadversarial ATE") + 
    geom_abline(slope=0) + theme_bw() + theme(text = element_text(size = 15)) + 
    theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1))
  plot(plot8) 
  
  res$graph.type <- factor(res$graph.type, levels=c("SBM", "small-world", "scale-free", "forest-fire"))
  
  plot5 <- ggplot(subset(res, size.of.dom==FALSE), aes(pt.adversaries, adversary.influence, color=method)) + 
    geom_smooth() + xlab("Adversarial fraction of network") + ylab("Adversary influence") + 
    facet_wrap(~graph.type, scales="free_x") + theme_bw() + theme(text = element_text(size = 15)) + 
    theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1))
  plot(plot5)
  
  plot6 <- ggplot(subset(res, size.of.dom==FALSE), aes(pt.adversaries, adversary.influence, color=method)) + 
    geom_smooth() + xlab("Adversarial fraction of network") + ylab("Adversary influence") + 
    facet_wrap(~graph.type, scales="free_x") + theme_bw() + theme(text = element_text(size = 15)) + 
    theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(fill=NA))) + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1))
  plot(plot6)
  
}