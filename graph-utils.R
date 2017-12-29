library(igraph)

build.graph.params <- function(configs, i) { 
  if(as.character(configs[i,]$graph.type) == "barabasi-albert") graph.params <- barabasi.params(configs[i,]$size, configs[i,]$power)
  if(as.character(configs[i,]$graph.type) == "small-world") graph.params <- sw.params(configs[i,]$size, configs[i,]$degree, configs[i,]$p)
  if(as.character(configs[i,]$graph.type) == "sbm") graph.params <- sbm.params(configs[i,]$size, configs[i,]$mu)
  if(as.character(configs[i,]$graph.type) %in% c("polyblogs", "facebook")) { 
    graph.params <- list()
    graph.params$graph.type <- as.character(configs[i,]$graph.type)
  }  
  
  return(graph.params)
}

barabasi.params <- function(n, power) { 
  graph.params <- list()
  graph.params$graph.type <- "barabasi-albert"
  graph.params$n <- n
  graph.params$power <- power
  
  return(graph.params)
}

sw.params <- function(n, degree, p) { 
  graph.params <- list()
  graph.params$graph.type <- "small-world"
  graph.params$degree <- 5
  graph.params$n <- n
  graph.params$p <- p
  
  return(graph.params)
}

sbm.params <- function(n, mu) { 
  graph.params <- list() 
  graph.params$graph.type <- "sbm"
  graph.params$n <- n
  graph.params$mu <- mu
  
  return(graph.params)
}

generate.graph <- function(graph.params) { 
  graph.type <- graph.params$graph.type
  
  if(graph.type == "small-world") {
    g <- watts.strogatz.game(1, graph.params$n, graph.params$degree, graph.params$p)
  }
  
  if(graph.type == "barabasi-albert") { 
    g <- barabasi.game(graph.params$n, graph.params$power, directed=FALSE)    
  }
  
  if(graph.type == "sbm") { 
    edg <- read.csv(paste0("binary_networks/nets/network-", graph.params$ind, "-", substr(as.character(graph.params$mu), 2, 3), ".dat"), sep="\t")
    adj <- matrix(0, graph.params$n, graph.params$n)
    for(i in 1:dim(edg)[1]) adj[edg[i,1],edg[i,2]] <- 1
    g <- graph_from_adjacency_matrix(adj, mode="undirected")
    
  }
  
  if(graph.type == "polyblogs") { 
    g <- read_graph(file = "graph-data/polblogs/polblogs.gml", format = "gml")   
    cl <- clusters(g)
    g <- induced_subgraph(g, which(cl$membership == which.max(cl$csize)))
  }
  
  # if(graph.type == "facebook") { 
  #   g <- read_graph(file = "graph-data/", format = "gml")   
  #   cl <- clusters(g)
  #   g <- induced_subgraph(g, which(cl$membership == which.max(cl$csize)))
  # }
  return(g)
}

check.dominating.set <- function(g, adversaries) { 
  adj <- get.adjacency(g, sparse=FALSE)
  adv <- as.matrix(adversaries)
  adj.adversaries <- (adj %*% t(adv)) + t(adv)
  return(sum(rowSums(adj.adversaries) > 0) == dim(adj)[1])
}

generate.clusters <- function(g, clustering) { 
  return(cluster_infomap(g)$membership)
}

dominate.greedy <- function(g,weight=NULL,proportion=1.0) {
  A <- get.adjacency(g, sparse=FALSE)
  od <- degree(g,mode="out")+1
  S <- NULL
  diag(A) <- 0
  n <- nrow(A)
  covered <- rep(0,n)
  
  while(sum(covered)<n*proportion){
    i <- which.max(od)
    cands <- which(od==od[i])
    i <- cands[sample(length(cands), 1)]
    
    covered[A[i,]>0] <- 1
    covered[i] <- 1
    S <- c(S,i)
    A[,covered>0] <- 0
    h <- graph.adjacency(A,mode="directed")
    od <- degree(h,mode="out")+1-covered
  }
  S
}

dominate.greedy.inf <- function(g,weight=NULL,proportion=1.0) {
  A <- get.adjacency(g, sparse=FALSE)
  od <- degree(g,mode="out")
  degree.inv <- diag(ifelse(od > 0, 1/od, 0))
  od <- colSums(degree.inv %*% A)
  
  S <- NULL
  diag(A) <- 0
  n <- nrow(A)
  covered <- rep(0,n)
  while(sum(covered)<n*proportion){
    i <- which.max(od)
    cands <- which(od==od[i])
    i <- cands[sample(length(cands), 1)]
    
    covered[A[i,]>0] <- 1
    covered[i] <- 1
    S <- c(S,i)
    
    od <- degree(g,mode="out")
    degree.inv <- diag(ifelse(od > 0, 1/od, 0))
    trans <- degree.inv %*% A
    trans[,S] <- 0
    trans[S,] <- 0
    od <- colSums(trans)
  }
  S
}