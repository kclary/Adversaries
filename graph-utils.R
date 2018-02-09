library(igraph)

build.graph.params <- function(configs, i) { 
  if(as.character(configs[i,]$graph.type) == "barabasi-albert") graph.params <- barabasi.params(configs[i,]$size, configs[i,]$power)
  if(as.character(configs[i,]$graph.type) == "small-world") graph.params <- sw.params(configs[i,]$size, configs[i,]$degree, configs[i,]$p)
  if(as.character(configs[i,]$graph.type) == "sbm") graph.params <- sbm.params(configs[i,]$size, configs[i,]$mu)
  if(as.character(configs[i,]$graph.type) == "forest-fire") graph.params <- ff.params(configs[i,]$size, configs[i,]$fw, configs[i,]$bw)
  if(as.character(configs[i,]$graph.type) %in% c("polyblogs", "facebook")) { 
    graph.params <- list()
    graph.params$graph.type <- as.character(configs[i,]$graph.type)
  }  
  
  return(graph.params)
}

get.graph.properties <- function(g) { 
  graph.properties <- list()
  graph.properties$g <- g
  graph.properties$degrees <- degree(g)
  graph.properties$n <- length(V(g))
  graph.properties$adj <- get.adjacency(g, sparse=FALSE)
  
  graph.properties$degree.inv <- solve(diag(graph.properties$degrees))
  graph.properties$transition <-  degree.inv %*% adj 
  
  return(graph.properties)
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
  graph.params$degree <- degree
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

ff.params <- function(n, fw, bw) { 
  graph.params <- list() 
  graph.params$graph.type <- "forest-fire"
  graph.params$n <- n
  graph.params$forward.prob <- fw
  graph.params$backward.prob <- bw
  return(graph.params)
}

generate.graph <- function(graph.params) { 
  graph.type <- graph.params$graph.type
  
  if(graph.type == "small-world") {
    g <- watts.strogatz.game(1, graph.params$n, graph.params$degree, graph.params$p)
  }
  
  if(graph.type == "barabasi-albert") { 
    if(graph.params$n == 500) add.edges <- rbinom(graph.params$n, 1, 0.1)+4
    if(graph.params$n == 1000) add.edges <- rbinom(graph.params$n, 1, 0.1)+8
    if(graph.params$n == 5000) add.edges <- rbinom(graph.params$n, 1, 0.1)+37
    
    g <- barabasi.game(graph.params$n, graph.params$power, out.seq=add.edges, directed=FALSE)    
  }
  
  if(graph.type == "sbm") { 
    edg <- read.csv(paste0("binary_networks/nets/sbm-", graph.params$n, "-", graph.params$mu, "-", graph.params$ind, "-adj.txt"), sep="\t", header=FALSE)
    edg <- as.matrix(edg)
    g <- graph_from_adjacency_matrix(edg, mode="undirected")
    
  }
  
  if(graph.type == "forest-fire") { 
    bw.factor <- graph.params$backward.prob/graph.params$forward.prob
    
    med_edges <- 0
    if(graph.params$n == 500) med_edges <-2020
    if(graph.params$n == 1000) med_edges <- 7650
    if(graph.params$n == 5000) med_edges <- 186825
    
    edge_count <- 0
    while(edge_count < (med_edges-med_edges*.1) | edge_count > (med_edges+med_edges*.1)) { 
      g <- forest.fire.game(n=graph.params$n, fw.prob = graph.params$forward.prob, bw.factor = bw.factor, directed = FALSE)  
      edge_count <- sum(get.adjacency(g))/2
    }
  }
  
  if(graph.type == "polyblogs") { 
    g <- read_graph(file = "graph-data/polblogs/polblogs.gml", format = "gml")   
    cl <- clusters(g)
    g <- induced_subgraph(g, which(cl$membership == which.max(cl$csize)))
  }
  
  if(graph.type == "facebook") { 
    edg <- read.csv(paste0("graph-data/facebook/facebook_combined.txt"), sep=" ", header=FALSE)
    edg <- as.matrix(edg)+1
    
    graph.params$n <- max(edg)
    adj <- matrix(0, graph.params$n, graph.params$n)
    for(i in 1:dim(edg)[1]) adj[edg[i,1],edg[i,2]] <- 1
    g <- graph_from_adjacency_matrix(adj)
  }
  
  return(g)
}

check.dominating.set <- function(graph.properties, adversaries) { 
  adv <- as.matrix(adversaries)
  adj.adversaries <- (graph.properties$adj %*% t(adv)) + t(adv)
  return(sum(rowSums(adj.adversaries) > 0) == dim(graph.properties$adj)[1])
}

generate.clusters <- function(g, clustering) { 
  return(cluster_infomap(g)$membership)
}

dominate.greedy <- function(graph.properties,weight=NULL,proportion=1.0) {
  A <- graph.properties$adj
  od <- degree(graph.properties$g,mode="out")+1
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

test.sw <- function() { 
  sw.params(1000, 3, 0.02)
}

test.sbm <- function() {
  sbm.params(1000, 0.2)
}

test.barabasi <- function() { 
  barabasi.params(1000, 0.3)
}

test.ff <- function() { 
  ff.params(1000, 0.37, 0.25)  
}

test.graph.properties <- function(graph.params) { 
  num.edges <- lapply(1:100, function(x) { 
    graph.params$ind <- x
    g <- generate.graph(graph.params)  
    return(sum(get.adjacency(g))/2)
  })
  print(median(unlist(num.edges)))
}

test.sbm.edges <- function() { 
  for(i in c(500,1000,5000)) { 
    for(j in c(0.1, 0.2, 0.3)) { 
      graph.params <- sbm.params(i,j)  
      cat(paste("sbm", i, j))
      test.graph.properties(graph.params)
    }  
  }
}

test.barabasi.edges <- function() { 
  for(i in c(500,1000,5000)) { 
    for(j in c(0.1, 0.3, 0.5)) { 
      graph.params <- barabasi.params(i,j)  
      cat(paste("barabasi", i, j))
      test.graph.properties(graph.params)
    }  
  }
}

test.sw.edges <- function() { 
  for(i in c(500, 1000, 5000)) { 
    z <- 0
    if(i == 500) z <- 4
    if(i == 1000) z <- 8
    if(i == 5000) z <- 37
    for(j in c(0.03, 0.05, 0.1)) { 
      graph.params <- sw.params(i, z, j)  
      cat(paste("sw", i, z, j))
      test.graph.properties(graph.params)
    }  
  }
}

test.ff.edges <- function() { 
  for(i in c(500, 1000, 5000)) { 
    fw <- 0
    bw <- 0
    
    if(i == 500) { 
      fw <- 0.36 
      bw <- 0.34
    }  
    if(i == 1000) { 
      fw <- 0.36 
      bw <- 0.365
    } 
    if(i == 5000) { 
      fw <- 0.36 
      bw <- 0.365
    } 
    graph.params <- ff.params(i, fw, bw)  
    cat(paste("ff", i, fw, bw))
    test.graph.properties(graph.params)
  }
}
