using Graphs

function create_graphs(i)
	graph_params = readdlm("all_graph_configurations.csv", ',')

	n = graph_params[i,9]
	graph_type = graph_params[i,2]

	if graph_type == "sbm"
		for trial=1:100
			mu = graph_params[i,6]

			label = "$graph_type-$n-$mu-$trial"
			netfile = "nets/$label-adj.txt"
			if !isfile(netfile)
				println(label)
				create_network(n, mu, netfile)
			end
			while !isfile(netfile) end
		end
	end
end

#=Generate and return network with community ground truth
int n: number of nodes"
float mu: inter-community mixing parameter=#
function create_network(n, mu, netfile)
	"generate network with SBM benchmark suite from Lancichinetti, Fortunato (2009)"
	avgd = ceil(10/1000*n)
	maxd = (100/1000*n)
	readall(`./benchmark -N $n -k $avgd -maxk $maxd -mu $mu`)

	adj = round(Int64, open(readdlm, "network.dat"))
	n = maximum(adj)
	num_edges = size(adj,1)

	g = simple_graph(n, is_directed=false)

	for i=1:num_edges
		neighbors = adj[i,:]
		if neighbors[1] != neighbors[2]
			add_edge!(g, neighbors[1], neighbors[2])
		end
	end

	adj = adjacency_matrix(g) + 0
	writedlm(netfile, adj, "\t")
	readall(`rm network.dat`)
	readall(`rm community.dat`)
end


t = readdlm("all_graph_configurations.csv", ',')
for i=1:length(t[:,1])
	create_graphs(i)
end