using Graphs

"Generate stochastic block models under various configurations"
function run_sbm_experiments()
	graph_sizes = [500,1000,5000]
	community_sizes = [0,1]
	community_sizes = [1]
	mus = collect(0.1:0.1:0.3)
	num_samples = 100

	for i=1:size(graph_sizes)[1]
		for j=1:size(community_sizes)[1]
			for m=size(mus)[1]
				for x=1:num_samples
					n = graph_sizes[i]
					community_sizes[j] == 0 ? (minc = 10; maxc = 50;) : (minc = 20; maxc = 100;) 
					mu = mus[m]

					netfile = 
					network = create_network(n, mu, minc, maxc, x)
					G = network[1]
					ground_truth = network[2]
					k = maximum(ground_truth[:,2])

					# get assignments from graph kernel clustering
					# calculate nmi between assignments, ground truth 
				end
			end
		end
	end
end

# "Generate and return network with n nodes from file"
# function create_network(n, mu, minc, maxc)
# 	run(`./benchmark -N $n -k 20 -maxk 50 -mu $mu -minc $minc -maxc $maxc`)

# 	adj = int(open(readdlm, "network.dat"))
# 	n = maximum(adj[:,1])
# 	num_edges = size(adj)[1]
	
# 	g = simple_inclist(n)

# 	for i=1:num_edges
# 		neighbors = adj[i,:]
# 		add_edge!(g, neighbors[1], neighbors[2])
# 	end

# 	ground_truth = int(open(readdlm, "community.dat"))

# 	return g, ground_truth
# end


#=Generate and return network with community ground truth
int n: number of nodes"
float mu: inter-community mixing parameter
int minc: minimum community size
int maxc: maximum community size=#
function create_network(n, mu, minc, maxc, offset, netfile)
	"generate network with SBM benchmark suite from Lancichinetti, Fortunato (2009)"
	readall(`binary_networks/./benchmark -N $n -k 20 -maxk 50 -mu $mu -minc $minc -maxc $maxc -netfile $netfile -comfile $comfile`)
		new_f = string("binary_networks/", netfile)
		"move files to generating directory to keep things clean"
		run(`mv $netfile $new_f`)
	end
end
