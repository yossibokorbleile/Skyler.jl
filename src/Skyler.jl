#= License
Copyright 2019, 2020 (c) Yossi Bokor

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#

__precompile__()


module Skyler

#### Requirements ####
using CSV
using Hungarian
using DataFrames
using LinearAlgebra
using Distances
using NearestNeighbors
using PlotlyJS
using Optim
using SimpleGraphs
using Random
using Distributions
using QuadGK
using Plots

#### Exports ####
export 	skyler,
		unittest,
		Determine_Dimension_Coordinates,
		Plot_Partition,
		Plot_Graph,
		Determine_Dimension_Coordinates,
		Load_Example

#### Background Functions ####

function Spherical_Shell(ball_tree, point, inner, outer) # this obtains all the samples in the spherical shell around a chosen sample 'point'

	ball_1 = inrange(ball_tree, point, inner, true)
	ball_2 = inrange(ball_tree, point, outer, true)
	annulus = setdiff(ball_2, ball_1)
	return annulus
end


function Determine_Dimension_Coordinates(coordinates, ball_tree, index, sample_epsilon, angle_condition, inner_radius, outer_radius) # determine whether the `index'th point in the sample looks 0 or 1 dimensional given the coordinates of the other points, a pre-generated ball tree, and the correct algorithm parameters of 'angle_conditon', 'inner_radius', and 'outer_radius'
	point = coordinates[:,index]
	distance_matrix = pairwise(Euclidean(), coordinates, dims=2)
	node_neighbours = Spherical_Shell(ball_tree, point, 0, outer_radius)
	node_neighbours = push!(node_neighbours, index)
	G_point = IntGraph()

	for i in node_neighbours
		add!(G_point,i)
	end

	for i in vlist(G_point)
		for j in vlist(G_point)
			if distance_matrix[i,j] <= 2*sample_epsilon
				add!(G_point, i, j)
			end
		end
	end

	nc_1 = num_components(G_point)
	if nc_1 != 1 #if the graph on the points within a ball of the chosen points is disconnected, then the point is not close to a vertex
		return 1
	else # if the graph on the points in a ball is connected, then we need to check if 0 or 1 dimensional structure
		node_neighbours_1 = Spherical_Shell(ball_tree, point, inner_radius, outer_radius)
		G_point_1 = IntGraph() # construct graph on points in a spherical shell
		for i in node_neighbours_1
			add!(G_point_1, i)
		end
		for i in vlist(G_point_1)
			for j in vlist(G_point_1)
				if distance_matrix[i,j] <= 3*sample_epsilon
					add!(G_point_1, i, j)
				end
			end
		end
		nc_2 =num_components(G_point_1)
		if nc_2 == 2 #if the number of connected components is 2, we need to check if they are roughly opposite each other or not
			cc = collect(SimpleGraphs.components(G_point_1))
			mid_point_1 = []
			for l in 1:size(coordinates)[1]
				mid_point_1_l = 0
				cc_1 = collect(cc[1])
				for p in cc_1
					mid_point_1_l += coordinates[l, p]
				end
			append!(mid_point_1, mid_point_1_l/(length(cc_1)))
			end
			mid_point_2= []
			for l in 1:size(coordinates)[1]
				mid_point_2_l = 0
				cc_2 = collect(cc[2])
				for p in cc_2
					mid_point_2_l += coordinates[l,p]
				end
			append!(mid_point_2, mid_point_2_l/(length(cc_2)))
			end
			mid_point_1 = mid_point_1 .- coordinates[:,index]
			mid_point_2 = mid_point_2 .- coordinates[:,index]
			cos_of_angle = dot(mid_point_1, mid_point_2)/(norm(mid_point_1)*norm(mid_point_2))
			if cos_of_angle >= cos(2*acos(angle_condition))
				return 0
			else
				return 1
			end
		else #if the number of connected components is not 2, return dimension 0
			return 0
		end
	end
end

function Dimension_Split(coordinates, sample_epsilon, angle_condition, inner_radius, outer_radius) # for each point in the sample, determine if it looks 0 or 1 dimensional

	number_of_points = size(coordinates)[2]

	ball_tree = BallTree(coordinates)
	dimension_split = Dict{Integer,Set}(0 => Set(), 1 =>Set())
	for i in 1:number_of_points
		dim = Determine_Dimension_Coordinates(coordinates, ball_tree, i, sample_epsilon, angle_condition, inner_radius, outer_radius)
		dimension_split[dim]=push!(dimension_split[dim], i)
	end
	
	return dimension_split
end

function Find_Groups(distance_matrix, dictionary_of_points, dimension, connection_threshold) # obtain the clusters of points which look 'dimensional'  dimensional
	points = collect(dictionary_of_points[dimension]) #create a set of points
	G = SimpleGraph()
	for i in points
		add!(G, i)
	end
	for i in points
		for j in points
			if distance_matrix[i,j] <= connection_threshold #connect points below the connection threshold to cluster
				add!(G, i, j)
			end
		end
	end
	
	components = []
	
	vertices = vlist(G)
	n_v = length(vertices)
	checked = Set()
	
	discovered=[]
	function DFS(Graph, vertex)
		append!(discovered,vertex)
		for u in neighbors(Graph, vertex)
			if u in(discovered)
			else
				DFS(Graph, u)
			end
		end
	end
	
	for i in vertices
		if i in checked
		else
			push!(checked, i)
			discovered =[]
			DFS(G, i)
			append!(components, [discovered])
			for n in discovered
				push!(checked, n)
			end
		end
	end
		
	return components
end



function Generate_Clusters(distance_matrix, dimension_split, v_threshold, e_threshold) # generate the clusters of vertices and edges using the above functions

	v_comps = Find_Groups(distance_matrix, dimension_split, 0, v_threshold)
	v_clusters = [collect(v_comps[i]) for i in 1:length(v_comps)]
	number_of_vertices = length(v_clusters)
	e_comps = Find_Groups(distance_matrix, dimension_split, 1, e_threshold)
	e_clusters = [collect(e_comps[i]) for i in 1:length(e_comps)]
	number_of_edges = length(e_clusters)

	return v_clusters, number_of_vertices, e_clusters, number_of_edges
end


function Generate_Graph_Structure(distance_matrix, coordinates, clusters, epsilon) # obtain the abstract structure of the graph
	
	vertex_clusters = clusters[1]
	edge_clusters = clusters[3]
	number_edges=size(edge_clusters)[1]
	number_vertices = size(vertex_clusters)[1]
	boundaries = []
	for i in 1:number_vertices
		append!(boundaries, [i])
	end
	ball_tree = BallTree(coordinates)

	for i in 1:number_edges #for each edge, we need to find the two vertex clusters it goes between
		i_boundary = []
		if length(i_boundary) !=2
			for j in edge_clusters[i]
				near = inrange(ball_tree, coordinates[:,j], 3*epsilon, true)
				if issubset(near, edge_clusters[i]) == false
					vertex_candidates = setdiff(near, edge_clusters[i])
					for l in 1:number_vertices
						if l in i_boundary
						elseif issubset(vertex_candidates, vertex_clusters[l]) == true
							append!(i_boundary, l)
						end
					end
				end
			end
		end

		append!(boundaries, [i_boundary])
	end
	return boundaries
end


function Create_Partition_DataFrame(coordinates, clusters) #create a dataframe we can save the partition in
	
	dimension = size(coordinates,1)
	vertex_clusters = clusters[1]
	v_n =  clusters[2]
	edge_clusters = clusters[3]
	e_n = clusters[4]
	partition = Array{Float64}(undef, 0, dimension+2)
	
	for i in 1:v_n
		for j in vertex_clusters[i]
			line = [i 0]
			for k in 1:dimension
				line = hcat(line, [coordinates[k,j]])
			end
			partition= vcat(partition, line)
		end
	end
	
	for i in 1:e_n
		for j in edge_clusters[i]
			ind = i+v_n
			line = [ind 1]
			for k in 1:dimension
				line = hcat(line, [coordinates[k,j]])
			end
			partition= vcat(partition, line)
		end
	end
	
	return convert(DataFrame, partition)
end


function Generate_Partition(v_clusters, e_clusters) #this was used for modelling the embedding, unsure about future uses
	lengths = []
	partitions = []

	for i in 1:length(v_clusters)
		append!(lengths, length(v_clusters[i]))
		append!(partitions, v_clusters[i])
	end
	
	for i in 1:length(e_clusters)
		append!(lengths, length(e_clusters[i]))
		append!(partitions, e_clusters[i])
	end

	return partitions, lengths
end


function Plot_Partition(coordinates, split, clusters) #plot the partiton with nice colours
	
	vertex_clusters = clusters[1]
	v_n = clusters[2]
	edge_clusters = clusters[3]
	e_n = clusters[4]
	to_plot = DataFrame(x=[], y=[] , index=[])

	for i in 1:v_n
		#vc_i = collect(vertex_clusters[i])
		for j in vertex_clusters[i]
			push!(to_plot, [coordinates[1, j], coordinates[2, j], i])
		end
	end

	for i in 1:e_n
		for j in edge_clusters[i]
			push!(to_plot, [coordinates[1,j], coordinates[2,j], i+v_n])
		end
	end
	p=scatter(to_plot; x=:x, y=:y, mode="markers", marker_size=3, group=:index)
end


function Plot_Graph(number_of_vertices, Edges, vertex_locations) #plot a graph (either exact or modelled)
	
	G = IntGraph(number_of_vertices)
	locs = Dict(i => [out[1,i],out[2,i]] for i in 1:size(vertex_locations)[2])#1 =>[out[1,1], out[2,1]], 2=> [out[1,2], out[2,2]], 3=> [out[1,3], out[2,3]], 4=> [out[1,4], out[2,4]], 5=> [out[1,5], out[2,5]] )

	for i in 1:length(Edges)
		add!(G, Edges[i][1], Edges[i][2])
	end
	
	embed(G, locs)

	draw(G)
end

#### Functions for modelling the underlying structure ####

#=
To update:
-keeps sigma fixed at the moment, this can be made optional once sped up
-can make faster through analytical gradient by autodiff
-should be able to rid of numerical integration by use of error function
=#


#=
notns:
n		---	Number of data points
d		---	Dimension of data
N		---	Number of strata pieces
N_0		---	Number of 0-dim strata pieces
N_1		---	Number of 1-dim strata pieces
data	---	n by d array		---	Sample from embedded graph
mu		---	N_0 by d dim array	---	Iniitial vertex sample locations
sigma	---	N dim array			---	Initial error on each strata piece
E		---	n by N array		---	Initial assignment for each data point
S		---	N dim str array		---	Abstract graph vertex/edge names in that order
=#
function model_fit(data, mu, sigma, E, S, EM_it)

	##pre-processing
	#get size of different dimensional strata
	println(mu)
	N_0 = size(mu)[1]

	N_1 = size(S)[1] - N_0

	N = N_0 + N_1

	n = size(data)

	d = size(mu)[2]
	
	#load parameters
	Pi = sum(E, dims = 1)./n
	MU		=	Dict{String,Array}(S[i] => mu[i,:] for i = 1:N_0)
	SIGMA	=	Dict{Any,Float64}(S[i] => sigma[i] for i = 1:N)
	PI		=	Dict{Any,Float64}(S[i] => Pi[i] for i = 1:N)


	function reload_param(mu_t,sigma_t,Pi_t; S_t = S, N_0_t = N_0, N_t = N)
		MU_t	=  	Dict{String,Array}(S_t[i] => mu_t[i,:] for i = 1:N_0_t)
		SIGMA_t	=  	Dict{Any,Float64}(S_t[i] => sigma_t[i] for i = 1:N_t)
		PI_t	=  	Dict{Any,Float64}(S_t[i] => Pi_t[i] for i = 1:N_t)
		return MU_t, SIGMA_t, PI_t
	end

	##density functions
	#funciton for gaussian
	function gaussian(x,v,s)
		linear = x - v
		return exp(-LinearAlgebra.dot(linear,linear)/(2*s))
	end

	#function for gaussian line
	function gaussian_line(x,t,v1,v2,s)
		linear = x - (t*v1 + (1-t)*v2)
		return exp(-LinearAlgebra.dot(linear,linear)/(2*s))
	end


	#function for vertex density
	function p_e(x,v,s)
		d = size(x)[1]
		return gaussian(x,v,s)/((2*pi*s)^(d/2))
	end

	#funciton for line density
	function p_l(x,v1,v2,s)
		d = size(x)[1]
		return quadgk(t -> gaussian_line(x,t,v1,v2,s), 0, 1)[1]/((2*pi*s)^(d/2))
	end

	#log function for vertex densities
	function log_p_e(x,v,s)
		d = size(x)[1]
		linear = x - v
		return -(d/2)*log(2*pi*s) - LinearAlgebra.dot(linear,linear)/(2*s)
	end

	#log function for line density
	function log_p_l(x,v1,v2,s)
		return log(p_l(x,v1,v2,s))
	end

	##Expectation Calculation
	#make density functions
	function densities(MU_t,SIGMA_t; S_t = S, N_0 = N_0 , N_1 = N_1, N = N)
		zero_dim = [x -> p_e(x,MU_t[strat],SIGMA_t[strat]) for strat = S_t[1:N_0]]
		one_dim  = [x -> p_l(x,MU_t[strat[1]],MU_t[strat[2]],SIGMA_t[strat]) for strat = S_t[N_0+1:N]]
		return vcat(zero_dim,one_dim)
	end

	function log_densities(mu_t,SIGMA_t; S_t = S, N_0 = N_0 , N_1 = N_1, N = N)
		MU_t    =   Dict{String,Array}(S_t[i] => mu_t[i,:] for i = 1:N_0)
		zero_dim = [x -> log_p_e(x,MU_t[strat],SIGMA_t[strat]) for strat = S_t[1:N_0]]
		one_dim  = [x -> log_p_l(x,MU_t[strat[1]],MU_t[strat[2]],SIGMA_t[strat]) for strat = S_t[N_0+1:N]]
		return vcat(zero_dim,one_dim)
	end


	#form density matrix
	function density_matrix(data,density_functions)
		return hcat([ (broadcast(rho, data))  for rho in density_functions ]...)
	end


	function expectation(data,MU_t,SIGMA_t,Pi_t)
		DM = density_matrix(data,densities(MU_t,SIGMA_t))
		E = Pi_t.*DM
		E = E./sum(E, dims = 2)
		return E
	end

	##Cost function
	function Cost(data,mu_t,SIGMA_t, Pi, E)
		log_DM = density_matrix(data,log_densities(mu_t,SIGMA_t))
		return sum(E.*(log_DM.+log.(Pi)))
	end
	
	##EM algorithm
	for it = 1:EM_it
		#E step
		E = expectation(data,MU,SIGMA,Pi)
		#M step
		Pi = sum(E, dims = 1)./n
		C = t -> -Cost(data,t,SIGMA, Pi, E)
		results = Optim.optimize(C,mu)
		mu = Optim.minimizer(results)
		println("Current vertices:")
		println(mu)
		println("Current Cost:")
		println(C(mu))
		#Update parameters
		MU, SIGMA, PI = reload_param(mu ,sigma ,Pi ; S_t = S, N_0_t = N_0, N_t = N)
	end
	return MU, SIGMA, PI
end



#### Main functions and wrappers ####
function skyler(points, sample_epsilon, angle_condition, inner_radius, outer_radius, vertex_threshold, edge_threshold; out = "model", EM_it = 20) #wrapper for obtaining the partition
	dimension = size(points,1)
	distances = pairwise(Euclidean(), points, dims=2)
	dim_split = Dimension_Split(points, sample_epsilon, angle_condition, inner_radius, outer_radius) #partition the sample
	clusters = Generate_Clusters(distances, dim_split, vertex_threshold, edge_threshold) #find the clusters in each dimensions
	boundaries = Generate_Graph_Structure(distances, points, clusters, sample_epsilon) #obtain the boundary relations
	partition = Create_Partition_DataFrame(points, clusters) #create the data frame which contains the partition of the points
	if out == "struct"
		return partition, boundaries
	elseif out == "test"
		return partition, boundaries, clusters
	elseif out == "model"
		number_vertices = clusters[2]
		vertex_guess = Array{Float64}(undef, dimension,0)
		for i in 1:number_vertices
			println("vertx $i has the following ", length(clusters[1][i])," points: ", clusters[1][i])
			v_i = zeros(dimension, 1)
			for j in clusters[1][i]
			println(j)
				v_i .+= points[:, j]
				println(v_i)
			end
			v_i = v_i/(length(clusters[1][i]))
			println(v_i)
			vertex_guess = hcat(vertex_guess, v_i)
		end
		println("vertex guess is:",vertex_guess)
		E = zeros(size(points,2), length(boundaries))
		
		for i in 1:number_vertices
			for j in clusters[1][i]
				E[j,i] =1
			end
		end
		
		for i in (number_vertices+1):length(boundaries)
			for j in clusters[3][i-number_vertices]
				E[j,i] =1
			end
		end
		
		sigma = [sample_epsilon for i in 1:length(boundaries)]
		
		S = []
		for i in 1:length(boundaries)
			if typeof(boundaries[i]) == Int64
				j = boundaries[i]
				push!(S, "$j")
			else
				J = []
				for k in 1:length(boundaries[i])
						push!(J, "$k")
					end
				push!(S, tuple(J...,))
			end
		end
		
		data = vcat([points[:,i] for i in 1:size(points,2)])
		
		v_guess = Array{Float64}(undef, size(vertex_guess,2), size(vertex_guess,1))
		
		for i in 1:size(vertex_guess,2)
			for j in 1:size(vertex_guess,1)
				v_guess[i,j] = vertex_guess[j,i]
			end
		end
		
		MU, SIGMA, PI = model_fit(data, v_guess, sigma, E, S, EM_it)
		mu_ans = collect(values(MU))
		mu_ans = vcat([datum' for datum in mu_ans]...)
		return mu_ans', boundaries
	else
		println("Error: output option 'out' selected is not valid, please choose from 'model' to obtain a model of the underlying structure, or 'struct' to receive a DataFrame indicating which strata a point as been associated with and the dimension of this strata, and the boundary relations.")
		return []
	end
end
	
#### Functions for examples ####
function Load_Example(i)
	dir = pwd()
	cd(@__DIR__)
	cd("..")
	cd("Examples")
	points = CSV.read("Sample-$i.csv", header=false)
	cd(dir)
	return convert(Array,points)
end

#### unittest ####

function test_1()
dir = pwd()
	cd(@__DIR__)
	cd("..")
	cd("Examples")
	
	points= CSV.read("Sample-1.csv", header=false)
	structure = skyler(points, 0.01, 0.08, 0.1, 1/4, 0.03, 0.09, out="test")
	saved_partition = CSV.read("Partition-1.csv", header=false)
	
	for i in 1:4
		if structure[1][!,i] != saved_partition[!,i]
			println("Error: test_1, $i, structure[1] = ")
			return structure[1]
		end
	end
	cd(dir)
	return []
end


function test_2()

	points = [0.0 1 2; 0 1 2]
	
	structure = skyler(points,0.01, 0.08, 0.1, 1/4, 0.03, 0.09, out = "test")
	
	if partition[3][1] == [[1],[2],[3]]
		return []
	else
		println("Error: test_2, structure[3][1] = ")
			return structure[3][1]
	end

end



function unittest()
	x = Array{Any}(undef, 2)
	
	x[1] = test_1() #expected answer empty
	x[2] = test_2() #expected answer empty
	
	for p 	= 	1:length(x)
		if !isempty(x[p])
			println(p)
			return x
		end
	end
	return []
end



end #module
