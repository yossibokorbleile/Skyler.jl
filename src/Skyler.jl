#= License
Copyright 2019, 2020, 2021 (c) Yossi Bokor

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#

__precompile__()

module Skyler

#### Requirements ####
using CSV
using DataFrames
using LinearAlgebra
using NearestNeighbors

#=
using Distributed
@everywhere using CSV
@everywhere using DataFrames
@everywhere using LinearAlgebra
@everywhere using NearestNeighbors
=#

using Distances
using Optim
using SimpleGraphs
using Random
using Distributions
using QuadGK
using Plots
using JLD
using SpecialFunctions
using SharedArrays
using Calculus
using NNlib
using FiniteDiff
using ProgressMeter

#### Exports ####

export 	skyler,
		unittest,
		Determine_Dimension,
		Plot_PPC,
		Plot_Graph,
		Load_Example,
		PartitionedPointCloud,
		StratifiedSpace,
		EmbeddedStratifiedSpace,
		Dimension_Split,
		Create_Partition_DataFrame
		
#### Structures ####

mutable struct PartitionedPointCloud
	Points::AbstractArray
	Strata::Dict{Int, Set}
	Dimensions::Dict{Int, Set}
	Boundaries::AbstractArray
end

mutable struct StratifiedSpace
	StrataDimensions::Dict{Int, Set}
	Boundaries::AbstractArray
end

mutable struct EmbeddedStratifiedSpace
	StratSpace::StratifiedSpace
	Equations::Dict{Int, Array}
end

#### Background Functions ####

function Spherical_Shell(ball_tree::BallTree, point::Array, inner::T, outer::U) where {T<: Number, U<: Number} #this obtains all the samples in the spherical shell around a chosen sample 'point'
	ball_1 = inrange(ball_tree, point, inner, true)
	ball_2 = inrange(ball_tree, point, outer, true)
	annulus = setdiff(ball_2, ball_1)
	return annulus
end


function Determine_Dimension(sample::A, ball_tree::BallTree, index::Int, sample_epsilon::P, radius::R) where {P<:Number,R<:Number, A<:Union{AbstractArray, LinearAlgebra.Adjoint}}# determine whether the `index'th point in the sample looks 0 or 1 dimensional given the coordinates of the other points, a pre-generated ball tree, and the correct algorithm parameters of 'angle_conditon', 'inner_radius', and 'outer_radius'
	point = sample[:,index]
	distance_matrix = pairwise(Euclidean(), sample, dims=2)
	node_neighbours = Spherical_Shell(ball_tree, point, 0., radius+sample_epsilon)
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
	n_p = size(sample,2) #number of samples
	nc_1 = num_components(G_point)
	if nc_1 != 1 #if the graph on the points within a ball of the chosen points is disconnected, then the point is not close to a vertex
		return 1
	else # if the graph on the points in a ball is connected, then we need to check if 0 or 1 dimensional structure
		node_neighbours_1 = Spherical_Shell(ball_tree, point, radius-1sample_epsilon, radius+sample_epsilon)
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
			for l in 1:size(sample,1)
				mid_point_1_l = 0
				cc_1 = collect(cc[1])
				for p in cc_1
					mid_point_1_l += sample[l, p]
				end
			append!(mid_point_1, mid_point_1_l/(length(cc_1)))
			end
			mid_point_2= []
			for l in 1:size(sample,1)
				mid_point_2_l = 0
				cc_2 = collect(cc[2])
				for p in cc_2
					mid_point_2_l += sample[l,p]
				end
			append!(mid_point_2, mid_point_2_l/(length(cc_2)))
			end
			mid_point_1 = mid_point_1 .- sample[:,index]
			mid_point_2 = mid_point_2 .- sample[:,index]
			#cos_of_angle = dot(mid_point_1, mid_point_2)/(norm(mid_point_1)*norm(mid_point_2))
			if dot(mid_point_1, mid_point_2) <= -radius^2 +2*radius*sample_epsilon +7*sample_epsilon^2
				return 1
			else
				return 0
			end
		else #if the number of connected components is not 2, return dimension 0
			return 0
		end
	end
end


function Dimension_Split(sample::A, sample_epsilon::S, radius::U) where {S<: Number, U<:Number, A<:AbstractArray}# for each point in the sample, determine if it looks 0 or 1 dimensional
	#if nprocs() == 1
		n_p = size(sample,2)
		ball_tree = BallTree(sample)
		dimension_split = Dict{Int,Set}(0 => Set(), 1 =>Set())
		println("Calculating local structure for each sample: ")
		@showprogress for i in 1:n_p
			dim = Determine_Dimension(sample, ball_tree, i, sample_epsilon,radius)
			dimension_split[dim]=push!(dimension_split[dim], i)
		end
		return dimension_split
	#else
	#	n_p = size(sample,2)
	#	@everywhere ball_tree=BallTree(sample)
	#end
end

function Find_Groups(distance_matrix::A, dictionary_of_points::Dict{Int, Set}, dimension::T, connection_threshold::W) where {T <: Number, W<: Number, A<:AbstractArray} # obtain the clusters of points which look 'dimensional'  dimensional
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



function Generate_Strata_Assignments(distance_matrix::Array, dimension_split::Dict{Int, Set}, v_threshold::T, e_threshold::W) where {T <: Number, W<:Number} # generate the clusters of vertices and edges using the above functions
	v_comps = Find_Groups(distance_matrix, dimension_split, 0, v_threshold)
	number_of_vertices = length(v_comps)
	e_comps = Find_Groups(distance_matrix, dimension_split, 1, e_threshold)
	number_of_edges = length(e_comps)
	change = []
	for i in 1:number_of_edges
		if size(e_comps[i]) == 1
			append!(change, e_comps[i][1])
		end
	end
	if length(change) !=0
		d_s = copy(dimension_split)
		for i in 1:length(change)
			d_s[0] = push!(d_s[0], change[i])
			d_s[1] = setdiff(d_s[1], Set([change[i]]))
		end
		v_comps = Find_Groups(distance_matrix, d_s, 0, v_threshold)
		number_of_vertices = length(v_comps)
		e_comps = Find_Groups(distance_matrix, d_s, 1, e_threshold)
		number_of_edges = length(e_comps)
	end
	strata_assignments = Dict{Int, Set}(1=>Set())
	strata_dimensions = Dict{Int,Set}(0=> Set(), 1=>Set())
	#println("Detected ", number_of_vertices," vertices and ", number_of_edges, " edges.")
	for i in 1:number_of_vertices
		strata_assignments[i] =Set(v_comps[i])
		strata_dimensions[0] =push!(strata_dimensions[0],i)
	end
	for i in 1:number_of_edges
		strata_assignments[number_of_vertices+i] =Set(e_comps[i])
		strata_dimensions[1] =push!(strata_dimensions[1],number_of_vertices+i)
	end
	return strata_assignments, strata_dimensions
end


function Generate_Abstract_Structure(points::A, dimension_split::Dict{Int, Set}, sample_epsilon::T, vertex_threshold::T, edge_threshold::U) where {U<:Number, T<:Number, A<:AbstractArray} # obtain the abstract structure of the graph
	dists = pairwise(Euclidean(), points, dims=2)
	corrected = false
	while corrected == false
		strata_assignments, strata_dimensions = Generate_Strata_Assignments(dists, dimension_split, vertex_threshold, edge_threshold)
		max_dim = maximum(keys(strata_dimensions))
		n_s = maximum(keys(strata_assignments))
		boundaries = zeros(n_s, n_s)
		@showprogress for i in 0:max_dim-1
			for j in strata_dimensions[i]
				s_j = points[:,collect(strata_assignments[j])]
				for k in strata_dimensions[i+1]
					s_k = points[:,collect(strata_assignments[k])]
					d_j_k = Distances.pairwise(Euclidean(), s_j, s_k, dims=2)
					if minimum(d_j_k) <= 3*sample_epsilon
						boundaries[j,k]=1
					end
				end
			end
		end
		println("\r The strata are ", strata_dimensions, " and the boundary matrix is  ", boundaries,".")
		change = []
		if isempty(strata_dimensions[1])
			corrected = true
			StratSpace = StratifiedSpace(strata_dimensions, boundaries)
			return StratSpace, strata_assignments
			break
		else
			for i in strata_dimensions[1]
				if sum(boundaries[:,i]) !=2
					append!(change, i)
				end
			end
		end
		if length(change) == 0
			corrected = true
			StratSpace = StratifiedSpace(strata_dimensions, boundaries)
			return StratSpace, strata_assignments
			break
		else
			println("Stray edge detected, addressing it.")
			for k in 1:length(change)
				dimension_split[0] = union(dimension_split[0], strata_assignments[change[k]])
				dimension_split[1] = setdiff(dimension_split[1], strata_assignments[change[k]])
			end
		end
	end
end



function Partition_Point_Cloud(points::A, sample_epsilon::P, radius::W, vertex_threshold::T, edge_threshold::U)::PartitionedPointCloud where {P<:Number, Q<:Number, R<:Number, W<:Number, T <: Number, U<:Number, A<:AbstractArray}
	dimension = size(points,1)
	distances = pairwise(Euclidean(), points, dims=2)
	dim_split = Dimension_Split(points, sample_epsilon, radius) #partition the sample into 0 or 1 dimensional using a dictionary
	#strata_assignments, strata_dimensions = Generate_Strata_Assignments(distances, dim_split, vertex_threshold, edge_threshold) #find the strata for the partioned point cloud
	StratSpace, strata_assignments = Generate_Abstract_Structure(points, dim_split, sample_epsilon, vertex_threshold, edge_threshold) #obtain the boundary relations
	change = []
	ppc = PartitionedPointCloud(points, strata_assignments, StratSpace.StrataDimensions, StratSpace.Boundaries)
	return ppc
end
#=
Deprecated this was used for modelling the embedding, unsure about future uses
function Generate_Partition(v_clusters, e_clusters)
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
=#

function Create_Partition_DataFrame(ppc::PartitionedPointCloud) #create a dataframe we can save the partitioned point cloud in
	dimension = size(ppc.Points,1)
	df = DataFrame(ppc.Points[:,1]')
	insert!(df, 1, [1], :index)
	insert!(df, 2, [1], :dimension)
	deleterows!(df, 1)
	for i in 0:maximum(keys(ppc.Dimensions))
		for j in ppc.Dimensions[i]
			for k in ppc.Strata[j]
				l = [j i]
				for n in 1:dimension
					l=hcat(l, ppc.Points[n,k])
				end
				push!(df, l)
			end
		end
	end
	return sort(df)
end

function Plot_PPC(ppc::PartitionedPointCloud) #plot the partiton with nice colours
	Plots.pyplot()
	if size(ppc.Points,1) == 2
		df = DataFrame(index=Int[], dimension=Int[], c1=[], c2=[])
		for i in 0:maximum(keys(ppc.Dimensions))
			for j in ppc.Dimensions[i]
				for k in ppc.Strata[j]
					push!(df, [j,i,ppc.Points[1,k], ppc.Points[2,k]])
				end
			end
		end
		Plots.scatter( df[!,:c1],df[!,:c2], group=df[!,:index])
	else
		println("Points live in 3 or more dimensions, projecting onto the first 3.")
		df = DataFrame(index=Int[], dimension=Int[], c1=[], c2=[], c3=[])
		for i in 0:maximum(keys(ppc.Dimensions))
			for j in ppc.Dimensions[i]
				for k in ppc.Strata[j]
					push!(df, [j,i,ppc.Points[1,k], ppc.Points[2,k], ppc.Points[3,k]])
				end
			end
		end
		Plots.scatter3d( df[!,:c1],df[!,:c2], df[!,:c3], group=df[!,:index])
	end
end

function Plot_PPC(partition::DataFrame)
	Plots.pyplot()
	dimension = size(partition[!,:point][1],1)
	if dimension == 2
		points = Array{Float64}(undef, 0,2)
		for i in 1:size(partition,2)
			points = hcat(points, partition[i,:point])
		end
		Plots.scatter( points[:,1], points[:,2], group=partition[!,:index])
	else
		println("Points live in 3 or more dimensions, projecting onto the first 3.")
		points = Array{Float64}(undef, 0,3)
		for i in 1:size(partition,1)
			points = vcat(points, partition[i,:point][1:3]')
		end
		Plots.scatter3d( points[:,1], points[:,2], points[:,3], group=partition[!,:index])
	end
end

function Plot_PPC(partition::DataFrame, c_1::Int, c_2::Int)
	Plots.pyplot()
	println("Projecting on to the following coordinates: $c_1 and $c_2.")
	points = Array{Float64}(undef, 0,2)
	for i in 1:size(partition,1)
		points = vcat(points, [partition[i,:point][c_1] partition[i,:point][c_2]])
	end
	Plots.scatter( points[:,1], points[:,2], group=partition[!,:index])
end

function Plot_PPC(partition::DataFrame, c_1::Int, c_2::Int, c_3::Int)
	Plots.pyplot()
	println("Projecting on to the following coordinates: $c_1, $c_2 and $c_3.")
	points = Array{Float64}(undef, 0,3)
	for i in 1:size(partition,1)
		points = vcat(points, [partition[i,:point][c_1] partition[i,:point][c_2] partition[i,:point][c_3]])
	end
	Plots.scatter3d( points[:,1], points[:,2], points[:,3], group=partition[!,:index])
end

function Plot_Graph(vertex_locations, edges, sample) #plot a graph (either exact or modelled)
	Plots.pyplot()
	if length(sample) == 2
		scatter(vertex_locations[sample[1], :], vertex_locations[sample[2], :], s=1)
		for i in 1:(size(edges,1)-1)
			plot!([vertex_locations[sample[1], edges[i,1]], vertex_locations[sample[1],edges[i,2]]], [vertex_locations[sample[2], edges[i,1]], vertex_locations[sample[2], edges[i,2]]], linestyle = :solid)
		end
		i=size(edges,1)
		plot!([vertex_locations[sample[1], edges[i,1]], vertex_locations[sample[1],edges[i,2]]], [vertex_locations[sample[2], edges[i,1]], vertex_locations[sample[2], edges[i,2]]], linestyle = :solid)
	elseif length(sample) == 3
		scatter3d(vertex_locations[sample[1], :], vertex_locations[sample[2], :], vertex_locations[sample[3], :],s=1)
		for i in 1:(size(edges,1)-1)
			plot!([vertex_locations[sample[1], edges[i,1]], vertex_locations[sample[1], edges[i,2]]], [vertex_locations[sample[2], edges[i,1]], vertex_locations[sample[2], edges[i,2]]], [vertex_locations[sample[3], edges[i,1]], vertex_locations[sample[3], edges[i,2]]],  linestyle =:solid)
		end
		i=size(edges,1)
		plot!([vertex_locations[sample[1], edges[i,1]], vertex_locations[sample[1], edges[i,2]]], [vertex_locations[sample[2], edges[i,1]], vertex_locations[sample[2], edges[i,2]]], [vertex_locations[sample[3], edges[i,1]], vertex_locations[sample[3], edges[i,2]]],  linestyle =:solid)
	else
		println("Please specify 2 or 3 sample to project onto.")
		return []
	end
end

function Plot_Graph(vertex_locations, edges) #plot a graph (either exact or modelled)
	Plots.pyplot()
	if size(vertex_locations,1) == 2
		scatter(vertex_locations[sample[1], :], vertex_locations[sample[2], :], s=1)
		for i in 1:(size(edges,1)-1)
			plot!([vertex_locations[1, edges[i,1]], vertex_locations[1,edges[i,2]]], [vertex_locations[2, edges[i,1]], vertex_locations[2, edges[i,2]]], linestyle = :solid)
		end
		i=size(edges,1)
		plot!([vertex_locations[1, edges[i,1]], vertex_locations[1,edges[i,2]]], [vertex_locations[2, edges[i,1]], vertex_locations[2, edges[i,2]]], linestyle = :solid)
	elseif size(vertex_locations,1) == 3
		scatter3d(vertex_locations[1, :], vertex_locations[2, :], vertex_locations[3, :],s=1)
		for i in 1:(size(edges,1)-1)
			plot!([vertex_locations[1, edges[i,1]], vertex_locations[1, edges[i,2]]], [vertex_locations[2, edges[i,1]], vertex_locations[2, edges[i,2]]], [vertex_locations[3, edges[i,1]], vertex_locations[3, edges[i,2]]],  linestyle =:solid,  c=:red)
		end
		i=size(edges,1)
		plot!([vertex_locations[1, edges[i,1]], vertex_locations[1, edges[i,2]]], [vertex_locations[2, edges[i,1]], vertex_locations[2, edges[i,2]]], [vertex_locations[3, edges[i,1]], vertex_locations[3, edges[i,2]]],  linestyle =:solid,  c=:red)
	else
		println("Projecting onto the first 3 coordinates.")
		scatter3d(vertex_locations[1, :], vertex_locations[2, :], vertex_locations[3, :],s=1)
		for i in 1:(size(edges,1)-1)
			plot!([vertex_locations[1, edges[i,1]], vertex_locations[1, edges[i,2]]], [vertex_locations[2, edges[i,1]], vertex_locations[2, edges[i,2]]], [vertex_locations[3, edges[i,1]], vertex_locations[3, edges[i,2]]],  linestyle =:solid)
		end
		i=size(edges,1)
		plot!([vertex_locations[1, edges[i,1]], vertex_locations[1, edges[i,2]]], [vertex_locations[2, edges[i,1]], vertex_locations[2, edges[i,2]]], [vertex_locations[3, edges[i,1]], vertex_locations[3, edges[i,2]]],  linestyle =:solid)
	end
end

function model_fit(data, muT, sigma, E, S; EM_it = 3)

	#####################################################################################

	#### Functions for modelling the underlying structure ####

	#=
	To update:
	-keeps sigma fixed at the moment, this can be made optional
	=#

	#=
	MAKE SURE S AND MUT ARE IN THE SAME ORDER
	=#

	#=
	notns:
	n		---	Number of data points
	d		---	Dimension of data
	N		---	Number of strata pieces
	N_0		---	Number of 0-dim strata pieces
	N_1		---	Number of 1-dim strata pieces
	data	---	n dim array of d arrays		---	Sample from embedded graph
	mu		---	N_0 by d dim array	---	Iniitial vertex sample locations
	sigma	---	N dim array			---	Initial error on each strata piece
	E		---	n by N array		---	Initial assignment for each data point
	S		---	N dim str array		---	Abstract graph vertex/edge names in that order
	=#

	# data is a n dim array of d dimensional vectors

	#####################################################################################

	####pre-processing

	#get size of different dimensional strata

	#println("S starts as $S")

	#println("mu starts as $muT")


	N_0 	= size(muT,1)

	mu 		= [muT[i,:] for i in 1:N_0]

	N_1 	= size(S ,1) - N_0

	N 		= N_0 + N_1

	n  		= size(data,1)

	d 		= size(data[1],1)

	#load parameters
	
	E = [E[i,:] for i = 1:n]

	Pi      = 	sum(E)./n

	
	#println("sigma in EM is", sigma)
		
	MU		=	Dict{Any,Array}(S[i] => mu[i] for i = 1:N_0)

	SIGMA	=	Dict{Any,Float64}(S[i] => sigma[i] for i = 1:N)

	PI		=	Dict{Any,Float64}(S[i] => Pi[i] for i =1:N)
 
	#println("MU: ", MU)

	#println("SIGMA: ", SIGMA)

	#println("PI: ", PI)

	function reload_param(mu_t,sigma_t,Pi_t; S_t = S, N_0_t = N_0, N_t = N)

		#updates dictionary values for:
		#mu_t 		N_0 by d float array
		#SIGMA_t 	N float array
		#Pi 		N_float array

		MU_t	=  	Dict{Any,Array}(S_t[i] => mu_t[i] for i = 1:N_0_t)

		SIGMA_t	=  	Dict{Any,Float64}(S_t[i] => sigma_t[i] for i = 1:N_t)

		PI_t	=  	Dict{Any,Float64}(S_t[i] => Pi_t[i] for i = 1:N_t)

		return MU_t, SIGMA_t, PI_t

	end

	####density functions

	#0-dim strata

	function gaussian(x,v,s)

		linear = x - v

		return exp(-dot(linear, linear)/(2*s))

	end

	function p_e(x,v,s)

		d = size(x)[1]

		den = gaussian(x,v,s)/((2*pi*s)^(d/2))

		return den

	end

	function log_p_e(x,v,s)

		d = size(x)[1]
		
		linear = x - v

		den= -(d/2)*log(2*pi*s^2) - dot(linear,linear)/(2*s^2)

		return den

	end

	#1-dim strata

	#with clipping

	function p_l(xj, v1, v2, s)

		d = size(xj,1)

		xj = - xj

		den = (exp((4*dot(v1 - v2, (v1 + v2)/2 + xj)^2 - 4*norm(v1 - v2)^2*norm((v1 + v2)/2 + xj)^2)/(8*s^2*norm(v1 - v2)^2))*pi^(1/2 - d/2)*(-erf((2*dot(v1 - v2, (v1 + v2)/2 + xj) - norm(v1 - v2)^2)/(2*sqrt(2)*s*norm(v1 - v2))) + erf((2*dot(v1 - v2, (v1 + v2)/2 + xj) + norm(v1 - v2)^2)/(2*sqrt(2)*s*norm(v1 - v2)))))/norm(v1 - v2)

		if den <= -6

			den = -6

		end

  		return den

  	end

  	function log_p_l(Xj, V1, V2, s)

		d = size(Xj,1)

		Xj = - Xj

		den = log(2^((-1 - d)/2)*pi^(1/2 - d/2)*s^(1 - d)) +
				log(-erf((2*dot(V1 - V2, (V1 + V2)/2 + Xj) - norm(V1 - V2)^2)/(2*sqrt(2)*s*norm(V1 - V2))) +
			   	erf((2*dot(V1 - V2, (V1 + V2)/2 + Xj) + norm(V1 - V2)^2)/(2*sqrt(2)*s*norm(V1 - V2)))) -
			 	log(norm(V1 - V2)) + (4*dot(V1 - V2, (V1 + V2)/2 + Xj)^2 - 4*norm(V1 - V2)^2*norm((V1 + V2)/2 + Xj)^2)/
			  	(8*s^2*norm(V1 - V2)^2)

		if den <= -6

			den = -6

		end

  		return den

  	end

  	#density evaluation matrix

	function densities(MU_t,SIGMA_t; S_t = S, N_0 = N_0 , N_1 = N_1, N = N)

		zero_dim 	= [x -> p_e(x,MU_t[strat],SIGMA_t[strat]) for strat = S_t[1:N_0]]

		one_dim  	= [x -> p_l(x,MU_t[strat[2]],MU_t[strat[1]],SIGMA_t[strat]) for strat = S_t[N_0+1:N]]

		return one_dim # vcat(zero_dim,one_dim)

	end

	function log_densities(mu_t,SIGMA_t; S_t = S, N_0 = N_0 , N_1 = N_1, N = N)

		MU_t    	=   Dict{Any,Array}(S_t[i] => mu_t[i] for i = 1:N_0)

		zero_dim 	=   [x -> log_p_e(x,MU_t[strat],SIGMA_t[strat]) for strat = S_t[1:N_0]]

		one_dim  	=   [x -> log_p_l(x,MU_t[strat[1]],MU_t[strat[2]],SIGMA_t[strat]) for strat = S_t[N_0+1:N]]

		return  vcat(zero_dim,one_dim)

	end

	function density_matrix(data,density_functions)
		
		dm 			= 	[ (broadcast(rho, data))  for rho in density_functions ]

		dm 			=  	hcat(dm...)

		dm 			= 	[dm[i,:] for i = 1:n]
		 
		return dm

	end

	####cost construction

	#expectation matrix for current parameters

	function expectation(data,MU_t,SIGMA_t,Pi_t)

		DM = density_matrix(data,log_densities(MU_t,SIGMA_t))
		
		logpi = log.(Pi_t)

		E = broadcast(t -> logpi .+ t, DM)

		E = broadcast(t -> softmax(t), E)

		return E

	end

	#cost for current vertex values

	function Cost(data,mu0,SIGMA_t, Pi, E)

		mu_t = [mu0[i,:] for i in 1:N_0]

		log_DM = density_matrix(data,log_densities(mu_t,SIGMA_t))

		#not needed in optimisation but can be checked for EM
		#log_pi = log.(Pi)
		#log_sum = broadcast(t -> t + log_pi, log_DM)

		w_sum = broadcast((e,t) -> sum(e.*t), E, log_DM)

		return sum(w_sum)/n

	end

	####analytic gradient construction

	#gradient values for cost function

	function grad_log_p_l(k, Xj, V1, V2, s)

		xj = Xj[k]

		v1 = V1[k]

		v2 = V2[k]

		if dot(V2 - V1, Xj) > dot(V1, V2 - V1)

			xj = xj - ((dot(V2 - V1, Xj) - dot(V1, V2 - V1) )/(norm(V2 - V1)^2))*(v2 - v1)

			Xj = Xj - ((dot(V2 - V1, Xj) - dot(V1, V2 - V1) )/(norm(V2 - V1)^2)).*(V2 - V1)


		elseif dot(V2 - V1, Xj) < dot(V2, V1 - V2)

			xj = xj - ((dot(V2 - V1, Xj) -  dot(V2, V1 - V2) )/(norm(V2 - V1)^2))*(v2 - v1)

			Xj = Xj - ((dot(V2 - V1, Xj) -  dot(V2, V1 - V2) )/(norm(V2 - V1)^2)).*(V2 - V1)

		end


		g_v1 = -(v1 + v2 + 2*xj)/(4*s^2) + (2*(-1 + exp(dot(V1 - V2, (V1 + V2)/2 + Xj)/s^2))*(v1 - v2)*
			dot(V1 - V2, (V1 + V2)/2 + Xj) + (3*v1 - v2 + 2*xj - exp(dot(V1 - V2, (V1 + V2)/2 + Xj)/s^2)*
		  	(v1 + v2 + 2*xj))*norm(V1 - V2)^2)/
			(exp((2*dot(V1 - V2, (V1 + V2)/2 + Xj) + norm(V1 - V2)^2)^2/(8*s^2*norm(V1 - V2)^2))*sqrt(2*pi)*s*
			(erf((-2*dot(V1 - V2, (V1 + V2)/2 + Xj) + norm(V1 - V2)^2)/(2*sqrt(2)*s*norm(V1 - V2))) +
			erf((2*dot(V1 - V2, (V1 + V2)/2 + Xj) + norm(V1 - V2)^2)/(2*sqrt(2)*s*norm(V1 - V2))))*
			norm(V1 - V2)^3) + ((-v1 + v2)*dot(V1 - V2, (V1 + V2)/2 + Xj)^2 +
			(s^2*(-v1 + v2) + (v1 + xj)*dot(V1 - V2, (V1 + V2)/2 + Xj))*norm(V1 - V2)^2)/(s^2*norm(V1 - V2)^4)


		g_v2 = -(v1 + v2 + 2*xj)/(4*s^2) + (-2*(-1 + exp(dot(V1 - V2, (V1 + V2)/2 + Xj)/s^2))*(v1 - v2)*
			dot(V1 - V2, (V1 + V2)/2 + Xj) + (-v1 - v2 - 2*xj + exp(dot(V1 - V2, (V1 + V2)/2 + Xj)/s^2)*
	 		(-v1 + 3*v2 + 2*xj))*norm(V1 - V2)^2)/
			(exp((2*dot(V1 - V2, (V1 + V2)/2 + Xj) + norm(V1 - V2)^2)^2/(8*s^2*norm(V1 - V2)^2))*sqrt(2*pi)*s*
			(erf((-2*dot(V1 - V2, (V1 + V2)/2 + Xj) + norm(V1 - V2)^2)/(2*sqrt(2)*s*norm(V1 - V2))) +
			erf((2*dot(V1 - V2, (V1 + V2)/2 + Xj) + norm(V1 - V2)^2)/(2*sqrt(2)*s*norm(V1 - V2))))*
			norm(V1 - V2)^3) + ((v1 - v2)*dot(V1 - V2, (V1 + V2)/2 + Xj)^2 -
			(s^2*(-v1 + v2) + (v2 + xj)*dot(V1 - V2, (V1 + V2)/2 + Xj))*norm(V1 - V2)^2)/(s^2*norm(V1 - V2)^4)

	   	return g_v1, g_v2

	end

 	function grad_log_p_l_n(k, Xj, V1, V2, s)

		xj = Xj[k]

		v1 = V1[k]

		v2 = V2[k]

                if dot(V2 - V1, Xj) > dot(V1, V2 - V1)

                        xj = xj - ((dot(V2 - V1, Xj) - dot(V1, V2 - V1) )/(norm(V2 - V1)^2))*(v2 - v1)

                        Xj = Xj - ((dot(V2 - V1, Xj) - dot(V1, V2 - V1) )/(norm(V2 - V1)^2)).*(V2 - V1)


                elseif dot(V2 - V1, Xj) < dot(V2, V1 - V2)

                        xj = xj - ((dot(V2 - V1, Xj) -  dot(V2, V1 - V2) )/(norm(V2 - V1)^2))*(v2 - v1)

                        Xj = Xj - ((dot(V2 - V1, Xj) -  dot(V2, V1 - V2) )/(norm(V2 - V1)^2)).*(V2 - V1)

                end

        function V_r(v_r, V)

        	V[k] = v_r

        	return V

        end

		log_p_l_proj = v_r -> log_p_l(Xj, V_r(v_r[1], V1), V_r(v_r[2], V2), s)

		g = Calculus.gradient(log_p_l_proj)

		return g([v1, v2])

	end

	function grad_log_cost!(Grad, data, mu_t, S_t, E_list, SIGMA_t)

		n  = size(data,1)

		d  = size(data[1],1)

		grad = Dict{Any,Array}(S_t[i] => zeros(Float64, 1 , d) for i = 1:N_0)

		Mu_t = Dict{Any,Array}(S_t[i] => mu_t[i,:] for i = 1:N_0)

		E_array = vcat(transpose(E_list)...)

		E_t = Dict{Any,Array}(S_t[i] => E_array[:,i] for i = 1:N)

		for s in S_t

			if s isa Integer

				grad[s] += transpose(sum(E_t[s].*broadcast( x -> ((MU[s] - x)./SIGMA[s]), data))/n)

			elseif s isa Tuple

				for k in 1:d

						for j in 1:n

							g_v1, g_v2 = grad_log_p_l_n(k, data[j], Mu_t[s[1]], Mu_t[s[2]], SIGMA_t[s])
							
							grad[s[1]][k] += E_t[s[1]][j]*g_v1/n

							grad[s[2]][k] += E_t[s[2]][j]*g_v2/n
							
							#println(grad[s[1]][k],grad[s[2]][k])
						end

				end

			end
			
		end


		for i in 1:N_0

			Grad[i,:] = - grad[S_t[i]]

		end

	end
	
	####EM algorithm

	mu00 = vcat([m for m in mu]'...)

	for it = 1:EM_it
		#E step
		
		E = expectation(data,MU,SIGMA,Pi)

		Pi      = 	sum(E)./n

		#M step

		mu0 = vcat([m for m in mu]'...)
		
		C = t -> -Cost(data,t,SIGMA, Pi ,E)

		#swap if don't want gradient clipping or set tol = Inf
		#G! = (g,t) -> FiniteDiff.finite_difference_gradient!(g,C,t)

		function G!(storage,t; tol = 0.1)

			FiniteDiff.finite_difference_gradient!(storage,C,t)

			storage[abs.(storage) .> tol] = sign.(storage[abs.(storage) .> tol])*tol

			return

		end

		results = Optim.optimize(C,G!,mu0,GradientDescent(),Optim.Options(g_tol = 1e-2))

		mu1 = Optim.minimizer(results)

		mu = [mu1[i,:] for i in 1:N_0]

		#println("Current Pi:")
		#println(Pi)

		#println("Current Sigma:")
		#println(sigma)

		#println("Current vertices:")
		#println(mu)

		print("For iteration: ", it," cost is: ", -C(mu1) + sum(log.(Pi)) )

		#Update parameters
		
		MU, SIGMA, PI = reload_param(mu, sigma, Pi ; S_t = S, N_0_t = N_0, N_t = N)

		#println(MU)
	

	end

	return MU, SIGMA, PI, E

end


function Save_PPC(ppc::PartitionedPointCloud, path::String)
	jldopen(path, "w") do file
		addrequire(file, Skyler)
		addrequire(file, LinearAlgebra)
		write(file, "PartitionedPointCloud", ppc)
	end
end # Save_PPC

#### Main functions and wrappers ####
function skyler(points::A, sample_epsilon::P, radius::R; out = "model", EM_it = 3, sig = sample_epsilon/2) where {P<:Number, Q<:Number, R<:Number, A<:Union{LinearAlgebra.Adjoint,AbstractArray}} #wrapper for obtaining the partition

	dimension = size(points,1)
	ppc = Partition_Point_Cloud(points, sample_epsilon, radius, (3*radius)/2 + 2*sample_epsilon, 3*sample_epsilon)
	if out == "PPC"
		return ppc
	elseif out == "DF"
			return Create_Partition_DataFrame(ppc)
	elseif out == "model"
		vertices = size(collect(ppc.Dimensions[0]),1)
		vertex_guess = Array{Float64}(undef, dimension,0)
		for i in 1:vertices
			v_i = zeros(dimension, 1)
			for j in collect(ppc.Strata[i])
				v_i .+= points[:, j]
			end
			v_i = v_i/(length(collect(ppc.Strata[i])))
			vertex_guess = hcat(vertex_guess, v_i)
		end
		#println("Vertex guess is: ",vertex_guess)
		E = zeros(size(points,2), size(ppc.Boundaries,1))
		
		for i in 1:vertices
			for j in collect(ppc.Strata[i])
				E[j,i] =1
			end
		end
		for i in collect(ppc.Dimensions[1])
			for j in collect(ppc.Strata[i])
				E[j,i] =1
			end
		end
		
		sigma = [sig for i in 1:size(ppc.Boundaries,1)]
	
		S = []
		for i in collect(ppc.Dimensions[0])
			push!(S, i)
		end
		for k in 1:maximum(keys(ppc.Dimensions))
			for i in collect(ppc.Dimensions[k])
				t =tuple(findall(x->x==1, ppc.Boundaries[:,i])...,)
				push!(S, t)
			end
		end

		data = vcat([points[:,i] for i in 1:size(points,2)])
		
		v_guess = Array{Float64}(undef, size(vertex_guess,2), size(vertex_guess,1))
		
		for i in 1:size(vertex_guess,2)
			for j in 1:size(vertex_guess,1)
				v_guess[i,j] = vertex_guess[j,i]
			end
		end
		#println("S is $S")
		MU, SIGMA, PI, E = model_fit(data, v_guess, sigma, E, S, EM_it= EM_it)
		mu_ans = collect(values(MU))
		mu_ans = vcat([datum' for datum in mu_ans]...)
		return ppc, mu_ans', vertex_guess
	else
		println("Error: output option 'out' selected is not valid, please choose from 'model' to obtain a model of the underlying structure, or 'struct' to receive a DataFrame indicating which strata a point as been associated with and the dimension of this strata, and the boundary relations.")
		return []
	end
end

function skyler(ppc::PartitionedPointCloud, sample_epsilon::P; EM_it=3, sig=sample_epsilon/2) where {P<:Number} #wrapper for modelling a point cloud
	points=ppc.Points
	dimension = size(points,1)
	vertices = size(collect(ppc.Dimensions[0]),1)
	vertex_guess = Array{Float64}(undef, dimension,0)
	for i in 1:vertices
		v_i = zeros(dimension, 1)
		for j in collect(ppc.Strata[i])
			v_i .+= points[:, j]
		end
		v_i = v_i/(length(collect(ppc.Strata[i])))
		vertex_guess = hcat(vertex_guess, v_i)
	end
	println("Vertex guess is: ",vertex_guess)
	E = zeros(size(points,2), size(ppc.Boundaries,1))
	
	for i in 1:vertices
		for j in collect(ppc.Strata[i])
			E[j,i] =1
		end
	end
	for i in collect(ppc.Dimensions[1])
		for j in collect(ppc.Strata[i])
			E[j,i] =1
		end
	end
	
	sigma = [sig for i in 1:size(ppc.Boundaries,1)]

	S = []
	for i in collect(ppc.Dimensions[0])
		push!(S, i)
	end
	for k in 1:maximum(keys(ppc.Dimensions))
		for i in collect(ppc.Dimensions[k])
			t =tuple(findall(x->x==1, ppc.Boundaries[:,i])...,)
			push!(S, t)
		end
	end
	

	data = vcat([points[:,i] for i in 1:size(points,2)])
		
	v_guess = Array{Float64}(undef, size(vertex_guess,2),size(vertex_guess,1))
		
	for i in 1:size(vertex_guess,2)
		for j in 1:size(vertex_guess,1)
			v_guess[i,j] = vertex_guess[j,i]
		end
	end
	MU, SIGMA, PI, E = model_fit(data, v_guess, sigma, E, S, EM_it=EM_it)
	mu_ans = collect(values(MU))
	mu_ans = vcat([datum' for datum in mu_ans]...)
	return mu_ans', vertex_guess
end #skyler for modelling a partioned point cloud
	
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
	
	points= convert(Array,CSV.read("Sample-1.csv", header=false))
	ppc = skyler(points, 0.01, 0.12, out="PPC")
	par = Create_Partition_DataFrame(ppc)
	saved_partition = CSV.read("Partition-1.csv")
	
	
	if par == saved_partition
		return []
	else
		println("Error: test_1 not passed.")
		return par
	end
end


function test_2()
	dir = pwd()
	cd(@__DIR__)
	cd("..")
	cd("Examples")
	
	points= convert(Array,CSV.read("Sample-2.csv",header=false))'
	ppc = skyler(points, 0.1, 1.2, out="PPC")
	par = Create_Partition_DataFrame(ppc)
	saved_partition = CSV.read("Partition-2.csv")
	
	if par == saved_partition
		return []
	else
		println("Error: test_2 not passed.")
		return par
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
