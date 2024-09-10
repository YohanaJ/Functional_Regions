using Distributed
using CSV
using DataFrames
using Distances
using DelimitedFiles
using Combinatorics
using SharedArrays


# number of cores
addprocs(64)

# set the directory
@everywhere cd("data/julia")
@everywhere include("code/functions_beta_dispersion.jl")


# read in the data
comm_coor = DataFrame(Arrow.Table(string("community_mat_7.arrow")))
cnames = names(comm_coor)
comm_coor = SharedArray(Matrix{Int64}(comm_coor))
grid_look = DataFrame(Arrow.Table(string("grid_lookup_7.arrow")))

# # sort the grid id
grid_look = grid_look[sortperm(grid_look.grid_num), :]

# get the plots and number
plots = sort(unique(comm_coor[:, 1]))
n = size(plots, 1)

# read in the files, both individual traits and overall distance
trfiles = string.("traits/", readdir("traits"))

for ii=1:size(trfiles, 1)
    print(ii,"\n")

    # create the distance matrix
    dist_mat = SharedArray(Matrix{Float32}(DataFrame(Arrow.Table(trfiles[ii]))))

    # run in parallel. 3 is just a placeholder since we're not doing weighted
    res_flat = pmap(x -> get_dist(x, n, 3, comm_coor, dist_mat), 1:n)

    # flatten the output
    res_flat = transpose(hcat(res_flat...))

    # sort by the first column to make sure plot ids align
    res_flat = res_flat[sortperm(res_flat[:, 1]), :]

    # make summetric
    res_flat[:, 2:size(res_flat, 2)] = res_flat[:, 2:size(res_flat, 2)] + transpose(res_flat[:, 2:size(res_flat, 2)])

    #convert to a data frame
    res_flat = DataFrame(res_flat, vcat(["new_grid_id"], [ string(round(Int32, x)) for x in grid_look.new_grid_id])) 
    res_flat[:, 1] = round.(Int64, grid_look.new_grid_id)

    # write the output
    Arrow.write(replace(trfiles[ii], "_mat_" => "_matrix_"), res_flat)

    # empty out the matrix
    res_flat = "hello"
end


### calculate distance between all lat/longs
ll = CSV.read(string("lat_long_7.csv"), DataFrame)
lld = pairwise(Haversine(), transpose(hcat(ll.new_lon, ll.new_lat)), transpose(hcat(ll.new_lon, ll.new_lat)))
lld2 = DataFrame(hcat(ll.new_grid_id, lld), vcat(["new_grid_id"], [string(round(Int32, x)) for x in ll.new_grid_id])) 
Arrow.write(string("/home/dsenn/Git/cluster_forests/data/julia/distance_matrix_Haversine_non_natives7.arrow"), lld2)
