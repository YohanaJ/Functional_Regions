using Distributed
using CSV
using DataFrames
using Arrow
using DelimitedFiles
using Combinatorics
using SharedArrays


# this function calculates beta dispersion for all pairs of plots
function get_dist(plot_ind::Int64, num_plots::Int64, abund_col::Int64, comm_coor::SharedArray, dist_mat::SharedArray) 

    # initialize the community size integers
    n1 = zeros(Int64, 1)
    n2 = zeros(Int64, 1)

    myres = zeros(Float64, num_plots+1)

    # myres[x, 1] = plot_ind
    myres[1] = plot_ind  

    # cycle through the chunks
    for x = (plot_ind+1):num_plots
        
        
        print(x)

        # get the current species
        @views s1 = comm_coor[comm_coor[:, 1] .== plot_ind, :]
        @views s2 = comm_coor[comm_coor[:, 1] .== x, :]
    
        n1[1] = size(s1, 1)
        n2[1] = size(s2, 1)

        # calculate the distancs
        @views mydst = dist_mat[s1[:,2], s2[:,2]]

        # get the nearest neighbor for each taxon in each community
        m1 = minimum(mydst, dims = 2)
        m2 = minimum(mydst, dims = 1)

        # get the average (i.e., beta dispersion)
        myres[x+1] = (sum(m1) + sum(m2)) / (n1[1] + n2[1]) 

    end

    return myres

end

