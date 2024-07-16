using Pkg
#Pkg.activate("IE")
using DataFrames
using Plots
using XLSX
using ColorSchemes
using Statistics
using LinearAlgebra
using Dates
using CSV

include("dataprocessing.jl")
include("output_for_TASEP.jl")
include("ratecalculation.jl")
include("IE_plotting.jl")


