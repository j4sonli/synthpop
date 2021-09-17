using Gurobi, JuMP, CSV, DataFrames, Random, RLEVectors, DelimitedFiles

Random.seed!(1033)

@time begin

const CURR_PUMA = "2500506"
const HH_POOL_FILE = "Optim_Data/hh_pool-20210915-gq.csv"
const INDV_POOL_FILE = "Optim_Data/indv_pool-20210915-gq.csv"
const SYN_HHS_NOTRACT_FILE = "Optim_Data/Output/syn_hhs_notract-20210915-gq.csv"
const SYN_INDVS_NOTRACT_FILE = "Optim_Data/Output/syn_indvs_notract-20210915-gq.csv"
const PUMA_MARG_FILES = ["puma_hhtype-gq.csv"]
const TRACT_MARG_FILES = ["tract_tenur_hhinc-gq.csv", "tract_nwork-gq.csv"]
const BLKGP_MARG_FILES = ["blkgp_tenur_hhsiz-gq.csv"]
const MARG_COLS = [8, 10, 4, 11]

## Get tracts
puma_tract_equiv = CSV.read("Optim_Data/2010_Census_Tract_to_2010_PUMA.csv", DataFrame)
curr_geo = filter(row -> row.STATEFP == parse(Float64, CURR_PUMA[1:2]) && row.PUMA5CE == parse(Float64, CURR_PUMA[3:7]), puma_tract_equiv)
const TRACT_GEOIDS = [string(lpad(string(sct[1]), 2, "0"), lpad(string(sct[2]), 3, "0"), lpad(string(sct[3]), 6, "0"))
                        for sct in zip(curr_geo.STATEFP, curr_geo.COUNTYFP, curr_geo.TRACTCE)]

## Process marginals
marginals = []
levelss = []
for marg_file in vcat(PUMA_MARG_FILES, TRACT_MARG_FILES, BLKGP_MARG_FILES)
        marg_df = CSV.read("Optim_Data/ACS_Marg/$marg_file", DataFrame)
        if marg_file in PUMA_MARG_FILES
                marg_df = filter(row -> string(row.puma_id) == CURR_PUMA, marg_df)
                select!(marg_df, Not(:puma_id))
        elseif marg_file in TRACT_MARG_FILES
                marg_df = filter(row -> string(row.geoid) in TRACT_GEOIDS, marg_df)
                select!(marg_df, Not(:geoid))
        elseif marg_file in BLKGP_MARG_FILES
                marg_df = filter(row -> string(row.geoid)[1:11] in TRACT_GEOIDS, marg_df)
                select!(marg_df, Not(:geoid))
        end
        push!(marginals, sum.(eachcol(marg_df)))
        push!(levelss, names(marg_df))
end

const puma_pop = sum(marginals[1])

## Process household pool
hh_pool_df = CSV.read(HH_POOL_FILE, DataFrame)
puma_df = filter(row -> string(row.puma) == CURR_PUMA, hh_pool_df)
pool = [row[:] for row in eachrow(Matrix(puma_df))]

const n = length(pool)

println("Household Pool Size: ", n)
println("Marginals: ", vcat(PUMA_MARG_FILES, TRACT_MARG_FILES, BLKGP_MARG_FILES))

#######################################################
## Define optimization model
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "LogFile", "gurobi_log.txt")

## Variables
@variable(model, x[1:n*2]) # define variable x
@variable(model, y[1:sum(length(lvls) for lvls in levelss)]) # define absolute value helper variable y

## Objective
@objective(model, Min, sum(y))

## Constraints
@constraint(model, con_lb[i=1:n*2], 0 <= x[i])
@constraint(model, con_ub[i=1:n*2], x[i] <= 1)
@constraint(model, con_rowsum[i=1:n], x[2*(i-1)+1] + x[2*(i-1)+2] == 1)
@constraint(model, con_col1sum, sum(x[2*(i-1)+1] for i in 1:n) == puma_pop)
@constraint(model, con_col2sum, sum(x[2*(i-1)+2] for i in 1:n) == n - puma_pop)

# add marginal cell constraints to enforce absolute value
marg_cell_i = 1
for mml in zip(MARG_COLS, marginals, levelss)
        cell_weight = length(mml[3]) * 2
        for llvl in enumerate(mml[3])
                @constraint(model, y[marg_cell_i] >= (sum(x[2*(ihh[1]-1)+1] for ihh in enumerate(pool) if ihh[2][mml[1]] == llvl[2]) - mml[2][llvl[1]])/(1+mml[2][llvl[1]])/cell_weight)
                @constraint(model, y[marg_cell_i] >= -(sum(x[2*(ihh[1]-1)+1] for ihh in enumerate(pool) if ihh[2][mml[1]] == llvl[2]) - mml[2][llvl[1]])/(1+mml[2][llvl[1]])/cell_weight)
                global marg_cell_i += 1
        end
end

## Optimize
status=optimize!(model)

println("******************************************************")
println("optimal objective value is = ", objective_value(model))

## extract selected households, and reset IDs
sol_matrix = transpose(reshape(value.(x), (2,n)))
synpop_hhids = [pool[irow[1]][1] for irow in enumerate(eachrow(sol_matrix)) if abs(irow[2][1] - 1) < 1e-3]

syn_hhs_notract = filter(row -> (string(row.puma) == CURR_PUMA) && (row.HHID in synpop_hhids), puma_df)
syn_hhs_notract.HHID = 1:DataFrames.nrow(syn_hhs_notract)
CSV.write(SYN_HHS_NOTRACT_FILE, syn_hhs_notract)

indv_pool_df = CSV.read(INDV_POOL_FILE, DataFrame)
syn_indvs_notract = filter(row -> (string(row.puma) == CURR_PUMA) && (row.HHID in synpop_hhids), indv_pool_df)
hhsizs = map(hhsiz->parse(Int, hhsiz[1]), syn_hhs_notract.hhsiz)
syn_indvs_notract.HHID = collect(rep(collect(1:DataFrames.nrow(syn_hhs_notract)), each = hhsizs))
syn_indvs_notract.INDVID = 1:DataFrames.nrow(syn_indvs_notract)
CSV.write(SYN_INDVS_NOTRACT_FILE, syn_indvs_notract)

end

valuex = value.(x)
writedlm("Optime_Results/selection_weights-20210915-gq.csv", transpose(valuex), "\t")
result = DataFrame(group = [lvl for marg in levelss for lvl in marg],
                  synpop = [round(sum(valuex[2*(ihh[1]-1)+1] for ihh in enumerate(pool) if ihh[2][mml[1]] == llvl[2]), digits=2) for mml in zip(MARG_COLS, marginals, levelss) for llvl in enumerate(mml[3])],
                  marginal = [mml[2][llvl[1]] for mml in zip(MARG_COLS, marginals, levelss) for llvl in enumerate(mml[3])])
result.error = result.synpop - result.marginal
result.perror = result.error ./ result.marginal * 100

CSV.write("Optim_Results/selection-20210915-gq.csv", result)
