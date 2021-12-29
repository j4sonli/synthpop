ENV["GUROBI_HOME"] = ARGS[9]
ENV["GRB_LICENSE_FILE"] = ARGS[10]

using Gurobi, JuMP, CSV, DataFrames, Random, RLEVectors, DelimitedFiles, StatsBase

Random.seed!(1033)

@time begin

const CURR_PUMA = ARGS[1]
const HH_POOL_FILE = ARGS[2]
const INDV_POOL_FILE = ARGS[3]
const SYN_HHS_FILE = ARGS[4]
const SYN_INDVS_FILE = ARGS[5]
const PUMA_TRACT_EQUIV_FILE = ARGS[6]
const MARG_DIR = ARGS[7]
const GUROBI_LOGFILE = ARGS[8]
const PUMA_MARG_FILES = ["puma_hhtype.csv"]
const TRACT_MARG_FILES = ["tract_tenur_hhinc.csv", "tract_nwork.csv"]
const BLKGP_MARG_FILES = ["blkgp_tenur_hhsiz.csv"]
const MARG_COLS = [11, 21, 4, 22]
const INDV_TRACT_MARG_FILES = ["tract_i_sex_i_age.csv"]
const INDV_BLKGP_MARG_FILES = ["blkgp_emply.csv"]
const INDV_MARG_COLS = [19, 20]

## Get tracts
puma_tract_equiv = CSV.read(PUMA_TRACT_EQUIV_FILE, DataFrame)
curr_geo = filter(row -> row.STATEFP == parse(Float64, CURR_PUMA[1:2]) && row.PUMA5CE == parse(Float64, CURR_PUMA[3:7]), puma_tract_equiv)
const TRACT_GEOIDS = [string(lpad(string(sct[1]), 2, "0"), lpad(string(sct[2]), 3, "0"), lpad(string(sct[3]), 6, "0"))
                        for sct in zip(curr_geo.STATEFP, curr_geo.COUNTYFP, curr_geo.TRACTCE)]

## Process marginals
indv_index = length(vcat(PUMA_MARG_FILES, TRACT_MARG_FILES, BLKGP_MARG_FILES))+1
marginals = []
levelss = []
for marg_file in vcat(PUMA_MARG_FILES, TRACT_MARG_FILES, BLKGP_MARG_FILES, INDV_TRACT_MARG_FILES, INDV_BLKGP_MARG_FILES)
        marg_df = CSV.read("$MARG_DIR$marg_file", DataFrame)
        if marg_file in PUMA_MARG_FILES
                marg_df = filter(row -> lpad(string(row.geoid), 7, "0") == CURR_PUMA, marg_df)
        elseif marg_file in TRACT_MARG_FILES || marg_file in INDV_TRACT_MARG_FILES
                marg_df = filter(row -> lpad(string(row.geoid), 11, "0") in TRACT_GEOIDS, marg_df)
        elseif marg_file in BLKGP_MARG_FILES || marg_file in INDV_BLKGP_MARG_FILES
                marg_df = filter(row -> lpad(string(row.geoid), 12, "0")[1:11] in TRACT_GEOIDS, marg_df)
        end
        select!(marg_df, Not(:geoid))
        push!(marginals, sum.(eachcol(marg_df)))
        push!(levelss, names(marg_df))
end

const puma_hh_pop = sum(marginals[1])

## Process household  and individual pools
hh_pool_df = CSV.read(HH_POOL_FILE, DataFrame)
hh_pool = [row[:] for row in eachrow(Matrix(hh_pool_df))]

const n = length(hh_pool)

indv_pool_df = CSV.read(INDV_POOL_FILE, DataFrame)
hh_pool_hhsizs = rle(indv_pool_df.HHID)[2]
# precompute individual variable level counts for each household; factors for
#     household-level variables can all be 1 or whatever
indv_factors = Any[]
for lvls in levelss[1:(indv_index-1)]
        push!(indv_factors, [[1 for _ in 1:n] for _ in 1:length(lvls)])
end
for (indv_marg_col, lvls) in zip(INDV_MARG_COLS, levelss[indv_index:length(levelss)])
        with_counts = combine(groupby(indv_pool_df, :HHID), indv_marg_col => hh_values -> Tuple(count(val->val==l, hh_values) for l in lvls))
        counts_mat = [[v for v in vals] for vals in with_counts[!,2]]
        counts_mat = [[row[col] for row in counts_mat] for col in 1:size(counts_mat[1])[1]]
        push!(indv_factors, counts_mat)
end

println("Household Pool Size: ", n)
println("Marginals: ", vcat(PUMA_MARG_FILES, TRACT_MARG_FILES, BLKGP_MARG_FILES, INDV_TRACT_MARG_FILES, INDV_BLKGP_MARG_FILES))

#######################################################
## Define optimization model
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "LogFile", GUROBI_LOGFILE)

## Variables
@variable(model, x[1:n*2]) # define variable x
@variable(model, y[1:sum(length(lvls) for lvls in levelss)]) # define absolute value helper variable y

## Objective
@objective(model, Min, sum(y))

## Constraints
@constraint(model, con_lb[i=1:n*2], 0 <= x[i])
@constraint(model, con_ub[i=1:n*2], x[i] <= 1)
@constraint(model, con_rowsum[i=1:n], x[2*(i-1)+1] + x[2*(i-1)+2] == 1)
@constraint(model, con_col1sum, sum(x[2*(i-1)+1] for i in 1:n) == puma_hh_pop)
@constraint(model, con_col2sum, sum(x[2*(i-1)+2] for i in 1:n) == n - puma_hh_pop)
@constraint(model, con_col1indvsum, sum(hh_pool_hhsizs[i] * x[2*(i-1)+1] for i in 1:n) == sum(marginals[indv_index]))

# add marginal cell constraints to enforce absolute value
marg_cell_i = 1
for (marg_i, (marg_col, marg, lvls, indv_facs)) in enumerate(zip(vcat(MARG_COLS, INDV_MARG_COLS), marginals, levelss, indv_factors))
        cell_weight = length(lvls) * 2
        for (lvl_i, (lvl, ifs)) in enumerate(zip(lvls, indv_facs))
                @constraint(model, y[marg_cell_i] >= ((sum(ifs[i] * x[2*(i-1)+1] for (i,hh) in enumerate(hh_pool) if marg_i>=indv_index || hh[marg_col]==lvl)) - marg[lvl_i])/(1+marg[lvl_i])/cell_weight)
                @constraint(model, y[marg_cell_i] >= -((sum(ifs[i] * x[2*(i-1)+1] for (i,hh) in enumerate(hh_pool) if marg_i>=indv_index || hh[marg_col]==lvl)) - marg[lvl_i])/(1+marg[lvl_i])/cell_weight)
                global marg_cell_i += 1
        end
end

## Optimize
status=optimize!(model)

println("******************************************************")
println("optimal objective value is = ", objective_value(model))

## extract selected households, and reset IDs
sol_matrix = transpose(reshape(value.(x), (2,n)))
sorted_hhids = sort([[abs(irow[2][1] - 1), hh_pool[irow[1]][1]] for irow in enumerate(eachrow(sol_matrix))])
synpop_hhids = [diff_hhid[2] for diff_hhid in sorted_hhids[1:puma_hh_pop]]

syn_hhs = filter(row -> (row.HHID in synpop_hhids), hh_pool_df)
syn_hhs.HHID = 1:DataFrames.nrow(syn_hhs)
CSV.write(SYN_HHS_FILE, syn_hhs)

syn_indvs = filter(row -> (row.HHID in synpop_hhids), indv_pool_df)
true_hhsizs = rle(syn_indvs.HHID)[2]
syn_indvs.HHID = collect(rep(collect(1:DataFrames.nrow(syn_hhs)), each = true_hhsizs))
syn_indvs.INDVID = 1:DataFrames.nrow(syn_indvs)
CSV.write(SYN_INDVS_FILE, syn_indvs)

end
