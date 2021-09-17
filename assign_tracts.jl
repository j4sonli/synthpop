using Gurobi, JuMP, CSV, DataFrames, Random

Random.seed!(1033)

@time begin

const CURR_PUMA = "2500506"
const SYN_HHS_NOTRACT_FILE = "Optim_Data/Output/syn_hhs_notract-20210915-gq.csv"
const SYN_INDVS_NOTRACT_FILE = "Optim_Data/Output/syn_indvs_notract-20210915-gq.csv"
const SYN_HHS_FILE = "Optim_Data/Output/syn_hhs-20210915-gq.csv"
const SYN_INDVS_FILE = "Optim_Data/Output/syn_indvs-20210915-gq.csv"
const TRACT_MARG_FILES = ["tract_tenur_hhinc.csv", "tract_nwork.csv", "tract_hhtype.csv"]
const BLKGP_MARG_FILES = ["blkgp_tenur_hhsiz.csv"]
const MARG_COLS = [10, 4, 8, 11]

## Get tracts and block groups
puma_tract_equiv = CSV.read("Optim_Data/2010_Census_Tract_to_2010_PUMA.csv", DataFrame)
curr_geo = filter(row -> row.STATEFP == parse(Float64, CURR_PUMA[1:2]) && row.PUMA5CE == parse(Float64, CURR_PUMA[3:7]), puma_tract_equiv)
const TRACT_GEOIDS = [string(lpad(string(sct[1]), 2, "0"), lpad(string(sct[2]), 3, "0"), lpad(string(sct[3]), 6, "0"))
                        for sct in zip(curr_geo.STATEFP, curr_geo.COUNTYFP, curr_geo.TRACTCE)]

blkgp_marg_df = CSV.read("Optim_Data/ACS_Marg/$(BLKGP_MARG_FILES[1])", DataFrame)
blkgp_marg_df = filter(row -> string(row.geoid)[1:11] in TRACT_GEOIDS, blkgp_marg_df)
const BLKGP_GEOIDS = string.(blkgp_marg_df.geoid)

# create tract-block group mapping
const tract_blkgp = [findall(tract -> tract == tract_geoid, SubString.(BLKGP_GEOIDS, 1, 11)) for tract_geoid in TRACT_GEOIDS]

## Process marginals
const m = length(BLKGP_GEOIDS)
marginals = []
levelss = []
const MARG_FILES = vcat(TRACT_MARG_FILES, BLKGP_MARG_FILES)
for marg_file in MARG_FILES
        marg_df = CSV.read("Optim_Data/ACS_Marg/$marg_file", DataFrame)
        marg_df = filter(row -> string(row.geoid) in (marg_file in TRACT_MARG_FILES ? TRACT_GEOIDS : BLKGP_GEOIDS), marg_df)
        select!(marg_df, Not(:geoid))
        # convert df into array of arrays
        marg = [col[:] for col in eachcol(Matrix(marg_df))]

        push!(marginals, marg)
        push!(levelss, names(marg_df))
end

const tract_pops = sum(marginals[1], dims=1)[1]
const blkgp_pops = sum(marginals[length(MARG_FILES)], dims=1)[1]

## Process population
syn_hhs_notract = CSV.read(SYN_HHS_NOTRACT_FILE, DataFrame)
puma_df = filter(row -> string(row.puma) == CURR_PUMA, syn_hhs_notract)
pop = [row[:] for row in eachrow(Matrix(puma_df))]

const n = length(pop)

println("Block groups: ", m)
println("Tracts: ", length(tract_pops))
println("Households: ", n)
println("Marginals: ", MARG_FILES)

#######################################################
## Define optimization model
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "LogFile", "gurobi_log.txt")

## Variables
@variable(model, x[1:n*m]) # define variable x
n_cells = sum(sum(length(lvl) for lvl in marg) for marg in marginals)
@variable(model, y[1:n_cells]) # define absolute value helper variable y

## Objective
@objective(model, Min, sum(y))

## Constraints
@constraint(model, con_lb[i=1:n*m], 0 <= x[i]) # lower bound constraint
@constraint(model, con_ub[i=1:n*m], x[i] <= 1) # upper bound constraint
@constraint(model, con_rowsum[i=1:n], sum(x[m*(i-1)+j] for j in 1:m) == 1)
@constraint(model, con_colsum[j=1:m], sum(x[m*(i-1)+j] for i in 1:n) == blkgp_pops[j])

# add marginal cell constraints to enforce absolute value
marg_cell_i = 1
for mml in zip(MARG_COLS, marginals, levelss)
        cell_weight = length(mml[3]) * length(mml[2][1])
        for llvl in enumerate(mml[3])
                if length(mml[2][1]) == m # block group marginal
                        for j in 1:m
                                @constraint(model, y[marg_cell_i] >= (sum(x[m*(ihh[1]-1)+j] for ihh in enumerate(pop) if ihh[2][mml[1]] == llvl[2]) - mml[2][llvl[1]][j])/(1+mml[2][llvl[1]][j])/cell_weight*blkgp_pops[j])
                                @constraint(model, y[marg_cell_i] >= -(sum(x[m*(ihh[1]-1)+j] for ihh in enumerate(pop) if ihh[2][mml[1]] == llvl[2]) - mml[2][llvl[1]][j])/(1+mml[2][llvl[1]][j])/cell_weight*blkgp_pops[j])
                                global marg_cell_i += 1
                        end
                else # tract marginal
                        for j in 1:length(tract_pops)
                                @constraint(model, y[marg_cell_i] >= (sum(sum(x[m*(ihh[1]-1)+k] for k in tract_blkgp[j]) for ihh in enumerate(pop) if ihh[2][mml[1]] == llvl[2]) - mml[2][llvl[1]][j])/(1+mml[2][llvl[1]][j])/cell_weight*tract_pops[j])
                                @constraint(model, y[marg_cell_i] >= -(sum(sum(x[m*(ihh[1]-1)+k] for k in tract_blkgp[j]) for ihh in enumerate(pop) if ihh[2][mml[1]] == llvl[2]) - mml[2][llvl[1]][j])/(1+mml[2][llvl[1]][j])/cell_weight*tract_pops[j])
                                global marg_cell_i += 1
                        end
                end
        end
end

## Optimize
status=optimize!(model)

println("******************************************************")
println("optimal objective value is = ", objective_value(model))

## assign households to optimized tracts
sol_matrix = transpose(reshape(value.(x), (m,n)))
blkgp_assignments = [BLKGP_GEOIDS[argmax(row)] for row in eachrow(sol_matrix)]
blkgp_by_hhid = DataFrame(HHID = [row[1] for row in pop],
                          geoid = blkgp_assignments)

syn_hhs = leftjoin(syn_hhs_notract, blkgp_by_hhid, on=:HHID)
select!(syn_hhs, Not(:tenur_hhinc_prox))
select!(syn_hhs, Not(:tenur_hhsiz_prox))
CSV.write(SYN_HHS_FILE, syn_hhs)

syn_indvs_notract = CSV.read(SYN_INDVS_NOTRACT_FILE, DataFrame)
syn_indvs = leftjoin(syn_indvs_notract, blkgp_by_hhid, on=:HHID)
CSV.write(SYN_INDVS_FILE, syn_indvs)

end


valuex = value.(x)
result = DataFrame(group = [string(lvl,"_",TRACT_GEOIDS[j]) for marg in levelss for lvl in marg for j in 1:m],
                   synpop = [round(sum(valuex[m*(ihh[1]-1)+j] for ihh in enumerate(pop) if ihh[2][mml[1]] == llvl[2]), digits=2) for mml in zip(MARG_COLS, marginals, levelss) for llvl in enumerate(mml[3]) for j in 1:m],
                   marginal = [mml[2][llvl[1]][j] for mml in zip(MARG_COLS, marginals, levelss) for llvl in enumerate(mml[3]) for j in 1:m])
result.error = result.synpop - result.marginal
result.perror = result.error ./ result.marginal * 100

CSV.write("Optim_Results/assignment-20210915-gq.csv", result)
