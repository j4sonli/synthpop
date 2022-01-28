# synthpop
Synthetic population generation in the United States at scale.

SETUP
------------

### Machine requirements

We recommend 16-32 Gb of memory for most PUMAs; some PUMAs with a large number of tracts or block groups may require 64 Gb or more. 

### Set up R

1. Download R. Please follow the instructions on https://cran.r-project.org/. The scripts have been tested with version 4.1.0. 
2. We recommended downloading RStudio. Please follow the instructions on https://www.rstudio.com/products/rstudio/download/.  
3. Download this repository and save as an R project (in RStudio, File > New Project > Existing Directory). 
4. Once inside the desired R environment for the project (renv if you use it), install the package dependencies by running 
   `install.packages(c("dplyr","tidyr","stringr","R.utils","tidycensus","bnlearn","foreach","doSNOW","ggplot2"))` in the R console. 

In order to use `tidycensus` to download PUMS and ACS data, you will need to US Census API key. To obtain one, visit https://api.census.gov/data/key_signup.html. 

### Set up Julia and Gurobi

1. Download Julia from https://julialang.org/downloads/. The scripts have been tested with version 1.6.1. 
2. Download Gurobi with an academic license following the instructions on https://www.gurobi.com/academia/academic-program-and-licenses/ under "Individual Academic Licenses". You will need a free Gurobi account. The scripts have been tested with Gurobi version 9.1.0. 
3. Open the Julia console and install the dependencies by running the following:

`import Pkg`

`Pkg.add.(["Gurobi", "JuMP", "CSV", "DataFrames", "Random", "RLEVectors", "DelimitedFiles", "StatsBase"])`

5. Note the paths to both your Gurobi installation directory and your activated license file (.lic). 

SETTINGS
------------
Most settings are found in synthpop.R. The ones you will have to configure for your setup before running are described here:
- Rscript_path: Path to the R executable (usually Rscript.exe)
- julia_path: Path to the Julia executable (usually julia.exe)
- my_census_api_key: Your Census API key (as a string)
- gurobi_home: Path to your Gurobi installation directory
- gurobi_license: Path to your Gurobi license file (usually gurobi.lic)

Other settings may optionally be edited:
- hh_pool_file: Output file of BN sampling for sampled household pool
- indv_pool_file: Output file of BN sampling for sampled individual pool
- syn_hhs_file: Output file of synpop selection optimization for households
- syn_indvs_file: Output file of synpop selection optimization for individuals
- syn_hhs_spatial_file: Output file of spatial assignment optimization for households
- syn_indvs_spatial_file: Output file of spatial assignment optimization for individuals
- gurobi_logfile: Path to desired logfile for Gurobi warnings and output

The set of marginals used can be changed in the synpop_selection.jl and spatial_assignment.jl scripts. 

USAGE
------------
The entire synpop generation process can be run with a single command in the terminal:

`Rscript synthpop.R -puma=<XXYYYYY> [-use_prev_marg] [-marg_dir=<"some/path/">] [-samp_rate=<num>] [-suffix=<"abcdefg">] [-ncores=<num>]`

- `puma` is the only required argument, and must be a 7-digit numeric code XXYYYYY where XX is the FIPS state code and YYYYY is the PUMA code. (e.g., 2500506 for Cambridge, MA)
- `use_prev_marg` indicates that the marginals should not be re-downloaded (use when running for a PUMA that has already been run before)
- `marg_dir` specifies the path to the directory in which to store and lookup marginals (default: "synthpop_data/acs_marginals/<puma_code>/")
- `samp_rate` indicates the oversampling rate for households; `-samp_rate=400` means that for a population of 10 households we will sample 40 (default: 400)
- `suffix` is an optional string appended to the output files (default: "")
- `ncores` specifies the number of cores for R to use when sampling (default: 4) 

Occasionally, an error might be encountered early on in process_pums.R that reads 

`Error: Your API call has errors. The API message returned is There was an error while running your query. We've logged the error and we'll correct it ASAP. Sorry for the inconvenience..`

The cause of the error is unknown, but it can usually be resolved by simply running the command again. 

Additionally, note that Gurobi may throw a misleading error that reads "LoadError: result index of attribute mathoptinterace.objectivevalue(1) out of bounds. There are currently 0 solution(s) in the model." if the available memory is insufficient. 