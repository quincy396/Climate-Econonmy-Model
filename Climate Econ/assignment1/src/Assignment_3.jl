# Set you working directory to the folder where you have the 'data' and 'src' folders with cd("path to that folder"), check with pwd()
# The first time you have to run: import Pkg, then Pkg.add(["CSV", "DataFrames", "Roots", "Plots"])

using CSV, DataFrames, Roots, Plots, Interpolations

include(joinpath(@__DIR__, "helper_functions.jl"))
############################################################################################################
# MODEL ANALYSIS AND EMISSIONS SCENARIO PARAMETERS TO CHANGE (will add more options as semester progresses).
############################################################################################################

    # Start and end year for the model
    start_year = 2015
    end_year   = 2300

    # Set RCP Scenario (current options = "rcp85").
    
    rcp_scenario = "rcp85"

#####################################################################################################
# READ IN MODEL INPUT DATA
#####################################################################################################

    #Read in radiative forcing and emissions data for RCP Scenario.
    raw_radforc8 = CSV.File(normpath(@__DIR__,"..","data", (rcp_scenario * "_MIDYEAR_RADFORCING.csv")), skipto=60, header = 59) |> DataFrame
    raw_emiss8   = CSV.File(normpath(@__DIR__,"..","data", (rcp_scenario * "_EMISSIONS.csv")), skipto=38, header = 37) |> DataFrame
    raw_conc8    = CSV.File(normpath(@__DIR__,"..","data", (rcp_scenario * "_MIDYEAR_CONCENTRATIONS.csv")), skipto=39, header = 38) |> DataFrame

    #Isolate data for the years model will be run.
    radforc8 = raw_radforc8[(raw_radforc8[!,"v YEARS/GAS >"].>=start_year) .& (raw_radforc8[!,"v YEARS/GAS >"].<=end_year),:]
    emiss8 = raw_emiss8[(raw_emiss8[!,"v YEARS/GAS >"].>=start_year) .& (raw_emiss8[!,"v YEARS/GAS >"].<=end_year),:]
    conc8 = raw_conc8[(raw_conc8[!,"v YEARS/GAS >"].>=start_year) .& (raw_conc8[!,"v YEARS/GAS >"].<=end_year),:]

    #Subtract CO2 RF from Total Anthropogenic RF to avoid double counting.
    exogenous_rf8 = radforc8[!,"TOTAL_INCLVOLCANIC_RF"] - radforc8[!,"CO2_RF"]

    #Add fossil fuel and land use change + other sources CO2 emissions together.
    co2_emissions8 = emiss8[!, "FossilCO2"] + emiss8[!, "OtherCO2"]

    #Get N2O concentrations (used in CO2 radiative forcing calculations).
    N2O_conc8 = conc8[!, "N2O"]
    #N2O_conc[3]

    co2_emissions = co2_emissions8
   

#######################################################################################################
# SET MODEL PARAMETER VALUES
########################################################################################################

    # Number of timesteps to run model.
    n_steps = length(start_year:end_year)

    #------------------------------
    # Carbon Cycle
    #------------------------------

    # ----Parameters----#
    CO2_0   = 278.0                             # Pre-industrial atmospheric concentration of CO2.
    r0      = 32.40                               # Pre-industrial iIRF100.
    rC      = 0.019                               # Increase in iIRF100 with cumulative carbon uptake (yr/GtC).
    rT      = 4.165                             # Increase in iIRF100 with warming (yr/C).
    a       = [0.2173,  0.2240, 0.2824, 0.2763]                  # Fraction of emissions entering each carbon pool (geological reabsorption[1], deep ocean invasion/equilibration[2], biospheric uptake/ocean thermocline invasion[3], rapid biospheric uptake/ocean mixed-layer invasion[4]).
    tau     = [10^6, 394.4, 36.54, 4.304]                  # Decay time constants for each carbon pool in 'a'.
    ppm2gtc = 2.123                             # Conversion factor between ppm and GtC (with 1 ppm = 2.123 GtC).
    
    #------------------------------
    # Climate Dynamics
    #------------------------------

    # ----Parameters----#
    a1    = -2.410*10^-7             # CO2 raditive forcing constant.
    b1    = 7.2*10^-4             # CO2 raditive forcing constant.
    c1    = -2.1*10^-4             # CO2 raditive forcing constant.
    N2O_0 = 270             # Pre-industrial atmospheric concentration of N2O.
    q     = [ 0.33, 0.41]        # q1 (thermal equilibration of deep ocean) & q2 (thermal adjustment of upper ocean) in KW-1 m2.
    d     = [230, 4.1]        # d1 (thermal equilibration of deep ocean) & d2 (thermal adjustment of upper ocean) in years.


#####################################################################################################
#MODEL
#####################################################################################################
# Turn model into function to run multiple scenarios for policy questions.
# Argument is 'co2_emissions,' a vector of total CO2 emissions (GtCO2).

# run_model = function(CO2_emiss){
function run_model(CO2_emiss)
    N2O_conc = N2O_conc8
    exogenous_rf = exogenous_rf8
    



    # Initialize dataframe to store results for each time step (each row is a year, each column is a new variable).
    # First set a variable to generate the appropriate number of missing/blank values
    filldf = fill(0.::Float64, end_year-start_year+1)
    # Create the dataframe called 'output', assign model years to a column (just for convenience)
    output = DataFrame(years = [start_year:1:end_year;], #and
    # Add columns for Carbon Cycle results, and
    alpha = filldf, R1 = filldf, R2 = filldf, R3 = filldf, R4 = filldf, CO2 = filldf, Cacc = filldf,
    # add for the Climate/Temperature results.
    temp_j1 = filldf, temp_j2 = filldf, CO2_rf = filldf, total_rf = filldf, temperature = filldf)

    #Loop through the model
    for t in 1:n_steps

        #############################################
        # Set Initial Conditions for All Components
        #############################################
        if (t==1)

            #----------------------------#
            #----INITIAL CARBON CYCLE----#
            #----------------------------#

            # Initialise the carbon pools to be correct for first timestep in numerical method (this is not in the paper, but taken from original FAIR source code).
            output[t, "R1"] = a[1] * CO2_emiss[t] / ppm2gtc * 0.5
            output[t, "R2"] = a[2] * CO2_emiss[t] / ppm2gtc * 0.5
            output[t, "R3"] = a[3] * CO2_emiss[t] / ppm2gtc * 0.5
            output[t, "R4"] = a[4] * CO2_emiss[t] / ppm2gtc * 0.5

            # Initial state-dependent scaling factor.
            output[t,"alpha"] = 1e-10

            # Initial atmospheric CO2 concentration.
            output[t,"CO2"] = CO2_0

            # Initial carbon stock perturbation.
            output[t,"Cacc"] = CO2_emiss[t]

            #--------------------------------#
            #----INITIAL CLIMATE DYNAMICS----#
            #--------------------------------#
            output[t,"CO2_rf"] = 0.0
            output[t,"total_rf"] = 0.0

            # Initial temperature change for two reponse times.
            output[t,"temp_j1"] = 0.0
            output[t,"temp_j2"] = 0.0

            # Set initial surface temperature anomaly to 0.0.
            output[t,"temperature"] = 0.0

        else
            #--------------------------------#
            #---- CARBON CYCLE EQUATIONS ----#
            #--------------------------------#
            # Equation 8 in FAIR
            # Hint: You will need to use "Cacc" and "temperature" output from period 't-1'.
            iIRFT100 = r0 +rC*output[t-1,"Cacc"] + rT*output[t-1,"temperature"]

            # Set an upper bound to avoid unstable/non-physical results.  Not in paper, taken from original FAIR source code.
            if iIRFT100 >= 97.0 
                # Set alpha to it's limiting value.
                output[t,"alpha"] = 113.7930278
            else 
                # Solve for alpha, given current state of carbon and climate systems.
                # Use the previous value of alpha as initial guess because alpha varies relatively smoothly over time.
                output[t,"alpha"] = find_alpha(iIRFT100, a, tau, output[t-1,"alpha"])
            end

            #Calculate updated carbon cycle time constants and CO2 concentrations in 4 carbon pools.
            for i in 1:4
                output[t, "R"*string(i)] = output[t-1, "R"*string(i)] * exp((-1.0/(tau[i]* output[t,"alpha"]))) + 0.5 * a[i] * (CO2_emiss[t] + CO2_emiss[t-1]) / ppm2gtc
            end

            # Calculate the change in CO2 concentrations across all pools and the current atmospheric concentration.
            output[t,"CO2"] = CO2_0 + output[t,"R1"] + output[t,"R2"] + output[t,"R3"] + output[t,"R4"]
            #println(output[t,"CO2"] - CO2_0)


            # Calculate accumulated perturbation of carbon stock.
            # See top of right column on page 7216
            # Hint, use Cacc[t-1] and just add on the result of that equation for the current period
            # Hint2, Cacc is in GtC, but CO2 concentrations are in ppm, so will need to convert them to GtC (we have a parameter for this).
            output[t,"Cacc"] = output[t-1,"Cacc"] + (CO2_emiss[t] - (output[t,"CO2"] - output[t-1,"CO2"])) * ppm2gtc
            #println(output[t,"Cacc"])
            #println((sum(CO2_emiss[1:t]) - (output[t,"CO2"] - CO2_0)))

            #------------------------------------#
            #---- CLIMATE DYNAMICS EQUATIONS ----#
            #------------------------------------#
# Calculate N_hat term and difference in CO2 concentrations as temporary variables (for convenience).
# HINT: Use "N_hat" and "CO2_diff" when calculatint "CO2_rf" just to make your code more compact & easier to read.
            N_hat = 0.5 * (N2O_conc[t] + N2O_0)
            CO2_diff = output[t,"CO2"] - CO2_0

           # Calculate CO2 radiative forcing
# HINT: See Table 1 in the the paper, "Radiative forcing of carbon dioxide, methane, and nitrous oxide: a significant revision of the methane radiative forcing."
            output[t,"CO2_rf"] = (a1*(CO2_diff)^2 +b1*abs(CO2_diff) + c1*N_hat +5.36)*log(output[t,"CO2"]/CO2_0)
            # Calculate total radiative forcing.
# HINT: Don't forget about exogenous sources of radiative forcing that we are not explicitly calculating.
            output[t,"total_rf"] = output[t,"CO2_rf"] + exogenous_rf[t]

            # Calculate temperature change for the two different thermal response times.
            output[t,"temp_j1"] = output[t-1,"temp_j1"] * exp((-1.0)/d[1]) + 0.5 * q[1] * (output[t-1,"total_rf"] + output[t,"total_rf"]) * (1 - exp((-1.0)/d[1]))
            output[t,"temp_j2"] = output[t-1,"temp_j2"] * exp((-1.0)/d[2]) + 0.5 * q[2] * (output[t-1,"total_rf"] + output[t,"total_rf"]) * (1 - exp((-1.0)/d[2]))

            #Calculate global mean surface temperature anomaly
#LESS USEFUL HINT: the equation is in the FAIR paper.
            output[t,"temperature"] = output[t,"temp_j1"] + output[t,"temp_j2"]
            #println(output[t,"temperature"])
        
            end
        end

    return(output)
end

my_results = run_model(co2_emissions)

x = my_results[1:172,"years"]
y = my_results[1:172,"temperature"]
plot(x,y,  title = "Global av Temperature above pre-industrial", label = "Our model", legend=:topleft, ylab="degrees C")

my_results[156,"years"]

ssp = CSV.File(normpath(@__DIR__,"..","data", ("iamc_db_total.csv")), skipto=6, header = 4) |> DataFrame
ssp_itp = copy(ssp)

itp = LinearInterpolation(ssp_itp[:,1], ssp_itp[:,3])
itp(2006)