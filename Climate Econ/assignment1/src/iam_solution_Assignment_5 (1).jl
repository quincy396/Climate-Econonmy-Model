# Set you working directory to the folder where you have the 'data' and 'src' folders with cd("path to that folder"), check with pwd()
# The first time you have to run: import Pkg, then Pkg.add(["CSV", "DataFrames", "Roots", "Plots"])
using CSV, DataFrames, Roots, Plots, Statistics

include(joinpath(@__DIR__, "helper_functions.jl"))
############################################################################################################
# MODEL ANALYSIS AND EMISSIONS SCENARIO PARAMETERS TO CHANGE (will add more options as semester progresses).
############################################################################################################

    # Start and end year for the model
    start_year = 2015
    end_year   = 2300
    # Number of timesteps to run model = number of years.
    n_steps = length(start_year:end_year)

    # Set RCP Scenario (current options = "rcp85").
    rcp_scenario = "rcp85"

#####################################################################################################
# SET CLIMATE POLICY
#####################################################################################################
# Base policy is no abatement
co2_policy = fill(0.::Float64, n_steps)


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
   

    #Read in single row csv with data for 2015, initial conditions
    # For some reason it reads a second row without telling it not to with 'footerskip'
    Init = DataFrame(CSV.File(normpath(@__DIR__,"..","data", ("Initial_"*string(start_year)*".csv")), footerskip =1, header = 1))

    #Subtract CO2 RF from Total Anthropogenic RF to avoid double counting.
    exogenous_rf = radforc[:,"TOTAL_INCLVOLCANIC_RF"] .- radforc[:,"CO2_RF"]

    #Get N2O concentrations (used in CO2 radiative forcing calculations).
    N2O_conc = conc[:,"N2O"]


#######################################################################################################
# SET MODEL PARAMETER VALUES
########################################################################################################

    #------------------------------
    # Carbon Cycle
    #------------------------------

    # ----Parameters----#
    CO2_0   = 278.0                             # Pre-industrial atmospheric concentration of CO2.
    r0      = 32.4                              # Pre-industrial iIRF100.
    rC      = 0.019                             # Increase in iIRF100 with cumulative carbon uptake (yr/GtC).
    rT      = 4.165                             # Increase in iIRF100 with warming (yr/C).
    a       = [0.2173, 0.2240, 0.2824, 0.2763] # Fraction of emissions entering each carbon pool (geological reabsorption[1], deep ocean invasion/equilibration[2], biospheric uptake/ocean thermocline invasion[3], rapid biospheric uptake/ocean mixed-layer invasion[4]).
    tau     = [10.0^6, 394.4, 36.54, 4.304]    # Decay time constants for each carbon pool in 'a'.
    ppm2gtc = 2.123                             # Conversion factor between ppm and GtC (with 1 ppm = 2.123 GtC).

    #------------------------------
    # Climate Dynamics
    #------------------------------

    # ----Parameters----#
    a1    = -2.4e-7         # CO2 raditive forcing constant.
    b1    = 7.2e-4          # CO2 raditive forcing constant.
    c1    = -2.1e-4         # CO2 raditive forcing constant.
    N2O_0 = 270       # Pre-industrial atmospheric concentration of N2O.
    q     = [0.33, 0.41]   # q1 (thermal equilibration of deep ocean) & q2 (thermal adjustment of upper ocean) in KW-1 m2.
    d     = [239.0, 4.1]   # d1 (thermal equilibration of deep ocean) & d2 (thermal adjustment of upper ocean) in years.

    #------------------------------
    # Economic & Impact Component (Solow Growth Model)
    #------------------------------
     # ----Parameters----#
    savings  = 0.22         # Savings rate.
    k_share  = 0.3          # Capital share in production function.
    l_share  = 1 - k_share  # Labor share in production function.
    k_deprec = 0.1          # Capital depreciation rate.
    k_0      = 223000      # Initial capital stock ($billions). From DICE

    # generate the appropriate number of missing/blank values
    

    ## Calculate tfp, techical factor productivity (fudge factor) from gross GDP, K, and L as the residual
    Solow = DataFrame(YEARS = [start_year:1:end_year;], K = filldf, GDP = ssp[!,"GDP/Pop"].*ssp[!,"Population"], #/10^3 # $billion/million people * million people = $Billion$
    pop = ssp.Population, tfp = filldf)
    #Set the first row, ready for the loop to fill the rest
    Solow[1,"K"] = k_0; 
    Solow[1,"tfp"] = Solow[1,"GDP"]/(Solow[1,"K"]^k_share*Solow[1,"pop"]^l_share)

    for t in 2:n_steps
        # Calculate K from previous year (+ investment - depreciation)
        Solow[t,"K"] =  Solow[t-1,"K"]*(1-k_deprec) + Solow[t-1,"GDP"]*savings   
        # Solve Solow growth equation for tfp (calculate the tfp for each year)
        Solow[t,"tfp"] = Solow[t,"GDP"]/(Solow[t,"K"]^k_share*Solow[t,"pop"]^l_share)
     end 

#  output a file with the tfp so we can use it later without recalculating
CSV.write(normpath(@__DIR__,"..","data","Solow_Julia.csv"), Solow)

#####################################################################################################
#MODEL
#####################################################################################################
#Turn model into function to run multiple scenarios for policy questions.
#Variable is 'co2_emissions,' a vector of total CO2 emissions (GtCO2).

# run_model = function(CO2_emiss){
function run_model(co2_policy)
    filldf = fill(0.::Float64, end_year-start_year+1)
    # Initialize dataframe to store results for each time step (each row is a year, each column is a new variable).
    # Create the dataframe called 'output', assign model years to a column (just for convenience) and
    output = DataFrame(years = [start_year:1:end_year;], alpha = filldf,
    # Add columns for Carbon Cycle results, and
    R1 = filldf, R2 = filldf, R3 = filldf, R4 = filldf, CO2 = filldf, Cacc = filldf,
    # add for the Climate/Temperature results.
    temp_j1 = filldf, temp_j2 = filldf, CO2_rf = filldf, total_rf = filldf, temperature = filldf,
    # add columns for the new emissions that we will calculate in our model
    Anthrop_co2_emissions = filldf,
    # Abatement cost elements
    net_emiss = filldf, sigma = filldf, policy_cost_coefficient = filldf, policy_costs_share = filldf, policy_costs_dollar = filldf, 
    # Economic factors
    gross_output = filldf, capital = filldf, investment = filldf, net_output = filldf, consumption = filldf, pc_consumption = filldf)

    #Loop through the model
    for t in 1:n_steps

        #############################################
        # Set Initial Conditions for All Components
        #############################################
        if (t==1)

            #-----------------------------#
            #--EMISSIONS INTO THE OUTPUT--#
            #-----------------------------#
            #Initial CO2 emission values (assume no policy in period 1).
            output[t, "Anthrop_co2_emissions"]   = ssp[t,"Emissions/CO2"]
            output[t, "net_emiss"]  = output[t, "Anthrop_co2_emissions"] * (1-co2_policy[t])

            #----------------------------#
            #----INITIAL CARBON CYCLE----#
            #----------------------------#

            # Initialise the carbon pools to be correct for first timestep in numerical method (this is not in the paper, but taken from original FAIR source code).
            output[t, "R1"] = Init[t,"R1"]
            output[t, "R2"] = Init[t,"R2"]
            output[t, "R3"] = Init[t,"R3"]
            output[t, "R4"] = Init[t,"R4"]

            # Initial state-dependent scaling factor.
            output[t,"alpha"] = Init[t,"alpha"]

            # Initical atmospheric CO2 concentration.
            output[t,"CO2"] = Init[t,"CO2"]

            # Initial carbon stock perturbation.
            output[t,"Cacc"] = Init[t,"Cacc"]

            #--------------------------------#
            #----INITIAL CLIMATE DYNAMICS----#
            #--------------------------------#
            output[t,"CO2_rf"] = Init[t,"CO2_rf"]
            output[t,"total_rf"] = Init[t,"total_rf"]

            # Initial temperature change for two response times.
            output[t,"temp_j1"] = Init[t,"temp_j1"]
            output[t,"temp_j2"] = Init[t,"temp_j2"]

            # Set initial surface temperature anomaly.
            output[t,"temperature"] = Init[t,"temperature"]

                        #-----------------------------------#
            #----INITIAL ECONOMIC CONDITIONS----#
            #-----------------------------------#

            # Initial capital stock.
            output[t, "capital"] = k_0

            # Initial gross economic output.
            output[t, "gross_output"] = Solow[t,"tfp"] * Solow[t, "K"]^k_share * Solow[t, "pop"]^l_share

            output[t, "sigma"] = (ssp[t, "Emiss/Ener"] * ssp[t, "Ener/GDP"]) * (44/12) #* 1000
            
            # Calculate abatement cost coefficients (following DICE equations, i.e. that's where the 1000 comes from).
            output[t, "policy_cost_coefficient"] = dice_backstop[t] * output[t, "sigma"] / 2.6 # /1000
            
            # Calculate abatement costs as a share of gross economic output.
            output[t, "policy_costs_share"] = output[t, "policy_cost_coefficient"] * co2_policy[t] ^ 2.6
            
            # Calculate abatement costs in billions.
            output[t, "policy_costs_dollar"] = output[t, "gross_output"] * output[t, "policy_costs_share"]
            
            output[t,"net_output"] = output[t,"gross_output"] - output[t, "policy_costs_dollar"]
            # Initial investment level.
            output[t, "investment"] = output[t, "net_output"] * savings

            # Initial aggregate consumption.
            output[t, "consumption"] = output[t, "net_output"] - output[t, "investment"]

            # Iniitial per capita consumption values (dollars / person).
            output[t, "pc_consumption"] = output[t, "consumption"] / ssp[t, "Population"] * 1000


        else
            #--------------------------------#
            #---- GROSS OUTPUT EQUATIONS ----#
            #--------------------------------#

            # Capital stock ($billions).
            output[t, "capital"] = (1 - k_deprec) * output[t-1, "capital"] + output[t-1, "investment"]

            # Gross output ($billions).
            output[t, "gross_output"] = Solow[t,"tfp"] * output[t, "capital"]^k_share * ssp[t, "Population"]^l_share
            
            #--------------------------------#
            #---- CO2 EMISSION EQUATIONS ----#
            #--------------------------------#

            # Add Total CO2 emissions based on Kaya factors and K from the model!, assuming that OtherCO2 is also included in the Kaya emissions.
            output[t, "Anthrop_co2_emissions"] = output[t, "gross_output"] * ssp[t, "Ener/GDP"] * ssp[t, "Emiss/Ener"]

            # Now we use the net emissions as out emissions in the model
            output[t,"net_emiss"] = output[t, "Anthrop_co2_emissions"]*(1-co2_policy[t])

            #--------------------------------#
            #---- CARBON CYCLE EQUATIONS ----#
            #--------------------------------#
            # Equation 8 in FAIR
            iIRFT100 = r0 + rC * output[t-1,"Cacc"] + rT * output[t-1,"temperature"]

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
                    output[t, "R"*string(i)] = output[t-1, "R"*string(i)] * exp((-1.0/(tau[i]* output[t,"alpha"]))) + 0.5 * a[i] * (output[t, "net_emiss"] + output[t-1, "net_emiss"]) / ppm2gtc
                end

            # Calculate the change in CO2 concentrations across all pools and the current atmospheric concentration.
            output[t,"CO2"] = output[t, "R1"] + output[t, "R2"] + output[t, "R3"] + output[t, "R4"] + CO2_0

            # Calculate accumulated perturbation of carbon stock.
            # See top of right column on page 7216
            # Hint, use Cacc[t-1] and just add on the result of that equation for the current period
            # Hint2, Cacc is in GtC, but CO2 concentrations are in ppm, so will need to convert them to GtC (we have a parameter for this).
            output[t,"Cacc"] = output[t-1,"Cacc"] + output[t, "net_emiss"] - (output[t,"CO2"] - output[t-1,"CO2"]) * ppm2gtc

            #------------------------------------#
            #---- CLIMATE DYNAMICS EQUATIONS ----#
            #------------------------------------#

            N_hat = 0.5 * (N2O_conc[t] + N2O_0)

            #COâ Radiative Forcing
            CO2_diff = output[t,"CO2"] - CO2_0
            output[t,"CO2_rf"] = (a1*CO2_diff^2 + b1*abs(CO2_diff) + c1*N_hat + 5.36) * log(output[t,"CO2"] / CO2_0)
            output[t,"total_rf"] = output[t,"CO2_rf"] + exogenous_rf[t]

            # Calculate temperature change for the two different thermal response times.
            output[t,"temp_j1"] = output[t-1,"temp_j1"] * exp((-1.0)/d[1]) + 0.5 * q[1] * (output[t-1,"total_rf"] + output[t,"total_rf"]) * (1 - exp((-1.0)/d[1]))
            output[t,"temp_j2"] = output[t-1,"temp_j2"] * exp((-1.0)/d[2]) + 0.5 * q[2] * (output[t-1,"total_rf"] + output[t,"total_rf"]) * (1 - exp((-1.0)/d[2]))

            #Calculate global mean surface temperature anomaly
            output[t,"temperature"] = output[t,"temp_j1"] + output[t, "temp_j2"]

            #------------------------------#
            #---- CO2 MITIGATION COSTS ----#
            #------------------------------#

            # Calculate sigma term for CO2 mitigation costs (scaled to match DICE units).
            # Sigma is Emissions-output ratio; MtCO2 per $1000, emissions per $GDP
            # Note on Units: gigatonne = 10^9 tonnes, trillion = 10^12. SSP data = GtC / Billion.
            output[t, "sigma"] = ssp[t, "Emiss/Ener"] * ssp[t, "Ener/GDP"]  * (44/12) #* 1000
            # Use sigma and the backstop price data to calculate the Policy/Abatement cost coefficients
            # Calculate abatement cost coefficients (following DICE equations, i.e. that's where the 1000 comes from).
            output[t, "policy_cost_coefficient"] = ssp[t, "Backstop"] * output[t, "sigma"] / 2.6 #/ 1000
            
            # Calculate abatement costs as a share of gross economic output.
            output[t, "policy_costs_share"] = output[t, "policy_cost_coefficient"] * co2_policy[t] ^ 2.6
            
            # Calculate abatement costs in billions.
            output[t, "policy_costs_dollar"] = output[t, "gross_output"] * output[t, "policy_costs_share"]

            #-----------------------------#
            #---- NET ECONOMIC OUTPUT ----#
            #-----------------------------#

            # Net economic output ($billions).
            output[t, "net_output"] = output[t, "gross_output"] - output[t, "policy_costs_dollar"]

            # Investment based on net output and a constant savings rate ($billions).
            output[t, "investment"] = output[t, "net_output"] * savings

            # Consumption ($billions).
            output[t, "consumption"] = output[t, "net_output"] - output[t, "investment"]

            # Per capita consumption (x1000 converts units to dollars/person).
            output[t, "pc_consumption"] = output[t, "consumption"] / ssp[t, "Population"] * 1000
        end
    end
    return(output)
end

# Run the model with no mitigation as a baseline

MyResults = run_model(co2_policy)
MyResults[:,"net_output"]
## Create mitigation policy vector index 1 = 2015
co2_policy30 = copy(co2_policy); co2_policy30[6] = .3 
MyResultsPolicy = run_model((co2_policy30))

plot(MyResults[:,"years"],  MyResults[:,"temperature"], label = "No Policy", legend=:topleft)
plot!(MyResultsPolicy[:,"years"],  MyResultsPolicy[:,"temperature"], label = "30% mitigation in 2021")

plot(MyResultsPolicy[:,"years"],  MyResults[:,"policy_costs_dollar"], main = "Cost of policy", label= "no cost, as check",legend=:topleft)
plot!(MyResultsPolicy[:,"years"],  MyResultsPolicy[:,"policy_costs_dollar"], label = "Policy")

plot(MyResults[:,"years"],  MyResults[:,"net_output"], main = "Net Economic Output", label= "Without policy",legend=:topleft)
plot!(MyResultsPolicy[:,"years"],  MyResultsPolicy[:,"net_output"], label= "With policy")

# Substract the net output with policy from the net output without policy to see the difference which is the cost of the policy itself
abatementcost = DataFrame(years = MyResultsPolicy[:,"years"], Abatement_Cost = MyResults[:,"net_output"] .- MyResultsPolicy[:,"net_output"])
# Calculate the year when the abatement cost is < 1 (billion) and after 2019 by filtering for both conditions (gets a vector of years), and then taking the first value
# & Do the calcluation as part of a string to use as the title of the plot
title = "Year when abatement costs \n will be < 1 billion\$\$: "*string(abatementcost[(abatementcost.years .>2019) .& (abatementcost.Abatement_Cost.<1),"years"][1])


plot(abatementcost.years, abatementcost.Abatement_Cost)

# Plot of abatement cost leaving out 2020 (for better y-range to see the difference)
plot(abatementcost.years[7:286], abatementcost.Abatement_Cost[7:286], xlab = "years", ylab = "billion\$\$", title = title, label = "cost from 2021")
