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
    n_steps = end_year-start_year+1

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
# SET KAYA AND SOLOW MODEL 
########################################################################################################
raw_ssp = DataFrame(CSV.File(normpath(@__DIR__,"..", "data", "KayaBackStopTFP.csv"), skipto=2, header = 1))
raw_backstop = DataFrame(YEARS = raw_ssp[!,"YEARS"], Backstop_Price = raw_ssp[!,"Backstop"])
ssp  = raw_ssp[(raw_ssp[!,"YEARS"].>= start_year) .& (raw_ssp[!,"YEARS"] .<= end_year), :]
dice_backstop = raw_backstop[raw_backstop[!,"YEARS"] .>= start_year .&& raw_backstop[!,"YEARS"] .<= end_year, "Backstop_Price"]


#DICE abatement cost
participation = 1
exp_control = 2.6


#Solow
Y_original = ssp[!, "GDP/Pop"].*ssp[!, "Population"]
L_original = ssp[!, "Population"]
Solow_a = 0.3
Solow_A = ssp[!, "Solowtfp"]



Init = DataFrame(CSV.File(normpath(@__DIR__,"..","data", ("Initial_2020.csv")), skipto =2, header = 1))


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


#####################################################################################################
#MODEL
#####################################################################################################
# Turn model into function to run multiple scenarios for policy questions.
# Argument is 'co2_emissions,' a vector of total CO2 emissions (GtCO2).

# run_model = function(CO2_emiss){

function run_model(CO2_emissions, co2_add)



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
    temp_j1 = filldf, temp_j2 = filldf, CO2_rf = filldf, total_rf = filldf, temperature = filldf, 
    # gdp and abatement cost results
    GDP = filldf, abate_cost = filldf, net_GDP = filldf, K = filldf,
    #kaya
    CO2_emiss = filldf, damages = filldf,investment = filldf,consumption = filldf)

    output[:,"CO2_emiss"] = CO2_emissions

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
           

        else
            #--------------------------------#
            #------- GDP EQUATIONS ----------#
            #--------------------------------#

            # Set GDP and K


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
                    output[t, "R"*string(i)] = output[t-1, "R"*string(i)] * exp((-1.0/(tau[i]* output[t,"alpha"]))) + 0.5 * a[i] * (output[t, "CO2_emiss"] + output[t-1, "CO2_emiss"]) / ppm2gtc
                end

            # Calculate the change in CO2 concentrations across all pools and the current atmospheric concentration.
            output[t,"CO2"] = output[t, "R1"] + output[t, "R2"] + output[t, "R3"] + output[t, "R4"] + CO2_0

            # Calculate accumulated perturbation of carbon stock.
            # See top of right column on page 7216
            # Hint, use Cacc[t-1] and just add on the result of that equation for the current period
            # Hint2, Cacc is in GtC, but CO2 concentrations are in ppm, so will need to convert them to GtC (we have a parameter for this).
            output[t,"Cacc"] = output[t-1,"Cacc"] + output[t, "CO2_emiss"] - (output[t,"CO2"] - output[t-1,"CO2"]) * ppm2gtc

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
            
        
            end
            
           

        end

    return(output)
end



# my_results = run_model(fill(0.0,n_steps),0,fill(0.0,n_steps))
# x = my_results[:,"years"]
# y = my_results[:,"temperature"]
# g = my_results[:,"net_GDP"]
# g = my_results[:,"consumption"]
# myd = my_results[:,"damages"]
# plot(x,y,  title = "Global av Temperature above 2015", label = "Our model", legend=:topleft, ylab="degrees C")

# my_results[:,"CO2_emiss"]
# my_results[:,"GDP"]


############
#Co2 emissions
#savings rate
savings_rate = DataFrame(CSV.File(normpath(@__DIR__,"..", "data", "savings_rate.csv"), skipto=2, header = 1))
gdp = DataFrame(CSV.File(normpath(@__DIR__,"..", "data", "gdp.csv"), skipto=2, header = 1))
population = DataFrame(CSV.File(normpath(@__DIR__,"..", "data", "population.csv"), skipto=2, header = 1))

rice = DataFrame(CSV.File(normpath(@__DIR__,"..", "data", "KayaBackStopTFP.csv"), skipto=2, header = 1))



##########
# Damage to GDP function   
#DICE


damage_coefficient = 0.00236
savings = 0.22
real_temp = output[t,"temperature"]
income_elasticity = (output[t,"GDP"]/output[1,"GDP"])^elasticity
output[t,"damages"] = damage_coefficient*real_temp^2 * income_elasticity*output[t,"GDP"]
#println(output[t,"damages"])
output[t,"net_GDP"] = output[t,"net_GDP"] - output[t,"damages"]
output[t,"investment"] = output[t,"net_GDP"]*savings
output[t,"consumption"] = output[t,"net_GDP"] - output[t,"investment"]
percap_cons = output[t,"consumption"]/population







########
# Social Cost of Carbon

add2020 = fill(0.0,n_steps)
add2020[6] = 1
add2020

q1 = run_model(fill(0.0,n_steps),0,add2020)
q1_d = q1[:,"damages"]
d_diff = q1_d[6:end]-myd[6:end]

c = [0.025,0.03,0.05]
costs = []

for i in 1:3
    discounted = []
    for k in 1:length(d_diff)
        discounted = append!(discounted,d_diff[k]/((1+c[i])^(k-1)))
    end
    costs = append!(costs, sum(discounted))
end
costs
co2_costs = costs.*(12/44)



