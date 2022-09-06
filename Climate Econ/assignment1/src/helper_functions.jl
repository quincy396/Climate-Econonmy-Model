using Roots
# library(rootSolve)

# This file contains various functions to help run and analyze an integrated assessment model (IAM).

##############################################################################################
# CALCULATE CARBON UPTAKE SCALING FACTOR
##############################################################################################
#
# This function solves for the roots of an equation to find the state-dependent
# carbon uptake scaling factor from the Finite Amplitude Impulse Response (FAIR)
# carbon cycle & climate model.
#
#---------------------
# Function Inputs:
#---------------------
#           iIRFT100 = 100 year integrated carbon impulse response.
#           a        = Vector with fraction of emissions entering different carbon pools:
#                           (1) geological reabsorption,
#                           (2) deep ocean invasion/equilibration,
#                           (3) biospheric uptake/ocean thermocline invasion,
#                           (4) rapid biospheric uptake/ocean mixed-layer invasion.
#           tau       = Vector with decay time constants for each carbon pool in 'a'.
#           alpha_0   = Initial guess for the value of alpha (to give the solver an idea
#                       of where to start searcing).
#---------------------
# Function Outputs:
#---------------------
#           alpha     = State-dependent carbon uptake scaling factor.
#
#------------------------------------------------------------------------------------------------

function find_alpha(iIRFT100::Float64, a::Vector{Float64}, tau::Vector{Float64}, alpha_0::Float64)
    # Create a function in terms of alpha to pass into solver using Equation 7 from FAIR paper.
    function f(alpha) #????????????????
        # f =
         -iIRFT100 + sum(alpha .* a .* tau .* (1 .- exp.(-100/alpha./tau)))
        # return f
    end
    # Pass f into solve to find alpha.
    alpha = find_zero(f, alpha_0)

     # Return result.
     return alpha
end


# find_alpha = function(iIRFT100, a, tau, alpha_0){

#     # Create a function in terms of alpha to pass into solver using Equation 7 from FAIR paper.
#     f = function(alpha){
#         -iIRFT100 + sum(alpha * a * tau * (1 - exp(-100/alpha/tau)))
#     }

#     # Pass f into solve to find alpha.
#     alpha = multiroot(f, start=alpha_0)$root

#     # Return result.
#     return(alpha)
# }
# #------------------------------------------------------------------------------------------------
