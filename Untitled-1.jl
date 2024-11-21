#%% open data 
#%% Open files 
using XLSX, DataFrames

# Set the working directory to the folder where your files are located
cd("C:/Users/gealy/OneDrive/Documents/git/great dsge/great_dsge")

# Load "3-MonthTreasurey.xlsx"
try
    global t3_MonthTreasurey = XLSX.readtable("3-MonthTreasurey.xlsx", 1) |> DataFrame
    println("Successfully loaded: 3-MonthTreasurey.xlsx")
catch e
    println("Error loading 3-MonthTreasurey.xlsx: ", e)
end

# Load "30_YearFixedRateMortgageAverage.xlsx"
try
    global t30_YearFixedRateMortgageAverage = XLSX.readtable("30_YearFixedRateMortgageAverage.xlsx", 1) |> DataFrame
    println("Successfully loaded: 30_YearFixedRateMortgageAverage.xlsx")
catch e
    println("Error loading 30_YearFixedRateMortgageAverage.xlsx: ", e)
end

# Load "Depositshousehold.xlsx"
try
    global Depositshousehold = XLSX.readtable("Depositshousehold.xlsx", 1) |> DataFrame
    println("Successfully loaded: Depositshousehold.xlsx")
catch e
    println("Error loading Depositshousehold.xlsx: ", e)
end

# Load "PCE.xlsx"
try
    global PCE = XLSX.readtable("PCE.xlsx", 1) |> DataFrame
    println("Successfully loaded: PCE.xlsx")
catch e
    println("Error loading PCE.xlsx: ", e)
end

# Load "PNFI.xlsx"
try
    global PNFI = XLSX.readtable("PNFI.xlsx", 1) |> DataFrame
    println("Successfully loaded: PNFI.xlsx")
catch e
    println("Error loading PNFI.xlsx: ", e)
end

# Load "ResidentialMortgages.xlsx"
try
    global ResidentialMortgages = XLSX.readtable("ResidentialMortgages.xlsx", 1) |> DataFrame
    println("Successfully loaded: ResidentialMortgages.xlsx")
catch e
    println("Error loading ResidentialMortgages.xlsx: ", e)
end

# Load "comph.xlsx"
try
    global comph = XLSX.readtable("comph.xlsx", 1) |> DataFrame
    println("Successfully loaded: comph.xlsx")
catch e
    println("Error loading comph.xlsx: ", e)
end

# Load "corporatelendingxlsx.xlsx"
try
    global corporatelendingxlsx = XLSX.readtable("corporatelendingxlsx.xlsx", 1) |> DataFrame
    println("Successfully loaded: corporatelendingxlsx.xlsx")
catch e
    println("Error loading corporatelendingxlsx.xlsx: ", e)
end

# Load "housepriceindex.xlsx"
try
    global housepriceindex = XLSX.readtable("housepriceindex.xlsx", 1) |> DataFrame
    println("Successfully loaded: housepriceindex.xlsx")
catch e
    println("Error loading housepriceindex.xlsx: ", e)
end

# Load "inflation.xlsx"
try
    global inflation = XLSX.readtable("inflation.xlsx", 1) |> DataFrame
    println("Successfully loaded: inflation.xlsx")
catch e
    println("Error loading inflation.xlsx: ", e)
end

# Load "moodyCorporateBondYield.xlsx"
try
    global moodyCorporateBondYield = XLSX.readtable("moodyCorporateBondYield.xlsx", 1) |> DataFrame
    println("Successfully loaded: moodyCorporateBondYield.xlsx")
catch e
    println("Error loading moodyCorporateBondYield.xlsx: ", e)
end
#%% Remove spaces before and after 
dfs = [
    t3_MonthTreasurey,
    t30_YearFixedRateMortgageAverage,
    Depositshousehold,
    PCE,
    PNFI,
    ResidentialMortgages,
    comph,
    corporatelendingxlsx,
    housepriceindex,
    inflation,
    moodyCorporateBondYield
]
# Loop through each DataFrame in your list and remove spaces from the second column name
for df in dfs
    try
        # Get the current column names
        col_names = names(df)

        # Check if there is a second column and remove spaces if necessary
        if length(col_names) >= 2
            new_name = strip(col_names[2])  # Remove leading and trailing spaces
            rename!(df, col_names[2] => new_name)
            println("Renamed second column to: '$new_name'")
        end
    catch e
        println("Error processing DataFrame: ", e)
        continue
    end
end

# Check the column names after renaming
for df in dfs
    println(names(df))
end

#%% remove all spaces in between  
using DataFrames

# List of DataFrame variables to process
dfs = [
    t3_MonthTreasurey,
    t30_YearFixedRateMortgageAverage,
    Depositshousehold,
    PCE,
    PNFI,
    ResidentialMortgages,
    comph,
    corporatelendingxlsx,
    housepriceindex,
    inflation,
    moodyCorporateBondYield
]
 
# Loop through each DataFrame in your list and remove spaces from the second column name
for df in dfs
    try
        # Get the current column names
        col_names = names(df)

        # Check if there is a second column and remove all spaces
        if length(col_names) >= 2
            new_name = replace(col_names[2], " " => "")  # Remove all spaces
            rename!(df, col_names[2] => new_name)
            println("Renamed second column to: '$new_name'")
        end
    catch e
        println("Error processing DataFrame: ", e)
        continue
    end
end
#%%
# Verify the updated column names
for df in dfs
    println(names(df))
end 

#%%data handling calculate home price growth rate 
# PCE Multiply only the second column by a billion
PCE[!, 2] .= PCE[!, 2] .* 1_000_000_000
#PNFI
PNFI[!, 2] .= PNFI[!, 2] .* 1_000_000_000
# residential mortgages 
ResidentialMortgages[!,2]  .= ResidentialMortgages[!,2] .* 1_000_000_000 
# Depositshousehold 
Depositshousehold[!,2] .=  Depositshousehold[!,2] .* 1_000_000_000 
#%% figure out the peak housing market 
using DataFrames, ShiftedArrays

# Calculate the growth rate for the 'House_Price_Index' column
lagged_values = ShiftedArrays.lag(housepriceindex[!, "HousePriceIndex"], 1)
growth_rates = [missing; diff(housepriceindex[!, "HousePriceIndex"]) ./ lagged_values[2:end] .* 100]

# Add the growth rate column to the DataFrame
housepriceindex[:, :growthrate] = growth_rates

#%%make sure all dates aline 

using Dates

# List of DataFrame variables to process
dfs = [
    t3_MonthTreasurey,
    t30_YearFixedRateMortgageAverage,
    Depositshousehold,
    PCE,
    PNFI,
    ResidentialMortgages,
    comph,
    corporatelendingxlsx,
    housepriceindex,
    inflation,
    moodyCorporateBondYield
]

 
#%% rename variable to match the model 
using DataFrames

# Define the mapping of original dataset names to model variable names
rename_map = Dict(
    "PCE" => "C",                          # Personal Consumption Expenditures
    "PNFI" => "I",                         # Private Nonresidential Fixed Investment
    "inflation" => "π",                    # Inflation (Gross Domestic Product: Implicit Price Deflator)
    "comph" => "W",                        # Nonfarm Business Sector: Real Compensation Per Hour
    "Depositshousehold" => "h",            # Housing Services / Loans to Households (HHMSDODNS)
    "Depositshousehold" => "b",            # Loans to Households (HHMSDODNS)
    "corporatelendingxlsx" => "L",         # Corporate Lending (MLBSNNCB)
    "Depositshousehold" => "D",            # Deposits (DABSHNO)
    "t3_MonthTreasurey" => "r",            # Policy Rate (TB3MS)
    "t30_YearFixedRateMortgageAverage" => "rh", # Mortgage Rates (MORTG)
    "moodyCorporateBondYield" => "rc",     # Corporate Lending Rates (BAA)
    "housepriceindex" => "qh",             # House Price Index (USSTHPI)
    "ResidentialMortgages" => "ht" # Housing Services (ht)
)


# Function to rename both the second column and the DataFrame variable
function rename_dataframe!(df_name::String, new_name::String)
    try
        # Get the DataFrame using its current name
        df = getfield(Main, Symbol(df_name))

        # Rename the second column (assume the first column is the date)
        second_col = names(df)[2]
        rename!(df, second_col => new_name)

        # Reassign the DataFrame to a new global variable with the new name
        @eval global $(Symbol(new_name)) = df

        println("Renamed DataFrame '$df_name' to '$new_name' and column to '$new_name'")
    catch e
        println("Error processing DataFrame '$df_name': ", e)
    end
end

# Loop through the rename map and apply the renaming function
for (old_name, new_name) in rename_map
    rename_dataframe!(old_name, new_name)
end

#%% make the time frequency match 
using DataFrames, Dates



# Step 2: Define target quarterly dates (January 1st, April 1st, July 1st, October 1st)
function get_quarterly_data(df)
    quarterly_dates = [Date(year, month, 1) for year in year(minimum(df.observation_date)):year(maximum(df.observation_date))
                                               for month in [1, 4, 7, 10]]
    
    # Step 3: Find closest date in the DataFrame for each quarterly date
    result = DataFrame()
    for target_date in quarterly_dates
        # Find the closest date in the DataFrame
        closest_row = argmin(abs.(df.observation_date .- target_date))
        row = df[closest_row, :]
        row[:observation_date] = target_date  # Assign the target quarterly date
        push!(result, row)
    end
    
    return result
end

# Step 4: Extract the quarterly data
t30_YearFixedRateMortgageAveragequarterly = get_quarterly_data(t30_YearFixedRateMortgageAverage)

# Step 5: Format the dates as requested (e.g., "1947-01-01")
t30_YearFixedRateMortgageAveragequarterly.observation_date = Dates.format.(t30_YearFixedRateMortgageAveragequarterly.observation_date, "yyyy-mm-dd")

#%% Merge dfs


using Chain, DataFrames

# Helper function to get the first column name as a Symbol
first_col(df) = Symbol(names(df)[1])

# Helper function to get the second column name as a Symbol
second_col(df) = Symbol(names(df)[2])

# Merge DataFrames using Chain and selecting only the second column of each DataFrame
mergeddf = @chain t3_MonthTreasurey begin
    select(first_col(t3_MonthTreasurey), second_col(t3_MonthTreasurey))
    innerjoin(select(t30_YearFixedRateMortgageAverage, first_col(t30_YearFixedRateMortgageAverage), second_col(t30_YearFixedRateMortgageAverage)), on = first_col(t30_YearFixedRateMortgageAverage))
    innerjoin(select(Depositshousehold, first_col(Depositshousehold), second_col(Depositshousehold)), on = first_col(Depositshousehold))
    innerjoin(select(PCE, first_col(PCE), second_col(PCE)), on = first_col(PCE))
    innerjoin(select(PNFI, first_col(PNFI), second_col(PNFI)), on = first_col(PNFI))
    innerjoin(select(ResidentialMortgages, first_col(ResidentialMortgages), second_col(ResidentialMortgages)), on = first_col(ResidentialMortgages))
    innerjoin(select(comph, first_col(comph), second_col(comph)), on = first_col(comph))
    innerjoin(select(corporatelendingxlsx, first_col(corporatelendingxlsx), second_col(corporatelendingxlsx)), on = first_col(corporatelendingxlsx))
    innerjoin(select(housepriceindex, first_col(housepriceindex), second_col(housepriceindex)), on = first_col(housepriceindex))
    innerjoin(select(inflation, first_col(inflation), second_col(inflation)), on = first_col(inflation))
    innerjoin(select(moodyCorporateBondYield, first_col(moodyCorporateBondYield), second_col(moodyCorporateBondYield)), on = first_col(moodyCorporateBondYield))
end
#%% solve contained 
#%% create subsets 
using Dates, DataFrames

# Define the date ranges
start_date1 = Date("1991-01-01")
end_date1 = Date("2006-01-01")

start_date2 = Date("1995-01-01")
end_date2 = Date("2008-01-01")

# Subset 1: From 1991-01-01 to 2006-01-01
target = filter(row -> row[:observation_date] >= start_date1 && row[:observation_date] <= end_date1, mergeddf)

# Subset 2: From 1995-01-01 to 2013-01-01
test = filter(row -> row[:observation_date] >= start_date2 && row[:observation_date] <= end_date2, mergeddf)

#%% contained maximizations focs to plug in to the dsge
#%% Patient household 
using Symbolics

# Define time variable
@variables t

# Define parameters and variables
@syms β_P a φ ε_z(t) ε_h(t)
@syms c_P(t) h_P(t) d_P(t) q_h(t) w_P(t) l_P(t) r(t) pi_var(t) λ_P

# Define the utility function U_P
U_P = β_P^t * ((1 - a) * ε_z(t) * log(c_P(t) - a * c_P(t - 1)) +
               ε_h(t) * log(h_P(t)) - 
               (h_P(t)^(1 + φ)) / (1 + φ))

# Define the budget constraint BC_P
BC_P = w_P(t) * l_P(t) + (1 + r(t - 1)) / pi_var(t) * d_P(t - 1) - 
       (c_P(t) + q_h(t) * h_P(t) + d_P(t))

# Define the Lagrangian L_P
@syms λ_P
L_P = U_P + λ_P * BC_P

# Compute the FOCs
FOC_c_P = expand_derivatives(Differential(c_P(t))(L_P)) ~ 0
FOC_h_P = expand_derivatives(Differential(h_P(t))(L_P)) ~ 0
FOC_d_P = expand_derivatives(Differential(d_P(t))(L_P)) ~ 0

# Solve for λ_P using the FOCs
λ_P_from_c_P = symbolic_linear_solve(FOC_c_P, λ_P)
λ_P_from_h_P = symbolic_linear_solve(FOC_h_P, λ_P)

# Substitute λ_P into the system to eliminate it
FOC_h_P_no_lambda = substitute(λ_P_from_h_P, λ_P => λ_P_from_c_P)

# Display Results
println("FOC with respect to c_P(t) (Patient households):")
display(FOC_c_P)

println("\nFOC with respect to h_P(t) (Patient households):")
display(FOC_h_P)

println("\nFOC with respect to d_P(t) (Patient households):")
display(FOC_d_P)

println("\nSimplified system without λ_P:")
display(FOC_h_P_no_lambda)

#%%linearize patient households
using Symbolics

# Define variables and parameters
@syms β_P ε_z(t) ε_h(t) φ a q_h(t) h_P(t) c_P(t) d_P(t)
@syms c_P_ss h_P_ss q_h_ss d_P_ss ε_z_ss ε_h_ss λ_P_ss  # Steady-state values
@syms Δc_P Δh_P Δq_h Δd_P  # Deviations from steady state

# Define FOCs (Simplified System)
FOC_c_P = ((1 - a) * (β_P^t) * ε_z(t)) / (c_P(t) - a * c_P(t-1)) - λ_P_ss
FOC_h_P = (β_P^t) * (ε_h(t) / h_P(t) - (h_P(t)^φ)) - q_h(t) * λ_P_ss
FOC_d_P = -λ_P_ss

# Steady-state conditions (formatted as substitutions)
steady_state_subs = Dict(
    ε_z(t) => ε_z_ss,
    ε_h(t) => ε_h_ss,
    c_P(t) => c_P_ss,
    h_P(t) => h_P_ss,
    q_h(t) => q_h_ss,
    d_P(t) => d_P_ss
)

# Define linearized variables (formatted as substitutions)
linearized_vars = Dict(
    c_P(t) => c_P_ss * (1 + Δc_P),  # Linearize c_P around steady state
    h_P(t) => h_P_ss * (1 + Δh_P),  # Linearize h_P around steady state
    q_h(t) => q_h_ss * (1 + Δq_h),  # Linearize q_h around steady state
    d_P(t) => d_P_ss * (1 + Δd_P)   # Linearize d_P around steady state
)

# Substitute steady state and deviations into FOCs
FOC_c_P_linearized = substitute(substitute(FOC_c_P, steady_state_subs), linearized_vars)
FOC_h_P_linearized = substitute(substitute(FOC_h_P, steady_state_subs), linearized_vars)
FOC_d_P_linearized = substitute(substitute(FOC_d_P, steady_state_subs), linearized_vars)

# Simplify the linearized equations
FOC_c_P_simplified = expand(FOC_c_P_linearized)
FOC_h_P_simplified = expand(FOC_h_P_linearized)
FOC_d_P_simplified = expand(FOC_d_P_linearized)

# Display Results
println("Linearized FOC with respect to c_P(t) (Patient households):")
display(FOC_c_P_simplified)

println("\nLinearized FOC with respect to h_P(t) (Patient households):")
display(FOC_h_P_simplified)

println("\nLinearized FOC with respect to d_P(t) (Patient households):")
display(FOC_d_P_simplified)

#%% impaticent household 
using Symbolics


# Define time variable
@variables t

# Define parameters and variables
@syms β_I a φ ε_z(t) ε_h(t)
@syms c_I(t) h_I(t) b_I(t) q_h(t) w_I(t) l_I(t) r_b(t) pi_var(t) m_I λ_I μ_I

# Define the utility function U_I
U_I = β_I^t * ((1 - a) * ε_z(t) * log(c_I(t) - a * c_I(t - 1)) +
               ε_h(t) * log(h_I(t)) - 
               (h_I(t)^(1 + φ)) / (1 + φ))

# Define the budget constraint BC_I
BC_I = w_I(t) * l_I(t) + b_I(t) - (c_I(t) + q_h(t) * h_I(t) + (1 + r_b(t)) / pi_var(t) * b_I(t - 1))

# Define the borrowing constraint Borrowing_C
Borrowing_C = (1 + r_b(t)) * b_I(t) - m_I * q_h(t+1) * h_I(t) * pi_var(t+1)

# Define the Lagrangian L_I
L_I = U_I + λ_I * BC_I + μ_I * Borrowing_C

# Compute the FOCs by differentiating the Lagrangian
FOC_c_I = expand_derivatives(Differential(c_I(t))(L_I)) ~ 0
FOC_h_I = expand_derivatives(Differential(h_I(t))(L_I)) ~ 0
FOC_b_I = expand_derivatives(Differential(b_I(t))(L_I)) ~ 0

# Solve for λ_I from FOC_c_I
λ_I_expr = symbolic_linear_solve(FOC_c_I, λ_I)

# Substitute λ_I into FOC_h_I to simplify
FOC_h_I_substituted = substitute(FOC_h_I, λ_I => λ_I_expr)

# Solve for μ_I from FOC_b_I
μ_I_expr = symbolic_linear_solve(FOC_b_I, μ_I)

# Substitute μ_I into the substituted FOC_h_I to eliminate both multipliers
FOC_h_I_no_lagrange = substitute(FOC_h_I_substituted, μ_I => μ_I_expr)

# Display Results
println("FOC with respect to c_I(t) (Impatient households):")
display(FOC_c_I)

println("\nFOC with respect to h_I(t) (Impatient households):")
display(FOC_h_I)

println("\nFOC with respect to b_I(t) (Impatient households):")
display(FOC_b_I)

println("\nSimplified FOC for h_I(t) without Lagrange multipliers:")
display(FOC_h_I_no_lagrange)

# Add interpretation for binding/non-binding constraint
println("\nWhen Borrowing Constraint is Binding (Active):")
println("The borrowing constraint is treated as an equality, and μ_I > 0.")

println("\nWhen Borrowing Constraint is Not Binding (Inactive):")
println("Set μ_I = 0 in the above equations to simplify.")

#%% impatient household linearize
using Symbolics

# Define variables and parameters
@syms β_I ε_z(t) ε_h(t) φ a q_h(t) h_I(t) c_I(t) b_I(t) r_b(t) λ_I μ_I m_I π_var(t)
@syms c_I_ss h_I_ss q_h_ss b_I_ss r_b_ss π_var_ss ε_z_ss ε_h_ss λ_I_ss μ_I_ss  # Steady-state values
@syms Δc_I Δh_I Δq_h Δb_I Δπ_var Δr_b  # Deviations from steady state

# Define FOCs (Impatient Households)
FOC_c_I = ((1 - a) * (β_I^t) * ε_z(t)) / (c_I(t) - a * c_I(t-1)) - λ_I ~ 0
FOC_h_I = (β_I^t) * (ε_h(t) / h_I(t) - (h_I(t)^φ)) - q_h(t) * λ_I - m_I * π_var(t+1) * q_h(t+1) * μ_I ~ 0
FOC_b_I = λ_I + (1 + r_b(t)) * μ_I ~ 0

# Steady-state conditions (formatted as substitutions)
steady_state_subs = Dict(
    ε_z(t) => ε_z_ss,
    ε_h(t) => ε_h_ss,
    c_I(t) => c_I_ss,
    h_I(t) => h_I_ss,
    q_h(t) => q_h_ss,
    b_I(t) => b_I_ss,
    π_var(t) => π_var_ss,
    r_b(t) => r_b_ss,
    λ_I => λ_I_ss,
    μ_I => μ_I_ss
)

# Define linearized variables (formatted as substitutions)
linearized_vars = Dict(
    c_I(t) => c_I_ss * (1 + Δc_I),          # Linearize c_I around steady state
    h_I(t) => h_I_ss * (1 + Δh_I),          # Linearize h_I around steady state
    q_h(t) => q_h_ss * (1 + Δq_h),          # Linearize q_h around steady state
    b_I(t) => b_I_ss * (1 + Δb_I),          # Linearize b_I around steady state
    r_b(t) => r_b_ss * (1 + Δr_b),          # Linearize r_b around steady state
    π_var(t) => π_var_ss * (1 + Δπ_var)     # Linearize π_var around steady state
)

# Substitute steady state and deviations into FOCs
FOC_c_I_linearized = substitute(FOC_c_I, steady_state_subs)
FOC_c_I_linearized = substitute(FOC_c_I_linearized, linearized_vars)

FOC_h_I_linearized = substitute(FOC_h_I, steady_state_subs)
FOC_h_I_linearized = substitute(FOC_h_I_linearized, linearized_vars)

FOC_b_I_linearized = substitute(FOC_b_I, steady_state_subs)
FOC_b_I_linearized = substitute(FOC_b_I_linearized, linearized_vars)

# Simplify the linearized equations
FOC_c_I_simplified = expand(FOC_c_I_linearized)
FOC_h_I_simplified = expand(FOC_h_I_linearized)
FOC_b_I_simplified = expand(FOC_b_I_linearized)

# Display Results
println("Linearized FOC with respect to c_I(t) (Impatient households):")
display(FOC_c_I_simplified)

println("\nLinearized FOC with respect to h_I(t) (Impatient households):")
display(FOC_h_I_simplified)

println("\nLinearized FOC with respect to b_I(t) (Impatient households):")
display(FOC_b_I_simplified)



#%% entrepreneurs 


using Symbolics

# Define time variable
@variables t

# Define parameters and variables
@syms β_E a φ δ ξ_1 ξ_2 m_E
@syms c_E(t) k_E(t) b_E(t) u(t) l_E(t) Y_E(t) A_E(t)
@syms q_k(t) r_be(t) pi_var(t) X_t λ_E μ_E

# Define utility function U_E
U_E = β_E^t * log(c_E(t) - a * c_E(t - 1))

# Define utilization cost F(u)
F_u = ξ_1 * (u(t) - 1) + 0.5 * ξ_2 * (u(t) - 1)^2

# Define production function Y_E
Y_E = A_E(t) * (u(t) * k_E(t-1))^φ * l_E(t)^(1 - φ)

# Define budget constraint BC_E
BC_E = Y_E / X_t + (1 + r_be(t-1)) / pi_var(t) * b_E(t-1) + q_k(t) * (1 - δ) * k_E(t-1) - (
    c_E(t) + q_k(t) * k_E(t) + b_E(t) + F_u * k_E(t-1)
)

# Define borrowing constraint Borrowing_C as equality (when binding)
Borrowing_C = (1 + r_be(t)) * b_E(t) - m_E * q_k(t+1) * (1 - δ) * k_E(t) * pi_var(t+1)

# Define the Lagrangian L_E
L_E = U_E + λ_E * BC_E + μ_E * Borrowing_C

# Compute the FOCs by differentiating the Lagrangian
FOC_c_E = expand_derivatives(Differential(c_E(t))(L_E)) ~ 0
FOC_k_E = expand_derivatives(Differential(k_E(t))(L_E)) ~ 0
FOC_u = expand_derivatives(Differential(u(t))(L_E)) ~ 0
FOC_b_E = expand_derivatives(Differential(b_E(t))(L_E)) ~ 0

# Solve for λ_E from FOC_c_E
λ_E_expr = symbolic_linear_solve(FOC_c_E, λ_E)

# Substitute λ_E into FOC_k_E
FOC_k_E_substituted = substitute(FOC_k_E, λ_E => λ_E_expr)

# Solve for μ_E from FOC_b_E
μ_E_expr = symbolic_linear_solve(FOC_b_E, μ_E)

# Substitute μ_E into FOC_k_E_substituted to eliminate both multipliers
FOC_k_E_no_lagrange = substitute(FOC_k_E_substituted, μ_E => μ_E_expr)

# Display Results
println("FOC with respect to c_E(t) (Entrepreneurs):")
display(FOC_c_E)

println("\nFOC with respect to k_E(t) (Entrepreneurs):")
display(FOC_k_E)

println("\nFOC with respect to u(t) (Entrepreneurs):")
display(FOC_u)

println("\nFOC with respect to b_E(t) (Entrepreneurs):")
display(FOC_b_E)

println("\nSimplified FOC for k_E(t) without Lagrange multipliers:")
display(FOC_k_E_no_lagrange)

# Add analysis for binding/non-binding constraint
println("\nWhen Borrowing Constraint is Binding (Active):")
println("The borrowing constraint is treated as an equality.")

println("\nWhen Borrowing Constraint is Not Binding (Inactive):")
println("Set μ_E = 0 in the above equations.")

#%% entrepreneurs linearize
using Symbolics

# Define variables and parameters
@syms β_E a φ δ ξ_1 ξ_2 q_k(t) u(t) k_E(t) l_E(t) A_E(t) X_t c_E(t) b_E(t) r_be(t) λ_E μ_E m_E π_var(t)
@syms c_E_ss k_E_ss u_ss q_k_ss b_E_ss l_E_ss r_be_ss π_var_ss λ_E_ss μ_E_ss A_E_ss X_ss  # Steady-state values
@syms Δc_E Δk_E Δu Δq_k Δb_E Δl_E Δr_be Δπ_var ΔA_E ΔX  # Deviations from steady state

# Define FOCs (Entrepreneurs)
FOC_c_E = (β_E^t) / (c_E(t) - a * c_E(t-1)) - λ_E ~ 0
FOC_k_E = -q_k(t) * λ_E + m_E * q_k(t+1) * π_var(t+1) * (-1 + δ) * μ_E ~ 0
FOC_u = ((k_E(t-1) * (l_E(t)^(1 - φ)) * A_E(t) * ((k_E(t-1) * u(t))^(-1 + φ)) * φ) / X_t -
         k_E(t-1) * (ξ_1 + (-1 + u(t)) * ξ_2)) * λ_E ~ 0
FOC_b_E = -λ_E + (1 + r_be(t)) * μ_E ~ 0

# Steady-state conditions (formatted as substitutions)
steady_state_subs = Dict(
    c_E(t) => c_E_ss,
    k_E(t) => k_E_ss,
    u(t) => u_ss,
    q_k(t) => q_k_ss,
    b_E(t) => b_E_ss,
    l_E(t) => l_E_ss,
    r_be(t) => r_be_ss,
    π_var(t) => π_var_ss,
    A_E(t) => A_E_ss,
    X_t => X_ss,
    λ_E => λ_E_ss,
    μ_E => μ_E_ss
)

# Define linearized variables (formatted as substitutions)
linearized_vars = Dict(
    c_E(t) => c_E_ss * (1 + Δc_E),
    k_E(t) => k_E_ss * (1 + Δk_E),
    u(t) => u_ss * (1 + Δu),
    q_k(t) => q_k_ss * (1 + Δq_k),
    b_E(t) => b_E_ss * (1 + Δb_E),
    l_E(t) => l_E_ss * (1 + Δl_E),
    r_be(t) => r_be_ss * (1 + Δr_be),
    π_var(t) => π_var_ss * (1 + Δπ_var),
    A_E(t) => A_E_ss * (1 + ΔA_E),
    X_t => X_ss * (1 + ΔX)
)

# Substitute steady state and deviations into FOCs
FOC_c_E_linearized = substitute(FOC_c_E, steady_state_subs)
FOC_c_E_linearized = substitute(FOC_c_E_linearized, linearized_vars)

FOC_k_E_linearized = substitute(FOC_k_E, steady_state_subs)
FOC_k_E_linearized = substitute(FOC_k_E_linearized, linearized_vars)

FOC_u_linearized = substitute(FOC_u, steady_state_subs)
FOC_u_linearized = substitute(FOC_u_linearized, linearized_vars)

FOC_b_E_linearized = substitute(FOC_b_E, steady_state_subs)
FOC_b_E_linearized = substitute(FOC_b_E_linearized, linearized_vars)

# Simplify the linearized equations
FOC_c_E_simplified = expand(FOC_c_E_linearized)
FOC_k_E_simplified = expand(FOC_k_E_linearized)
FOC_u_simplified = expand(FOC_u_linearized)
FOC_b_E_simplified = expand(FOC_b_E_linearized)

# Display Results
println("Linearized FOC with respect to c_E(t) (Entrepreneurs):")
display(FOC_c_E_simplified)

println("\nLinearized FOC with respect to k_E(t) (Entrepreneurs):")
display(FOC_k_E_simplified)

println("\nLinearized FOC with respect to u(t) (Entrepreneurs):")
display(FOC_u_simplified)

println("\nLinearized FOC with respect to b_E(t) (Entrepreneurs):")
display(FOC_b_E_simplified)




#%%  Wholesale bankbranch
using Symbolics

# Define time variable
@variables t

# Define parameters and variables
@syms β_R r k_b δ_b Ω v_b pi_var(t)
@syms B_t(t) D_t(t) K_b(t) K_b_next(t) R_b(t) R_d(t) Adj_k(t) ε_b(t)

# Define Lagrange multipliers
@syms λ_bs λ_cap

# Define adjustment cost for capital adequacy
Adj_k = (k_b / 2) * ((K_b(t) / B_t(t)) - v_b)^2 * K_b(t)

# Define the profit function for the wholesale branch
Profit = β_R^t * ((1 + R_b(t)) * B_t(t) - B_t(t+1) * pi_var(t+1) -
                  (1 + R_d(t)) * D_t(t) + D_t(t+1) * pi_var(t+1) -
                  K_b_next(t) * pi_var(t+1) + K_b(t) - Adj_k)

# Define the balance sheet constraint
Balance_Sheet_C = B_t(t) - (D_t(t) + K_b(t) + ε_b(t))

# Define the capital evolution constraint
Capital_Evolution_C = K_b_next(t) * pi_var(t+1) - ((1 - δ_b) * K_b(t) + Ω * ε_b(t))

# Define the Lagrangian L_W
L_W = Profit + λ_bs * Balance_Sheet_C + λ_cap * Capital_Evolution_C

# Compute the FOCs by differentiating the Lagrangian
FOC_B_t = expand_derivatives(Differential(B_t(t))(L_W)) ~ 0
FOC_D_t = expand_derivatives(Differential(D_t(t))(L_W)) ~ 0
FOC_K_b = expand_derivatives(Differential(K_b(t))(L_W)) ~ 0

# Solve for λ_bs from FOC_B_t
λ_bs_expr = symbolic_linear_solve(FOC_B_t, λ_bs)

# Solve for λ_cap from FOC_K_b
λ_cap_expr = symbolic_linear_solve(FOC_K_b, λ_cap)

# Substitute λ_bs into FOC_D_t
FOC_D_t_substituted_bs = substitute(FOC_D_t, λ_bs => λ_bs_expr)

# Substitute λ_cap into the result
FOC_D_t_no_lagrange = substitute(FOC_D_t_substituted_bs, λ_cap => λ_cap_expr)

# Display Results
println("FOC with respect to B_t (Wholesale Branch):")
display(FOC_B_t)

println("\nFOC with respect to D_t (Wholesale Branch):")
display(FOC_D_t)

println("\nFOC with respect to K_b (Wholesale Branch):")
display(FOC_K_b)

println("\nSimplified FOC for D_t without Lagrange multipliers:")
display(FOC_D_t_no_lagrange)

# Add interpretation for binding constraints
println("\nInterpretation:")
println("The balance sheet and capital evolution constraints are always binding.")

#%% wholesale branch linearize
using Symbolics

# Define variables and parameters
@syms β_R R_b(t) R_d(t) v_b K_b(t) B_t(t) λ_bs λ_cap δ_b k_b
@syms B_t_ss K_b_ss R_b_ss R_d_ss λ_bs_ss λ_cap_ss v_b_ss  # Steady-state values
@syms ΔB_t ΔK_b ΔR_b ΔR_d  # Deviations from steady state

# Define FOCs (Wholesale Branch)
FOC_B_t = λ_bs + (1 + R_b(t) + k_b * (-v_b + K_b(t) / B_t(t)) * K_b(t) * (K_b(t) / (B_t(t)^2))) * (β_R^t) ~ 0
FOC_D_t = -λ_bs + (-1 - R_d(t)) * (β_R^t) ~ 0
FOC_K_b = -λ_bs + (-1 + δ_b) * λ_cap + (1 + (-k_b * (-v_b + K_b(t) / B_t(t)) * K_b(t)) / B_t(t) -
             (1//2) * k_b * ((-v_b + K_b(t) / B_t(t))^2)) * (β_R^t) ~ 0

# Steady-state conditions (formatted as substitutions)
steady_state_subs = Dict(
    B_t(t) => B_t_ss,
    K_b(t) => K_b_ss,
    R_b(t) => R_b_ss,
    R_d(t) => R_d_ss,
    λ_bs => λ_bs_ss,
    λ_cap => λ_cap_ss
)

# Define linearized variables (formatted as substitutions)
linearized_vars = Dict(
    B_t(t) => B_t_ss * (1 + ΔB_t),          # Linearize B_t around steady state
    K_b(t) => K_b_ss * (1 + ΔK_b),          # Linearize K_b around steady state
    R_b(t) => R_b_ss * (1 + ΔR_b),          # Linearize R_b around steady state
    R_d(t) => R_d_ss * (1 + ΔR_d)           # Linearize R_d around steady state
)

# Substitute steady state and deviations into FOCs
FOC_B_t_linearized = substitute(FOC_B_t, steady_state_subs)
FOC_B_t_linearized = substitute(FOC_B_t_linearized, linearized_vars)

FOC_D_t_linearized = substitute(FOC_D_t, steady_state_subs)
FOC_D_t_linearized = substitute(FOC_D_t_linearized, linearized_vars)

FOC_K_b_linearized = substitute(FOC_K_b, steady_state_subs)
FOC_K_b_linearized = substitute(FOC_K_b_linearized, linearized_vars)

# Simplify the linearized equations
FOC_B_t_simplified = expand(FOC_B_t_linearized)
FOC_D_t_simplified = expand(FOC_D_t_linearized)
FOC_K_b_simplified = expand(FOC_K_b_linearized)

# Display Results
println("Linearized FOC with respect to B_t (Wholesale Branch):")
display(FOC_B_t_simplified)

println("\nLinearized FOC with respect to D_t (Wholesale Branch):")
display(FOC_D_t_simplified)

println("\nLinearized FOC with respect to K_b (Wholesale Branch):")
display(FOC_K_b_simplified)


#%% retail branch 
using Symbolics

# Define time variable
@variables t

# Define parameters and variables
@syms β_E k_r r_t b_t_n_star ε_b_n(t)
@syms r_b_n(t) b_t_n(t) r_b_e(t) b_t_e(t) R_b(t) B_t(t) Adj_r(t)

# Define demand function for loans
Loan_Demand = b_t_n(t) ~ (r_b_n(t) / r_t)^(-ε_b_n(t)) * b_t_n_star

# Define adjustment cost
Adj_r = (k_r / 2) * ((r_b_n(t) / r_t) - 1)^2 * r_b_n(t) * b_t_n(t)

# Define profit function
Profit = β_E^t * (r_b_n(t) * b_t_n(t) + r_b_e(t) * b_t_e(t) - R_b(t) * B_t(t) - Adj_r)

# Define the Lagrangian
@syms λ_loan
L_Retail = Profit + λ_loan * (b_t_n(t) - (r_b_n(t) / r_t)^(-ε_b_n(t)) * b_t_n_star)

# Compute FOCs by differentiating the Lagrangian
FOC_r_b_n = expand_derivatives(Differential(r_b_n(t))(L_Retail)) ~ 0
FOC_r_b_e = expand_derivatives(Differential(r_b_e(t))(L_Retail)) ~ 0

# Solve for loan demand λ_loan
λ_loan_expr = symbolic_linear_solve(FOC_r_b_n, λ_loan)

# Substitute λ_loan into FOC_r_b_e to eliminate the multiplier
FOC_r_b_e_no_lagrange = substitute(FOC_r_b_e, λ_loan => λ_loan_expr)

# Display Results
println("FOC with respect to r_b_n(t) (Retail Branch):")
display(FOC_r_b_n)

println("\nFOC with respect to r_b_e(t) (Retail Branch):")
display(FOC_r_b_e)

println("\nSimplified FOC for r_b_e(t) without Lagrange multipliers:")
display(FOC_r_b_e_no_lagrange)

#%% retail branch linearirize 
using Symbolics

# Define variables and parameters
@syms β_E r_b_n(t) r_b_e(t) r_t b_t_n_star b_t_n(t) b_t_e(t) ε_b_n(t) λ_loan k_r
@syms r_b_n_ss r_b_e_ss r_t_ss b_t_n_ss b_t_e_ss ε_b_n_ss λ_loan_ss k_r_ss  # Steady-state values
@syms Δr_b_n Δr_b_e Δr_t Δb_t_n Δb_t_e Δε_b_n Δλ_loan  # Deviations from steady state

# Define FOCs (Retail Branch)
FOC_r_b_n = (b_t_n_star * ((r_b_n(t) / r_t)^(-1 - ε_b_n(t))) * ε_b_n(t) * λ_loan) / r_t +
            ((-k_r * r_b_n(t) * b_t_n(t) * (-1 + r_b_n(t) / r_t)) / r_t +
             b_t_n(t) - (1//2) * k_r * b_t_n(t) * ((-1 + r_b_n(t) / r_t)^2)) * (β_E^t) ~ 0

FOC_r_b_e = (β_E^t) * b_t_e(t) ~ 0

# Steady-state conditions (formatted as substitutions)
steady_state_subs = Dict(
    r_b_n(t) => r_b_n_ss,
    r_b_e(t) => r_b_e_ss,
    r_t => r_t_ss,
    b_t_n(t) => b_t_n_ss,
    b_t_e(t) => b_t_e_ss,
    ε_b_n(t) => ε_b_n_ss,
    λ_loan => λ_loan_ss
)

# Define linearized variables (formatted as substitutions)
linearized_vars = Dict(
    r_b_n(t) => r_b_n_ss * (1 + Δr_b_n),        # Linearize r_b_n around steady state
    r_b_e(t) => r_b_e_ss * (1 + Δr_b_e),        # Linearize r_b_e around steady state
    r_t => r_t_ss * (1 + Δr_t),                 # Linearize r_t around steady state
    b_t_n(t) => b_t_n_ss * (1 + Δb_t_n),        # Linearize b_t_n around steady state
    b_t_e(t) => b_t_e_ss * (1 + Δb_t_e),        # Linearize b_t_e around steady state
    ε_b_n(t) => ε_b_n_ss * (1 + Δε_b_n),        # Linearize ε_b_n around steady state
    λ_loan => λ_loan_ss * (1 + Δλ_loan)         # Linearize λ_loan around steady state
)

# Substitute steady state and deviations into FOCs
FOC_r_b_n_linearized = substitute(FOC_r_b_n, steady_state_subs)
FOC_r_b_n_linearized = substitute(FOC_r_b_n_linearized, linearized_vars)

FOC_r_b_e_linearized = substitute(FOC_r_b_e, steady_state_subs)
FOC_r_b_e_linearized = substitute(FOC_r_b_e_linearized, linearized_vars)

# Simplify the linearized equations
FOC_r_b_n_simplified = expand(FOC_r_b_n_linearized)
FOC_r_b_e_simplified = expand(FOC_r_b_e_linearized)

# Display Results
println("Linearized FOC with respect to r_b_n(t) (Retail Branch):")
display(FOC_r_b_n_simplified)

println("\nLinearized FOC with respect to r_b_e(t) (Retail Branch):")
display(FOC_r_b_e_simplified)


#%% wage
using Symbolics


# Define time variable
@variables t

# Define parameters and variables
@syms β_L λ_L(t) k_w π_L(t) π_w(t) φ ε_L(t) l_t_star W_L(t) W_L_prev(t) P_t
@syms l_t(i, m) i m

# Define labor demand
labor_demand = l_t(i, m) ~ (W_L(t) / W_L_prev(t))^(-ε_L(t)) * l_t_star

# Define the profit function
profit_function = β_L^t * λ_L(t) * (
    (W_L(t) / P_t) * l_t(i, m) -
    (k_w / 2) * ((W_L(t) / W_L_prev(t)) - π_L(t-1) * π_w(t))^2 * (W_L(t) / P_t) -
    (l_t(i, m)^(1+φ)) / (1+φ)
)

# Define the Lagrangian
@syms λ_d
L_Labor = profit_function + λ_d * (l_t(i, m) - (W_L(t) / W_L_prev(t))^(-ε_L(t)) * l_t_star)

# Compute the FOC by differentiating the Lagrangian with respect to W_L(t)
FOC_W_L = expand_derivatives(Differential(W_L(t))(L_Labor))

# Simplify the FOC by substituting labor demand
FOC_W_L_simplified = substitute(FOC_W_L, l_t(i, m) => (W_L(t) / W_L_prev(t))^(-ε_L(t)) * l_t_star)

# Display Results
println("FOC with respect to W_L(t) (Nominal Wage):")
display(FOC_W_L)

println("\nSimplified FOC with labor demand substituted:")
display(FOC_W_L_simplified)

%## wages
using Symbolics

# Define variables and parameters
@syms β_L W_L(t) W_L_prev(t) ε_L(t) π_w(t) π_L(t) k_w l_t_star l_t(i, m) P_t λ_d λ_L(t)
@syms W_L_ss W_L_prev_ss ε_L_ss π_w_ss π_L_ss k_w_ss l_t_star_ss l_t_ss P_t_ss λ_d_ss λ_L_ss  # Steady-state values
@syms ΔW_L ΔW_L_prev Δε_L Δπ_w Δπ_L Δl_t ΔP_t Δλ_d Δλ_L  # Deviations from steady state

# Define FOC with respect to W_L(t)
FOC_W_L = (l_t_star * ((W_L(t) / W_L_prev(t))^(-1 - ε_L(t))) * ε_L(t) * λ_d) / W_L_prev(t) +
          (((-k_w * (W_L(t) / W_L_prev(t) - π_w(t) * π_L(t-1)) * W_L(t)) / W_L_prev(t) -
            (1//2) * k_w * ((W_L(t) / W_L_prev(t) - π_w(t) * π_L(t-1))^2)) / P_t +
           l_t(i, m) / P_t) * λ_L(t) * (β_L^t) ~ 0

# Steady-state conditions (formatted as substitutions)
steady_state_subs = Dict(
    W_L(t) => W_L_ss,
    W_L_prev(t) => W_L_prev_ss,
    ε_L(t) => ε_L_ss,
    π_w(t) => π_w_ss,
    π_L(t-1) => π_L_ss,
    l_t_star => l_t_star_ss,
    l_t(i, m) => l_t_ss,
    P_t => P_t_ss,
    λ_d => λ_d_ss,
    λ_L(t) => λ_L_ss
)

# Define linearized variables (formatted as substitutions)
linearized_vars = Dict(
    W_L(t) => W_L_ss * (1 + ΔW_L),
    W_L_prev(t) => W_L_prev_ss * (1 + ΔW_L_prev),
    ε_L(t) => ε_L_ss * (1 + Δε_L),
    π_w(t) => π_w_ss * (1 + Δπ_w),
    π_L(t-1) => π_L_ss * (1 + Δπ_L),
    l_t(i, m) => l_t_ss * (1 + Δl_t),
    P_t => P_t_ss * (1 + ΔP_t),
    λ_d => λ_d_ss * (1 + Δλ_d),
    λ_L(t) => λ_L_ss * (1 + Δλ_L)
)

# Substitute steady state and deviations into the FOC
FOC_W_L_linearized = substitute(FOC_W_L, steady_state_subs)
FOC_W_L_linearized = substitute(FOC_W_L_linearized, linearized_vars)

# Simplify the linearized equation
FOC_W_L_simplified = expand(FOC_W_L_linearized)

# Display Results
println("Linearized FOC with respect to W_L(t) (Nominal Wage):")
display(FOC_W_L_simplified)

    

#%% Capital 

using Symbolics

# Define time variable
@variables t

# Define parameters and variables
@syms β_E λ_E(t) k_i q_k(t) e_q(t) δ
@syms i_t i_prev k_t k_prev λ_E_next q_k_next e_q_next

# Utility function for capital producers
Utility_Capital = λ_E(t) * β_E^t * (
    q_k(t) * i_t - (k_i / 2) * ((i_t * e_q(t) / i_prev) - 1)^2 * i_t
)

# Price of capital equation
Price_Capital = q_k(t) * (1 - (k_i / 2) * ((i_t * e_q(t) / i_prev) - 1)^2) +
                k_i * (i_t * e_q(t) / i_prev) * ((i_t * e_q(t) / i_prev) - 1) +
                β_E * (λ_E_next / λ_E(t)) * (q_k_next * e_q_next / i_t) ~ 1

# Law of motion of capital
Capital_Law = k_t ~ (1 - δ) * k_prev + (1 - (k_i / 2) * ((i_t * e_q(t) / i_prev) - 1)^2) * i_t

# Define the Lagrangian for capital producers
@syms λ_C
L_Capital = Utility_Capital + λ_C * (k_t - ((1 - δ) * k_prev + (1 - (k_i / 2) * ((i_t * e_q(t) / i_prev) - 1)^2) * i_t))

# Compute first-order conditions (FOCs)
FOC_q_k = expand_derivatives(Differential(q_k(t))(L_Capital))  # FOC w.r.t. q_k(t)
FOC_i = expand_derivatives(Differential(i_t)(L_Capital))        # FOC w.r.t. i_t

# Display Results
println("Capital Producers' FOCs:")

println("\nFOC with respect to q_k(t) (Price of capital):")
display(FOC_q_k)

println("\nFOC with respect to i_t (Investment):")
display(FOC_i)

println("\nPrice of Capital Equation:")
display(Price_Capital)

println("\nLaw of Motion of Capital:")
display(Capital_Law)

#%% capital linearize it
using Symbolics

# Define variables and parameters
@syms β_E λ_E(t) λ_C q_k(t) i_t e_q(t) i_prev k_i
@syms q_k_ss i_t_ss e_q_ss i_prev_ss λ_E_ss λ_C_ss k_i_ss  # Steady-state values
@syms Δq_k Δi_t Δe_q Δi_prev Δλ_E Δλ_C  # Deviations from steady state

# Define FOCs
FOC_q_k = i_t * (β_E^t) * λ_E(t) ~ 0
FOC_i_t = (-1 + (i_t * k_i * (-1 + (i_t * e_q(t)) / i_prev) * e_q(t)) / i_prev +
           (1//2) * k_i * ((-1 + (i_t * e_q(t)) / i_prev)^2)) * λ_C +
          (q_k(t) + (-i_t * k_i * (-1 + (i_t * e_q(t)) / i_prev) * e_q(t)) / i_prev -
           (1//2) * k_i * ((-1 + (i_t * e_q(t)) / i_prev)^2)) * (β_E^t) * λ_E(t) ~ 0

# Steady-state conditions (formatted as substitutions)
steady_state_subs = Dict(
    q_k(t) => q_k_ss,
    i_t => i_t_ss,
    e_q(t) => e_q_ss,
    i_prev => i_prev_ss,
    λ_E(t) => λ_E_ss,
    λ_C => λ_C_ss
)

# Define linearized variables (formatted as substitutions)
linearized_vars = Dict(
    q_k(t) => q_k_ss * (1 + Δq_k),
    i_t => i_t_ss * (1 + Δi_t),
    e_q(t) => e_q_ss * (1 + Δe_q),
    i_prev => i_prev_ss * (1 + Δi_prev),
    λ_E(t) => λ_E_ss * (1 + Δλ_E),
    λ_C => λ_C_ss * (1 + Δλ_C)
)

# Substitute steady state and deviations into FOCs
FOC_q_k_linearized = substitute(FOC_q_k, steady_state_subs)
FOC_q_k_linearized = substitute(FOC_q_k_linearized, linearized_vars)

FOC_i_t_linearized = substitute(FOC_i_t, steady_state_subs)
FOC_i_t_linearized = substitute(FOC_i_t_linearized, linearized_vars)

# Simplify the linearized equations
FOC_q_k_simplified = expand(FOC_q_k_linearized)
FOC_i_t_simplified = expand(FOC_i_t_linearized)

# Display Results
println("Linearized FOC with respect to q_k(t) (Price of capital):")
display(FOC_q_k_simplified)

println("\nLinearized FOC with respect to i_t (Investment):")
display(FOC_i_t_simplified)

#%% Final goods market 
using Symbolics

# Define time variable
@variables t

# Define parameters and callable variables
@syms β_P λ_P(t) k_p π_t π_prev τ_p ν_p ε_y y_t y_t_j P_w
@syms P_t(j) P_prev(j)  # Declare P_t and P_prev as callable
@syms j  # Declare 'j' as a symbolic variable

# Define the demand function for final goods
Demand = y_t_j ~ (P_t(j) / P_t(t))^(-ε_y) * y_t

# Define the profit function
Profit = β_P^t * λ_P(t) * (
    P_t(j) * y_t_j - P_w * y_t_j -
    k_p * (P_t(j) / P_prev(j) - τ_p * π_prev^(1-ν_p))^2 * P_t(t) * y_t
)

# Define the Lagrangian
@syms λ_d
L_FinalGoods = Profit + λ_d * (y_t_j - (P_t(j) / P_t(t))^(-ε_y) * y_t)

# Compute the FOC by differentiating the Lagrangian with respect to P_t(j)
FOC_P_t_j = expand_derivatives(Differential(P_t(j))(L_FinalGoods))

# Simplify the FOC by substituting the demand function
FOC_P_t_j_simplified = substitute(FOC_P_t_j, y_t_j => (P_t(j) / P_t(t))^(-ε_y) * y_t)

# Display Results
println("FOC with respect to P_t(j) (Price of final good):")
display(FOC_P_t_j)

println("\nSimplified FOC with demand substituted:")
display(FOC_P_t_j_simplified)

#%%final goods 
using Symbolics

# Define variables and parameters
@syms β_P P_t(j) P_t(t) P_prev(j) y_t y_t_j ε_y λ_d λ_P(t) k_p π_prev τ_p ν_p
@syms P_t_ss P_prev_ss y_t_ss y_t_j_ss ε_y_ss λ_d_ss λ_P_ss π_prev_ss τ_p_ss k_p_ss ν_p_ss  # Steady-state values
@syms ΔP_t ΔP_prev Δy_t Δy_t_j Δλ_d Δλ_P Δπ_prev  # Deviations from steady state

# Define FOC with respect to P_t(j)
FOC_P_t_j = (y_t * ((P_t(j) / P_t(t))^(-1 - ε_y)) * ε_y * λ_d) / P_t(t) +
            (y_t_j + (-2 * k_p * y_t * P_t(t) * (P_t(j) / P_prev(j) - (π_prev^(1 - ν_p)) * τ_p)) / P_prev(j)) *
            (β_P^t) * λ_P(t) ~ 0

# Steady-state conditions (formatted as substitutions)
steady_state_subs = Dict(
    P_t(j) => P_t_ss,
    P_t(t) => P_t_ss,
    P_prev(j) => P_prev_ss,
    y_t => y_t_ss,
    y_t_j => y_t_j_ss,
    ε_y => ε_y_ss,
    λ_d => λ_d_ss,
    λ_P(t) => λ_P_ss,
    π_prev => π_prev_ss,
    τ_p => τ_p_ss,
    k_p => k_p_ss,
    ν_p => ν_p_ss
)

# Define linearized variables (formatted as substitutions)
linearized_vars = Dict(
    P_t(j) => P_t_ss * (1 + ΔP_t),
    P_t(t) => P_t_ss * (1 + ΔP_t),
    P_prev(j) => P_prev_ss * (1 + ΔP_prev),
    y_t => y_t_ss * (1 + Δy_t),
    y_t_j => y_t_j_ss * (1 + Δy_t_j),
    λ_d => λ_d_ss * (1 + Δλ_d),
    λ_P(t) => λ_P_ss * (1 + Δλ_P),
    π_prev => π_prev_ss * (1 + Δπ_prev)
)

# Substitute steady state and deviations into the FOC
FOC_P_t_j_linearized = substitute(FOC_P_t_j, steady_state_subs)
FOC_P_t_j_linearized = substitute(FOC_P_t_j_linearized, linearized_vars)

# Simplify the linearized equation
FOC_P_t_j_simplified = expand(FOC_P_t_j_linearized)

# Display Results
println("Linearized FOC with respect to P_t(j) (Price of final good):")
display(FOC_P_t_j_simplified)


#%%  dsge
using DSGE, Wavelets, CSV, DataFrames, Statistics, LinearAlgebra

# Assuming 'target' is your DataFrame already loaded in Julia
@assert exists(:target) "DataFrame 'target' not found."

# Data Preprocessing with MODWT Wavelets
data_columns = [:r, :rh, :D, :C, :I, :ht, :W, :L, :qh, :π, :rc]
data_matrix = Matrix(target[:, data_columns])

# Apply MODWT with levels corresponding to 16-32 quarters (4-8 years)
wavelet_filter = wavelet(Daubechies, 2)
levels = 4 # Adjust levels as needed

# Apply MODWT to each column
modwt_data = [modwt(data_matrix[:, i], wavelet_filter, levels) for i in 1:size(data_matrix, 2)]

# Extract approximation coefficients at the desired level
approximation = hcat([wt.approx[levels] for wt in modwt_data]...)

# Prepare Data for Estimation
observable_variables = approximation

# Model Specification
model = DSGE.Model()

# Calibrated Parameters
set_parameter!(model, :βP, 0.99250)   # Patient households discount factor
set_parameter!(model, :βI, 0.97500)   # Impatient households discount factor
set_parameter!(model, :βE, 0.97500)   # Entrepreneurs discount factor
set_parameter!(model, :α, 0.33000)    # Capital share
set_parameter!(model, :δ, 0.02500)    # Depreciation rate of physical capital
set_parameter!(model, :ϕ, 1.00000)    # Inverse of Frisch elasticity
set_parameter!(model, :mE, 0.50000)   # Loan-to-value for entrepreneurs
set_parameter!(model, :mI, 0.75000)   # Loan-to-value for impatient households
set_parameter!(model, :π̄, 1.00000)   # Net steady state inflation
set_parameter!(model, :v_b, 0.11000)  # Basel II capital requirement
set_parameter!(model, :Ω, 1.00000)    # Profits invested in new bank capital
set_parameter!(model, :ε_be, 1.85000) # Markup over entrepreneurial loans
set_parameter!(model, :ε_bh, 1.85000) # Markup over impatient households loans
set_parameter!(model, :ε_y, 11.0000)  # Markup in the goods market
set_parameter!(model, :ε_l, 6.00000)  # Markup in the labor market
set_parameter!(model, :δ_b, 0.05070)  # Degree of utilization of bank capital
set_parameter!(model, :ξ1, 0.05460)   # Coefficient for capital utilization
set_parameter!(model, :ξ2, 0.00546)   # Coefficient for capital utilization
set_parameter!(model, :j̄, 0.07500)   # Steady state value of housing consumption preference

# Estimated Parameters (Posterior Modes)
set_parameter!(model, :k_p, 135.5)    # Price stickiness
set_parameter!(model, :k_bh, 24.65)   # H. loans adjustment cost
set_parameter!(model, :k_be, 30.06)   # E. loans adjustment cost
set_parameter!(model, :k_i, 12.46)    # Investment adjustment cost
set_parameter!(model, :k_w, 101.4)    # Wage stickiness
set_parameter!(model, :k_kb, 1.919)   # Capital requirement adjustment cost
set_parameter!(model, :a, 0.762)      # Habit formation
set_parameter!(model, :ι_w, 0.476)    # Wage indexation
set_parameter!(model, :ι_p, 0.581)    # Price indexation
set_parameter!(model, :ϕ_R, 0.808)    # Policy rate smoothing
set_parameter!(model, :ϕ_π, 1.691)    # Inflation target
set_parameter!(model, :ϕ_y, 0.402)    # Output gap

# Additional Parameters (Need to initialize)
set_parameter!(model, :κ, 0.5)        # Investment adjustment cost parameter
set_parameter!(model, :ρ_A, 0.9)      # Persistence of technology shock
set_parameter!(model, :ρ_εz, 0.9)     # Persistence of consumption preference shock
set_parameter!(model, :ρ_εh, 0.9)     # Persistence of housing preference shock
set_parameter!(model, :ρ_R, 0.9)      # Persistence of monetary policy shock
set_parameter!(model, :λ_w, 0.5)      # Wage adjustment cost parameter
set_parameter!(model, :κ_p, 0.1)      # Price adjustment parameter

# Declare Variables
declare_variables!(model, [
    # Household variables
    :cH_t, :cH_t_1, :hH_t, :lH_t, :λ_t, :ε_z_t, :ε_h_t,
    # Patient household variables
    :cP_t, :hP_t, :dP_t, :w_t, :lP_t, :r_t_1, :π_t, :dP_t_1,
    # Impatient household variables
    :cI_t, :hI_t, :bI_t, :mI_t, :q_h_t1, :π_t1, :r_bh_t,
    # Entrepreneur variables
    :cE_t, :lE_P_t, :lE_I_t, :r_be_t_1, :bE_t_1,
    :q_k_t, :kE_t, :F_u_t, :kE_t_1, :yE_t, :x_t, :bE_t,
    :u_t, :A_t_E, :lE_t, :i_t, :κ, :k_t, :k_t_1,
    # Bank variables
    :B_t, :D_t, :K_b_t, :Adj_kb_t, :λP_t, :R_t, :R_d_t,
    # Policy variables
    :r_t, :r̄, :π̄, :y_t, :y_t_1, :ε_R_t,
    # Exogenous shocks
    :A_t, :A_t_1, :ε_A_t, :ρ_A,
    :ε_z_t, :ε_z_t_1, :ε_εz_t, :ρ_εz,
    :ε_h_t, :ε_h_t_1, :ε_εh_t, :ρ_εh,
    # Other variables
    :mct, :λ_w, :π_t1, :ι_w, :ι_p,
    # Total variables
    :c_t, :h_t, :G_t, :Adj_t_j, :L_t,
    # Marginal utilities
    :λI_t,
    # Include any other variables needed
])

# Equations
@equations model begin
# Define the model equations
@equations model begin
    #################################
    # 1. Labor Market Equations
    #################################
    
    # Wage Phillips Curve (Equation 3)
    kw * (π_ws_t - π_{w_{t-1}} * π_t^(1 - ω)) * π_ws_t ==
        βP * E_t[ (λ_{t+1}^s / λ_t^s) * kw * (π_{ws_{t+1}} - π_w_t * π_{t+1}^{1 - ω}) * π_{ws_{t+1}} * (l_{t+1}^s / l_t^s) ] +
        (1 - ε_t^l) * l_t^s + (ε_t^l * l_t^{1 + φ}) / (λ_t^s * w_t^s)
    
    # Definition of wage inflation (Equation 4)
    π_ws_t == (w_t^s / w_{t-1}^s) * π_t
    
    # Labor Demand (Equation 2)
    l_t^s(i, m) == ( (W_t^s(m) / W_{t-1}^s(m)) )^( - ε_t^l ) * l_t^s

    #################################
    # 2. Capital Producers Equations
    #################################
    
    # Price of capital equation (Equation 6)
    q_k_t * [ 1 - (ki / 2) * ( (i_t * e_k_t) / (i_{t-1} * e_{k_{t-1}}) - 1 )^2 ] +
    ki * ( (i_t * e_k_t) / (i_{t-1} * e_{k_{t-1}}) - 1 ) * ( (i_t * e_k_t) / (i_{t-1} * e_{k_{t-1}}) ) +
    βE * E_t[ (λ_{t+1}^E / λ_t^E) * ki * ( (i_{t+1} * e_{k_{t+1}}) / (i_t * e_k_t) - 1 ) *
              ( (i_{t+1} * e_{k_{t+1}}) / (i_t * e_k_t) ) * q_{k_{t+1}} ] == 1
    
    # Law of motion of capital (Equation 7)
    k_t == (1 - δ) * k_{t-1} + [ 1 - (ki / 2) * ( (i_t * e_k_t) / (i_{t-1} * e_{k_{t-1}}) - 1 )^2 ] * i_t
    
    #################################
    # 3. Final Goods Producers Equations
    #################################
    
    # New Keynesian Phillips Curve for price inflation (Equation 10)
    (1 - ε_t^p) + ε_t^p * (P_t^w / P_t) -
    kp * (π_t - π_{t-1}^p * π_t^{(1 - υ)}) * π_t +
    βP * E_t[ (λ_{t+1}^P / λ_t^P) * kp * (π_{t+1} - π_t^p * π_{t+1}^{(1 - υ)}) * π_{t+1}^2 * (y_{t+1} / y_t) ] == 0
    
    # Demand Function for Final Goods (Equation 9)
    y_t(j) == ( P_t(j) / P_t )^( - ε_t^p ) * y_t

    #################################
    # 4. Households Equations
    #################################
    
    # Patient Households' Budget Constraint (Equation 2)
    C_t^P + q_h_t * (H_t^P - H_{t-1}^P) + d_t^P ==
        W_t^P * L_t^P + (1 + r_{t-1}) * (d_{t-1}^P / π_t)
    
    # Impatient Households' Budget Constraint (Equation 3)
    C_t^I + q_h_t * (H_t^I - H_{t-1}^I) + (1 + r_{t-1}^{bH}) * (b_{t-1}^I / π_t) ==
        W_t^I * L_t^I + b_t^I
    
    # Borrowing Constraint for Impatient Households (Equation 4)
    (1 + r_t^{bH}) * b_t^I <= m_t^I * E_t[ q_h_{t+1} * H_{t+1}^I * π_{t+1} ]
    
    # Marginal Utility of Consumption for Patient Households
    λ_t^P == (1 - α) * ε_t^c / (C_t^P - α * C_{t-1}^P)
    
    # Marginal Utility of Consumption for Impatient Households
    λ_t^I == (1 - α) * ε_t^c / (C_t^I - α * C_{t-1}^I)
    
    #################################
    # 5. Entrepreneurs Equations
    #################################
    
    # Entrepreneurs' Budget Constraint (Equation 6)
    C_t^E + W_t^P * L_t^{P,E} + W_t^I * L_t^{I,E} + (1 + r_{t-1}^E) * (b_{t-1}^E / π_t) +
    q_k_t * k_t^E + F(u_t) * k_{t-1}^E ==
        (y_t^E / χ_t) + b_t^E + q_k_t * (1 - δ) * k_{t-1}^E
    
    # Entrepreneurs' Production Function (Equation 7)
    y_t^E == A_t^E * (k_{t-1}^E * u_t)^θ * (L_t^E)^(1 - α)
    
    # Borrowing Constraint for Entrepreneurs (Equation 8)
    (1 + r_t^E) * b_t^E <= m_t^E * E_t[ q_{k_{t+1}} * (1 - δ) * k_t^E * π_{t+1} ]
    
    # Marginal Utility of Consumption for Entrepreneurs
    λ_t^E == (1 - α) / (C_t^E - α * C_{t-1}^E)
    
    #################################
    # 6. Wholesale Branch Equations
    #################################
    
    # Capital Requirement Adjustment Cost (Equation 10)
    Adj_kb_t == (k_bb / 2) * ( (K_b_t / B_t) - ν_b )^2
    
    # Balance Sheet Constraint (Equation 11)
    B_t == D_t + K_b_t + ε_t^{kb}
    
    # Law of Motion for Bank Capital (Equation 12)
    K_b_t * π_t == (1 - δ_b) * K_{b_{t-1}} + Ω * D_{t-1}
    
    # Wholesale Interest Rate (Equation 13)
    R_t^B == r_t - k_bb * ( (K_b_t / B_t) - ν_b ) * (K_b_t / B_t)^2
    
    #################################
    # 7. Retail Branch Equations
    #################################
    
    # Demand Function for Loans (Equation 15)
    b_t^n == ( R_t^{bn}(j) / R_t^{bn} )^( - ε_t^n ) * b_t^n_star  # Assuming b_t^n_star is the benchmark loan demand
    
    # Adjustment Cost on Commercial Interest Rates (Equation 16)
    Adj_kn_t == (k_bn / 2) * ( (R_t^{bn} / R_{t-1}^{bn} - 1) )^2 * b_t^n
    
    # Interest Rate Setting Equation (Equation 17)
  # Interest Rate Setting Equation (Equation 17)
 1 - ε_t^n + ε_t^n * (R_t_bn / R_t^B) -
 k_bn * ( (R_t_bn / R_{t-1}_bn - 1) ) * (R_t_bn / R_{t-1}_bn ) +
 β * E_t[ (λI_{t+1} / λI_t) * k_bn * ( (R_{t+1}_bn / R_t_bn - 1) ) * (R_{t+1}_bn / R_t_bn )^2 * (b_{t+1}^n / b_t^n) ] == 0

    
    #################################
    # 8. Monetary Policy Equation
    #################################
    
    # Taylor Rule (Equation 18)
    (1 + r_t) == (1 + r̄)^(1 - φ_r) * (1 + r_{t-1})^(φ_r) *
                 (π_t / π̄)^(φ_π * (1 - φ_r)) *
                 (y_t / y_{t-1})^(φ_y * (1 - φ_r)) *
                 (1 + ε_t^R)
    
    #################################
    # 9. Resource Constraint and Aggregation
    #################################
    
    # Resource Constraint of the Economy (Equation 19)
    y_t == c_t +
           q_k_t * [ k_t - (1 - δ) * k_{t-1} ] +
           k_{t-1} * [ ξ1 * (u_t - 1) + (ξ2 / 2) * (u_t - 1)^2 ] +
           (δ_b * K_{b_{t-1}}) / π_t + G_t +
           Adj_kb_t + Adj_kn_t + Adj_p_t + Adj_w_t  # Include all adjustment costs
    
    # Aggregate Consumption (Equation 20)
    c_t == c_t^P + c_t^I + c_t^E + ε_t^c
    
    # Aggregate Housing (Equation 21)
    h == h_t^P + h_t^I   # Assuming total housing supply h is normalized to 1
    
    # Capital Utilization Cost Function
    F(u_t) == ξ1 * (u_t - 1) + (ξ2 / 2) * (u_t - 1)^2
    
    #################################
    # 10. Additional Equations and Definitions
    #################################
    
    # Marginal Cost
    mct == P_t^w / P_t
    
    # Price Level Evolution
    π_t == P_t / P_{t-1}
    
    # Aggregate Price Index (Definition)
    # Assuming Calvo pricing is not used, so P_t evolves according to the price-setting behavior
    # If needed, include an equation for P_t based on the model's specifics
    
    # Wage Index (Definition)
    # Similar to the price index, define W_t if required
    
    # Definitions of Real Wages
    w_t^P == W_t^P / P_t
    w_t^I == W_t^I / P_t
    w_t^s == W_t^s / P_t
    
    # Labor Market Clearing
    L_t^P + L_t^I + L_t^E == L_t^s
    
    # Capital Market Clearing
    k_t^E == k_t
    
    # Loan Market Clearing
    b_t^I + b_t^E == B_t
    
    # Housing Market Clearing
    H_t^P + H_t^I == h   # Total housing supply
    
    # Aggregate Output
    y_t == y_t^E   # Assuming entrepreneurs produce the aggregate output
    
    # Real Estate Price Dynamics
    # Include equations for q_h_t if applicable
    
    # Expectations Operator
    # Ensure that E_t[...] is properly defined in your modeling framework
    
    # Exogenous Processes
    # Preference Shock
    ln(ε_t^c) == ρ_c * ln(ε_{t-1}^c) + ε_{c_t}
    
    # Housing Demand Shock
    ln(ε_t^h) == ρ_h * ln(ε_{t-1}^h) + ε_{h_t}
    
    # Labor Supply Shock
    ln(ε_t^l) == ρ_l * ln(ε_{t-1}^l) + ε_{l_t}
    
    # Price Markup Shock
    ln(ε_t^p) == ρ_p * ln(ε_{t-1}^p) + ε_{p_t}
    
    # Monetary Policy Shock
    ε_t^R == ρ_R * ε_{t-1}^R + ε_{R_t}
    
    # Total Factor Productivity Shock
    ln(A_t^E) == ρ_A * ln(A_{t-1}^E) + ε_{A_t}
    
    # Loan-to-Value Ratios
    ln(m_t^I) == (1 - ρ_mI) * ln(m̄_I) + ρ_mI * ln(m_{t-1}^I) + ε_{mI_t}
    ln(m_t^E) == (1 - ρ_mE) * ln(m̄_E) + ρ_mE * ln(m_{t-1}^E) + ε_{mE_t}
    
    # Interest Rate Spread Shock
    ln(ε_t^{kb}) == ρ_kb * ln(ε_{t-1}^{kb}) + ε_{kb_t}
    
    # ... Include any other exogenous processes as needed
    
end

    
   
end

# Prior Distributions for Estimated Parameters
set_prior!(model, :k_p, Gamma(30, 20))
set_prior!(model, :k_bh, Gamma(6, 2.5))
set_prior!(model, :k_be, Gamma(3, 2.5))
set_prior!(model, :k_i, Gamma(4.5, 2))
set_prior!(model, :k_w, Gamma(80, 20))
set_prior!(model, :k_kb, Gamma(10, 5))
set_prior!(model, :a, Beta(0.5, 0.05))
set_prior!(model, :ι_w, Beta(0.5, 0.05))
set_prior!(model, :ι_p, Beta(0.5, 0.05))
set_prior!(model, :ϕ_R, Beta(0.75, 0.05))
set_prior!(model, :ϕ_π, Gamma(2.0, 0.1))
set_prior!(model, :ϕ_y, Gamma(0.25, 0.05))
# Set priors for new parameters if needed
set_prior!(model, :κ, Gamma(4, 1))
set_prior!(model, :ρ_A, Beta(0.5, 0.2))
set_prior!(model, :ρ_εz, Beta(0.5, 0.2))
set_prior!(model, :ρ_εh, Beta(0.5, 0.2))
set_prior!(model, :ρ_R, Beta(0.5, 0.2))
set_prior!(model, :λ_w, Gamma(2, 0.5))
set_prior!(model, :κ_p, Gamma(0.1, 0.05))

# Initial Values for Parameters (Posterior Modes)
# (Already set above)

# Estimation Settings
mh_settings = MetropolisHastingsSettings(
    num_draws = 500000,
    burn_in = 250000,
    num_chains = 2,
    target_acceptance_rate = 0.30  # Between 25% and 35%
)

# Compile the Model
compile_model!(model)

# Run Metropolis-Hastings Estimation
results = estimate(model, observable_variables; mh_settings)

# Access Estimation Results
chain1 = results.chains[1][250001:end, :]
chain2 = results.chains[2][250001:end, :]

# Combine Chains
combined_chain = vcat(chain1, chain2)

# Calculate Acceptance Rates
acceptance_rate_chain1 = results.acceptance_rates[1]
acceptance_rate_chain2 = results.acceptance_rates[2]

println("Acceptance Rate Chain 1: ", acceptance_rate_chain1)
println("Acceptance Rate Chain 2: ", acceptance_rate_chain2)

# End of Code






