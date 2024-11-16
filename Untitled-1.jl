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

#%%
using DSGE, Wavelets, CSV, DataFrames, Statistics, LinearAlgebra

# Assuming 'target' is your DataFrame already loaded in Julia
@assert exists(:target) "DataFrame 'target' not found."

# Data Preprocessing with MODWT Wavelets
data_columns = [:r, :rh, :D, :C, :I, :ht, :W, :L, :qh, :π, :rc]
data_matrix = Matrix(target[:, data_columns])

# Apply MODWT with levels corresponding to 16-32 quarters (4-8 years)
wavelet_filter = wavelet(Daubechies, 2)
levels = 5  # Adjust levels as needed

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
    # 1. Household Utility Function - First-Order Conditions
    # Patient Household Euler Equation with Habit Formation
    λ_t == (1 - a) * ε_z_t / (cP_t - a * cP_t_1) - a * βP * E[ (1 - a) * ε_z_{t+1} / (cP_{t+1} - a * cP_t) ]
    
    # Patient Household Labor Supply Condition
    (lP_t)^ϕ == λ_t * w_t
    
    # Patient Household Housing Demand Condition
    ε_h_t / hP_t == q_h_t * λ_t - βP * E[ q_h_{t+1} * λ_{t+1} ]
    
    # Impatient Household Euler Equation (similar structure)
    λI_t == (1 - a) * ε_z_t / (cI_t - a * cI_t_1) - a * βI * E[ (1 - a) * ε_z_{t+1} / (cI_{t+1} - a * cI_t) ]
    
    # Impatient Household Labor Supply Condition
    (lI_t)^ϕ == λI_t * w_t
    
    # Impatient Household Housing Demand Condition
    ε_h_t / hI_t == q_h_t * λI_t - βI * E[ q_h_{t+1} * λI_{t+1} ]
    
    # 2. Patient Households' Budget Constraint
    cP_t + q_h_t * hP_t + dP_t == w_t * lP_t + ((1 + r_t_1) / π_t) * dP_t_1
    
    # 3. Impatient Households' Budget Constraint
    cI_t + q_h_t * hI_t + (1 + r_bh_t) * bI_t == w_t * lI_t + bI_{t-1} / π_t
    
    # 4. Impatient Households' Borrowing Constraint
    (1 + r_bh_t) * bI_t <= mI_t * E[ q_h_t1 * hI_t * π_t1 ]
    
    # 5. Entrepreneur's Budget Constraint
    cE_t + w_t * lE_P_t + w_t * lE_I_t + (1 + r_be_t_1) * bE_t_1 + q_k_t * kE_t + F_u_t * kE_t_1 ==
        yE_t / x_t + bE_t + q_k_t * (1 - δ) * kE_t_1
    
    # Capital Utilization Cost Function
    F_u_t == ξ1 * (u_t - 1) + (ξ2 / 2) * (u_t - 1)^2
    
    # 6. Cobb-Douglas Production Function for Entrepreneurs
    yE_t == A_t_E * (kE_t_1 * u_t)^α * lE_t^(1 - α)
    
    # 7. Bank's Optimization - First-Order Conditions
    # Loan Rate Setting (Simplified)
    (1 + R_t) == (1 + r_t) * (1 + ε_R_t)
    
    # Capital Requirement Adjustment
    K_b_t == v_b * B_t
    
    # Adjustment Cost
    Adj_kb_t == (k_kb / 2) * ( (K_b_t / K_b_{t-1}) - 1 )^2
    
    # 8. Capital Accumulation with Adjustment Costs
    k_t == (1 - δ) * k_t_1 + i_t * (1 - (κ / 2) * ( (i_t / i_{t-1}) - 1 )^2 )
    
    # 9. Price Adjustment Equation (New Keynesian Phillips Curve)
    π_t == βP * E[ π_{t+1} ] + κ_p * (mct - 1)
    
    # 10. Wage Adjustment Equation
    w_t == (w_{t-1})^(ι_w) * (π_{t-1})^(ι_w) * (1 + λ_w * (lP_t)^ϕ)
    
    # 11. Market Clearing Conditions
    # Total Consumption
    c_t == cP_t + cI_t + cE_t
    
    # Total Housing
    h_t == hP_t + hI_t
    
    # Goods Market Clearing
    y_t == c_t + i_t + G_t + Adj_kb_t + F_u_t * kE_t_1 + δ_b * K_b_{t-1} * π_t + sum(Adj_t_j)
    
    # Labor Market Clearing
    lP_t + lI_t + lE_P_t + lE_I_t == L_t
    
    # Capital Market Clearing
    k_t == kE_t
    
    # 12. Taylor Rule for Monetary Policy
    (1 + r_t) == (1 + r̄)^(1 - ϕ_R) * (1 + r_t_1)^(ϕ_R) * (π_t / π̄)^(ϕ_π * (1 - ϕ_R)) *
                 (y_t / y_t_1)^(ϕ_y * (1 - ϕ_R)) * (1 + ε_R_t)
    
    # 13. Laws of Motion for Exogenous Shocks
    ln(A_t_E) == ρ_A * ln(A_t_1) + ε_A_t
    ln(ε_z_t) == ρ_εz * ln(ε_z_t_1) + ε_εz_t
    ln(ε_h_t) == ρ_εh * ln(ε_h_t_1) + ε_εh_t
    ε_R_t == ρ_R * ε_R_{t-1} + ε_εR_t
    
    # 14. Marginal Cost Definition
    mct == (w_t / A_t_E) * ( (1 - α) / α ) * ( kE_t_1 * u_t / lE_t )
    
    # 15. Resource Constraint (Redundant with Goods Market Clearing)
     c_t + i_t + G_t + F_u_t * kE_t_1 + Adj_kb_t == y_t
    
    # Include any additional equations needed for the model
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






