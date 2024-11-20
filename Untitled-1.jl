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


#%% unpaticent household 
using Symbolics

using Symbolics

# Define time variable
@variables t

# Define parameters and functions
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
L_P = U_P + λ_P * BC_P

# Compute the FOCs by differentiating the Lagrangian with respect to c_P(t), h_P(t), and d_P(t)
FOC_c_P = expand_derivatives(Differential(c_P(t))(L_P))
FOC_h_P = expand_derivatives(Differential(h_P(t))(L_P))
FOC_d_P = expand_derivatives(Differential(d_P(t))(L_P))

# Display the FOCs
println("FOC with respect to c_P(t) (Patient households):")
display(FOC_c_P)

println("\nFOC with respect to h_P(t) (Patient households):")
display(FOC_h_P)

println("\nFOC with respect to d_P(t) (Patient households):")
display(FOC_d_P)

# Define callable variables
@variables t
@syms β_I a φ ε_z(t) ε_h(t)
@syms c_I(t) h_I(t) b_I(t) q_h(t) w_I(t) l_I(t) r(t) pi_var(t) m_I(t) λ_I μ_I

# Impatient Household's Utility Function
U_I = β_I^t * ((1 - a) * ε_z(t) * log(c_I(t) - a * c_I(t-1)) +
               ε_h(t) * log(h_I(t)) - 
               (h_I(t)^(1 + φ)) / (1 + φ))

# Impatient Household's Budget Constraint
BC_I = w_I(t) * l_I(t) + b_I(t) - (c_I(t) + q_h(t) * h_I(t) + (1 + r(t-1)) / pi_var(t) * b_I(t-1))

# Impatient Household's Borrowing Constraint
Borrowing_C = (1 + r(t)) * b_I(t) - m_I(t) * q_h(t+1) * h_I(t) * pi_var(t+1)

# Define the Lagrangian for Impatient Household
L_I = U_I + λ_I * BC_I + μ_I * Borrowing_C

# Compute FOCs for Impatient Household
FOC_c_I = expand_derivatives(Differential(c_I(t))(L_I))
FOC_h_I = expand_derivatives(Differential(h_I(t))(L_I))
FOC_b_I = expand_derivatives(Differential(b_I(t))(L_I))

# Display the FOCs
println("Impatient Household's FOCs:")

println("\nFOC with respect to c_I(t) (Impatient households):")
display(FOC_c_I)

println("\nFOC with respect to h_I(t) (Impatient households):")
display(FOC_h_I)

println("\nFOC with respect to b_I(t) (Impatient households):")
display(FOC_b_I)

#%% entrepreneurs 

using Symbolics

# Define callable variables
@variables t
@syms β_E a φ ε_u(t) ε_ξ ξ_1 ξ_2 δ m_E(t)
@syms c_E(t) b_E(t) q_k(t) k_E(t) u(t) r_be(t) pi_var(t) Y_E(t) A_E(t) L_P(t) L_I(t) w_P(t) w_I(t)
@syms λ_E μ_E

# Utility function for Entrepreneurs
U_E = β_E^t * ((1 - a) * log(c_E(t) - a * c_E(t-1)))

# Budget constraint
BC_E = w_P(t) * L_P(t) + w_I(t) * L_I(t) + (1 + r_be(t-1)) / pi_var(t) * b_E(t-1) +
       q_k(t) * u(t) * k_E(t-1) + ε_ξ * (u(t) - 1) + (ξ_1 * (u(t) - 1) + ξ_2 * (u(t) - 1)^2) * k_E(t-1) -
       (c_E(t) + b_E(t) + q_k(t) * (1 - δ) * k_E(t-1))

# Production function
Prod_F = A_E(t) * (k_E(t-1) * u(t))^φ * L_I(t)^(1-φ)

# Borrowing constraint
Borrowing_C = (1 + r_be(t)) * b_E(t) - m_E(t) * q_k(t+1) * k_E(t) * pi_var(t+1)

# Define Lagrangian for Entrepreneurs
 L_E = U_E + λ_E * BC_E + μ_E * Borrowing_C

# Compute FOCs
FOC_c_E = expand_derivatives(Differential(c_E(t))(L_E))
FOC_k_E = expand_derivatives(Differential(k_E(t))(L_E))
FOC_u = expand_derivatives(Differential(u(t))(L_E))
FOC_b_E = expand_derivatives(Differential(b_E(t))(L_E))

# Display FOCs
println("\nFOC with respect to c_E(t) (Entrepreneurs' consumption):")
display(FOC_c_E)

println("\nFOC with respect to k_E(t) (Entrepreneurs' capital):")
display(FOC_k_E)

println("\nFOC with respect to u(t) (Utilization rate):")
display(FOC_u)

println("\nFOC with respect to b_E(t) (Entrepreneurs' borrowing):")
display(FOC_b_E)


#%% bankbranch
using Symbolics


# Define callable variables
@variables t
@syms β_R R_b(t) R_d(t) pi_var(t) b_T(t) D_T(t) K_b(t) Adj_kb k_b v_b δ_b ξ_kb Ω
@syms r(t) λ_bs λ_lom λ_basel

# Wholesale Branch Profit Function
Profit = β_R^t * (
    (1 + R_b(t)) * b_T(t) - b_T(t+1) * pi_var(t+1) -
    (1 + R_d(t)) * D_T(t) + D_T(t+1) * pi_var(t+1) -
    (K_b(t+1) * pi_var(t+1) - K_b(t)) -
    Adj_kb / 2 * ((k_b / b_T(t)) * K_b(t) - v_b)^2 * K_b(t)
)

# Adjustment cost term
Adj_kb_term = Adj_kb / 2 * ((k_b / b_T(t)) * K_b(t) - v_b)^2 * K_b(t)

# Capital/assets ratio constraint (Basel II)
Basel_II_C = Adj_kb_term

# Balance sheet constraint
Balance_Sheet_C = b_T(t) ~ D_T(t) + K_b(t) + ξ_kb

# Law of motion for bank capital
Capital_LoM = K_b(t+1) * pi_var(t+1) ~ (1 - δ_b) * K_b(t) + Ω * ξ_kb

# Wholesale Interest Rate Constraint
Interest_Rate_C = R_b(t) ~ r(t) - k_b * (k_b / b_T(t)) * (K_b(t) / b_T(t))^2

# Define the Lagrangian
L_W = Profit +
      λ_bs * (b_T(t) - (D_T(t) + K_b(t) + ξ_kb)) +
      λ_lom * ((1 - δ_b) * K_b(t) + Ω * ξ_kb - K_b(t+1) * pi_var(t+1)) +
      λ_basel * Basel_II_C

# Compute FOCs
FOC_b_T = expand_derivatives(Differential(b_T(t))(L_W))
FOC_D_T = expand_derivatives(Differential(D_T(t))(L_W))
FOC_K_b = expand_derivatives(Differential(K_b(t))(L_W))

# Display the FOCs
println("Wholesale Branch's FOCs:")

println("\nFOC with respect to b_T(t) (Total loans):")
display(FOC_b_T)

println("\nFOC with respect to D_T(t) (Deposits):")
display(FOC_D_T)

println("\nFOC with respect to K_b(t) (Bank capital):")
display(FOC_K_b)


# Define callable variables
@variables t
@syms β_R R_b(t) R_d(t) π(t) b_T(t) D_T(t) K_b(t) Adj_kb k_b v_b δ_b ξ_kb Ω
@syms r(t)

# Wholesale Branch Profit Function
Profit = β_R^t * ((1 + R_b(t)) * b_T(t) - b_T(t+1) * π(t+1) -
                  (1 + R_d(t)) * D_T(t) + D_T(t+1) * π(t+1) -
                  (K_b(t+1) * π(t+1) - K_b(t)) -
                  Adj_kb / 2 * ((k_b / b_T(t)) * K_b(t) - v_b)^2 * K_b(t))

# Adjustment cost term
Adj_kb_term = Adj_kb / 2 * ((k_b / b_T(t)) * K_b(t) - v_b)^2 * K_b(t)

# Capital/assets ratio constraint (Basel II)
Basel_II_C = Adj_kb_term

# Balance sheet constraint
Balance_Sheet_C = b_T(t) ~ D_T(t) + K_b(t) + ξ_kb

# Law of motion for bank capital
Capital_LoM = K_b(t+1) * π(t+1) ~ (1 - δ_b) * K_b(t) + Ω * ξ_kb

# Wholesale Interest Rate Constraint
Interest_Rate_C = R_b(t) ~ r(t) - k_b * (k_b / b_T(t)) * (K_b(t) / b_T(t))^2

# Define the Lagrangian
@syms λ_bs λ_lom λ_basel
L_W = Profit +
      λ_bs * (b_T(t) - (D_T(t) + K_b(t) + ξ_kb)) +
      λ_lom * ((1 - δ_b) * K_b(t) + Ω * ξ_kb - K_b(t+1) * π(t+1)) +
      λ_basel * Basel_II_C

# Compute FOCs
FOC_b_T = expand_derivatives(Differential(b_T(t))(L_W))
FOC_D_T = expand_derivatives(Differential(D_T(t))(L_W))
FOC_K_b = expand_derivatives(Differential(K_b(t))(L_W))

# Display the FOCs
println("Wholesale Branch's FOCs:")

println("\nFOC with respect to b_T(t) (Total loans):")
display(FOC_b_T)

println("\nFOC with respect to D_T(t) (Deposits):")
display(FOC_D_T)

println("\nFOC with respect to K_b(t) (Bank capital):")
display(FOC_K_b)

#%% retail branch 
using Symbolics

# Define callable variables
@variables t
@syms r_bn(t) r_bm(t) b_tn(t) b_tm(t) Adj_rbn(t) ε_bn(t) k_bn β_E λ_t λ_Loan

# Retail Branch Objective Function
Objective = r_bn(t) * b_tn(t) - r_bm(t) * b_tm(t) - Adj_rbn(t)

# Demand Function for Loans
Loan_Demand = b_tn(t) ~ (r_bn(t) / r_bm(t))^(-ε_bn(t)) * b_tm(t)

# Adjustment Cost on Retail Interest Rates
Adj_Cost = k_bn / 2 * ((r_bn(t) / r_bm(t) - 1)^2) * b_tm(t)

# Full Retail Branch Objective with Costs
Full_Objective = β_E^t * (r_bn(t) * b_tn(t) - r_bm(t) * b_tm(t) - Adj_Cost)

# Define Lagrangian
L_Retail = Full_Objective + λ_Loan * (b_tn(t) - (r_bn(t) / r_bm(t))^(-ε_bn(t)) * b_tm(t))

# Compute FOCs for Retail Branch
FOC_r_bn = expand_derivatives(Differential(r_bn(t))(L_Retail))
FOC_b_tn = expand_derivatives(Differential(b_tn(t))(L_Retail))
FOC_b_tm = expand_derivatives(Differential(b_tm(t))(L_Retail))

# Solve for λ_Loan from the loan demand constraint
λ_Loan_expr = symbolic_linear_solve(Loan_Demand, λ_Loan)

# Substitute λ_Loan into FOCs
FOC_r_bn_no_lambda = substitute(FOC_r_bn, λ_Loan => λ_Loan_expr) |> simplify
FOC_b_tn_no_lambda = substitute(FOC_b_tn, λ_Loan => λ_Loan_expr) |> simplify
FOC_b_tm_no_lambda = substitute(FOC_b_tm, λ_Loan => λ_Loan_expr) |> simplify

# Display Results
println("Retail Branch's FOCs:")

println("\nFOC with respect to r_bn(t) (Retail interest rate):")
display(FOC_r_bn_no_lambda)

println("\nFOC with respect to b_tn(t) (Total loans by retail branch):")
display(FOC_b_tn_no_lambda)

println("\nFOC with respect to b_tm(t) (Loan supply to retail branch):")
display(FOC_b_tm_no_lambda)


#%% wage
using Symbolics

# Define callable variables
using Symbolics

# Define variables and symbols
@variables t
@syms β_l W_L(t) P(t) l_i(t) k_w W_L_star(t) pi_var(t) pi_w_star(t) l_var(t)
@syms ε_l φ λ_l m λ_d

# Profit function for unions
Profit = β_l^t * (
    (W_L(t) / P(t)) * l_i(t) -
    k_w / 2 * ((W_L(t) / W_L_star(t-1)) - pi_var(t-1) * pi_w_star(t))^2 * W_L(t) / P(t) -
    (l_i(t)^(1+φ)) / (1+φ)
)

# Labor demand
Labor_Demand = l_i(t) ~ (W_L(t) / W_L_star(t-1))^(-ε_l) * l_var(t)

# New Keynesian Phillips Curve for Wage Inflation
NKPC = k_w * ((pi_w_star(t) / (pi_var(t-1) * pi_w_star(t))) - 1) * pi_w_star(t) ~
    β_l * λ_l * (k_w * ((pi_w_star(t+1) / (pi_var(t) * pi_w_star(t+1))) - 1)^2 / pi_var(t+1)) +
    (1 - ε_l) * l_i(t)^(1+φ) / λ_l

# Define Lagrangian
L_Labor = Profit + λ_d * (l_i(t) - (W_L(t) / W_L_star(t-1))^(-ε_l) * l_var(t))

# Compute FOCs for the labor market
FOC_W_L = expand_derivatives(Differential(W_L(t))(L_Labor))
FOC_l_i = expand_derivatives(Differential(l_i(t))(L_Labor))

# Solve for λ_d from the labor demand constraint
λ_d_expr = symbolic_linear_solve(Labor_Demand, λ_d)

# Substitute λ_d into FOCs
FOC_W_L_no_lambda = substitute(FOC_W_L, λ_d => λ_d_expr) |> simplify
FOC_l_i_no_lambda = substitute(FOC_l_i, λ_d => λ_d_expr) |> simplify

# Simplify NKPC for cleaner display
NKPC_simplified = simplify(NKPC)

# Display Results
println("Labor Market's FOCs:")

println("\nFOC with respect to W_L(t) (Nominal wage):")
display(FOC_W_L_no_lambda)

println("\nFOC with respect to l_i(t) (Labor supply):")
display(FOC_l_i_no_lambda)

println("\nNew Keynesian Phillips Curve for Wage Inflation:")
display(NKPC_simplified)

#%% Capital 
using Symbolics

# Define callable variables

using Symbolics

# Define variables
@variables t
@syms β_E λ_E q_k(t) i(t) i_E(t) k_i δ k(t) k_prev u(t) λ_C

# Define the utility function
Utility_Capital = β_E^t * λ_E * (q_k(t) * i(t) - k_i / 2 * (i_E(t) / k_prev - 1)^2 * i(t))

# Define the constraint (law of motion of capital)
Constraint_L = k(t) - (1 - δ) * k_prev - (1 - k_i / 2 * (i_E(t) / k_prev - 1)^2) * i(t) ~ 0

# Define the Lagrangian
L_Capital = Utility_Capital + λ_C * (k(t) - ((1 - δ) * k_prev + (1 - k_i / 2 * (i_E(t) / k_prev - 1)^2) * i(t)))

# Compute First-Order Conditions (FOCs)
FOC_q_k = expand_derivatives(Differential(q_k(t))(L_Capital))  # FOC with respect to q_k(t)
FOC_i = expand_derivatives(Differential(i(t))(L_Capital))      # FOC with respect to i(t)

# Solve for λ_C from the constraint
λ_C_expr = symbolic_linear_solve(Constraint_L, λ_C)

# Substitute λ_C into FOCs and other equations
FOC_q_k_no_lambda = substitute(FOC_q_k, λ_C => λ_C_expr)
FOC_i_no_lambda = substitute(FOC_i, λ_C => λ_C_expr)

# Simplify the expressions
FOC_q_k_no_lambda = simplify(FOC_q_k_no_lambda)
FOC_i_no_lambda = simplify(FOC_i_no_lambda)

# Price of Capital Equation
Price_Capital = simplify(q_k(t) * (1 - k_i / 2 * (i_E(t) / k_prev - 1)^2 + k_i * (i_E(t) / k_prev - 1) * (i_E(t) / k_prev)) +
                β_E * (λ_E / λ_E) * (u(t+1) * q_k(t+1) * (1 + k_i * (i_E(t+1) / i(t) - 1)) / u(t)) ~ 1)

# Law of Motion of Capital
Capital_Law = simplify(k(t) ~ (1 - δ) * k_prev + (1 - k_i / 2 * (i_E(t) / k_prev - 1)^2) * i(t))

# Display Results
println("Capital Producers' FOCs:")

println("\nFOC with respect to q_k(t) (Price of capital):")
display(FOC_q_k_no_lambda)

println("\nFOC with respect to i(t) (Investment):")
display(FOC_i_no_lambda)

println("\nPrice of Capital Equation:")
display(Price_Capital)

println("\nLaw of Motion of Capital:")
display(Capital_Law)


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






