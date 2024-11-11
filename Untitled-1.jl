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

#%% Cleaning data
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
 
 # Function to rename columns by removing leading/trailing spaces and replacing spaces with underscores
function remove_spaces!(df::DataFrame)
    for col in names(df)
        # Remove leading/trailing spaces and replace internal spaces with underscores
        cleaned_col = strip(replace(col, r"\s+" => "_"))
        if cleaned_col != col
            rename!(df, col => cleaned_col)
        end
    end
end

# Loop through each DataFrame and apply the renaming function
for df in dfs
    try
        remove_spaces!(df)
    catch e
        println("Error processing DataFrame: ", e)
        continue
    end
end

#%%data handling 
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
lagged_values = ShiftedArrays.lag(housepriceindex[!, "_House_Price_Index_"], 1)
growth_rates = [missing; diff(housepriceindex[!, "_House_Price_Index_"]) ./ lagged_values[2:end] .* 100]

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

# Function to filter rows based on the month of the date column
function filter_quarter_months!(df::DataFrame)
    date_col = names(df)[1]  # Assume the first column is the date column
    # Ensure the date column is parsed as `Date` type
    df[!, date_col] = Date.(df[!, date_col])  # Convert to Date type if necessary

    # Filter for rows where the month is in [01, 04, 07, 10]
    filtered_df = filter(row -> month(row[date_col]) in [1, 4, 7, 10], df)

    return filtered_df
end

# Loop through each DataFrame, apply the filter, and update the DataFrame
for i in 1:length(dfs)
    try
        dfs[i] = filter_quarter_months!(dfs[i])
        println("Successfully filtered DataFrame $i")
    catch e
        println("Error processing DataFrame $i: ", e)
    end
end
#%% rename variable to match the model 
using DataFrames

# Define the mapping of original dataset names to model variable names
rename_map = Dict(
    "PCE" => "C",                 # Consumption
    "PNFI" => "I",                # Investment
    "inflation" => "π",           # Inflation
    "comph" => "W",               # Wage
    "Depositshousehold" => "h",   # Housing Services / Loans to Households
    "corporatelendingxlsx" => "L",# Corporate Lending
    "housepriceindex" => "qh",    # House Price Index
    "moodyCorporateBondYield" => "rc", # Corporate Lending Rates
    "t3_MonthTreasurey" => "r",   # Policy Rate
    "t30_YearFixedRateMortgageAverage" => "rh" # Mortgage Rates
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

#%% define relevant time frame 
using DataFrames, Dates

# List of DataFrame variables to process
dfs = [
    :t3_MonthTreasurey,
    :t30_YearFixedRateMortgageAverage,
    :Depositshousehold,
    :PCE,
    :PNFI,
    :ResidentialMortgages,
    :comph,
    :corporatelendingxlsx,
    :housepriceindex,
    :inflation,
    :moodyCorporateBondYield
]

# Function to filter each DataFrame to a specific date range
function filter_date_range!(df_name::Symbol, start_date::Date, end_date::Date)
    try
        # Retrieve the DataFrame
        df = getfield(Main, df_name)

        # Ensure the first column is of Date type
        date_col = names(df)[1]
        if eltype(df[!, date_col]) != Date
            df[!, date_col] = Date.(df[!, date_col], dateformat"yyyy-mm-dd")
        end

        # Filter the DataFrame to the specified date range
        filtered_df = filter(row -> row[date_col] >= start_date && row[date_col] <= end_date, df)

        # Assign the filtered DataFrame to a new global variable with "_timeframe" suffix
        new_name = Symbol(string(df_name) * "_timeframe")
        @eval global $new_name = $filtered_df

        println("Filtered DataFrame '$df_name' to date range: $start_date to $end_date")
    catch e
        println("Error processing DataFrame '$df_name': ", e)
    end
end

# Define the start and end dates
start_date = Date("1991-01-01")
end_date = Date("2005-01-01")

# Loop through each DataFrame and apply the date range filter
for df_name in dfs
    filter_date_range!(df_name, start_date, end_date)
end

