"""
    importlabdata()

Import excel file with data from experiments with various conditions

"""
function importlabdata(filename::String)
    # load excel file into DataFrame
    raw_data = XLSX.readtable(filename, "Sheet1") |> DataFrame
    raw_data_matrix = Matrix(raw_data)
    
    raw_data_matrix = map(x -> tryparse(Float64, string(x)) !== nothing ? parse(Float64, string(x)) : NaN, raw_data_matrix)

    # create a table with condition number, replicate number, time (hr), mass (kg), NTP conc in Mx4
    processeddata = zeros(size(raw_data_matrix, 1), 8) 

    # fill processed data matrix
    processeddata[:,1] = raw_data_matrix[:,1] # condition number
    processeddata[:,2] = raw_data_matrix[:,2] # replicate number
    processeddata[:,3] = raw_data_matrix[:,5] ./60 # Min to hr
    processeddata[:,4] = raw_data_matrix[:,6] .* 0.001 # g to kg
    processeddata[:,5:8] = raw_data_matrix[:,7:10] #.* 0.001 # units from mM to M

    processed_data_df = DataFrame(processeddata, :auto)
    processed_data_df= DataFrame(processed_data_df, [:condition, :replicate, :time, :mass, :ATP, :UTP, :CTP, :GTP])

    processed_data_df = filter(col -> col[end] != 0, processed_data_df)
    processed_data_df = filter(row -> 0.0055 < row[:mass] < 0.0064, processed_data_df)    
    
    return processed_data_df
end





"""
    copydata()

function uses the excel datafile and makes a matrix that only has the condition number, dna imputs and T7 inputs.
It sorts it in the order such that it is used in the TASEP code

"""



function copydata(filename::String)
    # load excel file into DataFrame
    excel_data = XLSX.readtable(filename, "Sheet1") |> DataFrame
    excel_data_matrix = Matrix(excel_data)
    
    # convert each element to a Float64, else assign NaN
    excel_data_matrix = map(x -> tryparse(Float64, string(x)) !== nothing ? parse(Float64, string(x)) : NaN, excel_data_matrix)

    # create a table with condition, DNA input concetrations and T7 inputs for sorting
    condata = zeros(size(excel_data_matrix, 1), 3) 

    # fill processed data matrix
    condata[:,1] = excel_data_matrix[:,1] # condition
    condata[:,2] = excel_data_matrix[:,3] # DNA inputs
    condata[:,3] = excel_data_matrix[:,4] # T7 inputs
    
    con_data_df = DataFrame(condata, :auto)
    con_data_df= DataFrame(con_data_df, [:condition, :DNAinputs, :T7inputs])

    # combine and sort data within the function
    first_rows = combine(groupby(con_data_df, :condition)) do group
        group[1, 2:end]  # Extract columns excluding :condition for the first row
    end
    
    # Convert the first rows into a matrix
    CDT_matrix = Matrix(first_rows[:, 1:end]) 

    return CDT_matrix
end