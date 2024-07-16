"""
    generate_rate_uncertainty

    one function call that takes in the csv data file and ouputs the csv data file with DNA concentration, T7 concentration
        reaction rate and uncertainty, including titles
"""
function generate_rate_uncertainty(datafileivt::String, new_filename::AbstractString)
    
    sigma = 0.27 #set from looking at IE data up until 7/15, variable to change

    sorted_matrix = rate_for_TASEP(datafileivt, sigma) #generate the rates and uncertainties

    #add titles to the columns
    column_titles = ["DNA (ug/uL)", "T7 RNAP (U/uL)", "Rate (mAu*s/uM)", "Uncertainty (mAu*s/uM)"] 

    df = DataFrame(
        [sorted_matrix[:, i] for i in 1:size(sorted_matrix, 2)],
        Symbol.(column_titles)
    )

    #generate full path for output CSV file, such that it is stored in 'output' folder
    output_folder = "output"
    output_path = joinpath(pwd(), output_folder, "$new_filename.csv")

    # check if an existing file with this name exists, then rename it
    if isfile(output_path)
        # Rename existing file with timestamp
        timestamp = Dates.format(now(), "yy-mm-dd")
        renamed_filename = replace(output_path, ".csv" => "_$timestamp.csv")
        mv(output_path, renamed_filename, force=true)
        println("Existing file renamed to: $renamed_filename")
    end

    # Write the combined matrix to a new CSV file
    CSV.write(output_path, df)
end