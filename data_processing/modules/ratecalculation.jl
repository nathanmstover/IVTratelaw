"""
    rate_estimation(data::DataFrame, sigma::Float64)

Calculates the reaction rate for the provided experimental data per condition and the measurement error

"""

function rate_estimation(data::DataFrame, sigma::Float64)
    #data should include column with replicate number, time [hr], and preferred (ATP) concentration [M], including headers
    #sigma is the standard deviation for all observations in the given data set

    #separate data into df per replicate
    df_byrepl = groupby(data, :replicate)
 
    #find number of replicates
    n_repl = length(df_byrepl)

    #initialise storage df with rate and uncertainty
    rate_all_repl = zeros(n_repl) #column 1 = rate, column 2 = st dev
    
    #iterate over different replicates
    for (j, df_repl) in enumerate(df_byrepl)

        df_repl = DataFrame(df_repl)


        rate_error_intersect, variance, cov = rate_estimation_1repl_1cond(df_repl, sigma)

        rate_all_repl[j] = rate_error_intersect[1] #store the calculated rate
    end

    mean_rate = -mean(rate_all_repl.*4)
    variance_rate = var(rate_all_repl.*4)

    std_deviation = sqrt(variance_rate)

    return mean_rate, std_deviation
end  



""" 
    rate_estimation_1repl_1cond

Calculates the reaction rate for the provided experimental data per replicate per condition and the measurement error

"""

function rate_estimation_1repl_1cond(data::DataFrame, sigma::Float64)
    #data should include column with replicate number, time [hr], and preferred (ATP) concentration [M], including headers
    #sigma is the standard deviation for all observations in the given data set
 
    #initialise storage df with rate and uncertainty
    rate_error_intersect= zeros(1, 3) #column 1 = rate, column 2 = uncertainty, column 3 = intersection point
    
    n_times = size(data, 1)                                          # number of timepoints in that replicate
    X = [ones(n_times) data.time]                                      # construct matrix X = [1 timepoints]
    phi = sigma.^2*I(n_times)                                       #[sigma] = mM [phi] = mM^2
    Y = data.ATP    #[Y] =  mM                                      # actual results are the NTP concentrations   
        
    new_matrix = transpose(X)*inv(phi)*X

    if det(new_matrix) < 10^(-5)
            print("not invertable")
    else
    
        B = inv(transpose(X)*inv(phi)*X)*transpose(X)*inv(phi)*Y        # find B0 and B1 -> Y = B0 + B1*X
                                                                      #[B] = mM/hr 
        rate_error_intersect[1,1] = B[2]                          # store B1 (eg slope)
        rate_error_intersect[1,3] = B[1]                          # store intersection point
        
        cov = inv(transpose(X)*inv(phi)*X)                       # calculate covariance of b #[cov]=mM^2/hr^2
        times = data.time
        variance = [[x,1]'* cov*[x,1] for x in times]            # calculate the variance #[var]=mM^2/hr^2 -> check whether using right cov element
        variance_average = sum(variance) #mean(variance)   #average across time points
        rate_error_intersect[1,2] = sqrt(variance_average)        # store the st dev of the replicate
    
    end
    
    return rate_error_intersect, variance, cov
end  






"""
    rate_for_TASEP() - has to be updated with new outputs of rate estimation

estimate the rate for the entire data set and combines it with the pDNA and T7 concentrations
    
"""



function rate_for_TASEP(data::String, sigma::Float64)

    #prepare matrix with the right order of pDNA & T7 concentration 
    CDT_matrix = copydata(data)

    #call filtered data 
    filtered_data = importlabdata(data)

    #make dataframe of selected NTP
    filtereddata_NTP = select(filtered_data, [:condition, :replicate, :time, :ATP])

    #split up by different conditions, as function rate_estimation is called per condition
    filtereddata_NTP_bycond = groupby(filtereddata_NTP, :condition)

    # Initialize an empty array to store matrices
    output_matrices = []

    #iterate through different conditions
    for df_cond in filtereddata_NTP_bycond
        replicate_cond = df_cond[:,2]
        time_cond = df_cond[:,3]
        NTP_cond = df_cond[:,4]

        # create dataframe w heading
        data_cond = DataFrame(replicate = replicate_cond, time = time_cond, ATP = NTP_cond)

        mean_rate, std_deviation = rate_estimation(data_cond, sigma)
        # find the condition specific rate and variance

        # Store rate and uncertainty in a single matrix
        output_matrix = hcat(mean_rate, std_deviation)

        # Push the matrix into the output_matrices array
        push!(output_matrices, output_matrix)
    
    end

    # Combine all output matrices into one matrix
    combined_matrix = vcat(output_matrices...)
    
    # adding rate and uncertainty to CDT_matrix
    # Existing matrix
    existing_matrix = CDT_matrix

    # Matrix with new column(s) to add
    new_column_matrix = combined_matrix

    # Ensure dimensions match (number of rows in existing_matrix should match new_column_matrix)
    if size(existing_matrix, 1) == size(new_column_matrix, 1)
    # Concatenate existing_matrix and new_column_matrix horizontally (along columns)
        updated_matrix = hcat(existing_matrix, new_column_matrix)
    end

# Filter out rows where rate is negative
filtered_matrix = updated_matrix[(updated_matrix[:, 4] .>= 0), :]

###sorted and ready to input into TASEP
# Delete the first column
DT_matrix = filtered_matrix[:, 2:end]

# sorted by the DNA
sorted_matrix = DT_matrix[sortperm(DT_matrix[:, 1]), :] 

return sorted_matrix
end  