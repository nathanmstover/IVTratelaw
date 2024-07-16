"""
    plot_for_repl_cond

plot the raw data + estimated slope + confidence interval for a specific condition and replicate

"""

function plot_for_repl_cond(condition_n, replicate_n, filtered_data, sigma)

    #filter out data for specific condition and replicate number
    filtered_data_cond_repl = filter(row -> row[:condition] == condition_n && row[:replicate] == replicate_n, filtered_data)

    #rate estimation
    rate_error_intersect, variance, cov = rate_estimation_1repl_1cond(filtered_data_cond_repl, sigma)
  
    slope = rate_error_intersect[1,1]
    intersect = rate_error_intersect[1,3]

    x_data= filtered_data_cond_repl.time
    y_data = filtered_data_cond_repl.ATP

    p = scatter(x_data, y_data, label="Data", yerror = variance)

    y_pred = [slope*t + intersect for t in x_data]

    lower_bound = y_pred .- sqrt.(variance)
    upper_bound = y_pred .+ sqrt.(variance)

    plot!(p, x_data, y_pred, label = "prediction")

    plot!(p, x_data, lower_bound, fillrange = upper_bound, fillalpha = 0.3, label="CI", color="blue")

    xlabel!("Time")
    ylabel!("NTP Concentration")

end

"""
    plot_all_replicates

plot the raw data + estimated slope + confidence interval for all replicates within a specific condition

"""

function plot_all_replicates(condition_n, filtered_data, sigma)
    #filter out data for specific condition and replicate number
    filtered_data_cond = filter(row -> row[:condition] == condition_n, filtered_data)

     #split up by different conditions, as function rate_estimation is called per condition
     filtered_data_cond_repl = groupby(filtered_data_cond, :replicate)

     n_replicates = length(filtered_data_cond_repl)

     nrows = 1
     ncols = 3

     p = plot(layout = (nrows, ncols), size = (1000, 600))
     

     #iterate through different conditions
     for  (i, df_repl) in enumerate(filtered_data_cond_repl)

         df_repl = DataFrame(df_repl)
         #rate estimation
         rate_error_intersect, variance, cov = rate_estimation_1repl_1cond(df_repl, sigma)
  
         slope = rate_error_intersect[1,1]
         intersect = rate_error_intersect[1,3]

         x_data = df_repl[!, :time]
         y_data = df_repl[!, :ATP]

         scatter!(p, x_data, y_data, yerror = variance, subplot = i)

         y_pred = [slope*t + intersect for t in x_data]

         lower_bound = y_pred .- sqrt.(variance)
         upper_bound = y_pred .+ sqrt.(variance)

         plot!(p, x_data, y_pred, label = "Prediction", subplot = i)
         plot!(p, x_data, lower_bound, fillrange = upper_bound, fillalpha = 0.3, label="CI", color="blue", subplot = i)
         
         xlabel!("Time", subplot = i)
         ylabel!("NTP Concentration", subplot = i)
         ylims!(p, 1, 6)
 
     end
     return p

end


"""
    plot_all_replicates_all_conditions

plot the raw data + estimated slope + confidence interval for all conditions and replicates

"""


function plot_all_replicates_all_conditions(filtered_data, sigma)
    
    filtered_data_cond = groupby(filtered_data, :condition)

    n_conditions = length(filtered_data_cond)

    plot_list = []


    for (j, df_cond) in enumerate(filtered_data_cond)
        filtered_data_cond_repl = groupby(df_cond, :replicate)

        for (i, df_repl) in enumerate(filtered_data_cond_repl)

            df_repl = DataFrame(df_repl)
            #rate estimation
            rate_error_intersect, variance, cov = rate_estimation_1repl_1cond(df_repl, sigma)

            slope = rate_error_intersect[1,1]
            intersect = rate_error_intersect[1,3]

            x_data = df_repl[!, :time]
            y_data = df_repl[!, :ATP]

            subplot_index = (j, i)  

            p = scatter(x_data, y_data, yerror = variance)

            y_pred = [slope*t + intersect for t in x_data]

            lower_bound = y_pred .- sqrt.(variance)
            upper_bound = y_pred .+ sqrt.(variance)

            plot!(p, x_data, y_pred, label = "Prediction")
            plot!(p, x_data, lower_bound, fillrange = upper_bound, fillalpha = 0.3, label="CI", color="blue")
        
            xlabel!("Time")
            ylabel!("NTP Concentration")
            ylims!(2,5.5)

            push!(plot_list,p)
        end
    end

    # Create subplots
    ncols = 3
    nrows = n_conditions
    plot(plot_list..., layout=(nrows, ncols), size = (3000,2000))

end


"""
    bar_rates_all_repl

create bar plot of the rate and standard deviation for all replicates at specific condition

"""



function bar_rates_all_repl(condition_n, filtered_data, sigma)
    #filter out data for specific condition and replicate number
    filtered_data_cond = filter(row -> row[:condition] == condition_n, filtered_data)

     #split up by different conditions, as function rate_estimation is called per condition
     filtered_data_cond_repl = groupby(filtered_data_cond, :replicate)

     n_replicates = length(filtered_data_cond_repl)

     #p = plot(layout = (nrows, ncols), size = (1000, 600))
     slopes = zeros(n_replicates)
     replicate_n =  1:n_replicates
     st_dev = zeros(n_replicates)

     #iterate through different conditions
     for  (i, df_repl) in enumerate(filtered_data_cond_repl)

         df_repl = DataFrame(df_repl)
         #rate estimation
         rate_error_intersect, variance, cov = rate_estimation_1repl_1cond(df_repl, sigma)
  
         slope = rate_error_intersect[1,1]
         slopes[i] = -slope
         st_dev[i] = sqrt(cov[2,2])

     end
     
    # Filter out points where slope is zero
    non_zero_indices = findall(x -> x != 0, slopes)
    filtered_replicate_n = replicate_n[non_zero_indices]
    filtered_slopes = slopes[non_zero_indices]
    filtered_st_dev = st_dev[non_zero_indices]

    # Plot with filtered data
    scatter(filtered_replicate_n, filtered_slopes, yerror=filtered_st_dev, label="Slopes",
        xlabel="Replicate", ylabel="Slope", title="Bar Plot of Slopes")
end






"""
    function bar_rates_all_repl_all_cond
create bar plot of the rate and standard deviation for all conditions and replicates

"""


function bar_rates_all_repl_all_cond(filtered_data, sigma)
       
    filtered_data_cond = groupby(filtered_data, :condition)

    n_conditions = length(filtered_data_cond)

    plot_list = []

    #iterate through conditions
    for (j, df_cond) in enumerate(filtered_data_cond)
        filtered_data_cond_repl = groupby(df_cond, :replicate)

        slopes = zeros(3)
        replicate_n =  1:3
        st_dev = zeros(3)

        #iterate through replicates
        for (i, df_repl) in enumerate(filtered_data_cond_repl)

            df_repl = DataFrame(df_repl)
            #rate estimation
            rate_error_intersect, variance, cov = rate_estimation_1repl_1cond(df_repl, sigma)
  
            slope = rate_error_intersect[1,1]
            slopes[i] = -slope
            st_dev[i] = sqrt(cov[2,2])

           
        end
        title_str = "Condition $j "

        p = scatter(replicate_n, slopes, yerror = st_dev, legend = false, ylabel="Slope", title=title_str)
    
     push!(plot_list,p)
    end

    # Create subplots
    ncols = 1
    nrows = n_conditions

    plot(plot_list..., layout=(nrows, ncols), size = (500,2000))
    

end