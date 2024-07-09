function adjustExcel(filename::AbstractString, data::Matrix, start_row::Int, start_col::Int, new_filename::AbstractString)
    if isfile(filename)
        # Rename existing file with timestamp
        timestamp = Dates.format(now(), "yy-mm-dd")
        renamed_filename = replace(filename, ".xlsx" => "_$timestamp.xlsx")
        mv(filename, renamed_filename, force=true)
        println("Existing file renamed to: $renamed_filename")
    end
    
    # Open an existing workbook or create a new one
    wb = XLSX.openxlsx(filename, mode="w")
    
    # Select the first sheet in the workbook
    sheet = wb[1]
    
    # Write data to specific cells starting from start_row, start_col
    rows, cols = size(data) 
    for r in 1:rows
        for c in 1:cols
            XLSX.setdata!(sheet, start_row + r - 1, start_col + c - 1, data[r, c])
        end
    end
    
    # Save and close the workbook
    XLSX.writexlsx(filename, wb)
end