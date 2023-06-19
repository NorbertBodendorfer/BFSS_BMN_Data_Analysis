# The method ReadOutFileToDataTable is used to efficienlty read a number of output files and join the measurements in a single datatable. 
# Parallelization is possible but has not shown significant advantage in the usecases. 


printf("\n\nLoading read libraries. \n ")
library(data.table)
printf("\n")


# Returns dataframe containing columns file_list, header_size
GetOutputFileStructure <- function(file_list) {
	#print(file_list)

    # This is the header size under which the script will complain, as headers should be larger than this
    dangerously_low_header_size <- 30

	header_size_list <- c()

	for(f in file_list) {
		header_size <- as.numeric(try(system(p("awk '/#-----------/{ print NR; exit }' ", f), intern = TRUE))) 
		header_size_list <- c(header_size_list, header_size)
		if(length(header_size) == 0 || header_size < dangerously_low_header_size) {
			# We print header_size before stop because sprintf cannot handle if awk returned nothing to header_size. In that case, header_size reads "numeric(0)"
			printf("header_size read as ")
			print(header_size)
			stop(sprintf("Unusual header size of extracted from %s.", f))
		} 
	}

	structure <- data.frame(file_list, header_size_list)
	colnames(structure) <- c("file_name", "header_size")
	return(structure)
}










ReadOutFileToDataTable <- function(target_directory, file_pattern, itraj_start, itraj_stop, plot_labels, parallel) {

    ######## Objects / lists for the script


    ### Create list of columns to fread

    # Initialize col_vec with entry 1, which is for itraj, which we need to plot against and to remove doubles from icon starts. 
    col_vec <- c(1)

    # Define the column labels. These are all possible measurements extracted from the files
    column_labels <- c("itraj", "h_f-h_i", "num_cv", "n_b_CG", "CG_it", "energy", "Polyak", "s_trx2", "trx2(1)", "trx2(2)", "trx2(3)", "trx2(4)", "trx2(5)", "trx2(6)", "trx2(7)", "trx2(8)", "trx2(9)", "comm2", "Myers", "accept")

    # Loop over labels to plot and append the number w.r.t. column_labels to col_vec
    # col_vec is then passed to fread and reads only those columns
    for(x in plot_labels) {
        col_vec <- c(col_vec,match(x, column_labels))
    }

    # itraj should be included as a column in the data.table all_data
    plot_labels_pitraj <- c("itraj", plot_labels)


    # Get list of files to include in data
    file_pattern_to_read <- p("out", file_pattern, ".*")

    printf("Scanning target directory %s for regex %s:\n\n", target_directory, file_pattern_to_read)
    files_to_read <- list.files(path=target_directory, pattern=file_pattern_to_read, full.names=TRUE, recursive=FALSE)

    if(length(files_to_read)==0) {
        stop("No machting files found. Stopping.")
    }

    file_structure <- GetOutputFileStructure(files_to_read)


    if(parallel == 1) {

        printf("Parallel reading enabled.\n")

        library(foreach)
        library(doParallel)

        numCores <- detectCores()

        printf("Using %d cores.\n", numCores)
        registerDoParallel(numCores)

        # Loop through file list and add data to all_data frame
        all_data <- foreach(i=1:length(file_structure$file_name), .combine=rbind) %dopar% {
            f <- as.character(file_structure$file_name[i])
            header_size <- as.numeric(file_structure$header_size[i])

            printf("\tReading %s with header size %d.\n", f, header_size)	
            tryCatch({
                try(fread(f, skip=header_size, select = col_vec)) # The sep="" command messes up fread, needs to be dropped. 
            }, warning = function(w) {
                warning(w)
            }, error = function(e) {
                stop(e)
            }, finally = {
            })
            
        } 
        stopImplicitCluster()
    } else {
        library(foreach)
        
        # Loop through file list and add data to all_data frame
        all_data <- foreach(i=1:length(file_structure$file_name), .combine=rbind) %do% {
            f <- as.character(file_structure$file_name[i])
            header_size <- as.numeric(file_structure$header_size[i])

            printf("\tReading %s with header size %d.\n", f, header_size)	
            tryCatch({
                try(fread(f, skip=header_size, select = col_vec)) # The sep="" command messes up fread, needs to be dropped. 
            }, warning = function(w) {
                warning(w)
            }, error = function(e) {
                stop(e)
            }, finally = {
            })
            
        } 
    }



    colnames(all_data) <- plot_labels_pitraj

    printf("\nDone reading. %d total trajectories read.\n\n", nrow(all_data))

    if(nrow(all_data)==0) {
        stop("0 trajectories read.")
    }



    printf("\nApplying trajetory cuts from ReadOutFileToDataTable function:\n\n")

    printf("\tLower limit %d: ", itraj_start)
    before_low <- nrow(all_data)
    all_data <- all_data[itraj >= itraj_start]
    after_low <- nrow(all_data)

    printf("%d trajectories removed from beginning.\n", before_low-after_low)


    if( itraj_stop > 0) {
        printf("\tUpper limit %d: ", itraj_stop)
        before_high <- nrow(all_data)
        all_data <- all_data[itraj <= itraj_stop]
        after_high <- nrow(all_data)
        printf("%d trajectories removed from end.\n", before_high-after_high)
    } else {
        printf("\tNo upper limit\n")
    }

    # These are trajetories that have the same itraj, but different observations. 
    printf("\tRemoving doubles: ")
    before_double <- nrow(all_data)
    all_data <- unique(all_data, by = "itraj")
    after_double <- nrow(all_data)
    printf("%d doubles removed.\n", before_double-after_double)



    # Force all entries of all_data to be numeric. Some may be read as characters

    #printf("Turning data numeric.\n\n")
    #all_data <- mutate_all(all_data, function(x) as.numeric(as.character(x)))


    # Remove incomplete lines, where certain observables could not be evaluated. 
    printf("\tRemoving NAs: ")
    before_NA <- nrow(all_data)
    #all_data[is.na(all_data)] <- 0
    all_data <- na.omit(all_data)
    after_NA <- nrow(all_data)
    printf("%d trajectories removed.\n\n", before_NA-after_NA)




    return(all_data)

}







