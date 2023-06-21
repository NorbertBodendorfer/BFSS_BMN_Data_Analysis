#!/usr/bin/env Rscript





############# CHANGE THIS #################



# These are the objects that are read from the files and plotted
#observable_labels <- c("energy", "Polyak", "s_trx2", "Myers", "trx2(1)", "trx2(2)", "trx2(3)", "trx2(4)", "trx2(5)", "trx2(6)", "trx2(7)", "trx2(8)", "trx2(9)", "comm2","Myers")

observable_labels <- c("energy", "Polyak", "s_trx2", "Myers")



# Run script in parallel (0 or 1=parallel)
parallel = 0

############# STOP CHANGING ###############






######### Get command arguments
args = commandArgs(trailingOnly=TRUE)

target_directory <- args[1]
file_pattern <- args[2]
itraj_start <- as.numeric(args[3])
itraj_stop <- as.numeric(args[4])




######### Settings

source("r_GeneralSettings.R")


######### Function defines

source("r_GeneralFunctions.R")
source("r_ReadOutFilesToDataTable.R")
source("r_ErrorBars.R")




######### Start of program output

PrintScriptName("R observable extractor script")


printf("Arguments: target_directory file_pattern itraj_start itraj_stop \n\n")


printf("e.g.: ./r_observable_extractor.R prod_broad_output GN16S30T0.3M0.5D9_F12_ 2000 -1 \n\n")

printf("itraj_stop < 1: no upper limit on trajectory number.\n\n")



all_data <- ReadOutFileToDataTable(target_directory, file_pattern, itraj_start, itraj_stop, observable_labels, parallel)




printf("Saving matrix sizes to file...\n\n")

n_obs <- length(observable_labels)

S <- as.numeric(gsub(".*S(.+)T.*", "\\1", file_pattern))

N <- as.numeric(gsub(".*N(.+)S.*", "\\1", file_pattern))

T <- as.numeric(gsub(".*T(.+)M.*", "\\1", file_pattern))

M <- as.numeric(gsub(".*M(.+)D.*", "\\1", file_pattern))



#observable_frame <- data.frame(matrix(ncol=5, nrow=n_obs))
#colnames(observable_frame) <- c("observable", "value", "error", "S", "N")

observable_frame <- foreach(i=1:n_obs, .combine=rbind) %do% {
    print(observable_labels[i])
    data.table( observable = observable_labels[i], value = mean(all_data[, get(observable_labels[i])]), error = JackknifeDelta_WSingle(all_data[ , get(observable_labels[i])], AutocorrelationLength(all_data[ , get(observable_labels[i])]), 1), S = S, N = N, T = T, M = M)
}

print(observable_frame)


observable_file <- p("observables_", file_pattern, ".csv")

write.table(observable_frame, file = observable_file, qmethod = "double", row.names=FALSE, sep = ",")


printf("Observables written to %s.\n", observable_file)

