#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)





######## Objects / lists for the script


##### These are all reasonable things to plot against itraj, which is done in this script. 
####  c("h_f-h_i", "CG_it", "energy", "Polyak", "s_trx2", "trx2(1)", "trx2(2)", "trx2(3)", "trx2(4)", "trx2(5)", "trx2(6)", "trx2(7)", "trx2(8)", "trx2(9)", "comm2","Myers","accept")



############# CHANGE THIS #################

# These are the objects that are read and plotted
plot_labels <- c("Polyak", "s_trx2")



#thin out data to plot and process by this value, we keep only every n-th value
thin_out_by = 5


#parallel reading
parallel <- 0

############# STOP CHANGing ###############








######### Settings

source("r_GeneralSettings.R")


######### Function defines

source("r_GeneralFunctions.R")
source("r_ReadOutFilesToDataTable.R")






### Start of program output

printf("\n\n\n|------------------------|\n| R data analysis script |\n|------------------------|\n\n")

printf("Arguments: target_directory file_pattern itraj_start itraj_stop \n\n")


printf("e.g.: ./r_data_analysis.Rd BFSS_prod_output GBN48S48T0.872M0D25_F1_ 10000 50000 \n\n")

printf("Upper limit < 1: no upper limit.\n\n")



target_directory <- args[1]
file_pattern <- args[2]
itraj_start <- as.numeric(args[3])
itraj_stop <- as.numeric(args[4])


printf("\n\nLoading read libraries... \n ")
library(data.table)
printf("\n")



plot_labels_pitraj <- c("itraj", plot_labels)



all_data <- as.data.frame(ReadOutFileToDataTable(target_directory, file_pattern, itraj_start, itraj_stop, plot_labels, parallel))






printf("\n\nLoading plot library...  ")

invisible(library("ggplot2"))

invisible(library("dplyr"))


all_data <- mutate_all(all_data, function(x) as.numeric(as.character(x)))
all_data[is.na(all_data)] <- 0



printf("Done.\n\n")





printf("Making plots for")
print(plot_labels)


for(plot_loop in plot_labels){
	plot_file <- p("MChistory_", file_pattern, "_", plot_loop, ".pdf")
	
	pdf(plot_file)
	
	#thin out data
	all_data_MC_plot<-subset(all_data, all_data[,"itraj"] %% thin_out_by == 0)

	
	#dev.new(width=8, height=8, unit="cm") # for production plots in paper
	#Correct working version 
	plot(all_data_MC_plot[,"itraj"], all_data_MC_plot[,plot_loop], type="l", xlab="MC time", ylab=plot_loop, pch=4, cex=1.0, cex.axis=1.5, cex.lab=1.5)

    # other plot options for production
    # ylim=c(0, 0.7)


	printf("Saved to %s\n\n", plot_file)
	#abline(h = mean(all_data[,x]), col = "blue")
	dev.off()
}

printf("Done.\n")


