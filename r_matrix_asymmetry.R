#!/usr/bin/env Rscript





############# CHANGE THIS #################


# These are the objects that are read from the files and plotted
plot_labels <- c("s_trx2", "trx2(1)", "trx2(2)", "trx2(3)", "trx2(4)", "trx2(5)", "trx2(6)", "trx2(7)", "trx2(8)", "trx2(9)")
#plot_labels <- c("energy", "Polyak", "s_trx2", "Myers", "trx2(1)", "trx2(4)")


#thin out data to plot and process by this value
max_points_per_series = 1000


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

PrintScriptName("R matrix asymmetry script")


printf("Arguments: target_directory file_pattern itraj_start itraj_stop \n\n")


printf("e.g.: ./r_matrix_asymmetry_production_lowT.R prod_broad_output GN16S30T0.3M0.5D9_F12 2000 -1 \n\n")

printf("itraj_stop < 1: no upper limit on trajectory number.\n\n")



all_data <- ReadOutFileToDataTable(target_directory, file_pattern, itraj_start, itraj_stop, plot_labels, parallel)




printf("Saving matrix sizes to file...\n\n")

matrix_frame <- data.frame(matrix(ncol=4, nrow=10))
colnames(matrix_frame) <- c("matrix", "size", "error", "S")


for(i in 1:9) {
    printf("\nProcessing matrix %d\n", i)
    matrix_frame[i, "matrix"] <- i
    matrix_frame[i, "size"] <- mean(all_data[, get(p("trx2(", i, ")"))])

    acl <- AutocorrelationLength(all_data[ , get(p("trx2(", i, ")"))])
    error <- JackknifeDelta_WSingle(all_data[ , get(p("trx2(", i, ")"))], acl, 3)
    matrix_frame[i, "error"] <- error
}

printf("\nProcessing sum_trx2\n", i)

matrix_frame[10,"matrix"] <- 10
matrix_frame[10,"size"] <- mean(all_data[ ,get("s_trx2")])

acl <- AutocorrelationLength(all_data[ , get("s_trx2")])
error <- JackknifeDelta_WSingle(all_data[ , get("s_trx2")], acl, 3)
matrix_frame[10,"error"] <- error


S <- as.numeric(gsub(".*S(.+)T.*", "\\1", file_pattern))
for(i in 1:10) {
    matrix_frame[i, "S"] <- S
}

matrix_file <- p("matrix_sizes_", file_pattern, ".csv")

write.table(matrix_frame, file = matrix_file, qmethod = "double", row.names=FALSE, sep = ",")


N <- as.numeric(gsub(".*N(.+)S.*", "\\1", file_pattern))



printf("Loading plot library...\n\n")

invisible(library("ggplot2"))
#invisible(library("dplyr"))


if(nrow(all_data) > max_points_per_series) {
    thin_out_by <- nrow(all_data) %/% max_points_per_series
    printf("thin_out_by set to %d.\n\n", thin_out_by)
} else {
    thin_out_by <- 1
}




# Create Plots:
printf("Creating MC histories...\n\n")



plot_file <- paste0("plot_", args[1])
plot_file <- paste0(plot_file, "_")
plot_file <- paste0(plot_file, file_pattern)
plot_file <- paste0(plot_file, "_") 
plot_file <- paste0(plot_file, "Polyak")
plot_file <- paste0(plot_file, ".pdf")
#printf("%s\n", plot_file)
pdf(plot_file)

#color_define = c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2")

#color_define = c("#fff100", "#ff8C00", "#e81123", "#ec008C", "#68217a", "#00188f", "#00bcf2", "#00b294", "#009e49", "#bad80a")

#n=10

#color_define = rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)

color_define = c("darkblue", "deeppink2", "grey", "red", "chartreuse", "cyan", "yellow", "black", "green3", "darkorange")

all_data_MC_plot<-all_data[itraj %% thin_out_by == 0]

#ymax = max(2.0, 2*max(all_data_MC_plot[,4:9]))
    
#Correct working version 
#plot(all_data_MC_plot[, get("itraj")], all_data_MC_plot[,get("Polyak")], xlab="MC time", ylab="P", main=file_pattern, pch=4, ylim=c(0, ymax), yaxp  = c(0.0, 1.0, 10), col="blue")

plot(all_data_MC_plot[,get("itraj")], all_data_MC_plot[,get("s_trx2")]/10, xlab="MC time", ylab="Matrix size" , ylim=c(0, 2), yaxp  = c(0.0, 2, 20), pch=6, col=color_define[10], cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)


grid(nx=NULL, ny=NULL)	
abline(h=(1:40)*0.05, col="gray", lty=3)

#points(all_data_MC_plot[,get("itraj")], 10*all_data_MC_plot[,get("Myers")], pch=5, col="red")


for(i in 1:9) {
    #print(color_define[i])
    points(all_data_MC_plot[,get("itraj")], 2*all_data_MC_plot[,get(p("trx2(", i, ")"))], pch=10+i, col=color_define[i])
}

Ntext <- p("N=", N)

text(1, 1.95, Ntext, adj = c(0,0), cex = 1.2)


#points(all_data_MC_plot[,get("itraj")], all_data_MC_plot[,get("accept")], pch=20, col=20)

legend("bottom", legend=c(expression(2~Tr~(X[1]^2)), expression(2~Tr~(X[2]^2)), expression(2~Tr~(X[3]^2)), expression(2~Tr~(X[4]^2)), expression(2~Tr~(X[5]^2)), expression(2~Tr~(X[6]^2)), expression(2~Tr~(X[7]^2)), expression(2~Tr~(X[8]^2)), expression(2~Tr~(X[9]^2)), expression(R^2 / 10)), col=color_define, pch=c(11:19, 6), cex=0.8, ncol = 5)




dev.off()
printf(".")

