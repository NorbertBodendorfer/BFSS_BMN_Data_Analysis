#!/usr/bin/env Rscript





############# CHANGE THIS #################

##### These are all reasonable things to plot. 
####  c("h_f-h_i", "CG_it", "energy", "Polyak", "s_trx2", "trx2(1)", "trx2(2)", "trx2(3)", "trx2(4)", "trx2(5)", "trx2(6)", "trx2(7)", "trx2(8)", "trx2(9)", "comm2", "Myers", "accept")


# These are the objects that are read from the files and plotted
plot_labels <- c("Polyak", "energy", "s_trx2", "trx2(1)", "trx2(4)", "Myers")




#thin out data to plot and process by this value, we keep only every n-th value
thin_out_by = 5


# Run script in parallel (0 or 1=parallel)
parallel = 0

############# STOP CHANGING ###############



## Correlogram code from stackexchange 

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE, breaks = 24)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}




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




######### Start of program output

PrintScriptName("R correlogram script")


printf("Arguments: target_directory file_pattern itraj_start itraj_stop \n\n")


printf("e.g.: ./r_correlogram.R prod_broad_output GN24S12T0.543M2.0D9_F1 1 20000 \n\n")

printf("itraj_stop < 1: no upper limit on trajectory number.\n\n")



all_data <- ReadOutFileToDataTable(target_directory, file_pattern, itraj_start, itraj_stop, plot_labels, parallel)

ntraj <- length(all_data[,get("itraj")])


if(thin_out_by > 1) {
    printf("Thinning out by %d\n", thin_out_by)
    all_data<-all_data[itraj %% thin_out_by == 0]
}

# Remove itraj because we don't need it in correlogram
all_data <- all_data[,itraj:=NULL]



pairs_title <- ""
if(grepl("GBN", file_pattern)) {
    pairs_title <- p(pairs_title, "Gauged bosonic ")
} else if(grepl("GN", file_pattern)) {
    pairs_title <- p(pairs_title, "Gauged supersymmetric ")
} else if(grepl("UN", file_pattern)) {
    pairs_title <- p(pairs_title, "Ungauged supersymmetric ")
} else if(grepl("UBN", file_pattern)) {
    pairs_title <- p(pairs_title, "Ungauged supersymmetric ")
}

if(grepl("M0D", file_pattern)) {
    pairs_title <- p(pairs_title, "BFSS, ")
} else {
    pairs_title <- p(pairs_title, "BMN, ")
}

pairs_title <- p(pairs_title, "T=", getstr(file_pattern, 'T', 'M'), ", ")
pairs_title <- p(pairs_title, "mu=", getstr(file_pattern, 'M', 'D'), ", ")
pairs_title <- p(pairs_title, "N=", getstr(file_pattern, 'N', 'S'), ", ")
pairs_title <- p(pairs_title, "S=", getstr(file_pattern, 'S', 'T'), ", ")
pairs_title <- p(pairs_title, "d=", getstr(file_pattern, 'D', '_'), ", ")
pairs_title <- p(pairs_title, "#traj=", ntraj)



plot_file_name <- p("correlogram_", file_pattern, ".pdf")

pdf(plot_file_name)

pairs(all_data, diag.panel = panel.hist, lower.panel = panel.smooth, upper.panel = panel.cor, main=pairs_title)

dev.off()

printf("Correlogram done.\n\n")

