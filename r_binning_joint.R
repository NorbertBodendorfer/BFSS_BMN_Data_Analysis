#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 



############# CHANGE THIS #################

##### These are all reasonable things to plot. 
####  c("h_f-h_i", "CG_it", "energy", "Polyak", "s_trx2", "trx2(1)", "trx2(2)", "trx2(3)", "trx2(4)", "trx2(5)", "trx2(6)", "trx2(7)", "trx2(8)", "trx2(9)", "comm2","Myers","accept")


# These are the objects that are read from the files and plotted
plot_labels <- c("Polyak", "s_trx2", "energy", "Myers")

#Number of bins for the histogram
nbins=35


# Compute the autocorrelation time for jackknife error estimation in histograms. If not, n_jackbin is used for the number of bins
use_acl <- 1
n_jackbin <- 2


# If we have use_acl = 1, then w_histjack = mult_histjack * acl, i.e. we can take mult_histjack times wider bins as opposed to estimated using acl
mult_histjack <- 1



#Include error bars?
eb <- 1


# Measurement domain of P to be shown. Other measurement domains are taken from data 
min_measurement_P <- 0.0
max_measurement_P <- 0.8


# Number of standard deviations to cut off data from below for energy measurement for gauged bosonic runs. This is necessary due to outliers
nSD_GB <- 5

# Number of rows (temperatures) in the legend
n_row_templist <- 1

# Run script in parallel (0 or 1=parallel)
parallel <- 0

# Test normalisation of the histgrams and output result. Not really necessary as code seems to work reliably
test_normalisation <- 0

############# STOP CHANGING ###############








######### Settings

source("r_GeneralSettings.R")


######### Function defines

source("r_GeneralFunctions.R")
source("r_ReadOutFilesToDataTable.R")
source("r_ErrorBars.R")




## Read command arguments
nskip <- as.numeric(args[1])
target_directory <- args[2]
file_pattern <- args[3]
D=args[4]
flux=args[5]



# Start output
PrintScriptName("R joint histograms of various temperatures")

printf("Arguments: num_skip target_directory file_pattern D M T_1 T_2 ...\n\n\n")


printf("e.g.: ./r_binning_joint.R 1000 prod_broad_output GN24S12 9 2.0 0.542 0.543 0.544\n\n")


# We are expecting this many temperatures
ndata = length(args)-5
if(ndata < 1) {
    stop("Not enough arguemnts supplied, need at least 6, including 1 temperature at the end.")
}



library(dplyr, quietly = TRUE, warn.conflicts = FALSE)






### Vector of numbers corresponding to columns energy, Polyak, s_trx2


# Joint data frame to store the individual binning curves
df_joint<-data.frame(matrix(ncol = 6, nrow = 0))
colnames(df_joint)<-c("x", "y", "y_min", "y_max", "T", "Observable")



### Read data from all temperatures

library(foreach)

printf("\n\nStart reading data.\n\n")


for(i_data in 1:ndata) {

    temperature <- args[i_data+5]

    file_pattern_T <- p(file_pattern, "T", temperature, "M", flux, "D", D, "_")

    printf("Passing file pattern %s in folder %s to ReadOutFileToDataTable function.\n\n", file_pattern_T, target_directory)


    # Here we do not include any trajectory cuts as the script deletes simply the first nskip trajectories. This is useful here because we may have different MC chains forked from other runs at different MC times, so that specifying simply a global minimum itraj is not very useful
    all_data_T <- ReadOutFileToDataTable(target_directory, file_pattern_T, 1, -1, plot_labels, parallel)


    printf("\nApplying trajetory cuts from r_binning_joint script:\n\n")
    printf("\t%s: %f\n", "Max. |P| before cut", max(all_data_T[,get("Polyak")]))

    all_data_T <- all_data_T[itraj >= all_data_T[1, itraj]+nskip]

    printf("\tDeleting first %d rows.\n", nskip)
    printf("\t%s:  %f\n\n", "Max. |P| after cut", max(all_data_T[,get("Polyak")]))


    if(nrow(all_data_T)<1) {
        stop("all_data_T is empty.")
    }
        

    ndata <- nrow(all_data_T)
	
    # We now create the histogram data for all measurements categories in plot_loop and store the result in df_joint
	for(plot_loop in plot_labels){
        printf("\nCreating histogram for %s", plot_loop)

        min_measurement <- min(all_data_T[, get(plot_loop)], na.rm = TRUE) 
        max_measurement <- max(all_data_T[, get(plot_loop)], na.rm = TRUE) 

        if(plot_loop=="Polyak"){
         #   min_measurement <- min_measurement_P 
         #   max_measurement <- max_measurement_P
        }
        
        

        if(grepl('GB', file_pattern) & (plot_loop == "energy")) {
            # Delete wrong energy zero measurements from histogram
            min_measurement <- mean(all_data_T[, get(plot_loop)], na.rm = TRUE)-nSD_GB*sd(all_data_T[, get(plot_loop)], na.rm = TRUE)
            #min_measurement <- 5.8
            #max_measurement <- 7.1

            all_data_T <- all_data_T[plot_loop >= min_measurement]

            printf("\n\nWARNING: Minimum energy set to mean-5sd due to gauged bosonic energy zero measurements\n\n")
        }

        
        bin_width <- (max_measurement - min_measurement) / nbins 

        integer_bins <- c(0:(nbins)) 

        bins <- min_measurement + (max_measurement - min_measurement) / nbins * c(0:(nbins)) 		

        binning_result <- hist(all_data_T[, get(plot_loop)], breaks=bins, include.lowest = T, right=FALSE, plot = FALSE)


        # Should we compute error bars? (eb=1)
        if(eb==1) {
            if(use_acl == 1) {
                hje <- JackknifeHistogramACL(all_data_T[, get(plot_loop)], mult_histjack)
            } else {
                hje <- JackknifeHistogramFixedNBin(all_data_T[, get(plot_loop)], n_jackbin)
            }
        } else {
            # If no errors are computed, then the histogram jackknife errors hje are set to zero
            hje <- rep(0, nbins)
        }


        # - bin_width/2 here takes the left edge of the bin. Then we can use geom_step below
        df_add<-data.frame(x=(binning_result$mids-bin_width/2), y=binning_result$density, y_min=binning_result$density-hje, y_max=binning_result$density+hje, T=args[i_data+5], Observable=plot_loop)	
        
        # Add rightmost point for geom_step
        df_add_right<-data.frame(x=(binning_result$mids[nbins]+bin_width/2), y=binning_result$density[nbins], y_min=binning_result$density[nbins]-hje[nbins], y_max=binning_result$density[nbins]+hje[nbins], T=args[i_data+5], Observable=plot_loop)
        df_add<-rbind(df_add, df_add_right)


        #Test normalization
        if(test_normalisation == 1) {
            s=0
            for(i_bins in 1:nbins) {
                s=s+binning_result$density[i_bins]
            } 
            s=s* bin_width
            printf("Integrated distribution: %f\n",s)
        }    
        
        df_joint<-rbind(df_joint, df_add)
	
	}
	

	remove(all_data_T, df_add)

    printf("\n\n")

}

printf("\n\nLoading plot libraries...  ")

library("ggplot2")

printf("Done.\n\n\n")



#### Print binned result. 

printf("Creating binned histograms:\n")

for(plot_loop in plot_labels){

    printf("\t%s\n", plot_loop)
	
	df_plot <- df_joint[df_joint$Observable==plot_loop, ]

    bin_width <- df_plot$x[2]-df_plot$x[1] # restore bin_width for the current observable
	
	##### Redefinitions for nicer titles and axes
	plot_loop_print<-plot_loop

	if(plot_loop == "Polyak") { plot_loop_print<-"|P|"}
	if(plot_loop == "energy") { plot_loop_print<-expression(E/N^2)}
	if(plot_loop == "s_trx2") { plot_loop_print<-expression(R^2)}

    


    if(eb==1) {
        if(use_acl==1) {
            plot_file <- p("plot_jointbins_", plot_loop, "_", file_pattern, "M", flux, "D", D, "acl", mult_histjack, ".pdf")
        }
        else {
            plot_file <- p("plot_jointbins_", plot_loop, "_", file_pattern, "M", flux, "D", D, "nj", n_jackbin, ".pdf")
        }
    } else {
        plot_file <- p("plot_jointbins_", plot_loop, "_", file_pattern, "M", flux, "D", D, ".pdf")
    }

	plot_title <- p("Binned ", plot_loop_print, " of ", file_pattern, "M", flux, "D", D)
	
	

	hist_plot<-	ggplot(aes(x=x, y=y, color=T), data=df_plot) + 
        geom_step(size=1) + 
        xlab(plot_loop_print) + 
        ylab("frequency") + 
        theme_bw() + 
        theme(legend.position = "top") + 
        scale_x_continuous(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0)) + 
        theme(legend.justification = c(0, 1), legend.position = c(0.01, 0.99), legend.direction="horizontal", legend.box = "horizontal")+expand_limits(y = max(df_plot$y)*1.16) + 
        guides(col = guide_legend(nrow = n_row_templist)) + 
        theme(text = element_text(size=20)) + 
        theme(axis.text=element_text(size=20)) + 
        geom_errorbar(aes(x=df_plot[,"x"] + bin_width/2, ymin=df_plot[,"y_min"], ymax=df_plot[,"y_max"]), width=bin_width)


    ### Possible options

    ## add vertical line to plot (incidating e.g. P_gap)
    # + geom_vline(xintercept = 0.5, linetype="dashed")
	
	suppressMessages(ggsave(plot_file, hist_plot))
	
}


printf("Done.\n\n\n")



