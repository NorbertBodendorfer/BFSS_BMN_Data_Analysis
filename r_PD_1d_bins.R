#!/usr/bin/env Rscript





############# CHANGE THIS #################

##### These are all reasonable things to plot against itraj, which is done in this script. 
####  c("h_f-h_i", "CG_it", "energy", "Polyak", "s_trx2", "trx2(1)", "trx2(2)", "trx2(3)", "trx2(4)", "trx2(5)", "trx2(6)", "trx2(7)", "trx2(8)", "trx2(9)", "comm2","Myers","accept")


# These are the objects that are read from the files and plotted
plot_labels <- c("Polyak", "s_trx2", "energy", "Myers", "trx2(1)")
#plot_labels <- c("energy", "Polyak", "s_trx2", "Myers", "trx2(1)", "trx2(4)")




# Number of standard deviations to cut off data from below for energy measurement for gauged bosonic runs. 
nSD_GB <- 4

nbins <- 30


# Run script in parallel (0 or 1=parallel)
parallel = 0

############# STOP CHANGING ###############






######### Get command arguments
args = commandArgs(trailingOnly=TRUE)

target_directory <- args[1]
file_pattern <- args[2]
itraj_start <- as.numeric(args[3])
itraj_stop <- as.numeric(args[4])
P_min <- as.numeric(args[5])
P_gap <- as.numeric(args[6])
P_max_plot <- as.numeric(args[7])





######### Settings

source("r_GeneralSettings.R")


######### Function defines

source("r_GeneralFunctions.R")
source("r_ReadOutFilesToDataTable.R")




######### Start of program output

PrintScriptName("R 1d binned histograms")


printf("Arguments: target_directory file_pattern itraj_start itraj_stop P_min P_gap P_max_plot \n\n")


printf("e.g.: ./r_PD_1d_bins.R prod_broad_output GN24S12T0.543M2.0D9_F1 1 20000 0 0.5 0.7\n\n")

printf("itraj_stop < 1: no upper limit on trajectory number.\n\n")



all_data <- ReadOutFileToDataTable(target_directory, file_pattern, itraj_start, itraj_stop, plot_labels, parallel)


if(length(all_data[, get("itraj")])<1) {
    stop("all_data has zero rows")
}



###### BEGIN Main function define 

PD1dBinned <- function(x_label, y_label, do_nls) {

    printf("\n\nCreating 1d binned histogram for %s vs %s", y_label, x_label)
    if(do_nls==1) {
        printf(" with fitting.\n")
    } else {
        printf(" without fitting.\n")
    }

    all_data_subset <- all_data

	if(grepl('GB', file_pattern) & ((y_label == "energy") || (y_label == "energy") )) {
		# Delete extreme energy zero measurements from histogram
		min_measurement <- mean(all_data[, get("energy")], na.rm = TRUE)-nSD_GB*sd(all_data[, get("energy")], na.rm = TRUE)
        all_data_subset <- all_data_subset[energy >= min_measurement]
		printf("\nWARNING: Minimum energy set to mean - %d sd due to gauged bosonic energy zero measurements\n", nSD_GB)
	}

    if(x_label=="Polyak") {
        all_data_subset <- all_data_subset[Polyak <= P_max_plot]
    }

	
	##### Redefinitions for nicer titles and axes
    
    x_label_print = x_label
    y_label_print = y_label

	if(y_label == "Polyak") { y_label_print<-"|P|"}
	if(y_label == "energy") { y_label_print<-expression(E/N^2)}
	if(y_label == "s_trx2") { y_label_print<-expression(R^2)}


	if(x_label == "Polyak") { x_label_print<-"|P|"}
	if(x_label == "energy") { x_label_print<-expression(E/N^2)}
	if(x_label == "s_trx2") { x_label_print<-expression(R^2)}

		
	plot_file <- p("plot_1d_bins_", x_label, "_", y_label, "_", file_pattern, ".pdf")
	plot_title <- p("1d Binned ", y_label, " vs ", x_label, " of ", file_pattern, "     ntraj=", nrow(all_data))

	bin_table <- all_data_subset[, .(get(x_label), get(y_label))]

    x_min <- min(all_data_subset[, get(x_label)])
    x_max <- max(all_data_subset[, get(x_label)])


  
    x_diff <- x_max-x_min
    
    x_data <- c(x_min+(0:nbins)*x_diff/nbins)
    y_data <- replicate(nbins+1, 0.0)


    setkey(bin_table, V1)


    for(i in 1:nbins) {
        if(i < nbins ) {
            y_data[i] <- mean(bin_table[ V1 >= x_data[i] & V1 < x_data[i+1]][, V2])
        } else { # This is necessary to catch cases where there is only 1 data point in the right-most bin, which is otherwise missed as it is exactly the right boundary
            y_data[i] <- mean(bin_table[ V1 >= x_data[i] & V1 <= x_data[i+1]][, V2])
        }
        
        if(is.na(y_data[i])) { y_data[i] <- 0.0 }
    }

    # Add right-most point for geom_step
    y_data[nbins+1] <- y_data[nbins]



    plot_frame <- as.data.frame(data.frame(x=x_data, y=y_data))

    

    # Fitting only if do_nls parameter is set to 1
    if(do_nls==1) {	
        
        # Fit a power law + const
        trcFunc <- function(x,a,b,c){c * x^a+b}


        # We fit against fit_frame which takes as x-value the middle point of the bin
        fit_frame <- plot_frame

        if(x_label == "Polyak") {
            fit_frame <- subset(fit_frame, fit_frame$x <= P_gap)
            fit_frame <- subset(fit_frame, fit_frame$x >= P_min)
        }

        fit_frame$x <- fit_frame$x+x_diff/nbins*0.5

        # Last row was just for geom_step
        fit_frame <- fit_frame[-(nbins+1), ]

        #print(fit_frame)
      
        # Fit with x_label Polyakov only for P<=P_gap, where partial deconfinement should hold
        



        ZfitR <- nls(y~trcFunc(x,a,b,c), data=fit_frame, start=list(a=2.0,b=0.5,c=0.6))

        
        print(summary(ZfitR))
                    
        fitby <- p(formatC(signif(coef(ZfitR)[[3]],digits=4), digits=4,format="fg", flag="#"), "(", formatC(summary(ZfitR)$parameters[3,2], digits=1,format="fg", flag="#"), ")","* |P|^", formatC(signif(coef(ZfitR)[[1]],digits=4), digits=4,format="fg", flag="#"), "(", formatC(summary(ZfitR)$parameters[1,2], digits=1,format="fg", flag="#"), ")", " + ", formatC(signif(coef(ZfitR)[[2]],digits=4), digits=4,format="fg", flag="#"), "(", formatC(summary(ZfitR)$parameters[2,2], digits=1,format="fg", flag="#"), ")")
        
        res <- p("Resid. sum. ^2: ", formatC(signif(sum(resid(ZfitR)^2),digits=3), digits=3,format="fg", flag="#"))
        
        
        RMSE <- p("RMSE = ", formatC(signif(sqrt(sum(resid(ZfitR)^2)/nrow(all_data)),digits=2), digits=2,format="fg", flag="#"))


            
        yplace <- (coef(ZfitR)[[2]] + coef(ZfitR)[[3]] * 0.6^(coef(ZfitR)[[1]]))
        yplace_res <- (coef(ZfitR)[[2]] + coef(ZfitR)[[3]] * 0.5^(coef(ZfitR)[[1]]))


        myplot <- ggplot(aes(x=x, y=y), data=plot_frame) +
            geom_step(size=1, colour="blue") + 
            ylab(y_label_print) + 
            xlab(x_label_print) + 
        #    ggtitle(plot_title)  + 
            theme_bw() + 
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme(text = element_text(size=20)) + 
            theme(axis.text=element_text(size=15)) +
            stat_function(fun = function(x)coef(ZfitR)[[2]] + coef(ZfitR)[[3]] * (x^coef(ZfitR)[[1]])) + annotate("text", x = 0.25, y = yplace , label = fitby, size=5)
    } else {
        myplot <- ggplot(aes(x=x, y=y), data=plot_frame) +
            geom_step(size=1, colour="blue")  + 
            ylab(y_label_print)+xlab(x_label_print)+ggtitle(plot_title)  + 
            theme_bw() +scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0))
        #+ stat_function(fun = function(x)((x-fitparams[plot_loop, 2])/fitparams[plot_loop, 3])^(1/fitparams[plot_loop, 1]))+ylim(0, 1.0)
    }
	

	
	### Fit with fixed a
	#if(0==1) {
	#	ZfitR <- nls(y~trcFunc(x,2.5,b,c), data=binframe, start=list(b=2.2,c=0.3))
	#	fitby <- p("Fit (fixed exponent): R^2 = ", formatC(signif(coef(ZfitR)[[2]],digits=3), digits=3,format="fg", flag="#"), "* |P|^2.5 + ", formatC(signif(coef(ZfitR)[[1]],digits=3), digits=3,format="fg", flag="#"))
	#		res <- p("Resid. sum. ^2: ", formatC(signif(sum(resid(ZfitR)^2),digits=3), digits=3,format="fg", flag="#"))
    #
	#		myplot <- qplot(x,y,data=binframe, geom='bin2d', bins=100)+xlab("|P|")+ylab("s_trx2")+ggtitle(plot_title)+stat_function(fun = function(x)coef(ZfitR)[[1]] + coef(ZfitR)[[2]] * (x^2.5)) + annotate("text", x = 0.3, y = yplace , label = fitby)+ annotate("text", x = 0.2, y = yplace_res , label = res)
	#}
	
	suppressMessages(ggsave(myplot, file=plot_file))

    

}


###### END Main function define 





#### Print double binned result. 



printf("\n\nLoading plot library.\n\n")

library("ggplot2", quietly = TRUE, warn.conflicts = FALSE)

library(gplots, quietly = TRUE, warn.conflicts = FALSE)


printf("\n\nCreating 1d binned histograms\n\n")


PD1dBinned("Polyak", "s_trx2", 1) 



PD1dBinned("Polyak", "energy", 1) 


stop("Manual abort")


PD1dBinned("Polyak", "s_trx2", 0) 


PD1dBinned("Polyak", "energy", 0) 


PD1dBinned("energy", "s_trx2", 0) 


PD1dBinned("Polyak", "Myers", 0) 

