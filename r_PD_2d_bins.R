#!/usr/bin/env Rscript





############# CHANGE THIS #################

##### These are all reasonable things to plot. 
####  c("h_f-h_i", "CG_it", "energy", "Polyak", "s_trx2", "trx2(1)", "trx2(2)", "trx2(3)", "trx2(4)", "trx2(5)", "trx2(6)", "trx2(7)", "trx2(8)", "trx2(9)", "comm2","Myers","accept")


# These are the objects that are read from the files and plotted
plot_labels <- c("Polyak", "s_trx2", "energy", "Myers", "trx2(1)")
#plot_labels <- c("energy", "Polyak", "s_trx2", "Myers", "trx2(1)", "trx2(4)")




# Number of standard deviations to cut off data from below for energy measurement for gauged bosonic runs. 
nSD_GB <- 4

#thin out data to plot and process by this value
thin_out_by = 1


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




######### Start of program output

PrintScriptName("R 2d binned histograms")


printf("Arguments: target_directory file_pattern itraj_start itraj_stop \n\n")


printf("e.g.: prod_broad_output GN24S12T0.543M2.0D9_F1 1 20000 \n\n")

printf("itraj_stop < 1: no upper limit on trajectory number.\n\n")



all_data <- ReadOutFileToDataTable(target_directory, file_pattern, itraj_start, itraj_stop, plot_labels, parallel)


if(length(all_data[, get("itraj")])<1) {
    stop("all_data has zero rows")
}



###### BEGIN Main function define 

PD2dBinned <- function(x_label, y_label, do_nls) {

    printf("\n\nCreating 2d binned histogram for %s vs %s", y_label, x_label)
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

	
	##### Redefinitions for nicer titles and axes
    
    x_label_print = x_label
    y_label_print = y_label

	if(y_label == "Polyak") { y_label_print<-"|P|"}
	if(y_label == "energy") { y_label_print<-expression(E/N^2)}
	if(y_label == "s_trx2") { y_label_print<-expression(R^2)}


	if(x_label == "Polyak") { x_label_print<-"|P|"}
	if(x_label == "energy") { x_label_print<-expression(E/N^2)}
	if(x_label == "s_trx2") { x_label_print<-expression(R^2)}

		
	plot_file <- p("plot_2d_bins_", x_label, "_", y_label, "_", file_pattern, ".pdf")
	plot_title <- p("2d Binned ", y_label, " vs ", x_label, " of ", file_pattern, "     ntraj=", nrow(all_data))

	
	binframe <- data.frame(x=all_data_subset[, get(x_label)], y=all_data_subset[, get(y_label)])

    # Fitting only if do_nls parameter is set to 1
    if(do_nls==1) {	
        
        # Fit a power law + const
        trcFunc <- function(x,a,b,c){c * x^a+b}

        ZfitR <- nls(y~trcFunc(x,a,b,c), data=binframe, start=list(a=2.0,b=0,c=0.3))
        print(summary(ZfitR))
        
        print(summary(ZfitR)$parameters[1,2])
            
        fitby <- p(formatC(signif(coef(ZfitR)[[3]],digits=4), digits=4,format="fg", flag="#"), "(", formatC(summary(ZfitR)$parameters[3,2], digits=1,format="fg", flag="#"), ")","* |P|^", formatC(signif(coef(ZfitR)[[1]],digits=4), digits=4,format="fg", flag="#"), "(", formatC(summary(ZfitR)$parameters[1,2], digits=1,format="fg", flag="#"), ")", " + ", formatC(signif(coef(ZfitR)[[2]],digits=4), digits=4,format="fg", flag="#"), "(", formatC(summary(ZfitR)$parameters[2,2], digits=1,format="fg", flag="#"), ")")
        
        res <- p("Resid. sum. ^2: ", formatC(signif(sum(resid(ZfitR)^2),digits=3), digits=3,format="fg", flag="#"))
        
        
        RMSE <- p("RMSE = ", formatC(signif(sqrt(sum(resid(ZfitR)^2)/nrow(all_data)),digits=2), digits=2,format="fg", flag="#"))


            
        yplace <- (coef(ZfitR)[[2]] + coef(ZfitR)[[3]] * 0.6^(coef(ZfitR)[[1]]))
        yplace_res <- (coef(ZfitR)[[2]] + coef(ZfitR)[[3]] * 0.5^(coef(ZfitR)[[1]]))


        myplot <- qplot(x,y,data=binframe, geom='bin2d', bins=70)+ theme_bw()+xlab(x_label_print)+ylab(y_label_print) +scale_x_continuous(expand = c(0, 0))+scale_y_continuous(expand = c(0, 0))+theme(legend.justification = c(0, 1), legend.position = c(0.01, 0.99), legend.direction="vertical", legend.box = "vertical")+expand_limits(y = max(binframe$y)+(max(binframe$y)-min(binframe$y))*0.02)+theme(text = element_text(size=20))+theme(axis.text=element_text(size=20))+theme(legend.text=element_text(size=20))+stat_function(fun = function(x)coef(ZfitR)[[2]] + coef(ZfitR)[[3]] * (x^coef(ZfitR)[[1]])) + annotate("text", x = 0.3, y = yplace , label = fitby)+ annotate("text", x = 0.2, y = yplace_res , label = RMSE)
        
    } else {
        myplot <- qplot(x,y,data=binframe, geom='bin2d', bins=70)+ theme_bw()+xlab(x_label_print)+ylab(y_label_print) +scale_x_continuous(expand = c(0, 0))+scale_y_continuous(expand = c(0, 0))+theme(legend.justification = c(0, 1), legend.position = c(0.01, 0.99), legend.direction="vertical", legend.box = "vertical")+expand_limits(y = max(binframe$y)+(max(binframe$y)-min(binframe$y))*0.02)+theme(text = element_text(size=20))+theme(axis.text=element_text(size=20))+theme(legend.text=element_text(size=20))
    }
	
    
    ### Additional options:
    
    # Fix the x-limits in case of Polyakov loop
    #+scale_x_continuous(limits=c(-0.0001,0.7),expand = c(0, 0))

	### Version for publication
	#	myplot <- qplot(x,y,data=binframe, geom='bin2d', bins=70)+ theme_bw()+xlab("|P|")+ylab(expression(R^2)) +scale_x_continuous(limits=c(-0.0001,0.7),expand = c(0, 0))+scale_y_continuous(expand = c(0, 0))+theme(legend.justification = c(0, 1), legend.position = c(0.01, 0.99), legend.direction="vertical", legend.box = "vertical")+expand_limits(y = max(binframe$y)+(max(binframe$y)-min(binframe$y))*0.02)+theme(text = element_text(size=20))+theme(axis.text=element_text(size=20))+theme(legend.text=element_text(size=20))

	## Working version without publication tweaks:
	#myplot <- qplot(x,y,data=binframe, geom='bin2d', bins=100)+xlab("|P|")+ylab("s_trx2")+ggtitle(plot_title)+stat_function(fun = function(x)coef(ZfitR)[[2]] + coef(ZfitR)[[3]] * (x^coef(ZfitR)[[1]])) + annotate("text", x = 0.2, y = yplace , label = fitby)+ annotate("text", x = 0.2, y = yplace_res , label = RMSE)+xlim(0.0, 0.7)
	
	
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


printf("\n\nCreating 2d binned histograms\n\n")


PD2dBinned("Polyak", "s_trx2", 1) 


PD2dBinned("Polyak", "energy", 1) 


PD2dBinned("energy", "s_trx2", 0) 


PD2dBinned("Polyak", "Myers", 0) 

