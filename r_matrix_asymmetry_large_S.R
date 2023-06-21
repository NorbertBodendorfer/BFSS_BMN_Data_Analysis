#!/usr/bin/env Rscript


require(graphics)
library(RColorBrewer)
library("gplots")


############# CHANGE THIS #################


# Smallest S to use in linear extrapolation in 1/S
smallest_S <- 24
largest_S <- 72


divide_sum_by <- 9

#options(error=traceback)

############# STOP CHANGING ###############






######### Get command arguments
args = commandArgs(trailingOnly=TRUE)

UGB <- args[1]
file_pattern <- args[2]
N_fit_power <- as.numeric(args[3])




######### Settings

source("r_GeneralSettings.R")


######### Function defines

source("r_GeneralFunctions.R")
source("r_ReadPhaseFilesToDataTable.R")
source("r_ErrorBars.R")

library(geometry) 


######### Start of program output

PrintScriptName("R matrix asymmetry large S")

printf("Arguments: U/G(B)Nxx file_pattern_starting_at_T \n\n")

printf("e.g.: ./r_matrix_asymmetry_large_S_production_lowT.R GN16 T0.3M0.5D9_F \n\n")






file_pattern_read <- p("matrix_sizes_", UGB, "S*", file_pattern, "*.csv")


### Extract T
T <- as.numeric(gsub(".*T(.+)M.*", "\\1", file_pattern))
printf("T extracted as %f\n\n", T)



printf("\n\nLoading read libraries... \n ")
library(data.table)
printf("\n")






printf("Looking for files in current directory matching %s\n\n", file_pattern_read)

# Get list of files to include in data
#files_out <- list.files(pattern=file_pattern_read, full.names=TRUE, recursive=FALSE)



files_to_read <- Sys.glob(file_pattern_read)

printf("Found the following files:\n")
print(files_to_read)


printf("\n\nReading content into database:\n")

library(foreach)

all_data <- foreach(i=1:length(files_to_read), .combine=rbind) %do% {
        f <- as.character(files_to_read[i])
        printf("\tReading %s.\n", f)	
        tryCatch({
            try(fread(f, select = c(1,2,3,4))) # The sep="" command messes up fread, needs to be dropped. 
        }, warning = function(w) {
            warning(w)
        }, error = function(e) {
            stop(e)
        }, finally = {
        })
    }


printf("\nDone reading.\n\n")





S_list <- unique(all_data[ ,get("S")])
printf("S list extracted as\n")
#S_list <- as.numeric(S_list)
print(S_list)




S_ticks <- as.character(S_list)

S_ticks <- c(expression(infinity), S_ticks)


S_list <- as.numeric(S_list)





#Load plot libraries 
printf("\n\nLoading plot library...  ")
invisible(library("ggplot2"))
invisible(library("dplyr"))






# Add Sinv (inverse) column for linear fit

all_data <- as.data.frame(all_data)

all_data$Sinv <- 1 / all_data$S


# divide sum to fit on same plot
for(i in 1:(nrow(all_data)/10)) {
    all_data[i*10, "size"] <- all_data[i*10, "size"] / divide_sum_by
    all_data[i*10, "error"] <- all_data[i*10, "error"] / divide_sum_by
}

#print(all_data)



# vector of linear fits
pred <- c()

# Get vector of standard ggplot colors to give to the linear fits
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(10)




color_define = c("darkblue", "deeppink2", "grey", "red", "chartreuse", "cyan", "yellow", "black", "green3", "darkorange")

cols <- col2hex(color_define)
print(cols)



# linear function whose a,b values are set via args in statfunction
linear <- function(x, u.values){
    X <- cbind(1, x)
    return(X %*% u.values)
}


# create large S extrapolation
for(i in 1:10) {
	data_i <- subset(all_data, all_data[,"matrix"] == i)
    data_i <- subset(data_i, data_i$S>=smallest_S)
    data_i <- subset(data_i, data_i$S<=largest_S)


    num_S_fit <- length(S_list[S_list >= smallest_S])
    #print(num_S_fit)
    #print(data_i)

    printf("\n\n*********\n\n  Linear model with nls for matrix %d\n\n", i)

    linFunction <- function(Sinv,c,cS){c + cS * Sinv  }
    modelLin <- nls(size~linFunction(Sinv,c,cS), data=data_i, weights=1 / data_i$error^2, start=list(c=0,cS=0))
    print(summary(modelLin))


    extrapolated_data <- data.frame(matrix=i, size=summary(modelLin)$coefficients[1,1], error=summary(modelLin)$coefficients[1,2], S=999, Sinv=0)

    a <- summary(modelLin)$coefficients[1,1]
    b <- summary(modelLin)$coefficients[2,1]
    
    all_data <- rbind(all_data, extrapolated_data)

	pred <- c(pred, stat_function(size =1, fun = linear, color=cols[i], args=list(u.values =c(a, b))))
	
}




all_data$c=as.character(all_data$matrix)
#print(all_data)

for(i in 1:(nrow(all_data)/10)) {
    all_data[i*10, "c"] <- p("sum/", divide_sum_by)
}


x_label <- "L"

bar_width = 1/max(S_list)/10
#print(Nc_list)
#print(Nc_list^N_fit_power)

myplot <- ggplot(aes(x=Sinv, y=size, color=c), data=all_data) + 
            scale_x_continuous(expand = c(0, 0), limits = c(-0.005, max(1 / S_list)*1.1), breaks=c(0.0, S_list^(-1)), labels=S_ticks) + 
            theme_bw() + 
            geom_point(size=2, alpha=1.0) +
            xlab(x_label) + 
            ylab(p("Matrix size")) + 
            theme(text = element_text(size=20)) + 
            theme(axis.text=element_text(size=15)) + 
            pred + 
            geom_errorbar(aes(x=jitter(all_data[,"Sinv"]), ymin=(all_data[,"size"]-all_data[,"error"]), ymax=(all_data[,"size"]+all_data[,"error"])), width=bar_width) + 
            geom_vline(xintercept=smallest_S^(-1)*1.05, linetype="dashed", color = "black", size=0.5) + 
            geom_vline(xintercept=largest_S^(-1)*0.95, linetype="dashed", color = "black", size=0.5) + 
            theme(aspect.ratio=0.6) + 
            #labs(color=expression(paste("Tr X"^"2"))) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10)) +
            theme(legend.position="right") +
            theme(legend.title=element_blank()) +
            scale_color_manual(labels = c(expression(Tr~(X[1]^2)), expression(Tr~(X[2]^2)), expression(Tr~(X[3]^2)), expression(Tr~(X[4]^2)), expression(Tr~(X[5]^2)), expression(Tr~(X[6]^2)), expression(Tr~(X[7]^2)), expression(Tr~(X[8]^2)), expression(Tr~(X[9]^2)), expression(R^2 / 9)), values = cols)

            #+scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

### Additional options
#+geom_smooth(method=lm, fullrange=TRUE) # adds linear regression with error ranges from multiple data points

plot_file <- p("large_S_", UGB, file_pattern, ".pdf")

suppressMessages(ggsave(myplot, file=plot_file))



