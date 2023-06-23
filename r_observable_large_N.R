#!/usr/bin/env Rscript





############# CHANGE THIS #################


# Smallest S to use in linear extrapolation in 1/S
smallest_N <- 16

plot_labels <- c("Polyak", "energy", "s_trx2")

smallest_N_quadratic <- 10


divide_sum_by <- 10


linear_fit_nls <- TRUE
quadratic_fit <- TRUE



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
source("r_ErrorBars.R")

library(geometry) 
library(MASS) 


######### Start of program output

PrintScriptName("R observable large N extrapolation")

printf("Arguments: U/G(B) file_pattern_starting_at_S \n\n")

printf("e.g.: ./r_observable_large_N.R G S30T0.3M0.5D9_F \n\n")






file_pattern_read <- p("processed/observables_", UGB, "N*", file_pattern, "*.csv")


### Extract T
T <- as.numeric(gsub(".*T(.+)M.*", "\\1", file_pattern))
printf("T extracted as %f\n\n", T)

### Extract M
M <- as.numeric(gsub(".*M(.+)D.*", "\\1", file_pattern))
printf("M extracted as %f\n\n", M)

### Extract S
S <- as.numeric(gsub(".*S(.+)T.*", "\\1", file_pattern))
printf("S extracted as %f\n\n", S)


printf("\n\nLoading read libraries... \n ")
library(data.table)
printf("\n")






printf("Looking for files in current directory matching %s\n\n", file_pattern_read)

# Get list of files to include in data
#files_out <- list.files(pattern=file_pattern_read, full.names=TRUE, recursive=FALSE)



files_to_read <- Sys.glob(file_pattern_read)

# Remove already extrapolated data from previous run
files_to_read <- files_to_read[!grepl( "Noo", files_to_read, fixed = TRUE)]


printf("Found the following files:\n")
print(files_to_read)


printf("\n\nReading content into database:\n")

library(foreach)

all_data <- foreach(i=1:length(files_to_read), .combine=rbind) %do% {
        f <- as.character(files_to_read[i])
        printf("\tReading %s.\n", f)	
        tryCatch({
            try(fread(f, select = c(1,2,3,4,5,6,7))) # The sep="" command messes up fread, needs to be dropped. 
        }, warning = function(w) {
            warning(w)
        }, error = function(e) {
            stop(e)
        }, finally = {
        })
    }

#print(all_data)

printf("\nDone reading.\n\n")



# Remove all observables not contained in plot_labels
all_data <- all_data[observable %in% plot_labels]

# divide s_trx2 by divide_sum_by to fit into plot
all_data[observable %in% "s_trx2"][, 2] <- all_data[observable %in% "s_trx2"][, 2] / divide_sum_by
all_data[observable %in% "s_trx2"][, 3] <- all_data[observable %in% "s_trx2"][, 3] / divide_sum_by

n_obs <- length(plot_labels)

all_data <- subset(all_data, all_data$N>=smallest_N_quadratic)



print(all_data)

N_list <- unique(all_data[ ,get("N")])
printf("N list extracted as\n")
#N_list <- as.numeric(N_list)
print(N_list)




N_ticks <- as.character(N_list)

N_ticks <- c(expression(infinity), N_ticks)
#printf("S list for fitting:\n")
#print(N_ticks)
N_list <- as.numeric(N_list)




#Load plot libraries 
printf("\n\nLoading plot library...  ")
invisible(library("ggplot2"))
invisible(library("dplyr"))






# Add Ninv (inverse) column for linear fit

all_data <- as.data.frame(all_data)

all_data$Ninv <- 1 / all_data$N

#all_data$var <- all_data$error^2

all_data[all_data[, "observable"] == "s_trx2",]$Ninv <- (all_data[all_data[, "observable"] == "s_trx2",]$Ninv)^2

all_data[all_data[, "observable"] == "energy",]$Ninv <- (all_data[all_data[, "observable"] == "energy",]$Ninv)^2


# vector of linear fits
pred <- c()

# Get vector of standard ggplot colors to give to the linear fits
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(n_obs)



# linear function whose a,b values are set via args in statfunction
linear <- function(x, u.values){
    X <- cbind(1, x)
    return(X %*% u.values)
}



error_data = data.frame(matrix(nrow=0, ncol=7))

colnames(error_data) <- c("observable", "method", "value", "error", "S", "T", "M")

new_error_entry <- data.frame(observable="energy", method="Gravity cont.", value=-5.77*T^0.4 -3.5 * T^(11/5) , error=T^(11/5)*2, S=999, T=0.3, M=M)
error_data <- rbind(error_data, new_error_entry)

new_error_entry <- data.frame(observable="s_trx2", method="Gravity cont.", value=NA, error=NA, S=999, T=0.3, M=M)
error_data <- rbind(error_data, new_error_entry)

new_error_entry <- data.frame(observable="Polyak", method="Gravity cont.", value=NA, error=NA, S=999, T=0.3, M=M)
error_data <- rbind(error_data, new_error_entry)

printf("=========================================\n")
printf(" Ansatz: E(N) = c + cN / N^2 + cNN / N^4  \n")    
printf("=========================================\n\n")



# create large S extrapolation
for(i in 1:n_obs) {

    printf("===============================================================================\n\n")
    printf("Fitting for %s\n", plot_labels[i])

    data_i_full <- subset(all_data, all_data[,"observable"] == plot_labels[i])


    # data_i is reduced by those S that are not within the fitting range
	data_i <- subset(data_i_full, data_i_full$N>=smallest_N)




    ### Manual weighted linear regression with error being the error da induced on a via error propagation or the errors in the datapoints

    printf("\n\n*********\n\n  Manual weighted linear regression for %s\n\n", plot_labels[i])


    num_N_fit <- length(N_list[N_list >= smallest_N])


	# Calculate regression coefficients
    if(num_N_fit == 2) {
        b <- cov(data_i$value, data_i$Ninv)/var(data_i$Ninv)
        a <- mean(data_i$value)-b*mean(data_i$Ninv)
        #print(a)
        #print(b)
        xb <- mean(data_i$Ninv)
        dady <- 1/num_N_fit-xb*(data_i$Ninv-xb)/dot(data_i$Ninv-xb, data_i$Ninv-xb)
        da <- (dot(dady^2, (data_i$error)^2))^0.5
    } else if(num_N_fit > 2) {
        # We do weighted linear regression for 3 or more data points
        sum_inv_error2 <- sum(data_i$error^-2)
        gen_cov <- dot(data_i$Ninv/data_i$error, data_i$value/data_i$error)-sum(data_i$Ninv/data_i$error^2)*sum(data_i$value/data_i$error^2)/sum_inv_error2
        gen_var <- dot(data_i$Ninv/data_i$error, data_i$Ninv/data_i$error)-sum(data_i$Ninv/data_i$error^2)^2/sum_inv_error2
        b <- gen_cov / gen_var
        a <- sum(data_i$value/data_i$error^2-b*data_i$Ninv/data_i$error^2)/sum_inv_error2
        #print(a)
        #print(b)
        dbdy <- (data_i$Ninv-sum(data_i$Ninv/data_i$error^2)/sum_inv_error2) / gen_var / data_i$error^2
        dady <- 1/sum(data_i$Ninv/data_i$error^2)*(data_i$Ninv/data_i$error^2-sum(data_i$Ninv^2/data_i$error^2)*dbdy)
        da <- (dot(dady^2, data_i$error^2))^0.5
      #  db <- (/)^0.5    

    } else {
        stop("Less than two data poitns for linear fit.")
    }


    printf("\n%s\na: %f,     b: %f\n\n", plot_labels[i], a, b)

    if(plot_labels[i] == "energy") {
        a_N_E <- a
        b_N_E <- b
    }

    extrapolated_data <- data.frame(observable=plot_labels[i], value=a, error=da, N=999, Ninv=0, S=S, T=T, M=M)
    
    all_data <- rbind(all_data, extrapolated_data)

	pred <- c(pred, stat_function(fun = linear, color=cols[i], linetype="dashed", args=list(u.values =c(a, b))))


    if( linear_fit_nls ) {
        printf("\n\n*********\n\n  Linear model with nls for %s\n\n", plot_labels[i])

        linFunction <- function(Ninv,a,b){a + b * Ninv  }

        modelLin <- nls(value~linFunction(Ninv,a,b), data=data_i, weights=1 / data_i$error^2, start=list(a=0,b=0), control = nls.control(warnOnly = TRUE, maxiter=5000))
        print(summary(modelLin))


        new_error_entry <- data.frame(observable=plot_labels[i], method="Linear nls", value=summary(modelLin)$coefficients[2,1], error=summary(modelLin)$coefficients[2,2], S=S, T=T, M=M)


        error_data <- rbind(error_data, new_error_entry)

    }	

    # Fit quadratic function with all data


    if( quadratic_fit ) {
        printf("\n\n*********\n\n  Quadratic model with nls for %s\n\n", plot_labels[i])

        quadFunction <- function(Ninv,a,b,c){a + b * Ninv + c * Ninv * Ninv }

        modelQuad <- nls(value~quadFunction(Ninv,a,b,c), data=data_i_full, weights=1 / data_i_full$error^2, start=list(a=0,b=0,c=0.0), control = nls.control(warnOnly = TRUE))
        print(summary(modelQuad))


        new_error_entry <- data.frame(observable=plot_labels[i], method="Quadratic", value=summary(modelQuad)$coefficients[2,1], error=summary(modelQuad)$coefficients[2,2], S=S, T=T, M=M)

        error_data <- rbind(error_data, new_error_entry)

        extrapolated_data <- data.frame(observable=plot_labels[i], value=summary(modelQuad)$coefficients[1,1], error=summary(modelQuad)$coefficients[1,2], N=999, Ninv=0, S=S, T=T, M=M)

        all_data <- rbind(all_data, extrapolated_data)

        pred <- c(pred, stat_function(fun = quadFunction, linetype="solid", color=cols[i], args=summary(modelQuad)$coefficients[,1]))
    }	


    if(plot_labels[i] == "energy") {
        printf("\n\n*********\n\n  Quadratic model with nls for %s, fixed cN \n\n", plot_labels[i])

        quadFunction <- function(Ninv,c,cNN){c -3.81 * Ninv + cNN * Ninv * Ninv }

        modelQuad <- nls(value~quadFunction(Ninv,c,cNN), data=data_i_full, weights=1 / data_i_full$error^2, start=list(c=0,cNN=0.0), control = nls.control(warnOnly = TRUE))
        print(summary(modelQuad))


    }



	#print(pred)
	
}


printf("===============================================================================\n\n")
printf("Done fitting.\n\n")

printf("\nError data:\n")

print(error_data)

printf("\n\n")



all_data$c=as.character(all_data$observable)


all_data[all_data[, "observable"] == "s_trx2","c"] <- p("s_trx2/", divide_sum_by)



x_label <- "N"

bar_width = 1/max(N_list)^2/10

subheading <- sprintf("a_N = %.2f    b_N = %.2f", a_N_E, b_N_E)


myplot <- ggplot(aes(x=Ninv, y=value, color=c), data=all_data) +
            scale_x_continuous(expand = c(0, 0), limits = c(-2 * bar_width, max(1 / N_list^2)*1.1), breaks=c(0.0, N_list^(-2)), labels=N_ticks, minor_breaks = NULL) + 
          #  scale_y_continuous(expand = c(0, 0), limits = c(0.1, 0.7), breaks = ((0.1)+(0:6)*0.1) ) +
            theme_bw() + 
            geom_point(size=2, alpha=1.0) + 
            xlab(x_label) + 
            ylab(p("observable")) +
            theme(text = element_text(size=20)) +
            theme(axis.text=element_text(size=15)) + 
            scale_color_manual(labels = c("|P|", expression(E / N^2), expression(R^2 / 10)), values = cols) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
         #   theme(legend.spacing.x = unit(1.0, 'cm')) +
            pred +
            #geom_errorbar(aes(x=jitter(all_data[,"Ninv"]), ymin=(all_data[,"value"]-all_data[,"error"]), ymax=(all_data[,"value"]+all_data[,"error"])), width=bar_width) + 
            geom_errorbar(aes(x=all_data[,"Ninv"], ymin=(all_data[,"value"]-all_data[,"error"]), ymax=(all_data[,"value"]+all_data[,"error"])), width=bar_width) + 
            geom_vline(xintercept=smallest_N^(-2)*1.05, linetype="dashed", color = "black", size=0.5) +
            theme(aspect.ratio=1) + 
            theme(legend.position="bottom") +
            labs(subtitle = subheading) +
            labs(color=expression(paste(""))) 

### Additional options
#+geom_smooth(method=lm, fullrange=TRUE) # adds linear regression with error ranges from multiple data points

plot_file <- p("observable_large_N_", UGB, file_pattern, ".pdf")

suppressMessages(ggsave(myplot, file=plot_file))



data_file <- file_pattern_read <- p("processed/observables_", UGB, "Noo", file_pattern, ".csv")

data_file_error <- file_pattern_read <- p("processed/observables_error_", UGB, "Noo", file_pattern, ".csv")


large_N_data <- subset(all_data, all_data[,"Ninv"] == 0) 

large_N_data$Ninv <- NULL

large_N_data$c <- NULL

# Readjust strx2 to original value



large_N_data <- as.data.table(large_N_data)

large_N_data[observable %in% "s_trx2"][, 2] <- large_N_data[observable %in% "s_trx2"][, 2] * divide_sum_by
large_N_data[observable %in% "s_trx2"][, 3] <- large_N_data[observable %in% "s_trx2"][, 3] * divide_sum_by

printf("Large N data:\n")

print(large_N_data)


write.table(large_N_data, file = data_file, sep = ",", qmethod = "double", row.names=FALSE)


write.table(error_data, file = data_file_error, sep = ",", qmethod = "double", row.names=FALSE)






error_data_E <- error_data[error_data[,"observable"] == "energy", ]
error_data_P <- error_data[error_data[,"observable"] == "Polyak", ]
error_data_R <- error_data[error_data[,"observable"] == "s_trx2", ]


error_plot_E <- ggplot(aes(x=method, y=value, color=observable), data=error_data_E) +
            #scale_x_continuous(expand = c(0, 0), limits = c(-2 * bar_width, max(1 / S_list)*1.1), breaks=c(0.0, S_list^(-1)), labels=S_ticks) + 
            #scale_y_continuous(expand = c(0, 0), limits = c(0.1, 0.7), breaks = ((0.1)+(0:6)*0.1) ) +
            theme_bw() + 
            geom_point(size=2, alpha=1.0) + 
            xlab("Method") + 
            ylab(expression(E / N^2)) +
            theme(text = element_text(size=20)) +
            #theme(axis.text=element_text(size=10)) + 
            #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            #scale_color_manual(labels = c(expression(E / N^2), "|P|", expression(R^2 / 10)), values = cols) +
            scale_color_manual(values = cols[2]) +
            geom_errorbar(aes(x=error_data_E[,"method"], ymin=(error_data_E[,"value"]-error_data_E[,"error"]), ymax=(error_data_E[,"value"]+error_data_E[,"error"])), width=0.5)  +
            theme(aspect.ratio=1/4) +
            theme(axis.ticks.x=element_blank()) +
            theme(legend.position = "none")


error_plot_P <- ggplot(aes(x=method, y=value, color=observable), data=error_data_P) +
            #scale_x_continuous(expand = c(0, 0), limits = c(-2 * bar_width, max(1 / S_list)*1.1), breaks=c(0.0, S_list^(-1)), labels=S_ticks) + 
            #scale_y_continuous(expand = c(0, 0), limits = c(0.1, 0.7), breaks = ((0.1)+(0:6)*0.1) ) +
            theme_bw() + 
            geom_point(size=2, alpha=1.0) + 
            #xlab("Method") + 
            ylab("|P|") +
            theme(text = element_text(size=20)) +
            #theme(axis.text=element_text(size=10)) + 
            #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            scale_color_manual(values = cols[1]) +
            geom_errorbar(aes(x=error_data_P[,"method"], ymin=(error_data_P[,"value"]-error_data_P[,"error"]), ymax=(error_data_P[,"value"]+error_data_P[,"error"])), width=0.5)  +
            theme(aspect.ratio=1/4) +
            theme(axis.ticks.x=element_blank()) +
            theme(axis.text.x=element_blank()) +
            theme(axis.title.x = element_blank()) + 
            labs(subtitle = expression(1 / N^(2) ~ correction ~ coefficients)) +
            theme(legend.position = "none")


error_plot_R <- ggplot(aes(x=method, y=value, color=observable), data=error_data_R) +
            #scale_x_continuous(expand = c(0, 0), limits = c(-2 * bar_width, max(1 / S_list)*1.1), breaks=c(0.0, S_list^(-1)), labels=S_ticks) + 
            #scale_y_continuous(expand = c(0, 0), limits = c(0.1, 0.7), breaks = ((0.1)+(0:6)*0.1) ) +
            theme_bw() + 
            geom_point(size=2, alpha=1.0) + 
            #xlab("Method") + 
            ylab(expression(R^2 / 10)) +
            theme(text = element_text(size=20)) +
            scale_color_manual(values = cols[3]) +
            #theme(axis.text=element_text(size=10)) + 
            #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            #scale_color_manual(labels = c(expression(E / N^2), "|P|", expression(R^2 / 10)), values = cols) +
            geom_errorbar(aes(x=error_data_R[,"method"], ymin=(error_data_R[,"value"]-error_data_R[,"error"]), ymax=(error_data_R[,"value"]+error_data_R[,"error"])), width=0.5)  +
            theme(aspect.ratio=1/4) +
            theme(axis.ticks.x=element_blank()) +
            theme(axis.text.x = element_blank()) + 
            theme(axis.title.x = element_blank()) + 
            theme(legend.position = "none")


require(gridExtra)


plot_file_error <- p("observable_large_N_", UGB, file_pattern, "_errors.pdf")



pdf(plot_file_error)
grid.arrange(error_plot_P, error_plot_R, error_plot_E, ncol=1, nrow=3)
dev.off()
