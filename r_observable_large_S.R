#!/usr/bin/env Rscript





############# CHANGE THIS #################


# Smallest S to use in linear extrapolation in 1/S
smallest_S <- 33

largest_S <- 120


# Still a bug here with the colouring. The order here chooses the colouring. Needs to be P, E, s_trx2 for correct colouring at the moment. Colour passed right now to geom_point as observable. color as the c column doesn't work at the moment. 
plot_labels <- c("Polyak", "energy", "s_trx2")

divide_sum_by <- 10

quadratic_fit = TRUE

linear_fit_nls = TRUE

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

PrintScriptName("R observable large S extrapolation")

printf("Arguments: U/G(B)Nxx file_pattern_starting_at_T \n\n")

printf("e.g.: ./r_observable_large_S.R GN16 T0.3M0.5D9_F \n\n")






file_pattern_read <- p("processed/observables_", UGB, "S*", file_pattern, "*.csv")


### Extract T
T <- as.numeric(gsub(".*T(.+)M.*", "\\1", file_pattern))
printf("T extracted as %f\n\n", T)

### Extract M
M <- as.numeric(gsub(".*M(.+)D.*", "\\1", file_pattern))
printf("M extracted as %f\n\n", M)

### Extract N
N <- as.numeric(gsub(".*N(.+)", "\\1", UGB))
printf("N extracted as %f\n\n", N)



printf("\n\nLoading read libraries... \n ")
library(data.table)
printf("\n")


if(T >= 0.3) {
    if(N == 16) {
        smallest_S <- 48
    }
    else {
        smallest_S <- 48
    }
} else {
    smallest_S <- 36
}




printf("Looking for files in current directory matching %s\n\n", file_pattern_read)

# Get list of files to include in data
#files_out <- list.files(pattern=file_pattern_read, full.names=TRUE, recursive=FALSE)



files_to_read <- Sys.glob(file_pattern_read)

# Remove already extrapolated data from previous run
files_to_read <- files_to_read[!grepl( "Soo", files_to_read, fixed = TRUE)]


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


S_list <- unique(all_data[ ,get("S")])
printf("S list extracted as\n")
#S_list <- as.numeric(S_list)
print(S_list)

printf("\n\n")


S_ticks <- as.character(S_list)


S_ticks <- c(expression(infinity), S_ticks)
#printf("S list for fitting:\n")
#print(S_ticks)
S_list <- as.numeric(S_list)




error_data = data.frame(matrix(nrow=0, ncol=4))

colnames(error_data) <- c("observable", "method", "value", "error")



# Load theory values
E_theory <- read.csv("theory_values_E.csv")

E_theory <- E_theory[(E_theory[,"N"] == N),]
E_theory <- E_theory[(E_theory[,"M"] == M),]
E_theory <- E_theory[(E_theory[,"T"] == T),]


new_error_entry <- data.frame(observable="energy", method="Gravity", value=E_theory[1,"value"], error=E_theory[1,"error"])

error_data <- rbind(error_data, new_error_entry)

new_error_entry <- data.frame(observable="Polyak", method="Gravity", value=NA, error=NA)
error_data <- rbind(error_data, new_error_entry)

new_error_entry <- data.frame(observable="s_trx2", method="Gravity", value=NA, error=NA)
error_data <- rbind(error_data, new_error_entry)


theory_value = FALSE
if (dim(E_theory)[1] == 0) {
    printf("No theory value found for above N,M,T.\n\n")
} else {
    printf("Theory value found for above N,M,T.\n")
    print(E_theory)
    printf("\n")
    theory_value = TRUE
}



#Load plot libraries 
printf("\n\nLoading plot library...  ")
invisible(library("ggplot2"))
invisible(library("dplyr"))






# Add Sinv (inverse) column for linear fit

all_data <- as.data.frame(all_data)

all_data$Sinv <- 1 / all_data$S


print(all_data)



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






# create large S extrapolation
for(i in 1:n_obs) {


    printf("===============================================================================\n\n")
    printf("Fitting for %s\n", plot_labels[i])


    data_i_full <- subset(all_data, all_data[,"observable"] == plot_labels[i])


    # data_i is reduced by those S that are not within the fitting range
	data_i <- subset(data_i_full, data_i_full$S>=smallest_S)
    data_i <- subset(data_i, data_i$S<=largest_S)





    ### Manual weighted linear regression with error being the error da induced on a via error propagation or the errors in the datapoints

    printf("\n\n*********\n\n  Manual weighted linear regression for %s\n\n", plot_labels[i])

    num_S_fit <- length(S_list[S_list >= smallest_S])
    #print(num_S_fit)
    #print(data_i)

	# Calculate regression coefficients
    if(num_S_fit == 2) {
        b <- cov(data_i$value, data_i$Sinv)/var(data_i$Sinv)
        a <- mean(data_i$value)-b*mean(data_i$Sinv)
        #print(a)
        #print(b)
        xb <- mean(data_i$Sinv)
        dady <- 1/num_S_fit-xb*(data_i$Sinv-xb)/dot(data_i$Sinv-xb, data_i$Sinv-xb)
        da <- (dot(dady^2, (data_i$error)^2))^0.5
    } else if(num_S_fit > 2) {
        # We do weighted linear regression for 3 or more data points
        sum_inv_error2 <- sum(data_i$error^-2)
        gen_cov <- dot(data_i$Sinv/data_i$error, data_i$value/data_i$error)-sum(data_i$Sinv/data_i$error^2)*sum(data_i$value/data_i$error^2)/sum_inv_error2
        gen_var <- dot(data_i$Sinv/data_i$error, data_i$Sinv/data_i$error)-sum(data_i$Sinv/data_i$error^2)^2/sum_inv_error2
        b <- gen_cov / gen_var
        a <- sum(data_i$value/data_i$error^2-b*data_i$Sinv/data_i$error^2)/sum_inv_error2
        #print(a)
        #print(b)
        dbdy <- (data_i$Sinv-sum(data_i$Sinv/data_i$error^2)/sum_inv_error2) / gen_var / data_i$error^2
        dady <- 1/sum(data_i$Sinv/data_i$error^2)*(data_i$Sinv/data_i$error^2-sum(data_i$Sinv^2/data_i$error^2)*dbdy)
        da <- (dot(dady^2, data_i$error^2))^0.5
    } else {
        stop("Less than two data poitns for linear fit.")
    }

    printf("\n%s\na: %f,     da: %f,     b: %f\n\n", plot_labels[i], a, da, b)

    if(plot_labels[i] == "energy") {
        a_S_E <- a
        b_S_E <- b
    }

    extrapolated_data <- data.frame(observable=plot_labels[i], value=a, error=da, S=999, Sinv=0, N=N, T=T, M=M)
    
    all_data <- rbind(all_data, extrapolated_data)

	pred <- c(pred, stat_function(fun = linear, color=cols[i], linetype="dashed", args=list(u.values =c(a, b))))

    new_error_entry <- data.frame(observable=plot_labels[i], method="Linear da", value=a, error=da)

    error_data <- rbind(error_data, new_error_entry)



    if( linear_fit_nls ) {
        printf("\n\n*********\n\n  Linear model with nls for %s\n\n", plot_labels[i])

        linFunction <- function(Sinv,c,cS){c + cS * Sinv  }

        modelLin <- nls(value~linFunction(Sinv,c,cS), data=data_i, weights=1 / data_i$error^2, start=list(c=0,cS=0))
        print(summary(modelLin))

        new_error_entry <- data.frame(observable=plot_labels[i], method="Linear nls", value=summary(modelLin)$coefficients[1,1], error=summary(modelLin)$coefficients[1,2])

        error_data <- rbind(error_data, new_error_entry)

    }	

    # Fit quadratic function with all data


    if( quadratic_fit ) {
        printf("\n\n*********\n\n  Quadratic model with nls for %s\n\n", plot_labels[i])

        quadFunction <- function(Sinv,c,cS,cSS){c + cS * Sinv + cSS * Sinv * Sinv }

        data_i_quadratic <- subset(data_i_full, data_i_full$S<=largest_S)

        modelQuad <- nls(value~quadFunction(Sinv,c,cS,cSS), data=data_i_quadratic, weights=1 / data_i_quadratic$error^2, start=list(c=0.2,cS=-10,cSS=-30))
        print(summary(modelQuad))


        new_error_entry <- data.frame(observable=plot_labels[i], method="Quadratic", value=summary(modelQuad)$coefficients[1,1], error=summary(modelQuad)$coefficients[1,2])

        error_data <- rbind(error_data, new_error_entry)

        extrapolated_data <- data.frame(observable=plot_labels[i], value=summary(modelQuad)$coefficients[1,1], error=summary(modelQuad)$coefficients[1,2], S=999, Sinv=0, N=N, T=T, M=M)
        all_data <- rbind(all_data, extrapolated_data)

        pred <- c(pred, stat_function(fun = quadFunction, linetype="solid", color=cols[i], args=summary(modelQuad)$coefficients[,1]))
    }	

    

}

printf("===============================================================================\n\n")
printf("Done fitting.\n\n")

printf("\nError data:\n")

print(error_data)

printf("\n\n")

all_data$c=0

for(i in 1:n_obs) {
    all_data[all_data[, "observable"] == plot_labels[i], "c"] <- cols[i]
}



x_label <- "L"

bar_width = 1/max(S_list)/10

subheading <- sprintf("a_S_E = %.2f    b_S_E = %.2f", a_S_E, b_S_E)


data_plot <- ggplot(aes(x=Sinv, y=value, color=observable), data=all_data) +
            scale_x_continuous(expand = c(0, 0), limits = c(-2 * bar_width, max(1 / S_list)*1.1), breaks=c(0.0, S_list^(-1)), labels=S_ticks, minor_breaks = NULL) + 
            scale_y_continuous(expand = c(0, 0), limits = c(-0.1, 0.7), breaks = ((-0.1)+(0:8)*0.1) ) +
            theme_bw() + 
            geom_point(size=2, alpha=1.0) + 
            xlab(x_label) + 
            ylab(p("observable")) +
            theme(text = element_text(size=20)) +
            theme(axis.text=element_text(size=10)) + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            scale_color_manual(labels = c("|P|", expression(E / N^2), expression(R^2 / 10)), values = cols) +
          #  scale_fill_identity() +
            pred +
            #geom_errorbar(aes(x=jitter(all_data[,"Sinv"]), ymin=(all_data[,"value"]-all_data[,"error"]), ymax=(all_data[,"value"]+all_data[,"error"])), width=bar_width) + 
            geom_errorbar(aes(x=all_data[,"Sinv"], ymin=(all_data[,"value"]-all_data[,"error"]), ymax=(all_data[,"value"]+all_data[,"error"])), width=bar_width) + 
            geom_vline(xintercept=smallest_S^(-1)*1.05, linetype="dashed", color = "black", size=0.5) +
            geom_vline(xintercept=largest_S^(-1)*0.95, linetype="dashed", color = "black", size=0.5) +
            theme(aspect.ratio=1) + 
            theme(legend.position="bottom") +
          #  labs(subtitle = subheading) +
            {if(theory_value) geom_point(aes(x=0, y=E_theory[1,"value"]), colour="darkgreen", shape=3) }+
            {if(theory_value)  geom_errorbar(aes(x=0, ymin=(E_theory[1,"value"]-E_theory[1,"error"]), ymax=(E_theory[1,"value"]+E_theory[1,"error"])), width=bar_width, colour="darkgreen") }+
            labs(color=expression(paste(""))) 
            #+geom_smooth(method="nls", formula=y~a+b*x+c*x^2, start=c(a=0.0, b=0.0, c=0.0), se=FALSE)


plot_file <- p("observable_large_S_", UGB, file_pattern, ".pdf")

suppressMessages(ggsave(data_plot, file=plot_file))



data_file <- file_pattern_read <- p("processed/observables_", UGB, "Soo", file_pattern, ".csv")

data_file_error <- file_pattern_read <- p("processed/observables_errors_", UGB, "Soo", file_pattern, ".csv")


large_S_data <- subset(all_data, all_data[,"Sinv"] == 0) 

large_S_data$Sinv <- NULL

large_S_data$c <- NULL





# Readjust strx2 to original value

large_S_data <- as.data.table(large_S_data)

large_S_data[observable %in% "s_trx2"][, 2] <- large_S_data[observable %in% "s_trx2"][, 2] * divide_sum_by
large_S_data[observable %in% "s_trx2"][, 3] <- large_S_data[observable %in% "s_trx2"][, 3] * divide_sum_by

printf("Large S data:\n")

print(large_S_data)


write.table(large_S_data, file = data_file, sep = ",", qmethod = "double", row.names=FALSE)


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


plot_file_error <- p("observable_large_S_", UGB, file_pattern, "_errors.pdf")



pdf(plot_file_error)
grid.arrange(error_plot_P, error_plot_R, error_plot_E, ncol=1, nrow=3)
dev.off()

### Additional options
#+geom_smooth(method=lm, fullrange=TRUE) # adds linear regression with error ranges from multiple data points


#suppressMessages(ggsave(error_plot_E, file=plot_file_error))

