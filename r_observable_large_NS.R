#!/usr/bin/env Rscript





############# CHANGE THIS #################


# Smallest S to use in linear extrapolation in 1/S
smallest_N <- 10

largest_N <- 100


smallest_S <- 24

largest_S <- 192


plot_labels <- c("energy") #, "Polyak", "s_trx2")

divide_sum_by <- 10


set_tolerance <- 1e-10
set_minFactor <- 1e-12
set_maxiter <- 50000


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
library(MASS) 



######### Start of program output

PrintScriptName("R observable simulatenous large N large S extrapolation")

printf("Arguments: U/G(B) file_pattern_starting_at_T \n\n")

printf("e.g.: ./r_observable_large_NS.R G T0.3M0.5D9_F \n\n")




file_pattern_read <- p("observables_", UGB, "*", file_pattern, "*.csv")


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
            try(fread(f, select = c(1,2,3,4,5))) # The sep="" command messes up fread, needs to be dropped. 
        }, warning = function(w) {
            warning(w)
        }, error = function(e) {
            stop(e)
        }, finally = {
        })
    }

print(all_data)

printf("\nDone reading.\n\n")




all_data <- subset(all_data, all_data$S<999)
all_data <- subset(all_data, all_data$N<999)




# Remove all observables not contained in plot_labels
all_data <- all_data[observable %in% plot_labels]

# divide s_trx2 by divide_sum_by to fit into plot
all_data[observable %in% "s_trx2"][, 2] <- all_data[observable %in% "s_trx2"][, 2] / divide_sum_by

n_obs <- length(plot_labels)


print(all_data)

N_list <- unique(all_data[ ,get("N")])
S_list <- unique(all_data[ ,get("S")])
printf("N list extracted as\n")
print(N_list)


printf("S list extracted as\n")
print(S_list)



N_ticks <- as.character(N_list)

S_ticks <- as.character(S_list)



N_ticks <- c(expression(infinity), N_ticks)
S_ticks <- c(expression(infinity), S_ticks)



N_list <- as.numeric(N_list)
Ninv_list <- c(0, N_list^(-2))
Sinv_list <- c(0, S_list^(-1))

print("NInvList")
print(Sinv_list)
print("N_ticks")
print(S_ticks)
#




#Load plot libraries 
printf("\n\nLoading plot library...  ")
invisible(library("ggplot2"))
invisible(library("dplyr"))






# Add Ninv (inverse) column for linear fit

all_data <- as.data.frame(all_data)

all_data$Ninv <- 1 / all_data$N^2

all_data$Sinv <- 1 / all_data$S

all_data$weights <- 1 / all_data$error^2






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


all_data$predicted <- NA




# create large S extrapolation
for(i in 1:n_obs) {
	data_i <- subset(all_data, all_data[,"observable"] == plot_labels[i])
    data_i <- subset(data_i, data_i$N>=smallest_N)
    data_i <- subset(data_i, data_i$N<=largest_N)
    data_i <- subset(data_i, data_i$S>=smallest_S)
    data_i <- subset(data_i, data_i$S<=largest_S)

    #print(data_i)
    print(plot_labels[i])




    printf("==============================================================================\n")
    printf(" Ansatz: E(N,S) = c + cN / N^2 + cS / S + cNN / N^4 + cNS / N^2 S + cSS / S^2 \n")    
    printf("==============================================================================\n\n")


    a1 = - 10
    a2 = 5.8
    fun.BMN_mu05 <- function(T) 1/T^(36/5)*(1.97004*10^(-12)+T*(-6.37306*10^(-11)+T*(8.52492*10^(-10)+T*(-6.91845*10^(-9)+T*(9.70163*10^(-8)+T*(-1.02708*10^(-7)+T*(-0.0000878033+T*(-6.95932*10^(-7)+T*(0.0209683 + T*(-5.04764*10^(-9)+7.41*T))))))))))
    fun.BMN_mu05_all_corrections <- function(T)  ( fun.BMN_mu05(T) + a1 * T^(23/5) + a2 * T^(29/5))

    gravity_prediction <- fun.BMN_mu05_all_corrections(T)

    printf("Expected large N value for mu=0.5, T=%.2f is %.3f\n\n\n", T, gravity_prediction) 

    printf("Datapoints are weighted by their error bars, so Residual standard error should be about 1 for a good fit.\n\n\n")


    printf("\n\nModel with all corrections at order two in Ninv and Sinv\n")
    quadFunction <- function(Ninv,Sinv,c,cN,cS,cNN,cNS,cSS){c + cN * Ninv + cS * Sinv + cSS * Sinv * Sinv + cNS * Sinv * Ninv + cNN * Ninv * Ninv}
    modelQuadAll <- nls(value~quadFunction(Ninv,Sinv,c,cN,cS,cNN,cNS,cSS), data=data_i, weights=data_i$weights, start=list(c=0.0,cN=-3.81,cS=0.0,cSS=0.0,cNN=0.0,cNS=0.0), nls.control(warnOnly = TRUE, tol = set_tolerance, minFactor = set_minFactor, maxiter = set_maxiter))
    print(summary(modelQuadAll))

   # printf("\n\nModel with all corrections at order two in Ninv and Sinv, no 1/N^4 and cNS = 0\n")
   # quadFunction <- function(Ninv,Sinv,c,cN,cS,cSS){c + cN * Ninv + cS * Sinv + cSS * Sinv * Sinv}
   # modelQuadAll <- nls(value~quadFunction(Ninv,Sinv,c,cN,cS,cSS), data=data_i, weights=data_i$weights, start=list(c=0.0,cN=0.0,cS=0.0,cSS=0.0), control = nls.control(warnOnly = TRUE, tol = set_tolerance, minFactor = set_minFactor, maxiter = set_maxiter))
   # print(summary(modelQuadAll))


    # Store coefficients as starting values for subsequent fits
    s <- summary(modelQuadAll)$coefficients


    # Additional models to try
    if(1==0) {

        printf("\n\nModel with all corrections at order two in Ninv and Sinv, no 1/N^4\n")
        quadFunction <- function(Ninv,Sinv,c,cN,cS,cNS,cSS){c + cN * Ninv + cS * Sinv + cSS * Sinv * Sinv + cNS * Sinv * Ninv}
        modelQuad <- nls(value~quadFunction(Ninv,Sinv,c,cN,cS,cNS,cSS), data=data_i, weights=data_i$weights, start=list(c=0.0,cN=0.0,cS=0.0,cSS=0.0,cNS=0.0), control = nls.control(warnOnly = TRUE))
        print(summary(modelQuad))


        printf("\n\nModel with all corrections at order two in Ninv and Sinv, no 1/N^4 and cNS = 0\n")
        quadFunction <- function(Ninv,Sinv,c,cN,cS,cSS){c + cN * Ninv + cS * Sinv + cSS * Sinv * Sinv}
        modelQuad <- nls(value~quadFunction(Ninv,Sinv,c,cN,cS,cSS), data=data_i, weights=data_i$weights, start=list(c=0.0,cN=0.0,cS=0.0,cSS=0.0), control = nls.control(warnOnly = TRUE, tol = set_tolerance, minFactor = set_minFactor, maxiter = set_maxiter))
        print(summary(modelQuad))


        printf("\n\nModel with fixed finite N^2 corrections, only c, cS and cSS free, no 1/N^4 or 1/N^2 S\n")
        quadFunction <- function(Ninv,Sinv,c,cS,cSS){c - 3.81 * Ninv + cS * Sinv + cSS * Sinv * Sinv}
        modelQuad <- nls(value~quadFunction(Ninv,Sinv,c,cS,cSS), data=data_i, weights=data_i$weights, start=list(c=s[1],cS=s[3],cSS=s[6]), control = nls.control(warnOnly = TRUE))
        print(summary(modelQuad))

        printf("\n\nModel with fixed finite N^2 corrections, only c, cS and cSS free, no 1/N^4 or 1/N^2 S, cN = -3.48\n")
        quadFunction <- function(Ninv,Sinv,c,cS,cSS){c - 3.48 * Ninv + cS * Sinv + cSS * Sinv * Sinv}
        modelQuad <- nls(value~quadFunction(Ninv,Sinv,c,cS,cSS), data=data_i, weights=data_i$weights, start=list(c=s[1],cS=s[3],cSS=s[6]), control = nls.control(warnOnly = TRUE, tol = 1e-10, minFactor = 1e-12, maxiter = 50000))
        print(summary(modelQuad))

        printf("\n\nModel with all corrections at order two in Ninv and Sinv, but fixed cN = -3.81 and no 1/N^4\n")
        quadFunction <- function(Ninv,Sinv,c,cS,cNS,cSS){c - 3.81 * Ninv + cS * Sinv + cSS * Sinv * Sinv + cNS * Sinv * Ninv}
        modelQuad <- nls(value~quadFunction(Ninv,Sinv,c,cS,cNS,cSS), data=data_i, weights=data_i$weights, start=list(c=s[1],cS=s[3],cSS=s[6],cNS=s[5]), control = nls.control(warnOnly = TRUE))
        print(summary(modelQuad))

        printf("\n\nModel with all corrections at order two in Ninv and Sinv, but fixed cN = -3.48 and no 1/N^4\n")
        quadFunction <- function(Ninv,Sinv,c,cS,cNS,cSS){c - 3.48 * Ninv + cS * Sinv + cSS * Sinv * Sinv + cNS * Sinv * Ninv}
        modelQuad <- nls(value~quadFunction(Ninv,Sinv,c,cS,cNS,cSS), data=data_i, weights=data_i$weights, start=list(c=s[1],cS=s[3],cSS=s[6],cNS=s[5]), control = nls.control(warnOnly = TRUE, tol = set_tolerance, minFactor = set_minFactor, maxiter = set_maxiter))
        print(summary(modelQuad))

        printf("\n\nModel with all corrections at order two in Ninv and Sinv, but fixed cN = -3.81\n")
        quadFunction <- function(Ninv,Sinv,c,cS,cNN,cNS,cSS){c - 3.81 * Ninv + cS * Sinv + cSS * Sinv * Sinv + cNS * Sinv * Ninv + cNN * Ninv * Ninv}
        modelQuad <- nls(value~quadFunction(Ninv,Sinv,c,cS,cNN,cNS,cSS), data=data_i, weights=data_i$weights, start=list(c=s[1],cS=s[3],cSS=s[6],cNN=s[4],cNS=s[5]), control = nls.control(warnOnly = TRUE, tol = set_tolerance, minFactor = set_minFactor, maxiter = set_maxiter))
        print(summary(modelQuad))

        printf("Predicted continuum energy at N=16: %f\n\n",
            predict(modelQuad, newdata = data.frame(
            Ninv = 1/16^2,
            Sinv = 1/100000000
            ), type=response) )

        printf("\n\nModel with all corrections at order two in Ninv and Sinv, but fixed cN = -3.48\n")
        quadFunction <- function(Ninv,Sinv,c,cS,cNN,cNS,cSS){c - 3.48 * Ninv + cS * Sinv + cSS * Sinv * Sinv + cNS * Sinv * Ninv + cNN * Ninv * Ninv}
        modelQuad <- nls(value~quadFunction(Ninv,Sinv,c,cS,cNN,cNS,cSS), data=data_i, weights=data_i$weights, start=list(c=s[1],cS=s[3],cSS=s[6],cNN=s[4],cNS=s[5]), control = nls.control(warnOnly = TRUE, tol = set_tolerance, minFactor = set_minFactor, maxiter = set_maxiter))
        print(summary(modelQuad))

        printf("\n\nModel with all corrections at order two in Ninv and Sinv, but fixed cN = -3.81 and cNS = 0\n")
        quadFunction <- function(Ninv,Sinv,c,cS,cNN,cSS){c - 3.81 * Ninv + cS * Sinv + cSS * Sinv * Sinv  + cNN * Ninv * Ninv}
        modelQuad <- nls(value~quadFunction(Ninv,Sinv,c,cS,cNN,cSS), data=data_i, weights=data_i$weights, start=list(c=s[1],cS=s[3],cSS=s[6],cNN=s[4]), control = nls.control(warnOnly = TRUE))
        print(summary(modelQuad))

        printf("Predicted continuum energy at N=16: %f\n\n",
            predict(modelQuad, newdata = data.frame(
            Ninv = 1/16^2,
            Sinv = 1/100000000
            ), type=response) )


        printf("\n\nModel with all corrections at order two in Ninv and Sinv, but fixed cN = -3.81, and cNN = 550\n")
        quadFunction <- function(Ninv,Sinv,c,cS,cNS,cSS){c - 3.81 * Ninv + cS * Sinv + cSS * Sinv * Sinv + cNS * Sinv * Ninv + 550 * Ninv * Ninv}
        modelQuad <- nls(value~quadFunction(Ninv,Sinv,c,cS,cNS,cSS), data=data_i, weights=data_i$weights, start=list(c=s[1],cS=s[3],cSS=s[6],cNS=s[5]), control = nls.control(warnOnly = TRUE))
        print(summary(modelQuad))

        printf("\n\nModel with all corrections at order two in Ninv and Sinv, but fixed cN = -3.81, and cNN = 550, cS = 9.67, cSS = -36 \n")
        quadFunction <- function(Ninv,Sinv,c,cNS){c - 3.81 * Ninv + 9.67 * Sinv -36 * Sinv * Sinv + cNS * Sinv * Ninv + 550 * Ninv * Ninv}
        modelQuad <- nls(value~quadFunction(Ninv,Sinv,c,cNS), data=data_i, weights=data_i$weights, start=list(c=s[1],cNS=s[5]), control = nls.control(warnOnly = TRUE))
        print(summary(modelQuad))

        printf("\n\nModel with all corrections at order two in Ninv and Sinv, + 1 / S S NN\n")
        quadFunction <- function(Ninv,Sinv,c,cN,cS,cNN,cNS,cSS, cSSN){c + cN * Ninv + cS * Sinv + cSS * Sinv * Sinv + cNS * Sinv * Ninv + cNN * Ninv * Ninv +  cSSN * Sinv * Sinv * Ninv}
        modelQuad <- nls(value~quadFunction(Ninv,Sinv,c,cN,cS,cNN,cNS,cSS, cSSN), data=data_i, weights=data_i$weights, start=list(c=0.0,cN=0.0,cS=0.0,cSS=0.0,cNN=0.0,cNS=0.0, cSSN=0.0), control = nls.control(warnOnly = TRUE))
        print(summary(modelQuad))

        printf("\n\nModel with all corrections at order two in Ninv and Sinv, + 1 / S S NN, + 1 / S S NN NN\n")
        quadFunction <- function(Ninv,Sinv,c,cN,cS,cNN,cNS,cSS, cSSN, cSSNN){c + cN * Ninv + cS * Sinv + cSS * Sinv * Sinv + cNS * Sinv * Ninv + cNN * Ninv * Ninv + cSSN * Sinv * Sinv * Ninv + cSSNN * Sinv * Sinv * Ninv * Ninv}
        modelQuad <- nls(value~quadFunction(Ninv,Sinv,c,cN,cS,cNN,cNS,cSS, cSSN, cSSNN), data=data_i, weights=data_i$weights, start=list(c=0.0,cN=0.0,cS=0.0,cSS=0.0,cNN=0.0,cNS=0.0, cSSN=0.0, cSSNN=0.0), control = nls.control(warnOnly = TRUE))
        print(summary(modelQuad))

    }



    model <- modelQuadAll

    dof <- summary(model)$df[2]
    large_NS_value <- summary(model)$coefficients[1,1]
    large_NS_error <- summary(model)$coefficients[1,2]


    res_by_error <- summary(model)$residuals / data_i$error

    all_data[all_data$observable == plot_labels[i], ]$predicted <- predict(model, newdata = data.frame(
        Ninv = all_data[all_data$observable == plot_labels[i], ]$Ninv,
        Sinv = all_data[all_data$observable == plot_labels[i], ]$Sinv
        ), type=response)
    
    if(plot_labels[i] == "energy") {
        a <- model$coefficients[1]
        b_N <- model$coefficients[2]
        b_S <- model$coefficients[3]
        c_S <- model$coefficients[4]
    }
   
  
}


all_data$predicted_normalized <- (  all_data$predicted - all_data$value) / all_data$error

large_NS_data <- data.frame(observable="energy", value=large_NS_value, error=large_NS_error, S=999, N=999, Ninv=0, Sinv=0, weights=0, predicted=gravity_prediction, predicted_normalized=(large_NS_value-gravity_prediction)/large_NS_error)


all_data <- rbind(all_data, large_NS_data)



subheading <- sprintf("Gravity: %.3f    Model: %.3f \u00B1 %.3f", gravity_prediction, large_NS_value, large_NS_error)


title <- sprintf("2D fit of E/N^2 up to quadratic order for T=%.2f, mu=0.5", T)

data_plot <- ggplot(aes(x=Sinv, y=Ninv, color=predicted_normalized), data=subset(all_data, all_data[,"observable"] == "energy")) +
    #geom_point(aes(size = (predicted_normalized*predicted_normalized)^0.25)) +
    #scale_colour_gradient2() + 
    ggtitle(title) +
    scale_fill_gradient2(name=expression(frac(model - data, sigma[data])), low="blue", high="red", mid="white") +
    geom_point(aes(fill=predicted_normalized), colour="black", pch=21, size=5) +
    scale_x_continuous(breaks=Sinv_list, labels=S_ticks, minor_breaks = NULL) + 
    scale_y_continuous(breaks=Ninv_list, labels=N_ticks, minor_breaks = NULL) + 
    geom_hline(yintercept=smallest_N^(-2)*1.05, linetype="dashed", color = "black", size=0.5) +
    geom_vline(xintercept=smallest_S^(-1)*1.05, linetype="dashed", color = "black", size=0.5) +
    geom_hline(yintercept=largest_N^(-2)/1.05-0.0003, linetype="dashed", color = "black", size=0.5) +
    geom_vline(xintercept=largest_S^(-1)/1.05, linetype="dashed", color = "black", size=0.5) +
    theme_bw() +  
    xlab("L") + 
    ylab(p("N")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(subtitle = subheading) +
    theme(legend.position="bottom")



plot_file <- p("observable_large_NS_", UGB, file_pattern, ".pdf")

suppressMessages(ggsave(data_plot, file=plot_file))


printf("\n\nHistogram of normalized deviations\n\n")

max_sigma <- 3
bin_width_histogram <- 0.01
nbins <- 2* max_sigma / bin_width_histogram +1


bins <- -max_sigma - bin_width_histogram / 2 + bin_width_histogram * (0:(nbins))
#print(bins)
mids <- bins + bin_width_histogram / 2
mids <- mids[1:(length(mids)-1)]



hist_dev <- ggplot(data=subset(all_data, all_data[,"observable"] == "energy"), aes(x=predicted_normalized)) + 
  geom_histogram(binwidth=bin_width_histogram, aes(y=..count../sum(..count..)/bin_width_histogram)) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1)) +
  theme_bw() + 
  xlab("normalized deviation of model") + 
  ylab("density") + 
  theme(legend.position="bottom") 

suppressMessages(ggsave(hist_dev, file="hist.pdf"))


printf("\n\nIntegrated histogram of normalized deviations for KS test\n\n")



test_data <- subset(all_data, all_data[,"observable"] == "energy")$predicted_normalized



ks_data <- hist(test_data, breaks=bins, plot=FALSE)$density * bin_width_histogram


fgauss <- function(x, mu, sigma)
{
    return (exp(-0.5 * ((x - mu)/sigma)^2) / (2 * 3.14159265359)^0.5 / sigma)
}

error_data_fit = data.frame(mids, ks_data)

 printf("\n\nFit histogram of fitting errors with Gaussian: \n")
    modelError <- nls(ks_data~fgauss(mids, mu, sigma), data=error_data_fit, start=list(mu=0, sigma=1))
    print(summary(modelError))

printf("\n\nFit histogram of fitting errors with Gaussian, mu=0: \n")
    modelError <- nls(ks_data~fgauss(mids, 0, sigma), data=error_data_fit, start=list(sigma=1))
    print(summary(modelError))






for(i in 2:length(ks_data)) {
    ks_data[i] <- ks_data[i] + ks_data[i-1]
}

ks_data <- c(ks_data,ks_data[length(ks_data)])   


bin_width_ks <- 0.01

nbins_ks <- 2* max_sigma / bin_width_ks 

erf_data <- data.frame(matrix(ncol=2, nrow=nbins_ks+1))
erf_data[,1] <-  -max_sigma - bin_width_ks / 2 + bin_width_ks * (0:(nbins_ks))
erf_data[,2] <- pnorm(erf_data[,1])





d_max = 0

for(i in 1:(length(ks_data)-1)) {
    d_max = max(d_max, abs(ks_data[i]-pnorm(bins[i], mean = 0, sd = 1)))
    d_max = max(d_max, abs(ks_data[i]-pnorm(bins[i+1], mean = 0, sd = 1)))
}

vcrit_95 = 1.358 / (dof)^0.5



printf("Critical 95%% value for %d DOF: %f\n", as.integer(dof), vcrit_95)
printf("d_max from KS-Test: %f\n\n", d_max)

if(d_max < vcrit_95) {
    printf("KS-Test passed.\n\n")
} else {
    printf("KS-Test NOT passed.\n\n")
}




ks_frame <- data.frame(mids, ks_data[1:(length(ks_data)-1)])


colnames(ks_frame) <- c("mids", "ks_data")


printf("\n\nFit KS histogram with Error function: \n")
modelKS <- nls(ks_data~pnorm( (mids - mu) / sigma), data=ks_frame, start=list(mu=0, sigma=1))
print(summary(modelKS))

mu <- summary(modelKS)$coefficients[1,1]
sigma <- summary(modelKS)$coefficients[2,1]



erf_data_fit <- data.frame(matrix(ncol=1, nrow=nbins_ks+1))
erf_data_fit <- pnorm((erf_data[,1]-mu) / sigma)


d_max_fit = 0

for(i in 1:(length(ks_data)-1)) {
    d_max_fit = max(d_max_fit, abs(ks_data[i]-pnorm(bins[i], mean = mu, sd = sigma)))
    d_max_fit = max(d_max_fit, abs(ks_data[i]-pnorm(bins[i+1], mean = mu, sd = sigma)))
}

printf("d_max from KS-Test with fitted pnorm: %.2f\n\n", d_max_fit)



vcrit_95_fit = 0.895 / ( (dof)^0.5 - 0.01 + 0.83 /(dof)^0.5 )



printf("Critical 95%% value for %d DOF: %f\n", as.integer(dof), vcrit_95_fit)

if(d_max_fit < vcrit_95_fit) {
    printf("KS-Test with fitted pnorm passed.\n\n")
} else {
    printf("KS-Test with fitted pnorm NOT passed.\n\n")
}

sub = sprintf("d_max = %.2f,   95%% critical value for %d DOF = %.2f", d_max, dof, vcrit_95)


if(d_max < vcrit_95) {
    sub = sprintf("%s,   KS-Test passed!", sub)
} else {
    sub = sprintf("%s,   KS-Test NOT passed!", sub)
}

sub_fit = sprintf("d_max_fit = %.2f,   95%% critical value for %d DOF = %.2f", d_max_fit, dof, vcrit_95_fit)



hist_int_dev <- ggplot() + 
#  labs(title = "KS Test against standard (fitted) normal") +
#  annotate("text", x = -1, y = 1, label = sub_fit, color="blue") +
#  labs(subtitle = sub) +
  geom_step(aes(x=bins, y=ks_data)) +
  geom_line(aes(x=erf_data[,1], y=erf_data[,2]), color="blue") +
#  geom_line(aes(x=erf_data[,1], y=erf_data_fit), linetype="dashed", color = "blue", label="Fitted error function") +
  theme_bw() + 
 theme(text = element_text(size=20)) + 
  xlab("normalized deviation of model < x") + 
  ylab("density") + 
  theme(legend.position="right") 

hist_int_file <- p("observable_large_NS_KS_test_", UGB, file_pattern, ".pdf")

suppressMessages(ggsave(hist_int_dev, file=hist_int_file))
