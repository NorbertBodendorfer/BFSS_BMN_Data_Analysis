###### Error Bar functions for histograms etc.


# jackknife bin width = jack_mult * autocorrelation length

plot_acl <- 1


# Calculate the autocorrelation length of sample data

AutocorrelationLength <- function(hist_data) {

    ndata <- length(hist_data)

    # Compute autocorrelation length for jackknife method
    auto_correlation <- acf(hist_data, lag.max=ndata, plot = FALSE)
 	   
    if(plot_acl == 1) {
        pdf(file="acl_plot.pdf")
        plot(auto_correlation)
        dev.off()
    }

    # Find the first negative value of the acf, this will be w_jack
    i_ac = 1
    while( (auto_correlation$acf[i_ac]>0) & (i_ac < ndata) ) {
        i_ac <- i_ac+1
        #print(auto_correlation$acf[i_ac])
    }

    if(i_ac > ndata/2) {
        warning("Autocorrelation length > ndata/2.")
    } else {
        #printf("Autocorrelation length: %d.\n", i_ac)
    }

    return(i_ac)

}



# Calculates the Jackknife error on a data series with autocorrelation length acl

JackknifeDelta_WSingle <- function(data_series_jack, acl, jack_mult) {

    ndata <- length(data_series_jack)

    if((acl*jack_mult>=ndata / 3)) {
        warning("Autocorrelation length > ndata/2. Using 2 Jackknife bins.")
        n_jack <- 2
		w_jack <- ndata %/% n_jack
    } else {
        w_jack <- jack_mult * acl
		n_jack <- ndata %/% w_jack
        printf("Autocorrelation length: %d <=> %d jackknife bins with w = %d * acl.\n", acl, n_jack, jack_mult)
    }

    # Create Jackknife bin means
    jack_data <- rep(0, n_jack)

    for(k_jack in 1:n_jack) {
        jack_data[k_jack] <- mean(data_series_jack[((k_jack-1)*w_jack+1):(k_jack*w_jack)])
	}


    # Calculate variance of jackknife data and obtain Delta_w, the jackknife error
    jack_var <- var(jack_data)

    Delta_w = (jack_var / n_jack )^0.5

    printf("Jackknife with %d bins: w=%d, Delta=%f \n", n_jack, w_jack, Delta_w)

    return(Delta_w)
}




JackknifeHistogramACL <- function(hist_data, mult_histjack) {
   
    ndata <- length(hist_data)

    acl <- AutocorrelationLength(hist_data)
    
    
    # We need at least 2 bins
    if((acl*mult_histjack<=ndata / 2)) {
        n_histjack <- ndata %/% (acl*mult_histjack)
        w_histjack <- ndata %/% n_histjack
        printf("\nAutocorrelation length: %d  = %.0f%% of data <=> %d jackknife bins with w = %d * acl.\n", acl, acl/ndata*100, n_histjack, mult_histjack)
    } else {
        w_histjack <- ndata %/% 2
        n_histjack <- 2
        warning(sprintf("%d * Autocorrelation length > ndata / 2 for %s and T=%s.", mult_histjack, plot_loop, temperature))
        printf("\nAutocorrelation length: %d = %.0f%% of data <=> 2 jackknife bins with w = ndata/2. (see warning)\n", acl, acl/ndata*100)
        #print(ndata)
    }
    
    # Create histogram from jackknife bin
    binned_table <- foreach(k_histjack=1:n_histjack, .combine=rbind) %do% {
        h <- hist(hist_data[((k_histjack-1)*w_histjack+1):(k_histjack*w_histjack)], breaks=bins, include.lowest = T, right=FALSE, plot = FALSE)
        data.table((h$mids-bin_width/2), h$density)
    }

    # Initialize return object histogram jack errors (hje)
    hje <- rep(0, nbins)

    for(i_bin in 1:nbins) {
        bin_left <- (h$mids-bin_width/2)[i_bin]
        # print(binned_table[V1==bin_left])
        histjack_error <- sd(binned_table[V1==bin_left, V2]) / n_histjack^0.5
        hje[i_bin] <- histjack_error
    }
    return(hje)
}


JackknifeHistogramFixedNBin <- function(hist_data, n_jackbin) {
   
    ndata <- length(hist_data)    
    
    # We need at least 2 bins
    
    w_histjack <- ndata %/% n_jackbin
    n_histjack <- n_jackbin
    printf("\nUsing %d jackknife bins with w = ndata/%d in JackknifeHistogramFixedNBin.\n", n_jackbin, n_jackbin) 

    # Create histogram from jackknife bin
    binned_table <- foreach(k_histjack=1:n_histjack, .combine=rbind) %do% {
        h <- hist(hist_data[((k_histjack-1)*w_histjack+1):(k_histjack*w_histjack)], breaks=bins, include.lowest = T, right=FALSE, plot = FALSE)
        data.table((h$mids-bin_width/2), h$density)
    }

    # Initialize return object histogram jack errors (hje)
    hje <- rep(0, nbins)

    for(i_bin in 1:nbins) {
        bin_left <- (h$mids-bin_width/2)[i_bin]
        histjack_error <- sd(binned_table[V1==bin_left, V2]) / n_histjack^0.5
        hje[i_bin] <- histjack_error
    }
    return(hje)
}

