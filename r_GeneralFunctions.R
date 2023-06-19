#### Import printf function 
printf <- function(...) cat(sprintf(...))

#### Define paste function for multiple strings
p <- function(..., sep='') {
    paste(..., sep=sep, collapse=sep)
}



#### Print name of R script in format

#|------------------------|
#| R data analysis script |
#|------------------------|

PrintScriptName <- function(name) 
{
	name <- p("| ", name, " |")
	border <- p("|", strrep("-", (nchar(name)-2)), "|")
	printf("\n\n%s\n%s\n%s\n\n\n", border, name, border)
}




##### Extracting a String Between 2 Characters in R
##### By Eric Cai - The Chemical Statistician


# getstr() is my customized function
# it extracts a string between 2 characters in a string variable
getstr = function(mystring, initial.character, final.character)
{
 
     # check that all 3 inputs are character variables
     if (!is.character(mystring))
     {
          stop('The parent string must be a character variable.')
     }
 
     if (!is.character(initial.character))
     {
          stop('The initial character must be a character variable.')
     }
 
 
     if (!is.character(final.character))
     {
          stop('The final character must be a character variable.')
     }

  
 
     # pre-allocate a vector to store the extracted strings
     snippet = rep(0, length(mystring))
 
 
 
     for (i in 1:length(mystring))
     {
          # extract the initial position
          initial.position = gregexpr(initial.character, mystring[i])[[1]][1] + 1
  
          # extract the final position
          final.position = gregexpr(final.character, mystring[i])[[1]][1] - 1
 
          # extract the substring between the initial and final positions, inclusively
          snippet[i] = substr(mystring[i], initial.position, final.position)
     }
 
     return(snippet)
}