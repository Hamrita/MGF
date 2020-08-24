# This code allows you to calculate and evaluate the various moments and centralized moments of different Probability Distributions.
# Its purpose is to aid in the computation and understanding of these concepts so that students can better grap these concepts.

#####################################################################################################################
### These are the list of parameters for all of the distributions set them to the values
### you wish to examine in the MGF_evaluator function as desired.
### These are the default values for all parameters.

# For the Normal Distribution
#mu = 0 # Mean for the Normal Distribution.
#sigma = 1 # Standard Deviation for the Normal Distribution.

# For the Bernoulli, Geometric, Negative Binomial and Binomial Distribution
#p = 0.5 # The success probabiltiy for a Bernoulli, Geometric, Negative Binomial and Binomial Distribution.
#n = 10 # The number of trials for a Binomial Distribution.
#r = 5 # The number of successes for the Negative Binomial Distribution.

# For the Chi-square, Gamma, Exponential and Poisson Distributions.
#k = 2 # Degrees of freedom for the Chi-Square Distribution.
#lambda = 5 # for the Poisson, Exponential and Gamma Distribution.
#alpha = 1 # for the Gamma Distribution.
#####################################################################################################################

#' Moment Generating Function.
#'
#' mgf returns the moment generating function of the specified distribution
#'@param distribution One of the following distributions as a string: {"Bernoulli", "Binomial", "Geometric", "Negative Binomial", "Poisson",  "Exponential", "Gamma", "Normal", "Chi-Square"}.
#'@return The Moment Generating Function of the specified distribution.
#'@examples
#'mgf("Bernoulli")
#'mgf("Binomial")
#'mgf("Normal")
#'@export
mgf <- function(distribution){
  if (distribution == "Bernoulli"){
    return("1-p+p*exp(t)")
  }
  if (distribution == "Geometric"){
    return("((p*exp(t))/(1-(1-p)*exp(t)))")
  }
  if (distribution == "Negative Binomial"){
    return("((p*exp(t))/(1-(1-p)*exp(t)))^r")
  }
  if (distribution == "Binomial"){
    return("(1-p+p*exp(t))^n")
  }
  else if (distribution == "Normal"){
    return("exp(mu*t + 0.5*sigma^2*t^2)")
  }
  else if (distribution == "Exponential"){
    return("lambda/(lambda-t)")
  }
  else if (distribution == "Gamma"){
    return("(lambda/(lambda-t))^alpha")
  }
  else if (distribution == "Poisson"){
    return("exp(lambda*(exp(t)-1))")
  }
  else if (distribution == "Chi-Square" || distribution == "Chi Square"){
    return("(1-2*t)^(-k/2)")
  }
}
#'Cumulant Generating Function.
#'
#'cgf returns the cumulant generating function of the specified distribution
#'the derivatives of the cgf correspond to the centralized moments of the distribution centered about the mean.
#'@param distribution One of the following distributions as a string: {"Bernoulli", "Binomial", "Geometric", "Negative Binomial", "Poisson",  "Exponential", "Gamma", "Normal", "Chi-Square"}.
#'@return The Cumulant Generating Function of the specified distribution.
#'@examples
#'cgf("Bernoulli")
#'cgf("Binomial")
#'cgf("Normal")
#'@export
cgf <- function(distribution){
  if (distribution == "Bernoulli"){
    return("log(1-p+p*exp(t))")
  }
  if (distribution == "Binomial"){
    return("log((1-p+p*exp(t))^n)")
  }
  if (distribution == "Geometric"){
    return("log((p*exp(t))/(1-(1-p)*exp(t)))")
  }
  if (distribution == "Negative Binomial"){
    return("r*log((p*exp(t))/(1-(1-p)*exp(t)))")
  }
  else if (distribution == "Normal"){
    return("mu*t + 0.5*sigma^2*t^2")
  }
  else if (distribution == "Exponential"){
    return("log(lambda/(lambda-t))")
  }
  else if (distribution == "Gamma"){
    return("log((lambda/(lambda-t))^alpha)")
  }
  else if (distribution == "Poisson"){
    return("lambda*(exp(t)-1)")
  }
  else if (distribution == "Chi-Square" || distribution == "Chi Square"){
    return("log((1-2*t)^(-k/2))")
  }
}


#'Moment and Cumulant Generating Function Evaluator.
#'
#'This function evaluates the nth order derivative of each mgf and cgf for the selected probability distribution. It provdes both the formula as well as the value for each.
#'@param distribution One of the following distributions as a string: {"Bernoulli", "Binomial", "Geometric", "Negative Binomial", "Poisson",  "Exponential", "Gamma", "Normal", "Chi-Square"}.
#'@param t The value of t at which the derivatives are evaluated, default value is 0.
#'@param order_of_moment The order of the moments and centralized moments of the distribution that you are interested in, default = 1.
#'@param mu The mean of the Normal Distribution, default = 0.
#'@param sigma The Standard Deviation of the Normal Distribution default = 1.
#'@param n The number of trials for a Binomial Distribution, default = 10
#'@param p The success probabiltiy for a Bernoulli, Geometric, Negative Binomial and Binomial Distribution, default = 0.5.
#'@param r The number of successes for the Negative Binomial Distribution, default = 5.
#'@param lambda Rate parameter for the Poisson, Exponential and Gamma Distributions, default = 5.
#'@param alpha Shape parameter for the Gamma Distribution, default = 1.
#'@param k Degrees of freedom for the Chi-Square Distribution, default = 2.
#'@return
#' Returns the MGF, The formula of the nth order derivatives of the MGF and CGF and their values for the specified distribution and parameters.
#'
#'@examples
#'MGF_evaluator(distribution = "Binomial",order_of_moment = 2,n = 20)
#'MGF_evaluator(distribution = "Normal",t = 0,order_of_moment = 4)
#'MGF_evaluator(distribution = "Poisson",t = 0,order_of_moment = 2)
#'MGF_evaluator(distribution = "Exponential",t = 0,order_of_moment = 2)
#'@export
MGF_evaluator <- function(distribution,t = 0, order_of_moment =1,mu =0,sigma =1,n = 10,p = 0.5,r = 5,lambda = 5,alpha = 1,k = 2){

  # First check if the order of the moment is greater than 0 so that derivatives can be calculated appropriately
  if(order_of_moment > 0){
    # Use if statements to create the appropriate suffix for the order of the moment.
    if(order_of_moment == 1){
      suffix = "st"
    }
    else if(order_of_moment == 2){
      suffix = "nd"
    }
    else if(order_of_moment == 3){
      suffix = "rd"
    }
    else{
      suffix = "th"
    }
    # Print out the MGF of the selected probability distribution
    cat("The Moment Generating Function for the ",distribution,"distribution is: \n")
    print(mgf(distribution))

    # Parse the MGFs to make them into expression format so that their derivatives can be
    # calculated.
    mgf_for_derivative <- parse(text = mgf(distribution))
    cgf_for_derivative <- parse(text = cgf(distribution))

    # Take the derivatives of the MGF and CGF up to the desired term for calculating their moments.
    d <- D(mgf_for_derivative,"t")
    ln_d <- D(cgf_for_derivative,"t")
    if (order_of_moment != 1){
      for (i in 1:(order_of_moment-1)){
        d <- D(d,"t")
        ln_d <- D(ln_d,"t")
      }}

    # Print the moment generating function and the cumulant generating function
    cat("The formula for the ",order_of_moment,suffix, " moment is: \n",sep = "")
    print(d)
    cat("The value of the ",order_of_moment,suffix," moment is: \n",sep = "")
    print(eval(d))
    cat("The formula of the ",order_of_moment,suffix, " centralized moment is: \n",sep = "")
    print(ln_d)
    cat("The ",order_of_moment,suffix, " centralized moment's value is: \n", sep = "")
    print(eval(ln_d))
  }
  else {
    return(mgf(distribution))
  }
}

