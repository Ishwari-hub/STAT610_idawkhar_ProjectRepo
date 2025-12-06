

#Load data (csv file) 
setwd("/Users/ishwaridawkhar/Desktop/STATS Computing/Judicial") 
data <- read.csv("judicial.csv", stringsAsFactors = FALSE)
head(data)

#Load data (all citations text file) 
allcites <- read.table("allcites.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(allcites) <- c("caseid_citingcase", "caseid_citedcase")
head(allcites)

#To understand time range of the data
summary(data$year)

#adding the cited_year and citing_year to citation data. The year is needed for filtering citations later
case_years <- data$year
names(case_years) <- data$caseid

allcites$citing_year <- case_years[allcites$caseid_citingcase]
allcites$cited_year  <- case_years[allcites$caseid_citedcase]

head(allcites)

#test adj matrix logic on subset of citations data 
test_sub <- allcites[1:10, ]

#extracting unique case IDs for creating the adj matrix 
##cases_test <- sort(unique(c(test_sub$caseid_citingcase, test_sub$caseid_citedcase)))

#mapping case IDs to matrix rows and columns 
##case_index <- seq_along(cases_test)
##names(case_index) <- cases_test

#using index to assign row and column values 
##rows <- case_index[as.character(test_sub$caseid_citingcase)] 
##cols <- case_index[as.character(test_sub$caseid_citedcase)]

#create adjacency matrix 
##A_test <- matrix(0, nrow = length(cases_test), ncol = length(cases_test))
##A_test[cbind(rows, cols)] <- 1
##A_test

## building functions 

##adjacency matrix function 
adj_matrix_func <- function(test_sub){
  cases_test <- sort(unique(c(test_sub$caseid_citingcase, test_sub$caseid_citedcase)))
  
  #mapping case IDs to matrix rows and columns 
  case_index <- seq_along(cases_test)
  names(case_index) <- cases_test
  
  #using index to assign row and column values 
  rows <- case_index[as.character(test_sub$caseid_citingcase)] 
  cols <- case_index[as.character(test_sub$caseid_citedcase)]
  
  #create adjacency matrix 
  A_test <- matrix(0, nrow = length(cases_test), ncol = length(cases_test))
  A_test[cbind(rows, cols)] <- 1
  
  return (A_test)
}

## function to calculate hub scores and authority scores
score_calc_func <- function(A_test,max_iterations,tol){
  n <- nrow(A_test)
  h <- rep(1, n)   # initial hub scores
  a <- rep(1, n)   # initial authority scores
  
  for (iter in 1:max_iterations) {
    
    a_new <- normalize_func(t(A_test) %*% h)
    h_new <- normalize_func(A_test %*% a_new)
    
    # check convergence BEFORE updating
    if (max(abs(a_new - a)) < tol && max(abs(h_new - h)) < tol) {
      return(list(authority_score = a_new, hub_score = h_new, iterations = iter))
    }
    
    a <- a_new
    h <- h_new
  }
  
  return(list(authority_score = a, hub_score = h, iterations = max_iterations))
}
  

##function to normalize the calculated scores
normalize_func <- function(x){
  return(x / sqrt(sum(x^2)))
}

##function to validate scores 
val_func_subset <- function(subset_test_data, a, h) {
  subset_test_data$computed_auth <- a
  subset_test_data$computed_hub  <- h
  
  # correlations
  auth_cor <- cor(subset_test_data$auth, subset_test_data$computed_auth)
  hub_cor  <- cor(subset_test_data$hub, subset_test_data$computed_hub)
  
  list(compare = subset_test_data,
       authority_correlation = auth_cor,
       hub_correlation = hub_cor)
}


A_test <- adj_matrix_func(test_sub)
print(A_test)

test_run <- score_calc_func(A_test, max_iterations = 100, tol = 1e-6)
print(test_run$authority_score)
print(test_run$hub_score)

cases_test <- sort(unique(c(test_sub$caseid_citingcase, test_sub$caseid_citedcase)))

subset_test_data <- data.frame(
  caseid = cases_test,
  auth   = runif(length(cases_test)),  
  hub    = runif(length(cases_test))   
)

validation <- val_func(data, 
                       a = test_run$authority_score, 
                       h = test_run$hub_score, 
                       caseid = cases_test)

validation$authority_correlation
validation$hub_correlation


  

  
  
  


