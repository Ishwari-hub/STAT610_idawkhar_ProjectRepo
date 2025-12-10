

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


##building sparse matrix 

cases <- sort(unique(c(allcites$caseid_citingcase, allcites$caseid_citedcase)))
n_cases <- length(cases)
n_cases

#assigning index to the cases 
case_index <- seq_along(cases)
names(case_index) <- cases

## making rows and columns for the sparse matrix 
rows <- case_index[as.character(allcites$caseid_citingcase)]
columns <- case_index[as.character(allcites$caseid_citedcase)]

length(rows); length(columns)

## using rows and columns to build a sparse matrix
## creating a sparse adj matrix 
library(Matrix)

adj_matrix_func <- function(citation_subset) {
  # unique case IDs in the data subset
  cases <- sort(unique(c(citation_subset$citing_caseid, citation_subset$cited_caseid)))
  
  # Mapping case IDs to row/column indices
  case_index <- seq_along(cases)
  names(case_index) <- cases
  
  rows <- case_index[as.character(citation_subset$citing_caseid)]
  cols <- case_index[as.character(citation_subset$cited_caseid)]
  
  # sparse adjacency matrix
  A_sparse <- sparseMatrix(
    i = rows,
    j = cols,
    x = 1,
    dims = c(length(cases), length(cases)),
    dimnames = list(as.character(cases), as.character(cases))  
  )
  
  return(A_sparse)
}

## tweaking the HITS algorthim function 
score_calc_func <- function(A_sparse,max_iterations,tol){
  n <- nrow(A_sparse)
  h <- rep(1, n)   
  a <- rep(1, n)   
  
  for (iter in 1:max_iterations) {
    
    a_new <- normalize_func(t(A_sparse) %*% h)
    h_new <- normalize_func(A_sparse %*% a_new)
    
    # checking convergence before updating
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
validate_func_subset <- function(case_subset, a_vec, h_vec, cases_vector) {
  index_positions <- match(case_subset$caseid, cases_vector)
  
  # extracting computed scores safely
  computed_auth <- a_vec[index_positions]
  computed_hub  <- h_vec[index_positions]
  
  result <- data.frame(
    caseid        = case_subset$caseid,
    parties       = case_subset$parties,
    computed_auth = computed_auth,
    computed_hub  = computed_hub,
    original_auth = case_subset$auth,
    original_hub  = case_subset$hub
  )
  
  list(
    result_df = result
  )
}

## creating a subset of data to check the scores 

##testing scores for Roe v. Wade 

roe_cases <- data[grep("Roe v. Wade", data$parties, ignore.case = TRUE), ]
roe_cases$caseid  

roe_index <- which(cases == 25347)  #caseID

#calculating scores for the entire function 
full_scores <- score_calc_func(A_sparse, max_iterations = 100, tol = 1e-6)

# Extracting the computed scores of  HITS algorithm
authority_roe <- full_scores$authority_score[roe_index]
hub_roe       <- full_scores$hub_score[roe_index]

#calling validation function to confirm values for Roe V. Wade 
validation_roe <- validate_func_subset(
  roe_cases,
  full_scores$authority_score,
  full_scores$hub_score,
  cases
)
validation_roe


brown_cases <- data[grep("Brown v. Board of Education", data$parties, ignore.case = TRUE), ]
brown_cases$caseid  

#call function 
full_scores <- score_calc_func(A_sparse, max_iterations = 100, tol = 1e-6)

# Extracting the computed scores of HITS algorithm 
authority_brown <- full_scores$authority_score[brown_index]
hub_brown       <- full_scores$hub_score[brown_index]

#calling validation function to confirm values for Brown v. board of education
validation_brown <- validate_func_subset(
  brown_cases,
  full_scores$authority_score,
  full_scores$hub_score,
  cases
)
validation_brown


