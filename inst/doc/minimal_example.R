## ---- message = FALSE---------------------------------------------------------
# install.packages("spectre")

## ---- message = FALSE---------------------------------------------------------
# install.packages("devtools")
# devtools::install_github("r-spatialecology/spectre") #Uncomment when repo is public

## ---- message = FALSE---------------------------------------------------------
library("gdm")

# create random species composition
set.seed(42) 

nspecies <- 20
nsites <- 15
presence_prob <- 0.3 # probability for each species to be present at each site

get_species_list <- function(nspecies, nsites, presence_prob)
{
  # generate a siteXspecies list with a random set of species presences/absences
  # fill list with random species presences
  m <- matrix(nrow = nspecies, ncol = nsites, data = 0)
  
  for (row in 1:ncol(m)) {
    for (col in 1:nrow(m)) {
      if (runif(1) < presence_prob) {
        m[col, row] <- 1
      }
    }
  }
  
  mode(m) <- "integer"
  
  return(m)
}

# random species composition
species_list <- get_species_list(nspecies = nspecies, nsites = nsites, 
                                 presence_prob = presence_prob) 

# calculation of Bray-Curtis dissimilarity with the gdm package: 
# bioData is required by the gdm package, but does not affect 
# the observed Bray-Curtis dissimilarity we will use later
bioData <- data.frame(site_id = 1:nsites, x_coords = rep(13, nsites), 
                      y_coords = rep(10, nsites)) 

bioData <- cbind(bioData, t(species_list)) 

predData <- data.frame(site_id = 1:nsites, preds = runif(nsites))

sitepairs <- gdm::formatsitepair(bioData = bioData, bioFormat = 1, abundance = FALSE, 
                                 siteColumn = "site_id",
                                 XColumn = "x_coords", YColumn = "y_coords", 
                                 predData = predData)

gdm_result <- gdm::gdm(sitepairs, geo = TRUE) 

## ---- message = FALSE, warning = FALSE----------------------------------------
library("spectre")

# Calculate objective_matrix from (modelled) alpha-diversity and Bray-Curtis dissimilarity
alpha_list <- colSums(species_list) # alpha-diversity of random species community

objective_matrix <- spectre::generate_commonness_matrix_from_gdm(gdm_predictions = gdm_result$observed, 
                                                                 alpha_list = alpha_list)

# Solve composition 
res <- spectre::run_optimization_min_conf(alpha_list = alpha_list,
                                          total_gamma = nspecies,
                                          target = objective_matrix,
                                          max_iterations = 1000) # n iterations

## ---- message = FALSE---------------------------------------------------------
error_c <- spectre::calc_commonness_error(x = res, objective_matrix = objective_matrix)

## ---- echo = TRUE, include = TRUE, out.width="50%"----------------------------
# With an increasing number of iterations, the solution matrix improved
spectre::plot_error(x = res)

# Plot commonness error between objective matrix and solution matrix
spectre::plot_commonness(x = res, target = objective_matrix)

