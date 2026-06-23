# ============================================================
# Hands-on workshop: Sequence and cluster analysis
# ============================================================
#
# In this exercise, we will apply sequence analysis to characterize
# lifetime Earned Income Tax Credit (EITC) eligibility trajectories
# across adulthood.
#
# The data and workflow are adapted from:
#
# Pacca L, et al. Lifetime Patterns of Earned Income Tax Credit
# Eligibility and Cognition: A Sequence Analysis Approach.
#
# The broader research project had two aims:
#
#   1. Identify common patterns of EITC eligibility between
#           ages 25-45, focusing on differences in timing,
#           duration, and order of eligibility.
#
#   2. Evaluate whether these patterns are associated with
#           cognition later in life.
#
## In this workshop, we will focus on Aim 1 and use sequence and
# cluster analysis to identify trajectory groups.

# The goal is to:
#   1. Create individual-level EITC eligibility sequences
#   2. Visualize trajectory patterns
#   3. Compute dissimilarities between trajectories
#   4. Identify trajectory groups using cluster analysis
#
# These trajectory groups can then be used as exposures in
# conventional statistical analyses, such as regression models,
# to address Aim 2.
#
# Throughout the exercise, several methodological decisions are left
# open for discussion. There is no single correct solution.
#
# Contacts: 
#Lucia Pacca, pacca@wustl.edu
#Jilly Hebert, hebertj@wustl.edu
# ============================================================

# ============================================================
# Install and load required packages
# ============================================================

# Uncomment if needed:
# install.packages(c(
#   "cluster",
#   "WeightedCluster",
#   "TraMineR"
# ))

library(cluster)         # clustering algorithms
library(WeightedCluster) # cluster quality measures
library(TraMineR)        # sequence analysis
library(TraMineRextras)  # additional sequence visualizations
library(ggseqplot)       # sequence plots
library(BurStMisc)       # visualize data
library(tidyverse)       # cleaning

# ============================================================
# STEP 0: Import the cleaned sequence data
# ============================================================

# The dataset is already cleaned and contains one row per participant.
# Each age-specific column contains the participant's EITC eligibility status at that age.

EITC_wide <- as.data.frame(read_csv("EITC_Data.csv"))

View(EITC_wide) #Inspect the data

# ============================================================
# NOTE ON MISSING DATA
# ============================================================
#
# In the original NLSY79 data, some income variables used to determine EITC eligibility were missing due to non-response.
#
# Prior to creating the trajectories, missing income values were linearly interpolated between adjacent observations.
#
# This approach is commonly used for income variables because it preserves gradual changes in earnings over time and avoids
# creating abrupt shifts that may result from carrying values forward or backward.
#
# As a result, the EITC eligibility trajectories used in this exercise contain no missing values.
#
# In practice, handling missing data is an important decision in sequence analysis and may involve interpolation, carrying
# values forward/backward, multiple imputation, or treating missingness as a separate state.
# ============================================================

# ============================================================
# STEP 1a: Define the sequence states (categories)
# ============================================================

# State definitions:
# 1 = EITC ineligible, no income
# 2 = EITC eligible
# 3 = EITC ineligible, income above eligibility threshold

state_alphabet <- c("1", "2", "3")

state_labels <- c(
  "EITC ineligible, no income",
  "EITC eligible",
  "EITC ineligible, > threshold"
)

# Colors used throughout the workshop

No_income      <- "#999999"  # grey
EITC_eligible  <- "#E41A1C"  # red
Above_threshold <- "#377EB8" # blue

state_colors <- c(
  No_income,
  EITC_eligible,
  Above_threshold
)

# ============================================================
# STEP 1b: Create the individual trajectories (sequences)
# ============================================================

# Identify sequence columns
seq_cols <- grep("^eitc_eligibility", names(EITC_wide))

# Rename sequence columns to ages (for graphs)
names(EITC_wide)[seq_cols] <- 25:45

# Create the sequence object used by TraMineR

EITC_seq <- seqdef(
  
  EITC_wide,             # Dataset containing the sequences
  var = seq_cols,        # Columns that contain the age-specific states
  alphabet = state_alphabet, # Possible values/states in the sequences
  labels = state_labels, # Labels displayed in plots and tables
  start = 25,            # First position in the sequence corresponds to age 25
  xtstep = 5,            # Display x-axis labels every 5 years
  cpal = state_colors    # Color palette used in sequence plots
  
)

# Notice that the sequence object contains only the trajectory states used for sequence analysis; 
#non-sequence variables such as participant ids are not included.
View(EITC_seq)

# ============================================================
# DISCUSSION QUESTION
# ============================================================
# 1. Why do we distinguish between:
#      - EITC eligible
#      - EITC ineligible because income is too high
#      - EITC ineligible because the individual is not working?
#
# 2. Why did we choose ages 25-45 as the observation window?
#
# 3. Why do we use age (rather than calendar year) as the time scale?
#
# 4. Why do we measure EITC eligibility annually?
#    What alternative time units could be considered?
# ============================================================

# ============================================================
# STEP 1c: Visualize individual trajectories
# ============================================================

# Sequence index plots display one row per participant.
# Each color represents the participant's state at a given age.

seqIplot(
  EITC_seq,
  sortv = "from.start",     # Sort sequences according to their initial state
  with.legend = TRUE,
  cex.legend = 0.8,         # Legend size
  main = "EITC Eligibility Trajectories",
  xlab ="Age (years)"
  )

# ============================================================
#  STEP 1c DISCUSSION QUESTIONS
# ============================================================
#
# 1. What patterns do you notice in the index plot?
#
# 2. What sources of heterogeneity can you identify?
#    Consider:
#      - Timing of EITC eligibility
#      - Duration of eligibility
#      - Transitions between states
#
# 3. What makes it difficult to summarize these trajectories using traditional approaches?
##
# ============================================================

# ============================================================
# STEP 2: Measure trajectory dissimilarity
# ============================================================

# In this step, we compare each participant's trajectory to every
# other participant's trajectory.
#
# The result is a dissimilarity matrix: a square matrix where each cell represents how different two trajectories are.

#Analytical Decision 1: Cost
# For this exercise, we use the classic OM specification:
# substitution cost = 2
# insertion/deletion cost = 1 (one substitution=one insertion + one deletion)
#
# Discussion:
# - Higher indel costs emphasize timing differences
# - Lower indel costs allow sequences to "shift" in time
# - Alternative approaches include transition-rate costs

# For this exercise, choose ONE of the following approaches:

# ------------------------------------------------------------
# OPTION A: Optimal Matching (OM)
# ------------------------------------------------------------

# OM is useful when we are interested in differences in the
# duration of states and overall trajectory patterns.
#
# OM allows:
#   - substitutions
#   - insertions
#   - deletions

costmatrix <- seqsubm(
  EITC_seq,
  method = "CONSTANT",
  time.varying = FALSE
)

dist_om <- seqdist(
  EITC_seq,
  method = "OM",
  sm = costmatrix,
  indel = 1
)

corner(dist_om, n = 15)


# ------------------------------------------------------------
# OPTION B: Hamming Distance
# ------------------------------------------------------------

# Hamming distance is useful when we are interested in the
# timing of states and transitions.
#
# Hamming compares sequences position-by-position and does not
# allow insertions or deletions.

dist_ham <- seqdist(
  EITC_seq,
  method = "HAM",
  sm = costmatrix
)

corner(dist_ham, n = 15)

# ------------------------------------------------------------
# CHOOSE YOUR DISTANCE MATRIX
# ------------------------------------------------------------

# For the rest of the exercise, choose ONE distance matrix.
# Below, uncomment the line corresponding to your choice.

dist_seq <- dist_om   # Use this if choosing Optimal Matching
#dist_seq <- dist_ham  # Use this if choosing Hamming distance

# ------------------------------------------------------------
# STEP 2 DISCUSSION QUESTIONS
# ------------------------------------------------------------

#What do the values in the distance matrix represent?
  #What does a distance of 0 mean?
  #What do larger vs. smaller distances mean?
  
#How does the distance matrix help us move from visualization to clustering?

# ------------------------------------------------------------

#The distance matrix is the main output of sequence analysis. 
#It summarizes the similarity between every pair of trajectories and serves as the input for cluster analysis.

# ============================================================
# STEP 3: Cluster the sequences
# ============================================================

# The distance matrix summarizes how different each pair of trajectories is. 
#We now use that distance matrix to identify groups of participants with similar EITC eligibility trajectories.

# In this exercise, we compare two commonly used clustering approaches:
#
#   1. Hierarchical agglomerative clustering
#   2. Partitioning Around Medoids (PAM)
#
# We will compare cluster quality across both approaches before choosing a final solution.

# ------------------------------------------------------------
# STEP 3A: Hierarchical agglomerative clustering
# ------------------------------------------------------------

# Ward's method starts with each participant as their own cluster
# and progressively merges similar clusters.

cluster_hier <- hclust(
  as.dist(dist_seq),
  method = "ward.D2"
)

# Assess cluster quality across different numbers of clusters

hier_quality <- as.clustrange(
  cluster_hier,
  diss = as.dist(dist_seq)
)

summary(hier_quality, max.rank = 10)

# Plot cluster quality indicators

plot(
  hier_quality,
  stat = c("ASW", "HC", "CH"),
  norm = "zscore",
  col = c("#6666ff", "#cc0000", "#008000"),
  main = "Cluster Quality: Hierarchical Clustering"
)

# ------------------------------------------------------------
# STEP 3B: Partitioning Around Medoids (PAM)
# ------------------------------------------------------------

# PAM identifies representative trajectories, called medoids,
# and assigns participants to the closest medoid.
#
# Unlike hierarchical clustering, PAM requires us to evaluate
# a range of possible numbers of clusters.

pam_quality <- wcKMedRange(
  dist_seq,
  kvals = 2:20
)

summary(pam_quality, max.rank = 10)

# Plot cluster quality indicators

plot(
  pam_quality,
  stat = c("ASW", "HC", "CH"),
  norm = "zscore",
  lwd = 2,
  cex = 2,
  col = c("#6666ff", "#cc0000", "#008000"),
  legendpos = "topright",
  main = "Cluster Quality: PAM"
)

# ============================================================
# DISCUSSION QUESTIONS
# ============================================================
#
# 1. Which clustering algorithm appears to perform better?
#    Consider ASW, HC, and CH.
#
# 2. Do all cluster quality indicators point to the same number
#    of clusters?

# ============================================================
# STEP 3C: Choose and visualize a cluster solution
# ============================================================

# Based on the cluster quality indicators and substantive interpretability, choose a number of clusters.
# Cluster quality indicators may suggest relatively simple solutions (e.g., 3–5 clusters). 
# However, solutions with more clusters may capture additional meaningful heterogeneity in EITC eligibility trajectories.
#
# As in many sequence analyses, the final choice should balance statistical quality and substantive interpretability.

n_clusters <- 10   # Replace with your preferred solution

# Option 1: Hierarchical clustering

cluster_final_hier <- cutree(
  cluster_hier,
  k = n_clusters
)

# Option 2: PAM clustering

cluster_final_pam <- pam(
  as.dist(dist_seq),
  k = n_clusters,
  diss = TRUE
)$clustering

# Visualize the resulting cluster solution
# Replace cluster_final_pam with cluster_final_hier if preferred

seqIplot(
  EITC_seq,
  group = cluster_final_pam,
  sortv = "from.start",
  with.legend = TRUE,
  cex.legend = 0.8,
  main = paste("PAM solution:", n_clusters, "clusters"),
  xlab = "Age (years)"
)


# ============================================================
# DISCUSSION QUESTIONS
# ============================================================
#
# 1. Which number of clusters would you recommend?
#
# 2. Compare several cluster solutions (e.g., 5, 8, 10).
#    Which solution provides the best balance between
#    simplicity and meaningful heterogeneity?
#
# 3. What trajectory patterns do you observe in your preferred cluster solution?
#
# 4. Can you come up with cluster names that describe the different patterns?
# ============================================================

# ------------------------------------------------------------
# EXTENSION 1: Compare dissimilarity measures
# ------------------------------------------------------------

# In this exercise, we used Optimal Matching (OM), which captures
# overall trajectory patterns and differences in duration.
#
# To assess the sensitivity of your results, return to STEP 2
# and repeat the analysis using Hamming distance.

# Compare:
#   - Cluster quality indicators
#   - Preferred number of clusters
#   - Resulting trajectory groups

# DISCUSSION QUESTIONS:
#
# 1. How did the cluster solutions change?
# 2. Which trajectories were most affected?
# 3. Which measure seems most appropriate for our research question?

# ------------------------------------------------------------
# EXTENSION 2: Alternative cluster visualizations
# ------------------------------------------------------------

# Index plots show individual trajectories, but other visualizations
# can make cluster interpretation easier.

# State distribution plot:
# Shows the proportion of individuals in each state at each age.

seqdplot(
  EITC_seq,
  group = cluster_final_pam,
  with.legend = TRUE,
  main = "State Distribution Plots by Cluster",
  xlab = "Age (years)"
)

# Modal state plot:
# Shows the most common state at each age within each cluster.

seqmsplot(
  EITC_seq,
  group = cluster_final_pam,
  with.legend = TRUE,
  main = "Modal State Plots by Cluster",
  xlab = "Age (years)"
)

# DISCUSSION QUESTIONS:
#
# 1. What information is easier to see in the state distribution plots?
# 2. What information is lost in the modal state plots?
# 3. How do these plots help you name and interpret the clusters?
