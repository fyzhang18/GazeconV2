# Data Processing Script Example
  # This script uses a Monte Carlo approach that involves randomly permuting the observed data to create a null distribution. This approach has been described
  # in the EEG and MEG literature as a cluster-based permutation analysis. The specific example below involves doing one-sample t-tests against chance, but the 
  # approach can also been applied for paired- and independent-samples t-tests. This script is an example that is based upon some code we used to analyze for SRCD 2019,
  # so it's currently a little hard-coded.

# Summary
  # We process the data in several steps. First, we format our data so that it can be read into EyetrackingR, a custom package developed by Dink and Ferguson (2015). They
  # have some excellent online examples of how to prepare your data for formatting and even have some code they developed for doing cluster-based permutation. Unfortunately, 
  # their script doesn't allow you to do one-sample t-tests against chance in these permutations, but I was able to modify and elaborate upon some old matlab code that we have previously used to 
  # do these analyses. This is the code below. 

# Use EyetrackingR function to format data so it can be read into the package
Time_Course_Data_3 <- make_eyetrackingr_data(Time_Course_Data_2_flipped, participant_column = "subject_id", trial_column = "trial",
                                             trackloss_column = "Track_Loss", time_column = "time", aoi_columns = c("target","distractor"),
                                             treat_non_aoi_looks_as_missing = TRUE)

# Specify onset and offset time for analysis window
response_window <- subset_by_window(Time_Course_Data_3, 
                                    window_start_time = 0, 
                                    window_end_time = 2800, 
                                    rezero = FALSE)

# Additional data formatting. Here we specify our independent variable (IV). Since it has two levels, we will do our chance comparison separately for each level of the IV. 
response_time <- make_time_sequence_data(response_window, time_bin_size = 1, predictor_columns = "which_item", aois ="target", summarize_by = "subject_id")


# Permutations
# Initialize variables to automate the process and create functions that will be necessary
absmax <- function(x) {x[which.max( abs(x) )]}

# This indicates our variable of interest
current_item <- "not-to-be-changed"

# This is the predictor column if we are doing a within- or between-subject comparison
predictor_column <- "which_item"

# This is the number of iterations we want for the permutation
niters = 1000

# This is a variable that will contain each null value from the permutation. In this script, each value in this null distribution 
# will contain the sum of the largest cluster detected in our permutation.
null_distribution <- vector(mode = "numeric", length = niters)

# Make wide subject-weighted data
data_wide <- dcast(response_time, subject_id + which_item ~ TimeBin, value.var = "Prop")

# Subset out data from other trial type so we are only looking at our variable of interest (for one-sample t-tests against chance)
data_wide_current_item <- data_wide %>%
  filter(which_item == current_item)

# Then isolate that data
data_wide_current_item_1 <- data_wide_current_item %>%
  select(-subject_id,-which_item)

# Make wide trial-weighted data
time_sequence_data <- make_time_sequence_data(Time_Course_Data_3, time_bin_size = 1, predictor_columns = c("which_item","subject_id","trial"),
                                              aois = "Target",summarize_by = "subject_id")
time_sequence_data_1 <- dcast(time_sequence_data, subject_id + which_item + trial ~ TimeBin, value.var = "Prop")

# Subset out data from other trial type and then isolate that data
time_sequence_data_current_item <- time_sequence_data_1 %>%
  filter(which_item == current_item) %>%
  select(-subject_id,-which_item,-trial)

# Pull out the identifying information
identifying_information <- time_sequence_data_1 %>% 
  filter(which_item == current_item) %>%
  select(subject_id,which_item,trial)

# Select column indices for permutation. Our analysis window starts 1000 ms after the onset of the trial and ends at the end of the trial.
x <- which(colnames(time_sequence_data_current_item) == "1000")
y <- which(colnames(time_sequence_data_current_item) == "2783")

# Pull out data for analysis window
current_data_subset <- time_sequence_data_current_item[,x:y]
current_data <- cbind(identifying_information,current_data_subset)

# Prepare for loops and initialize variables for the permutation
# Calculate number of rows
n = nrow(time_sequence_data_current_item)
# Set some p-values
num_sub <- length(unique(identifying_information$subject_id))
threshold_t = qt(p = 1 - .05/2,
                 df = num_sub-1)

# Run the for loop
for(h in 1:niters){
  shuffled_data <- current_data
  z <- ncol(shuffled_data)
  
  # Flip a coin. If heads, we will flip the 0's and 1's for all rows in a given trial. If trails, they will remain the same
  shuffled_data$RandomNumbers <- sample(c(0,1), n, rep = T)
  
  # For loop for shuffling data. We flip the values contained in a given row (e.g., one trial) if the coin comes up heads, and keep the same if the coin comes up tails. 
  for(i in 1:n){
    if(shuffled_data$RandomNumbers[i] == 1){
      shuffled_data[i,4:z] <- replace(shuffled_data[i,4:z],shuffled_data[i,4:z] == 1, 4)
      shuffled_data[i,4:z] <- replace(shuffled_data[i,4:z],shuffled_data[i,4:z] == 0, 1)
      shuffled_data[i,4:z] <- replace(shuffled_data[i,4:z],shuffled_data[i,4:z] == 4, 0)}
    else if(shuffled_data$RandomNumbers[i] == 0){
      shuffled_data[i,4:z] <- shuffled_data[i,4:z]}
  }
  
  # Aggregate by subject
  subject_weighted_mean <- shuffled_data %>%
    group_by(subject_id) %>%
    summarize_at(c(4:z),mean,na.rm = TRUE)
  
  # T-Test on each column of data
  results_permuted <- lapply(subject_weighted_mean[,2:109],t.test,alternative = c("two.sided"), mu = .50)
  
  # Call these results up
  resultsmatrix_permuted <- do.call(cbind,results_permuted)
  
  # Call T-Test Results
  t_tests_clusters_permuted <- as.data.frame(resultsmatrix_permuted[c("statistic"),] > threshold_t)
  t_tests_clusters_1_permuted <- as.data.frame(resultsmatrix_permuted[c("statistic"),])
  t_tests_clusters_2_permuted <- as.data.frame(resultsmatrix_permuted[c("p.value"),] < .025) 
  
  # Create the results that we will use to calculate cluster size and clean them up
  new_data_permuted <- gather(t_tests_clusters_1_permuted, time, t_statistic)
  new_data_1_permuted <- gather(t_tests_clusters_2_permuted, time_2, p.value)
  new_data_1_permuted <- cbind(new_data_permuted,new_data_1_permuted)
  new_data_1_permuted$time = as.numeric(gsub("\\X","",new_data_1_permuted$time))
  
  new_data_2_permuted <- new_data_1_permuted %>%
    select(t_statistic,p.value,time)
  new_data_3_permuted <- cbind(new_data_2_permuted,t_tests_clusters_permuted)
  new_data_3_permuted$cluster_sum <- vector(mode = "numeric",length = length(new_data_3_permuted$time))
  new_data_4_permuted <- new_data_3_permuted %>%
    select(t_statistic,time,cluster_sum,p.value) %>%
    rename(cluster = p.value)
  
  # Calculate the cluster sizes for the permuted data
  for(i in 1:length(new_data_4_permuted$time)){
    if(i == 1 & new_data_4_permuted$cluster[i] == TRUE)
      new_data_4_permuted$cluster_sum[i] <- new_data_4_permuted$t_statistic[i]
    else if(i == 1 & new_data_4_permuted$cluster[i] != TRUE)
      new_data_4_permuted$cluster_sum[i] <- 0
    else if(new_data_4_permuted$cluster[i-1] == "TRUE" & new_data_4_permuted$cluster[i] == TRUE)
      new_data_4_permuted$cluster_sum[i] <- new_data_4_permuted$t_statistic[i] + new_data_4_permuted$cluster_sum[i-1]
    else if(new_data_4_permuted$cluster[i] == "TRUE" & new_data_4_permuted$cluster[i-1] == FALSE)
      new_data_4_permuted$cluster_sum[i] <- new_data_4_permuted$t_statistic[i]
    else if(new_data_4_permuted$cluster[i-1] == "FALSE")
      new_data_4_permuted$cluster_sum[i] <- 0
  }
  
  null_distribution[h] <- absmax(new_data_4_permuted$cluster_sum)
  print(sprintf("Iteration: %d", h))
}

write.csv(null_distribution, file = "non_lookers_null_distribution.csv")

# T-Test
results_permuted <- lapply(data_wide_current_item_1[,61:168],t.test,alternative = c("two.sided"), mu = .50)

# Call these results up
resultsmatrix_permuted <- do.call(cbind,results_permuted)

# Call T-Test Results
t_tests_clusters_permuted <- as.data.frame(resultsmatrix_permuted[c("statistic"),] > threshold_t)
t_tests_clusters_1_permuted <- as.data.frame(resultsmatrix_permuted[c("statistic"),])
t_tests_clusters_2_permuted <- as.data.frame(resultsmatrix_permuted[c("p.value"),] < .025) 

# Create the results that we will use to calculate cluster size and clean them up
new_data_permuted <- gather(t_tests_clusters_1_permuted, time, t_statistic)
new_data_1_permuted <- gather(t_tests_clusters_2_permuted, time_2, p.value)
new_data_1_permuted <- cbind(new_data_permuted,new_data_1_permuted)
new_data_1_permuted$time = as.numeric(gsub("\\X","",new_data_1_permuted$time))

new_data_2_permuted <- new_data_1_permuted %>%
  select(t_statistic,p.value,time)
new_data_3_permuted <- cbind(new_data_2_permuted,t_tests_clusters_permuted)
new_data_3_permuted$cluster_sum <- vector(mode = "numeric",length = length(new_data_3_permuted$time))
new_data_4_permuted <- new_data_3_permuted %>%
  select(t_statistic,time,cluster_sum,p.value) %>%
  rename(cluster = p.value)

# Calculate the cluster sizes for the permuted data
for(i in 1:length(new_data_4_permuted$time)){
  if(i == 1 & new_data_4_permuted$cluster[i] == TRUE)
    new_data_4_permuted$cluster_sum[i] <- new_data_4_permuted$t_statistic[i]
  else if(i == 1 & new_data_4_permuted$cluster[i] != TRUE)
    new_data_4_permuted$cluster_sum[i] <- 0
  else if(new_data_4_permuted$cluster[i-1] == "TRUE" & new_data_4_permuted$cluster[i] == TRUE)
    new_data_4_permuted$cluster_sum[i] <- new_data_4_permuted$t_statistic[i] + new_data_4_permuted$cluster_sum[i-1]
  else if(new_data_4_permuted$cluster[i] == "TRUE" & new_data_4_permuted$cluster[i-1] == FALSE)
    new_data_4_permuted$cluster_sum[i] <- new_data_4_permuted$t_statistic[i]
  else if(new_data_4_permuted$cluster[i-1] == "FALSE")
    new_data_4_permuted$cluster_sum[i] <- 0
}

# Calculate the cluster sizes for the permuted data
for(i in 1:length(new_data_4_permuted$time)){
  if(i == 1 & new_data_4_permuted$cluster[i] == TRUE)
    new_data_4_permuted$cluster_sum[i] <- new_data_4_permuted$t_statistic[i]
  else if(i == 1 & new_data_4_permuted$cluster[i] != TRUE)
    new_data_4_permuted$cluster_sum[i] <- 0
  else if(new_data_4_permuted$cluster[i-1] == "TRUE" & new_data_4_permuted$cluster[i] == TRUE)
    new_data_4_permuted$cluster_sum[i] <- new_data_4_permuted$t_statistic[i] + new_data_4_permuted$cluster_sum[i-1]
  else if(new_data_4_permuted$cluster[i] == "TRUE" & new_data_4_permuted$cluster[i-1] == FALSE)
    new_data_4_permuted$cluster_sum[i] <- new_data_4_permuted$t_statistic[i]
  else if(new_data_4_permuted$cluster[i-1] == "FALSE")
    new_data_4_permuted$cluster_sum[i] <- 0
}          

# Pull Max Cluster Sizes
for(i in 1:length(new_data_4_permuted$time)){
  if(new_data_4_permuted$cluster_sum[i] != 0 & new_data_4_permuted$cluster_sum[i+1] == 0)
    new_data_4_permuted$cluster_sizes[i] <- new_data_4_permuted$cluster_sum[i]
  else if(new_data_4_permuted$cluster_sum[i] != 0 & new_data_4_permuted$cluster_sum[i+1] != 0)
    new_data_4_permuted$cluster_sizes[i] <- "FALSE"
  else if(new_data_4_permuted$cluster_sum == 0)
    new_data_4_permuted$cluster_sizes[i] <- "FALSE"
  else if(i == length(new_data_4_permuted$time) & new_data_4_permuted$cluster_sum[i] != 0)
    new_data_4_permuted$cluster_sizes[i] <- new_data_4_permuted$cluster_sum[i]
}

new_data_4_permuted$cluster_sizes[i] <- new_data_4_permuted$cluster_sum[i]

Clusters <- new_data_4_permuted %>%
  select(cluster_sizes) %>%
  filter(cluster_sizes != FALSE) 

# Mark the trial onset
Count_Onset = 1
for(i in 1:length(new_data_4_permuted$time)){
  if(new_data_4_permuted$cluster[i] == TRUE & i == 1){
    Clusters$Onset[1] <- new_data_4_permuted$time[i]
    Count_Onset = Count_Onset + 1}
  else if(new_data_4_permuted$cluster[i] == TRUE & new_data_4_permuted$cluster[i+1] == FALSE & new_data_4_permuted$cluster[i-1] != TRUE ){
    Clusters$Onset[Count_Onset] <- new_data_4_permuted$time[i]
    Count_Onset = Count_Onset + 1}
  else if(new_data_4_permuted$cluster[i] == TRUE & new_data_4_permuted$cluster[i-1] == FALSE & i != 1){
    Clusters$Onset[Count_Onset] <- new_data_4_permuted$time[i]
    Count_Onset = Count_Onset + 1}
  else if(new_data_4_permuted$cluster[i] == TRUE & new_data_4_permuted$cluster[i+1] == TRUE & i != 1){
    new_data_4_permuted$time[i] <- new_data_4_permuted$time[i]}
}

# Mark the trial offset
Count_Offset = 1
for(i in 1:length(new_data_4_permuted$time)){
  if(new_data_4_permuted$cluster[i] == TRUE & i == 1){
    new_data_4_permuted$time[i] <- new_data_4_permuted$time[i]}
  else if( new_data_4_permuted$cluster[i-1] == FALSE & new_data_4_permuted$cluster[i] == TRUE & new_data_4_permuted$cluster[i+1] == FALSE){
    Clusters$Offset[Count_Offset] <- new_data_4_permuted$time[i] + 16.67
    Count_Offset = Count_Offset + 1}
  else if(new_data_4_permuted$cluster[i-1] == TRUE & new_data_4_permuted$cluster[i] == TRUE & new_data_4_permuted$cluster[i+1] == FALSE & i != 1){
    Clusters$Offset[Count_Offset] <- new_data_4_permuted$time[i] + 16.67
    Count_Offset = Count_Offset + 1}
  else if(new_data_4_permuted$cluster[i] == TRUE & new_data_4_permuted$cluster[i+1] == TRUE & i != 1){
    new_data_4_permuted$time[i] <- new_data_4_permuted$time[i]}
  else if(i == length(new_data_4_permuted$time) & new_data_4_permuted$cluster[i] == TRUE){
    Clusters$Offset[Count_Offset] <- new_data_4_permuted$time[i] + 16.67}
}
Clusters$Offset[Count_Offset] <- new_data_4_permuted$time[i] + 16.67


Clusters$cluster_sizes <-as.numeric(Clusters$cluster_sizes)
null_distribution <- read.csv(file = "non_lookers_null_distribution.csv")

k = 1
for(k in 1:length(Clusters$cluster_sizes)){
  doesnt_exceed_how_many_null_values <- sum(abs(Clusters$cluster_sizes[k]) < abs(null_distribution))
  Clusters$p.value[k] = (doesnt_exceed_how_many_null_values+1)/1001
  doesnt_exceed_how_many_null_values <- NULL
}

write.csv(file = "non_lookers_clusters",Clusters)

# Calculate confidence intervals, means, etc for time-course plot
lookers_CI <- as.data.frame(resultsmatrix_permuted[c("conf.int"),])
lookers_mean <- as.data.frame(resultsmatrix_permuted[c("estimate"),])
lookers_CI_and_mean <- rbind(lookers_CI,lookers_mean)
lookers_CI_and_mean$Measure <- c("Lower_CI","Upper_CI","Mean")
lookers_plot_data <- melt(lookers_CI_and_mean)
lookers_plot_data_2 <- dcast(lookers_plot_data, variable ~ Measure, value.var = "value")
lookers_plot_data_2$Time <- as.numeric(gsub("X","", lookers_plot_data_2$variable))
lookers_plot_data_2$Change <- "Yes"

# First time-course plot
one <-   ggplot(lookers_plot_data_2, aes(Time, group =1)) + 
  geom_line(aes(y=Mean), colour="blue") + 
  geom_ribbon(aes(ymin=Lower_CI, ymax=Upper_CI), alpha=0.2, fill= "blue")+
  annotate("rect",xmin = 1000, xmax = 1450, ymin = -1, ymax = 2, fill = "dark grey", alpha = .5)+
  annotate("rect",xmin = 2050, xmax = 2533,  ymin = -1, ymax = 2, fill = "dark grey", alpha = .5)+
  annotate("rect",xmin = 2717, xmax = 2800, ymin = -1, ymax = 2, fill = "dark grey", alpha = .5)+
  annotate("text",x = 1225, y = 0.9, label = "***", size = 12)+
  annotate("text",x = 2291.5, y = 0.9, label = "**", size = 12)+
  annotate("text",x = 2760, y = 0.9, label = "*", size = 12)+
  geom_hline(aes(yintercept = .50), color = "black")+      
  geom_vline(aes(xintercept = 800), color = "black", linetype="solid")+
  coord_cartesian(ylim = c(0,1), xlim = c(1000,2900))+xlab("Time (in ms)")+ylab("Preference for Overtly Attended Item")+
  scale_x_continuous(breaks = c(1000,1200,1400,1600,1800,2000,2200,2400,2600,2800),
                     labels = c(0,200,400,600,800,1000,1200,1400,1600,1800), expand = c(0,0))+
  scale_y_continuous(breaks = c(0,.10,.20,.30,.40,.50,.60,.70,.80,.90,1.0), expand = c(0,0))+
  geom_hline(aes(yintercept = .50), color = "black")+
  geom_vline(aes(xintercept = 1000), color = "black", linetype="solid")+
  theme_classic()+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))
one

# Save plot
ggsave(filename = "Change_Preference_Score_Test_Period_Lookers.pdf", plot = one, width =10, height= 6)
