################################# Diversity analysis script ##############################

#Libraries
library(vegan)
library(dplyr)
#Functions
alpha_diversity <- function(dataset){
  
  roles <- unique(dataset$Role)
  
  
  matrix<-c()
  dataset_month<-dataset
  
  
  
  for (role in roles){
    
    dataset_role <-dataset_month[dataset_month$Role==role,]
    
    skill_counts <-dataset_role[,4:length(colnames(dataset_role))]
    
    #final_role<-data.frame(colSums(skill_counts))
    final_role<-skill_counts
    #Sum_all
    N<-sum(final_role)
    ##Shannon Index
    H<-round(diversity(final_role),3)
    ##Simpson Index 1-D
    D<-round(diversity(final_role,"simpson"),3)
    ##Inverse Simpson index 1/D
    invD<-round(diversity(final_role,"invsimpson"),3)
    ##Species Richness
    S <- apply(final_role, 1, function(row) sum(row > 0))
    
    ## Pielous J
    J <- round(H/log(S),3)
    ## Margalef
    R <- round(S/log(N),3)
    
    df_indices <- data.frame(cbind(mean(S),mean(H),mean(D),mean(invD)))
    df_indices$Role<-role
    
    colnames(df_indices)<-c('Richness','Shannon','Simpson','Inverse_Simpson','Role')
    
    
    if (length(matrix)==0){
      matrix<-df_indices
    }else{
      matrix <-rbind(matrix,df_indices)
      
    }
    
    
    
    
  }
  return(matrix)
  
  
}

#Main code

###### load the binary datataset of Role-Month skill
load('./sample_dataset.Rda')
dataset_jobs<-df


## Matrix construction
months <- unique(dataset_jobs$Month)
roles <- unique(dataset_jobs$Role)
matrix <- c()
for(role in roles){
  
  data_role <-dataset_jobs[dataset_jobs$Role==role,]
  
  month_counts<-table(data_role$Month)
  
  no_selection<-min(month_counts)
  
  for(month in months){
    
    data_month <- data_role[data_role$Month==month,]
    
    data_month<-sample_n(data_month,no_selection)
    
    skills <-data_month[,3:ncol(data_month)]
    
    numeric_counts <- sapply(skills, function(column) sum(column))
    
    final_role <-data.frame(t(numeric_counts))
    
    colnames(final_role)<-colnames(skills)
    
    final_role$Role <- role
    
    final_role$Month <- month
    
    if (length(matrix)==0){
      matrix <- final_role
    }
    else{
      matrix<-rbind(matrix,final_role)
    }
  }
  
  
  
}
dataset_aggregated<- matrix[, c("Role", "Month",setdiff(names(matrix), c("Role", "Month")))]

jobs_dataset<-dataset_jobs[rowSums(dataset_jobs[3:ncol(dataset_jobs)])!=0,]
dataset_jobs<-jobs_dataset


## Sort the dataset
roles <-unique(dataset_aggregated$Role)
dataset_sorted <-dataset_aggregated[order(dataset_aggregated$Role),]

## Run alpha diversity
res_alpha <-alpha_diversity(dataset_sorted)
res_alpha<- res_alpha[, c("Role",setdiff(names(res_alpha), c("Role")))]


## Calculate the frequencies of the skills
dataset_sorted_freq <-dataset_sorted
row_sums <- rowSums(dataset_sorted_freq[3:ncol(dataset_sorted_freq)])
dataset_sorted_freq[3:ncol(dataset_sorted_freq)]<-dataset_sorted_freq[3:ncol(dataset_sorted_freq)]/row_sums



## Run beta diversity - ordination
bray_dist <- vegdist(dataset_sorted_freq[,3:ncol(dataset_sorted_freq)], method = "bray",na.rm=TRUE)
dataset_sorted_freq$Role<-as.factor(dataset_sorted_freq$Role)
dataset_sorted_freq$Month<-as.factor(dataset_sorted_freq$Month)
#check for main effect of Month (Time)
res_month<-  adonis2( bray_dist ~ Month, data = dataset_sorted_freq, permutations = 999)  
#check for main effects of Role
res_role <-adonis2( bray_dist ~ Role, data = dataset_sorted_freq, permutations = 999)
#Test for Role and Month if both significant 
res<-adonis2( bray_dist ~ Role + Month, data = dataset_sorted_freq, permutations = 999)




## Pairwise test
#devtools::install_github("jeffkimbrel/jakR")
library(jakR)

rownames(dataset_sorted_freq) <-1:nrow(dataset_sorted_freq)
res<-pairwiseAdonis(dataset_sorted_freq[3:ncol(dataset_sorted_freq)],factors=dataset_sorted_freq$Role,"bray")
res<-data.frame(res)


## PCoA
# Perform PCoA (Principal Coordinates Analysis) on Bray-Curtis dissimilarity
bray_dist <- vegdist(dataset_sorted_freq[,3:ncol(dataset_sorted_freq)], method = "bray")
pcoa_result <- cmdscale(bray_dist, k = 2, eig = TRUE)
pcoa_df <- data.frame(PC1 = pcoa_result$points[, 1],
                      PC2 = pcoa_result$points[, 2],
                      Role = dataset_sorted_freq$Role,
                      Month=dataset_sorted_freq$Month)

library(ggplot2)
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Role)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCoA of Ecosystems Based on Bray-Curtis Dissimilarity",
       x = "PC1", y = "PC2")




### ISA analysis
library(indicspecies)
groups <-as.factor(dataset_jobs$Role)
indval <- multipatt(dataset_jobs[3:ncol(dataset_jobs)], groups,duleg=TRUE,
                    control = how(nperm=99)) 
