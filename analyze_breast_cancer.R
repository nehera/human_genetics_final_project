# Load packages
library(tidyverse)

# Load breast cancer data
data <- load("Final_Project_TCGA_Breast_Data.rdata")

# Standardize gene expression
standardized_matrix <- scale(t(Expression.matrix), scale = TRUE, center = TRUE)

# Let's examine the inter-gene expression correlation
gene_corr_matrix <- cor(standardized_matrix)

# Creating a heatmap
ggplot(data = reshape2::melt(gene_corr_matrix)) +
  geom_tile(aes(x = Var1, y = Var2, fill = value)) +
  viridis::scale_fill_viridis(option = "magma", limits = c(-1, 1)) +
  theme_minimal() +
  labs(title = "Correlation Heatmap of Rows")

# Show correlation distribution
corr_vector <- gene_corr_matrix[lower.tri(gene_corr_matrix)]

# Creating a dataframe with the vector
data <- data.frame(correlation = corr_vector)

# Creating the density plot
ggplot(data, aes(x = correlation)) +
  geom_density() +
  labs(x = "Correlation", y = "Density", title = "Density Plot of Inter-Gene Correlation")

# Creating boxplot
ggplot(data, aes(x = correlation)) +
  geom_boxplot() +
  labs(x = "Correlation", y = "Density", title = "Boxplot of Inter-Gene Correlation")

# Get the quantiles of the correlations
corr_quantiles <- quantile(corr_vector) # 50% of all pairwise correlations are between -.1 and .16
corr_IQR <- IQR(corr_vector)
# corr=.4 or -.4 is along the way to the outlier cutpoint, so seemingly reasonable
corr_quantiles[4] + corr_IQR
