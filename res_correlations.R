# Load required libraries
library(R.matlab)
library(ggplot2)
library(ggpubr)


#find the average across n=4 subjects
for(subject in c("sub-0002", "sub-0004", "sub-0006", "sub-0007")) {
  # Extract variables
  vertex_res_map_old <- as.numeric(mat_data[[1]])
  vertex_resistance_map <- as.numeric(mat_data[[2]])
  mat_data <- readMat(paste0(subject, "_region_res_cell_column.mat"))
  # Create a data frame
  data <- data.frame(vertex_res_map_old = vertex_res_map_old, vertex_resistance_map = vertex_resistance_map)
  
  #filter zeroes
  data <- data[which(data$vertex_resistance_map != 0),]
  
  subavg_vert_res_old <- c(subavg_vert_res_old, mean(data$vertex_res_map_old))
  subavg_vert_res_new <- c(subavg_vert_res_new, mean(data$vertex_resistance_map))
  
}

###summary stats
mean(subavg_vert_res_new)
min(subavg_vert_res_new)
max(subavg_vert_res_new)

mean(subavg_vert_res_old)
min(subavg_vert_res_old)
max(subavg_vert_res_old)

###repeat all for sub-0002 for visualization purposes
# Load the .mat file
setwd("C:/Users/vik16/Documents/Baillet Lab Manuscript/")
mat_data <- readMat("sub-0002_region_res_cell_column.mat")

subavg_vert_res_old <- c()
subavg_vert_res_new <- c()
# Extract variables
vertex_res_map_old <- as.numeric(mat_data[[1]])
vertex_resistance_map <- as.numeric(mat_data[[2]])

# Create a data frame
data <- data.frame(vertex_res_map_old = vertex_res_map_old, vertex_resistance_map = vertex_resistance_map)

#filter zeroes
data <- data[which(data$vertex_resistance_map != 0),]
# Calculate Spearman correlation
correlation <- cor(data$vertex_res_map_old, data$vertex_resistance_map, method = "spearman")

# Create scatterplot with larger fonts and title
p <- ggplot(data, aes(x = vertex_res_map_old, y = vertex_resistance_map)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "R from GM resistivity in cortical column (Ohms)",
       y = "R from dendrite and neuron density (Ohms)",
       title = "Dendritic vs. Columnar Resistance at Each Vertex, sub-0002") +
  annotate("text", x = max(data$vertex_res_map_old) * 0.95, y = max(data$vertex_resistance_map) * 0.95, 
           label = paste("Spearman r =", round(correlation, 2)), 
           hjust = 1, vjust = 1, size = 6, color = "darkred") +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

# Display the plot
print(p)
