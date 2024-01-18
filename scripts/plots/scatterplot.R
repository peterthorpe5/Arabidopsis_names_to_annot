# Load the required libraries for plotting
library(ggplot2)
library(dplyr)
library(plotly)
library(htmlwidgets)



setwd("/Users/PThorpe001/Library/CloudStorage/OneDrive-UniversityofDundee/ggs_arabidopsis/demeth_mutants_DE/genes/scatter_plot")

# Specify the file paths
file_path1 <- "vir1_low_vs_vir1_high.GLM.edgeR.DE_results"
file_path2 <- "col0_low_vs_col0_high.GLM.edgeR.DE_results"



# Read the data from the first file
data1 <- read.table(file_path1, header = TRUE, sep = "\t", row.names = 1)

# Read the data from the second file
data2 <- read.table(file_path2, header = TRUE, sep = "\t", row.names = 1)


# Merge the datasets based on common gene identifiers using dplyr
merged_data <- inner_join(data1, data2, by = "Row.names", suffix = c("_vir1", "_col0"))


# Create a scatter plot with merged data, a red line of best fit, a density heatmap,
# and dashed black lines at +2 and -2 log-fold change on both axes
ggplot_obj <- ggplot(merged_data, aes(x = logFC_vir1, y = logFC_col0)) +
   geom_point(alpha = 0.6) +
   geom_smooth(method = "lm", se = FALSE, color = "red") +
   geom_bin2d(bins = 200, color = "white") +
   scale_fill_gradient(low = "blue", high = "red") +
   geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "black") +
   geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
   labs(title = "Scatter Plot with Density Heatmap",
        x = "logFC (vir1_high_vs_vir1_low)",
        y = "logFC (col0_low_vs_col0_high)") +
   theme_minimal()


ggsave("scatter_plot_density_heatmap.pdf", ggplot_obj, scale = 1, width = 9, height = 7, 
       dpi = 1200, limitsize = TRUE)



# make an interactive plot so the user can get the names of interest

# Create an interactive scatter plot with hover information using a different color scale
plotly_obj <- plot_ly(merged_data, x = ~logFC_vir1, y = ~logFC_col0, type = 'scatter', mode = 'markers',
                      marker = list(color = ~logFC_vir1, colorscale = 'Jet'))  # Change 'Jet' to another scale

# Add gene names as hover text
plotly_obj <- plotly_obj %>%
  add_trace(text = ~Row.names, hoverinfo = "text")

# Customize the color scale
color_scale <- 'Jet'
plotly_obj$marker$colorscale <- color_scale

# Add dashed lines at +- 2.0 log-fold change on both the x and y-axes using layout
plotly_obj <- plotly_obj %>%
  layout(shapes = list(
    list(type = 'line', x0 = -2, x1 = -2, y0 = -Inf, y1 = Inf, line = list(color = 'black', dash = 'dash')),
    list(type = 'line', x0 = 2, x1 = 2, y0 = -Inf, y1 = Inf, line = list(color = 'black', dash = 'dash')),
    list(type = 'line', x0 = -Inf, x1 = Inf, y0 = -2, y1 = -2, line = list(color = 'black', dash = 'dash')),
    list(type = 'line', x0 = -Inf, x1 = Inf, y0 = 2, y1 = 2, line = list(color = 'black', dash = 'dash'))
  ))

# Save the Plotly interactive plot as an HTML file
saveWidget(plotly_obj, file = "interactive_scatter_plot_density_heatmap.html", selfcontained = TRUE)


# Convert ggplot to a Plotly interactive plot
plotly_obj <- ggplotly(ggplot_obj)

# Save the Plotly interactive plot as an HTML file
saveWidget(plotly_obj, file = "ggplot_interactive_scatter_plot_density_heatmap.html")


#######################################


# Specify the file paths
file_path1 <- "vir1_low_vs_vir1_high.GLM.edgeR.DE_results"
file_path2 <- "col0_low_vs_col0_high.GLM.edgeR.DE_results"



# Read the data from the first file
data1 <- read.table(file_path1, header = TRUE, sep = "\t", row.names = 1)

# Read the data from the second file
data2 <- read.table(file_path2, header = TRUE, sep = "\t", row.names = 1)


# Merge the datasets based on common gene identifiers using dplyr
merged_data <- inner_join(data1, data2, by = "Row.names", suffix = c("_vir1", "_col0"))



# Create an interactive scatter plot with hover information using a different color scale
plotly_obj <- plot_ly(merged_data, x = ~logFC_vir1, y = ~logFC_col0, 
                      type = 'scatter', mode = 'markers',
                      marker = list(color = ~ifelse(logFC_vir1 <= -2 | logFC_vir1 >= 2 | logFC_col0 <= -2 | logFC_col0 >= 2, 'rgba(255, 0, 0, 0.7)', 'rgba(169,169,169,0.8)'),
                                    colorscale = 'Jet'))  # Change 'Jet' to another scale

# Add gene names as hover text
plotly_obj <- plotly_obj %>%
  add_trace(text = ~Row.names, hoverinfo = "text")

# Customize the color scale
color_scale <- 'Jet'
plotly_obj$marker$colorscale <- color_scale

# Add dashed lines at +- 2.0 log-fold change on both the x and y-axes using layout
plotly_obj <- plotly_obj %>%
  layout(shapes = list(
    list(type = 'line', x0 = -2, x1 = -2, y0 = -Inf, y1 = Inf, 
         line = list(color = 'grey', dash = 'dash')),
    list(type = 'line', x0 = 2, x1 = 2, y0 = -Inf, y1 = Inf, 
         line = list(color = 'grey', dash = 'dash')),
    list(type = 'line', x0 = -Inf, x1 = Inf, y0 = -2, y1 = -2, 
         line = list(color = 'grey', dash = 'dash')),
    list(type = 'line', x0 = -Inf, x1 = Inf, y0 = 2, y1 = 2, 
         line = list(color = 'grey', dash = 'dash'))
  ))

# Save the Plotly interactive plot as an HTML file
saveWidget(plotly_obj, file = "interactive_scatter_plot_density_heatmap.html", selfcontained = TRUE)


##########################################################################################
# add annotation


annot_data <- read.table("annot_modified.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


annot_data_short <- read.table("annot", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


#  'gene' is the gene identifier in your annotation data
merged_data <- merge(merged_data, annot_data, by.x = "Row.names", by.y = "gene", all.x = TRUE)

# Create an interactive scatter plot with hover information using a different color scale
plotly_obj <- plot_ly(merged_data, x = ~logFC_vir1, y = ~logFC_col0, type = 'scatter', mode = 'markers',
                      marker = list(color = ~ifelse(logFC_vir1 <= -2 | logFC_vir1 >= 2 | logFC_col0 <= -2 | logFC_col0 >= 2, 'rgba(255, 0, 0, 0.9)', 'rgba(230, 40, 40,0.99)'),
                                    colorscale = 'Jet',
                                    size = 3.5))  # Adjust size as needed

# Function to split a string into multiple lines of a specified length
split_text <- function(text, length_limit) {
  ifelse(nchar(text) <= length_limit, text, sapply(strwrap(text, width = length_limit, simplify = FALSE), paste, collapse = "<br>"))
}


# Add gene names and annotation information as hover text (substr(full_annot, 1, 250)) -  to truncate)
plotly_obj <- plotly_obj %>%
  add_trace(text = ~paste("Gene: ", Row.names, "<br>Annotation: ", split_text(full_annot, 100)) , 
            hoverinfo = "text",hovertemplate = "%{text}<extra></extra>")


# Customize the color scale
color_scale <- 'Jet'
plotly_obj$marker$colorscale <- color_scale


# Customize the layout to set the maximum width of hover text
plotly_obj <- plotly_obj %>%
  layout(
    hoverlabel = list(
      bgcolor = "white",  # Set background color for better readability
      font = list(size = 12),  # Set font size for better readability
      align = "centre",  # Align text to the left for multiline support
      bordercolor = "black",  # Set border color
      borderwidth = 1,  # Set border width
      namelength = -1  # Allow unlimited length for the text
    )
  )


# Add dashed lines at +- 2.0 log-fold change on both the x and y-axes using layout
plotly_obj <- plotly_obj %>%
  layout(shapes = list(
    list(type = 'line', x0 = -2, x1 = -2, y0 = -Inf, y1 = Inf, line = list(color = 'black', dash = 'dash')),
    list(type = 'line', x0 = 2, x1 = 2, y0 = -Inf, y1 = Inf, line = list(color = 'black', dash = 'dash')),
    list(type = 'line', x0 = -Inf, x1 = Inf, y0 = -2, y1 = -2, line = list(color = 'black', dash = 'dash')),
    list(type = 'line', x0 = -Inf, x1 = Inf, y0 = 2, y1 = 2, line = list(color = 'black', dash = 'dash'))
  ))



# Set axis labels
plotly_obj <- plotly_obj %>%
  layout(title = "Scatter Plot Col0 vs vir1 at 17 and 27 Degrees",
         xaxis = list(title = "logFC (vir1_low_vs_vir1_high)"),
         yaxis = list(title = "logFC (col0_low_vs_col0_high)"))


# Save the Plotly interactive plot as an HTML file
saveWidget(plotly_obj, file = "interactive_scatter_FULL_ANNOTATION_plot_density_heatmap.html", selfcontained = TRUE)





##########################################################################################
# add annotation


annot_data <- read.table("annot_modified.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


annot_data_short <- read.table("annot", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


#  'gene' is the gene identifier in your annotation data
merged_data <- merge(merged_data, annot_data, by.x = "Row.names", by.y = "gene", all.x = TRUE)

# Create an interactive scatter plot with hover information using a different color scale
plotly_obj <- plot_ly(merged_data, x = ~logFC_vir1, y = ~logFC_col0, type = 'scatter', mode = 'markers',
                      marker = list(color = ~ifelse(logFC_vir1 <= -2 | logFC_vir1 >= 2 | logFC_col0 <= -2 | logFC_col0 >= 2, 'rgba(255, 0, 0, 0.9)', 'rgba(230, 40, 40,0.99)'),
                                    colorscale = 'Jet',
                                    size = 3.5))  # Adjust size as needed

# Function to split a string into multiple lines of a specified length
split_text <- function(text, length_limit) {
  ifelse(nchar(text) <= length_limit, text, sapply(strwrap(text, width = length_limit, simplify = FALSE), paste, collapse = "<br>"))
}


# Add gene names and annotation information as hover text (substr(full_annot, 1, 250)) -  to truncate)
plotly_obj <- plotly_obj %>%
  add_trace(text = ~paste("Gene: ", Row.names, "<br>Annot: ", annot, "<br>Annotation: ", split_text(full_annot, 100)) , 
            hoverinfo = "text",hovertemplate = "%{text}<extra></extra>")


# Customize the color scale
color_scale <- 'Jet'
plotly_obj$marker$colorscale <- color_scale


# Customize the layout to set the maximum width of hover text
plotly_obj <- plotly_obj %>%
  layout(
    hoverlabel = list(
      bgcolor = "white",  # Set background color for better readability
      font = list(size = 12),  # Set font size for better readability
      align = "centre",  # Align text to the left for multiline support
      bordercolor = "black",  # Set border color
      borderwidth = 1,  # Set border width
      namelength = -1  # Allow unlimited length for the text
    )
  )


# Add dashed lines at +- 2.0 log-fold change on both the x and y-axes using layout
plotly_obj <- plotly_obj %>%
  layout(shapes = list(
    list(type = 'line', x0 = -2, x1 = -2, y0 = -Inf, y1 = Inf, line = list(color = 'black', dash = 'dash')),
    list(type = 'line', x0 = 2, x1 = 2, y0 = -Inf, y1 = Inf, line = list(color = 'black', dash = 'dash')),
    list(type = 'line', x0 = -Inf, x1 = Inf, y0 = -2, y1 = -2, line = list(color = 'black', dash = 'dash')),
    list(type = 'line', x0 = -Inf, x1 = Inf, y0 = 2, y1 = 2, line = list(color = 'black', dash = 'dash'))
  ))



# Set axis labels
plotly_obj <- plotly_obj %>%
  layout(title = "Scatter Plot Col0 vs vir1 at 17 and 27 Degrees",
         xaxis = list(title = "logFC (vir1_low_vs_vir1_high)"),
         yaxis = list(title = "logFC (col0_low_vs_col0_high)"))


# Save the Plotly interactive plot as an HTML file
saveWidget(plotly_obj, file = "interactive_scatter_FULL_ANNOTATION_plot.html", selfcontained = TRUE)

