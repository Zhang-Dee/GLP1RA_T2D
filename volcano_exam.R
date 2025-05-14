library(ggVolcano)

## basic example code
# load the data
data(deg_data)

# use the function -- add_regulate to add a regulate column 
# to the DEG result data. 
data <- add_regulate(deg_data, log2FC_name = "log2FoldChange",
                     fdr_name = "padj",log2FC = 1, fdr = 0.05)

# plot
ggvolcano(data, x = "log2FoldChange", y = "padj", pointSize = 1.5,
          label = "row", label_number = 0, output = FALSE)
