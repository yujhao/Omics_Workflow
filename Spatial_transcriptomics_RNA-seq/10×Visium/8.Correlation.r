library(ggplot2)
library(dplyr)

scoredata = read.delim("gene_module_scores.xls", row.names = 1)
scoredata = scoredata[, -c(1:3)]
g = c("FN1","COL4A1","COL1A2","COL6A2","THBS1")
exp = FetchData(data_ob, vars = g, slot = "scale.data")
df = cbind(data_ob@meta.data, exp, scoredata)
dir.create("corplot", showWarnings = FALSE)
setwd("corplot")

for(x in unique(df$group)){
for(y in unique(df$clusters)){
for(i in g){
for(j in colnames(scoredata)){
        d = df[df$group == x & df$clusters == y, ]
        r = cor.test(d[,i], d[,j], method = "spearman")
        p = ggplot(d, aes(x = .data[[i]], y = .data[[j]])) +
            geom_point(color = "#ff4040", alpha = 0.4, size = 0.8) +
            geom_smooth(method = "lm", se = FALSE, color = "black", linetype = 2, linewidth = 0.8) +
            annotate("text",
                    x = Inf, y = Inf,
                    label = paste0("R=", round(r$estimate, 3), "\nP=", formatC(r$p.value, format = "e", digits = 2)),
                    hjust = 1.1, vjust = 1.5) +
            xlab(i) + ylab(j) +
            ggtitle(paste(x, y, sep = "_")) +
            theme_bw() +
            theme(
            plot.title = element_text(hjust = 0.5),
            panel.grid = element_blank()
            ) +
            coord_cartesian(clip = "off")
        dir.create(file.path("plot", x), showWarnings = FALSE)
        ggsave(file.path("plot", x, paste0(y, "_", i, "_", j, ".pdf")), p, width = 4, height = 4)
        ggsave(file.path("plot", x, paste0(y, "_", i, "_", j, ".png")), p, width = 4, height = 4, dpi = 600)
}}}}

