suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(argparser))

p <- arg_parser("gene expression cumulative plot")
p <- add_argument(p, "--exp_file", help = "gene max expression file.")
p <- add_argument(p, "--out_prefix", help = "output prefix")
argv <- parse_args(p)

exp <- read.delim(argv$exp_file)

theme_onmath <- function(base_size = 14) {
  theme_bw() + theme(panel.background = element_blank(), 
                     panel.grid = element_blank(),
                     plot.background = element_blank(), 
                     panel.border = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(color = "black", face = "bold", size = rel(1)), 
                     axis.title = element_text(face = "bold",size = base_size), 
                     axis.title.x = element_text(vjust = -0.2, size = rel(1)),
                     axis.title.y = element_text(angle = 90, vjust = 2, size = rel(1)), 
                     plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5), 
                     legend.key = element_blank(), legend.title = element_text(face = "italic"),
                     strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"), 
                     strip.text = element_text(face = "bold"))
}

exp$log_tpm = log10(exp$tpm)
p <- ggplot(exp, aes(log_tpm, colour = type)) + 
  stat_ecdf() +
  scale_x_continuous(breaks = seq(-1, 4, 1), labels = c(0.1, 1, 10, 100, 1000, 10000)) +
  theme_onmath() + ylab('Cumulative fraction') + xlab('Max (TPM)')


ggsave(paste(argv$out_prefix, 'cumulative_distribution.png', sep='.'), plot = p, width = 8, height = 6,
       dpi = 300, type = "cairo")
ggsave(paste(argv$out_prefix, 'cumulative_distribution.pdf', sep='.'), plot = p, width = 8, height = 6,
       device = cairo_pdf)