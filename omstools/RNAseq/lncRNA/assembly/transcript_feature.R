suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(scales))
suppressMessages(library(ggthemes))
suppressMessages(library(dplyr))
suppressMessages(library(argparser))
suppressMessages(library(data.table))
options(stringsAsFactors = FALSE)

type_portion <- function(type_count) {
  colnames(type_count) <- c('Count', 'Type')
  type_count_df <- data.frame(table(type_count))
  type_count_p_tb <- type_count_df %>%
    group_by(Type) %>%
    mutate(Por = Freq / sum(Freq))
  type_count_p_df <- data.frame(type_count_p_tb)
  type_count_p_df$Count <- as.numeric(as.character(type_count_p_df$Count))
  return(type_count_p_df)
}


p <- arg_parser("transcript/gene feature plot")
p <- add_argument(p, "--feature_file_dir", help = "feature data directory")
p <- add_argument(p, "--detail", help = "plot by detail lncRNA type", flag=TRUE)
p <- add_argument(p, "--out_dir", help = "feature plot output directory", default='./')
argv <- parse_args(p)

tr_feature_file <- file.path(argv$feature_file_dir, 'Transcript_feature.txt')
gene_feature_file <- file.path(argv$feature_file_dir, 'Gene_feature.txt')
exon_intron_file <- file.path(argv$feature_file_dir, 'Exon_Intron_feature.txt')

# tr_feature_file <- 'Transcript_feature.txt'
# gene_feature_file <- 'Gene_feature.txt'
# exon_intron_file <- 'Exon_Intron_feature.txt'
# argv$detail <- TRUE

tr_feature <- read.delim(tr_feature_file)
gene_feautre <- read.delim(gene_feature_file)
exon_intron_feature <- fread(exon_intron_file)

set1_left <- c(brewer.pal(9, "Set1")[2:4], brewer.pal(9, "Set1")[6:9])
mypal <- colorRampPalette(set1_left)
colors <- c("#A50026", "orange")
plot_order <- c('protein_coding', 'TUCP')

tr_len_name = 'Processed_Transcript_Size_Distribution'
if ( ! argv$detail ) {
  tr_feature_plot <- mutate(tr_feature, 
                       Transcript_biotype = ifelse(
                         Transcript_biotype %in% c('protein_coding', 'TUCP'),
                         Transcript_biotype,"lncRNA"))
  exon_intron_feature_plot <- mutate(exon_intron_feature, 
                                     Transcript_biotype = ifelse(
                                       Transcript_biotype %in% c('protein_coding', 'TUCP'),
                                       Transcript_biotype, "lncRNA"))
  
  plot_order <- c(plot_order, 'lncRNA')
  colors <- c(colors, "#313695")
  names(colors) <- plot_order
  plot_width = 8
  plot_height = 8
  suffix = ''
} else {
  tr_feature_plot <- tr_feature
  exon_intron_feature_plot <- exon_intron_feature
  lnc_df <- filter(tr_feature, ! Transcript_biotype %in% c('protein_coding', 'TUCP'))
  lnc_type <- sort(unique(lnc_df$Transcript_biotype))
  plot_order <- c(plot_order, lnc_type)
  colors <- c(colors, mypal(length(lnc_type)))
  names(colors) <- plot_order
  plot_unit = sqrt(length(colors))
  plot_width = 3*plot_unit
  plot_height = 3*plot_unit
  suffix <- '_detail'
}

tr_feature_plot$Transcript_biotype <- factor(tr_feature_plot$Transcript_biotype, levels = plot_order)
len_max <- round(quantile(tr_feature_plot$Transcript_length,0.98),-1)


## length distribution
p <- ggplot(tr_feature_plot,aes(Transcript_length,fill = Transcript_biotype)) + 
  geom_density(alpha = 0.8) + xlim(0,len_max) + 
  scale_fill_manual(values=colors) +
  guides(fill=guide_legend(title=NULL)) +
  xlab('Size (bp)') + ylab('Density') + 
  labs(title="Processed Transcript Size Distribution")
if (argv$detail) {
  p <- p + facet_wrap(~Transcript_biotype, scales = "free")
}
ggsave(filename=paste(argv$out_dir, '/', tr_len_name, suffix,'.png', sep=''),
       type = 'cairo-png',plot = p ,width = plot_width, height = plot_width,dpi=300)
ggsave(filename=paste(argv$out_dir, '/', tr_len_name, suffix,'.pdf', sep=''),
       plot = p ,width = plot_width, height = plot_width)

## exon number

tr_exon_df <- tr_feature_plot[,c('Exon_number', 'Transcript_biotype')]
tr_exon_p_plot_df <- type_portion(tr_exon_df)
tr_exon_p_plot_df <- filter(tr_exon_p_plot_df, Count <=15)
exon_max <- max(tr_exon_p_plot_df$Por) + 0.05
tr_exon_p_plot_df$Type <- factor(tr_exon_p_plot_df$Type, levels = plot_order)
exon_plot_name <- 'Exon_per_transcript'

p <- ggplot(tr_exon_p_plot_df, aes(x=Count, y=Por, fill=Type)) +
  geom_bar(position="dodge",stat="identity") + 
  scale_x_continuous(breaks=seq(1,15)) +
  xlab('Number of exon(s)') + 
  ylab('Proportion of transcripts(%)') + 
  labs(title="Exon per transcript")
if (argv$detail) {
  p <- p + facet_wrap(~Type, scales = "free") + 
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values=colors, guide=FALSE)
    
} else {
  p <- p + scale_y_continuous(labels = percent, limits = c(0, exon_max), 
                              breaks = seq(0, exon_max, 0.1)) +
    scale_fill_manual(values=colors) +
    guides(fill=guide_legend(title=NULL)) 
}

ggsave(filename=paste(argv$out_dir, '/', exon_plot_name, suffix,'.png', sep=''),
       type = 'cairo-png', plot = p, width = plot_width, height = plot_width, dpi=300)
ggsave(filename=paste(argv$out_dir, '/', exon_plot_name, suffix,'.pdf', sep=''),
       plot = p, width = plot_width, height = plot_width)

# gene isoform
if ( ! argv$detail) {
  gene_isoform <- gene_feautre[,c('Transcript_number', 'Gene_biotype')]
  gene_isoform_p_df <- type_portion(gene_isoform)
  gene_isoform_p_df$Type <- as.character(gene_isoform_p_df$Type)
  isf_max <- max(gene_isoform_p_df$Por) + 0.05
  gene_isoform_f_p_df <- filter(gene_isoform_p_df, Count <= 10)
  gene_isoform_f_p_df$Type <- factor(gene_isoform_f_p_df$Type, levels = plot_order)
  p <- ggplot(gene_isoform_f_p_df, aes(x=Count, y=Por, fill=Type)) +
    geom_bar(position="dodge",stat="identity") + 
    scale_x_continuous(breaks=seq(1,10)) +
    scale_y_continuous(labels = percent, limits = c(0, isf_max), 
                       breaks = seq(0, isf_max, 0.1)) +
    scale_fill_manual(values=colors) +
    guides(fill=guide_legend(title=NULL)) +
    xlab('Number of isoform(s)') + 
    ylab('Proportion of Genes(%)') + 
    labs(title="Isoform per Gene locus")
  plot_name <- 'Isoform_per_gene_locus'
  ggsave(filename=paste(argv$out_dir, '/', plot_name, suffix,'.png', sep=''),
         type = 'cairo-png', plot = p, width = plot_width, height = plot_width, dpi=300)
  ggsave(filename=paste(argv$out_dir, '/', plot_name, suffix,'.pdf', sep=''),
         plot = p, width = plot_width, height = plot_width)
}

# exon/intron length
exon_intron_feature_plot$feature_id <- paste(exon_intron_feature_plot$Chr, exon_intron_feature_plot$Start,
                                             exon_intron_feature_plot$End, exon_intron_feature_plot$Strand,
                                             sep='_')
exon_intron_plot <- exon_intron_feature_plot[, c('Feature', 'Length', 'feature_id', 'Transcript_biotype')]
exon_intron_rd_plot <- unique( exon_intron_plot[ , 1:4 ] )
exon_intron_rd_plot$Transcript_biotype <- factor(exon_intron_rd_plot$Transcript_biotype,
                                                 levels = plot_order)
len_max <- round(quantile(exon_intron_rd_plot$Length,0.9),-1)
plot_name <- 'Exon_Intron_Size_Distribution'
type_colors <- c("#A50026", "#313695")
  
if ( ! argv$detail) {
  p <- ggplot(exon_intron_rd_plot,aes(Length,fill=Transcript_biotype)) +   
    geom_density(alpha=0.8) + xlim(0, len_max) + 
    scale_fill_manual(values=colors) +
    facet_wrap(~Feature)

} else {
  p <- ggplot(exon_intron_rd_plot,aes(Length,fill=Feature)) +
    geom_density(alpha=0.8) + xlim(0, len_max) + 
    scale_fill_manual(values=type_colors) +
    facet_wrap(~Transcript_biotype)
}
p <- p + guides(fill=guide_legend(title=NULL)) +
  xlab('Size (bp)') + ylab('Density') + 
  labs(title="Exon/Intron Size Distribution")

ggsave(filename=paste(argv$out_dir, '/', plot_name, suffix,'.png', sep=''),
       type = 'cairo-png', plot = p, width = plot_width, height = plot_width, dpi=300)
ggsave(filename=paste(argv$out_dir, '/', plot_name, suffix,'.pdf', sep=''),
       plot = p, width = plot_width, height = plot_width)  
