library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# Helper function to read enrichment results
read_tag <- function(path, label, direction) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    df$Category <- label
    df$Direction <- direction
    df
}

# ============================================================
# GO BARPLOTS
# ============================================================

go_files <- list(
    list(path = "enrichment_ORA/GO_BP_UP.csv", label = "BP", dir = "Up"),
    list(path = "enrichment_ORA/GO_BP_DOWN.csv", label = "BP", dir = "Down"),
    list(path = "enrichment_ORA/GO_CC_UP.csv", label = "CC", dir = "Up"),
    list(path = "enrichment_ORA/GO_CC_DOWN.csv", label = "CC", dir = "Down"),
    list(path = "enrichment_ORA/GO_MF_UP.csv", label = "MF", dir = "Up"),
    list(path = "enrichment_ORA/GO_MF_DOWN.csv", label = "MF", dir = "Down")
)

go_all <- bind_rows(lapply(go_files, function(x) read_tag(x$path, x$label, x$dir)))

top_n <- 10
go_top <- go_all %>%
    group_by(Category, Direction) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(
        neg_log10_fdr = -log10(p.adjust),
        Description = str_wrap(Description, width = 45)
    )

go_colors <- c(BP = "#1b9e77", CC = "#d95f02", MF = "#7570b3")

make_go_plot <- function(data, direction_label) {
    d <- data %>% filter(Direction == direction_label)
    
    ggplot(d, aes(x = neg_log10_fdr, 
                  y = reorder(Description, neg_log10_fdr), 
                  fill = Category)) +
        geom_col(width = 0.7) +
        scale_fill_manual(values = go_colors, name = "GO Category") +
        facet_wrap(~ Category, scales = "free_y", ncol = 1) +
        labs(title = paste0("GO Enrichment - ", direction_label, "regulated"),
             x = expression(-log[10]~"(FDR)"), y = NULL) +
        theme_bw(base_size = 11) +
        theme(strip.background = element_rect(fill = "grey90"),
              strip.text = element_text(face = "bold"),
              legend.position = "none")
}

p_go_up <- make_go_plot(go_top, "Up")
p_go_down <- make_go_plot(go_top, "Down")

ggsave("GO_barplot_UP.pdf", p_go_up, width = 9, height = 13)
ggsave("GO_barplot_DOWN.pdf", p_go_down, width = 9, height = 13)

# ============================================================
# KEGG DOTPLOTS
# ============================================================

kegg_all <- bind_rows(
    read_tag("enrichment_ORA/KEGG_UP.csv", "KEGG", "Up"),
    read_tag("enrichment_ORA/KEGG_DOWN.csv", "KEGG", "Down")
)

kegg_top <- kegg_all %>%
    group_by(Direction) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 15) %>%
    ungroup() %>%
    mutate(Description = str_wrap(Description, width = 45))

make_dot_plot <- function(data, direction_label) {
    d <- data %>% filter(Direction == direction_label)
    
    ggplot(d, aes(x = FoldEnrichment, 
                  y = reorder(Description, FoldEnrichment),
                  size = Count, color = p.adjust)) +
        geom_point() +
        scale_color_gradient(low = "red", high = "blue", name = "FDR") +
        scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
        labs(title = paste0("KEGG Enrichment - ", direction_label, "regulated"),
             x = "Fold Enrichment", y = NULL) +
        theme_bw(base_size = 11)
}

p_kegg_up <- make_dot_plot(kegg_top, "Up")
p_kegg_down <- make_dot_plot(kegg_top, "Down")

ggsave("KEGG_dotplot_UP.pdf", p_kegg_up, width = 7, height = 9)
ggsave("KEGG_dotplot_DOWN.pdf", p_kegg_down, width = 7, height = 9)

# ============================================================
# REACTOME DOTPLOTS
# ============================================================

react_all <- bind_rows(
    read_tag("enrichment_ORA/Reactome_UP.csv", "Reactome", "Up"),
    read_tag("enrichment_ORA/Reactome_DOWN.csv", "Reactome", "Down")
)

react_top <- react_all %>%
    group_by(Direction) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 15) %>%
    ungroup() %>%
    mutate(Description = str_wrap(Description, width = 45))

p_react_up <- make_dot_plot(react_top, "Up")
p_react_down <- make_dot_plot(react_top, "Down")

ggsave("Reactome_dotplot_UP.pdf", p_react_up, width = 7, height = 9)
ggsave("Reactome_dotplot_DOWN.pdf", p_react_down, width = 7, height = 9)

cat("✓ All visualization files generated.\n")