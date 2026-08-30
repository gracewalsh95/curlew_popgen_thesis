# plot_theme.R
#
# Source this file at the top of every R plotting script:
#   source("plot_theme.R")
#
# Defines: pop_colors, pop_order, set_pop(), theme_ch4, fill_pop

library(ggplot2)

pop_colors <- c("Ireland" = "#7CAE00",
                "Britain" = "#F8766D",
                "Museum"  = "#00BFC4",
                "Sweden"  = "#C77CFF")

pop_order <- c("Ireland", "Britain", "Museum", "Sweden")

# Apply factor order to a data frame's population column
set_pop <- function(df, col = "Pop") {
  df[[col]] <- factor(df[[col]], levels = pop_order)
  df
}

# Shared theme for all plots
theme_ch4 <- theme_classic(base_size = 12) +
  theme(
    axis.text.x     = element_text(hjust = 1),
    legend.position = "none",
    plot.tag        = element_text(face = "plain")
  )

fill_pop <- scale_fill_manual(values = pop_colors)
