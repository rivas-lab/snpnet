#' generate a 2d heatmap plot comparing the PRS (Z-score if geno_z == TRUE) vs phenotype.
#'
#' @param plot_df PRS and phenotype columns for (typically test) set of individuals
#' @param plot_bin2d_x size of the x-axis grid in 2d heatmap
#' @param plot_bin2d_y size of the y-axis grid in 2d heatmap
#' @param geno_score_col column name of the score in plot_df
#' @param phe_col column name of the phenotype data in plot_df
#' @param geno_z whether to take the Z-score transformation of the score
#'
plot_PRS_vs_phe <- function(plot_df, plot_bin2d_x=NULL, plot_bin2d_y=NULL, geno_score_col="geno_score", phe_col="phe", geno_z=TRUE){
    # rename the columns
    plot_df %>%
    rename(!!'phe' := all_of(phe_col)) %>%
    rename(!!'geno_score' := all_of(geno_score_col)) -> plot_dff

    if(geno_z){
        plot_dff %>% mutate(
            geno_score_z = (geno_score - mean(geno_score)) / sd(geno_score)
        ) -> plot_dff
        if(is.null(plot_bin2d_x)) plot_bin2d_x <- 0.05
    }else{
        plot_dff %>% mutate(geno_score_z = geno_score) -> plot_dff
        if(is.null(plot_bin2d_x)){
            # compute the bin size for x-axis
            plot_bin2d_x <- diff(quantile(plot_dff$geno_score, c(.4, .6))) / 4
        }
    }

    if(is.null(plot_bin2d_y)){
        # compute the bin size for y-axis
        plot_bin2d_y <- diff(quantile(plot_dff$phe, c(.4, .6))) / 4
    }

    # plot the bin2d plot for the 99.9% coverage
    plot_dff %>%
    filter(
        # 99.9% coverage
        quantile(plot_dff$phe, .0005) < phe,
        quantile(plot_dff$phe, .9995) > phe
    ) %>%
    ggplot(aes(x = geno_score_z, y = phe)) +
    geom_bin2d(binwidth = c(plot_bin2d_x, plot_bin2d_y)) +
    scale_fill_continuous(type = "viridis") +
    theme_bw(base_size=16) +
    labs(x = ifelse(geno_z, 'snpnet PRS (Z-score)', 'snpnet PRS'))
}


#' generate Violin plot comparing the distribution of PRSs (Z-score if geno_z == TRUE)
#' stratified by case/control status.
#' plot_df should contain PRS and phenotype columns for (typically test) set of individuals
#' the column names can be specified in geno_score_col and phe_col
plot_PRS_binomial <- function(plot_df, geno_score_col="geno_score", phe_col="phe", geno_z=TRUE){
    # rename the columns
    plot_df %>%
    rename(!!'phe' := all_of(phe_col)) %>%
    rename(!!'geno_score' := all_of(geno_score_col)) -> plot_dff

    if(geno_z){
        plot_dff %>% mutate(
            geno_score_z = (geno_score - mean(geno_score)) / sd(geno_score)
        ) -> plot_dff
    }else{
        plot_dff %>% mutate(geno_score_z = geno_score) -> plot_dff
    }

    plot_dff %>%
    left_join(
        data.frame(phe_str=c('Control', 'Case'), phe=c(1, 2), stringsAsFactors=F), by='phe'
    ) %>%
    ggplot(aes(x = reorder(as.factor(phe_str), phe), y = geno_score_z, color=as.factor(phe))) +
    geom_violin() +
    geom_boxplot(outlier.size = 0, outlier.stroke = 0, width = 0.2) +
    stat_summary(
        fun=mean, geom="errorbar",
        aes(ymax = ..y.., ymin = ..y..),
        width = 1.1, linetype = "dashed"
    ) +
    theme_bw(base_size=16)+
    theme(legend.position = "none") +
    labs(x = "phenotype", y = ifelse(geno_z, "snpnet PRS (Z-score)", "snpnet PRS"))
}

plot_PRS_bin_vs_phe <- function(summary_plot_df, horizontal_line){
    summary_plot_df %>%
    mutate(
        x_ticks_labels = if_else(
            startsWith(bin_str, "0%"),
            paste0("[", bin_str, "]"),
            paste0("(", bin_str, "]")
        )
    ) %>%
    ggplot(aes(x=reorder(x_ticks_labels, -u_bin), y=mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = l_err, ymax = u_err)) +
    geom_hline(yintercept = horizontal_line, color='gray')+
    theme_bw(base_size=16) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))+
    labs(x = "The percentile of snpnet PRS")
}

plot_PRS_bin_vs_OR <- function(summary_plot_df){
    summary_plot_df %>%
    mutate(mean = OR, l_err = l_OR, u_err = u_OR) %>%
    plot_PRS_bin_vs_phe(1) +
    labs(y = "Odds ratio [SE]")
}
