correct_genomewide_burden <- function(genes_df) {
    mu_genome = sum(genes_df$gamma_per_sd * sqrt(genes_df$burden_score))/sum(genes_df$burden_score)
    message(paste("Genome-wide mean effect (mu_genome) applied per-sd:", mu_genome*sqrt(sum(genes_df$burden_score))))
    genes_df$gamma_per_sd <- genes_df$gamma_per_sd - mu_genome * sqrt(genes_df$burden_score)
    return(genes_df)
}