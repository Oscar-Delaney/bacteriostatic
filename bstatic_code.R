source("stochastic.R")
library(patchwork)
library(cowplot)

# Check whether the target was hit
target_hit <- function(sol, target = 1e2, strains = c("N_A", "N_B")) {
    target_times <- sol %>%
      filter(variable %in% strains) %>%
      group_by(rep, time) %>%
      summarise(total = sum(value), .groups = "drop") %>%
      group_by(rep) %>%
      summarise(t = time[min(which(diff(sign(total - target)) != 0), Inf)]) %>%
      pull(t)
    !is.na(target_times)
}

# Run many simulations with prespecified parameters
run_sims <- function(summary, zeta_A = c(N_S = 1, N_A = 28, N_B = 1, N_AB = 28),
zeta_B = c(N_S = 1, N_A = 1, N_B = 28, N_AB = 28), delta = 0.25, rep = 1, dose_gap = 10,
influx = 3 * c(C_A = 1, C_B = 1), m_A = 1e-9, m_B = 1e-9, d_ = 0, init_A = 0,
init_B = 0, R0 = 1e8, i_A_B = 0, i_B_A = 0, data = FALSE) {
    get_mult <- function(var, i) {
        varname <- paste0("mult_", var)
        return(ifelse(is.null(summary[i, varname]), 1, summary[i, varname]))
    }
    adjust_zeta <- function(zeta, drug, resistance_multiplier) { # `drug` should be "A" or "B"
        zeta[paste0("N_", drug)] <- zeta[paste0("N_", drug)] * resistance_multiplier
        zeta["N_AB"] <- zeta["N_AB"] * resistance_multiplier
        return(zeta)
    }

    summary$bstatic_A <- 1 - summary$bcidal_A
    summary$bstatic_B <- 1 - summary$bcidal_B
    for (i in seq_len(nrow(summary))) {
        cycl <- ifelse(summary$therapy[i] == "Cycling", TRUE, FALSE)
        res <- switch(as.character(summary$resources[i]),
            "Abundant" = 3, "Intermediate" = 1.5, "Limiting" = 0)
        d <- d_ + ifelse(cycl, 0.1, 0.35)
        sol <- simulate(
            seed = i * rep,
            init = c(N_S = ifelse(cycl, 5e8, 1e10),
              N_A = init_A, N_B = init_B, N_AB = 0) * get_mult("init", i),
            R0 = R0 * 10 ^ res,
            k = 1e8 * get_mult("k", i),
            alpha = 1 * get_mult("delta", i),
            supply = 1e8,
            mu = 1 * get_mult("mu", i),
            bcidal_A = summary$bcidal_A[i],
            bcidal_B = summary$bcidal_B[i],
            bstatic_A = summary$bstatic_A[i],
            bstatic_B = summary$bstatic_B[i],
            zeta_A = adjust_zeta(zeta_A, "A", get_mult("resistance_zeta", i)),
            zeta_B = adjust_zeta(zeta_B, "B", get_mult("resistance_zeta", i)),
            delta = (delta + res * 0.05) * get_mult("delta", i),
            time = 60,
            tau = 1e4,
            rep = rep,
            dose_gap = dose_gap,
            influx =  influx * (1 + !cycl) * get_mult("influx", i),
            cycl = cycl,
            m_A = m_A * get_mult("mutation_rate", i),
            m_B = m_B * get_mult("mutation_rate", i),
            i_A_B = i_A_B, i_B_A = i_B_A,
            d_A = d, d_B = d
        )[[1]]
        wins <- 1 - target_hit(sol)
        summary[i, c("wins", "ymin", "ymax")] <- c(mean(wins),
            binom.test(sum(wins), length(wins))$conf.int)
        cat(paste0(i, "/", nrow(summary), "\n"))
    }
    if (data) {
        return(sol)
    } else {
        return(summary)
    }
}

# Create some helper plots for annotation
bottom_plot <- ggplot() +
  annotate("segment", x = 0, xend = 0.25, y = 0, yend = 0, arrow = arrow(type = "closed", ends = "first", length = unit(0.2, "inches")), linewidth = 0.5) +
  annotate("segment", x = 1, xend = 0.75, y = 0, yend = 0, arrow = arrow(type = "closed", ends = "first", length = unit(0.2, "inches")), linewidth = 0.5) +
  annotate("text", x = 0.15, y = 0, label = "bacteriostatic", vjust = 2, size = 8) +
  annotate("text", x = 0.85, y = 0, label = "bactericidal", vjust = 2, size = 8) +
  annotate("text", x = 0.5, y = 0, label = expression(theta["A"]), size = 15) +
  theme_void()

side_plot <- ggplot() +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 0.25, arrow = arrow(type = "closed", ends = "first", length = unit(0.2, "inches")), linewidth = 0.5) +
  annotate("segment", x = 0, xend = 0, y = 1, yend = 0.75, arrow = arrow(type = "closed", ends = "first", length = unit(0.2, "inches")), linewidth = 0.5) +
  annotate("text", x = 0, y = 0.15, label = "bacteriostatic", angle = 90, vjust = -1, size = 8) +
  annotate("text", x = 0, y = 0.85, label = "bactericidal", angle = 90, vjust = -1, size = 8) +
  annotate("text", x = 0, y = 0.5, label = expression(theta["B"]), angle = 90, size = 15) +
  theme_void()

blank_plot <- ggplot() + theme_void()

# Create a 2D tile plot of the results
main_plot <- function(summary, titles = TRUE) {
    labels <- expand.grid(bcidal_A = 0, bcidal_B = 1, therapy = unique(summary$therapy), resources = unique(summary$resources))
    labels$label <- LETTERS[1:nrow(labels)]

    p <- ggplot(summary, aes(x = bcidal_A, y = bcidal_B)) +
        geom_tile(aes(fill = wins)) +
        labs(fill = "P(extinct)") +
        theme_minimal() +
        theme(
            axis.title = element_blank(),
            axis.text = element_text(size = 25),
            legend.title = element_text(size = 25),
            legend.text = element_text(size = 20),
            legend.position = "bottom",
            legend.key.width = unit(2, "cm"),
            legend.spacing.x = unit(1, "cm"),
            strip.text = element_text(size = 25, face = "bold")
        )

    if (nrow(unique(summary[, c("therapy", "resources")])) > 1) {
        all_side_plots <- side_plot / side_plot / side_plot / blank_plot +
            plot_layout(heights = c(1, 1, 1, 0.22))

        all_bottom_plots <- bottom_plot + bottom_plot

        p <- p +
            facet_grid(rows = vars(resources), cols = vars(therapy)) +
            geom_text(data = labels, aes(label = label), vjust = 1, hjust = 0, size = 15, fontface = "bold") +
            scale_fill_gradient(low = "white", high = "blue")
        if (!titles) {p <- p + theme(strip.text = element_blank())}

        ratio <- 20
    } else {
        all_side_plots <- side_plot / blank_plot +
            plot_layout(heights = c(1, 0.15))

        all_bottom_plots <- bottom_plot

        p <- p + scale_fill_gradient(low = "white", high = "blue")

        ratio <- 10
    }
    
    data_plot <- p + theme(legend.position = "none")

    final_plot <- ((all_side_plots | ((data_plot / all_bottom_plots) +
    plot_layout(heights = c(ratio, 1)))) + plot_layout(widths = c(1, ratio))) /
              get_legend(p) + plot_layout(heights = c(ratio, 1))

    return(final_plot)
}

### Figure 1
summary <- expand.grid(bcidal_A = seq(0, 1, 0.05), bcidal_B = 0,
    therapy = "Cycling", resources = "Abundant")
sol <- run_sims(summary[nrow(summary), ], rep = 1e1,
    influx = c(C_A = 6, C_B = 0), dose_gap = 5, m_B = 0, data = TRUE)
dynamics <- log_plot(sol, use = c("N_S", "N_A", "R")) +
    annotate("text", x = 0, y = Inf, label = "A", hjust = 0.5, vjust = 1.5,
        size = 15, fontface = "bold")
mono_high_res <- run_sims(summary, rep = 1e3,
    influx = c(C_A = 6, C_B = 0), dose_gap = 5, m_B = 0)
mono <- ggplot(mono_high_res, aes(x = bcidal_A, y = wins)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_light() +
    labs(
        x = NULL,
        y = "P(extinct)"
    ) +
    theme(
        axis.title = element_text(size = 35),
        axis.text = element_text(size = 25),
        plot.margin = unit(c(0, 0, 0, 2), "cm")
    ) +
    annotate("text", x = 0, y = Inf, label = "B", hjust = 1, vjust = 1.5,
        size = 15, fontface = "bold")
right <- mono / bottom_plot + plot_layout(heights = c(9, 1))
left <- dynamics / blank_plot + plot_layout(heights = c(1, 0))

# print as a pdf
pdf("figs/fig1.pdf", width = 20, height = 10)
left | right
dev.off()

save(mono_high_res, file = "figs/fig1.rdata")

### Figure 2
summary <- expand.grid(bcidal_A = seq(0, 1, 0.05), bcidal_B = seq(0, 1, 0.05),
    therapy = c("Combination", "Cycling"), resources = c("Abundant", "Intermediate", "Limiting"))
multi <- run_sims(summary, rep = 1e3)

pdf("figs/fig2.pdf", width = 20, height = 25)
main_plot(multi)
dev.off()

save(multi, file = "figs/fig2.rdata")

### Parameter sensitivity test
summary <- expand.grid(
    bcidal_A = c(0, 1), # Outer limits of figure 2.
    bcidal_B = c(0, 1), # Outer limits of figure 2.
    therapy = c("Combination", "Cycling"), # Same as figure 2.
    resources = c("Abundant", "Intermediate", "Limiting"), # Same as figure 2.
    mult_init = c(0.1, 1, 10),
    mult_mutation_rate = c(0.1, 1, 10),
    mult_influx = c(0.5, 1, 2),
    #TODO: Do we also want to vary drug elimination rate?
    mult_resistance_zeta = c(0.5, 1, 2),
    mult_mu = c(0.9, 1, 1.1),
    mult_k = c(0.9, 1, 1.1),
    mult_delta = c(0.9, 1, 1.1),
    mult_alpha = c(0.9, 1, 1.1)
)
sensitivity <- run_sims(summary, rep = 1e2)

set.seed(42)
sensitivity_wide <- sensitivity |>
    tidyr::pivot_wider(
        id_cols = c(therapy, resources, names(sensitivity)[grep("mult_", names(sensitivity))]),
        names_from = c(bcidal_A, bcidal_B),
        names_glue = "{ifelse(bcidal_A == 1, 'bcidal', 'bstatic')}_{ifelse(bcidal_B == 1, 'bcidal', 'bstatic')}",
        values_from = wins
    )
sensitivity_wide$outcome <- dplyr::select(sensitivity_wide, tidyselect::matches("b[^_]*_b[^_]*")) |>
    apply(1, function(x) {
        if(all(x == 0)) {
            return("always_survives")
        }
        if(all(x == 1)) {
            return("always_extinct")
        }
        max_case <- seq_along(x)[x == max(x)]
        if(length(max_case) > 1) {
            max_case <- sample(max_case, 1)
        }
        return(names(x)[max_case])
    }) |>
    unlist()

create_sensitivity_plot <- function(dat) {
    ggplot(dat) +
        geom_bar(aes(x = paste(therapy, resources, sep = "\n"), y = proportion, fill = outcome), stat = "identity") +
        labs(x = "Therapy type & resource availability", y = "Proportion of parameter combinations", fill = "Best combination") +
        scale_fill_manual(
            labels = c(
                "always_survives" = "None (Always Survives)",
                "always_extinct" = "None (Always Extinct)",
                "bcidal_bcidal" = "Cidal/Cidal",
                "bcidal_bstatic" = "Cidal/Static",
                "bstatic_bcidal" = "Static/Cidal",
                "bstatic_bstatic" = "Static/Static"
            ),
            values = c(
                "always_survives" = "#888888",
                "always_extinct" = "#555555",
                "bcidal_bcidal" = "#482173",
                "bcidal_bstatic" = "#25858E",
                "bstatic_bcidal" = "#2BB07F",
                "bstatic_bstatic" = "#C2DF23"
            )
        ) +
        scale_y_continuous(expand = c(0, 0))
}

sensitivity_plot <- sensitivity_wide |>
    dplyr::group_by(therapy, resources, outcome) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop_last") |>
    dplyr::mutate(proportion = count / sum(count)) |>
    create_sensitivity_plot()
pdf("figs/sensitivity.pdf", width = 8, height = 7)
print(sensitivity_plot)
dev.off()

sensitivity_breakdown_plot <- do.call(rbind,
    lapply(names(sensitivity_wide)[grep("^mult_", names(sensitivity_wide))], function(mult_var_name) {
        sensitivity_wide |>
            dplyr::group_by(therapy, resources, !!rlang::sym(mult_var_name), outcome) |>
            dplyr::summarise(count = dplyr::n(), .groups = "drop_last") |>
            dplyr::mutate(proportion = count / sum(count)) |>
            dplyr::rename(mult_var_value = !!rlang::sym(mult_var_name)) |>
            dplyr::mutate(mult_var_name = mult_var_name)
    }
)) |>
    dplyr::mutate(mult_var_value = factor(ifelse(mult_var_value == 1, "Unchanged", ifelse(mult_var_value < 1, "Decreased", "Increased")), levels = c("Decreased", "Unchanged", "Increased"), ordered = TRUE)) |>
    dplyr::mutate(mult_var_name = factor(mult_var_name, levels = c("mult_init", "mult_mutation_rate", "mult_influx", "mult_resistance_zeta", "mult_mu", "mult_k", "mult_delta", "mult_alpha"), ordered = TRUE)) |>
    create_sensitivity_plot() +
    facet_grid(
        rows = vars(mult_var_name),
        cols = vars(mult_var_value),
        labeller = labeller(mult_var_name = as_labeller(c(
            mult_init = "N(0)",
            mult_mutation_rate = "m",
            mult_influx = "C(0)",
            mult_resistance_zeta = "z",
            mult_mu = "mu",
            mult_k = "K",
            mult_delta = "delta",
            mult_alpha = "alpha"
        ), default = label_parsed))
    ) +
    theme(panel.spacing.y = unit(1, "lines"))
pdf("figs/sensitivity_breakdown.pdf", width = 24, height = 56)
print(sensitivity_breakdown_plot)
dev.off()

save(sensitivity, file = "figs/sensitivity.rdata")

### Figure S1
summary <- subset(summary, therapy == "Cycling" & resources == "Abundant")
cs <- run_sims(summary, rep = 1e3, zeta_A = c(N_S = 1, N_A = 28, N_B = 0.5, N_AB = 28),
    zeta_B = c(N_S = 1, N_A = 0.5, N_B = 28, N_AB = 28), influx = 1.7 * c(C_A = 1, C_B = 1))

pdf("figs/figS1.pdf", width = 10, height = 10)
main_plot(cs)
dev.off()

save(cs, file = "figs/figS1.rdata")

### Figure S2
quick_degrade <- run_sims(summary, rep = 1e3, influx = 30 * c(C_A = 1, C_B = 1), d_ = 0.4)

pdf("figs/figS2.pdf", width = 10, height = 10)
main_plot(quick_degrade)
dev.off()

save(quick_degrade, file = "figs/figS2.rdata")

### Figure S3
pre_existing <- run_sims(summary, rep = 1e3, m_A = 0, m_B = 0, init_B = 5, influx = 2 * c(C_A = 1, C_B = 1))

pdf("figs/figS3.pdf", width = 10, height = 10)
main_plot(pre_existing)
dev.off()

save(pre_existing, file = "figs/figS3.rdata")

### Figure S4
synergism <- run_sims(summary, rep = 1e3, i_A_B = 3, i_B_A = 3, influx = 1.5 * c(C_A = 1, C_B = 1))

pdf("figs/figS4.pdf", width = 10, height = 10)
main_plot(synergism)
dev.off()

save(synergism, file = "figs/figS4.rdata")

### Figure S5
antagonism <- run_sims(summary, rep = 1e3, i_A_B = -3, i_B_A = -3, influx = 7 * c(C_A = 1, C_B = 1))

pdf("figs/figS5.pdf", width = 10, height = 10)
main_plot(antagonism)
dev.off()

save(antagonism, file = "figs/figS5.rdata")


### Combined Figure
cr <- run_sims(summary, rep = 1e3, zeta_A = c(N_S = 1, N_A = 28, N_B = 2, N_AB = 28),
    zeta_B = c(N_S = 1, N_A = 2, N_B = 28, N_AB = 28))
save(cr, file = "figs/cr.rdata")
main_plot(cr)
further <- rbind(
  cs %>% mutate(resources = "CS/CR", therapy = "1"),
  cr %>% mutate(resources = "CS/CR", therapy = "2"),
  synergism %>% mutate(resources = "interaction", therapy = "1"),
  antagonism %>% mutate(resources = "interaction", therapy = "2"),
  quick_degrade %>% mutate(resources = "other", therapy = "1"),
  pre_existing %>% mutate(resources = "other", therapy = "2")
)
pdf("figs/further.pdf", width = 20, height = 25)
main_plot(further, titles = FALSE)
dev.off()

### Figure dynamics
summary <- expand.grid(bcidal_A = 1, bcidal_B = 1,
    therapy = c("Cycling", "Combination"), resources = "Abundant")

comb_results <- run_sims(summary[2, ], rep = 1e1, data = TRUE)
comb_graph <- log_plot(comb_results, use = c("N_S", "N_A", "N_B", "N_AB", "R")) +
    annotate("text", x = 0, y = Inf, label = "A", hjust = 0.8, vjust = 1.5,
        size = 15, fontface = "bold") +
    theme(legend.position = "none")  # Remove legend from comb_graph

cycl_results <- run_sims(summary[1, ], rep = 1e1, data = TRUE)
cycl_graph <- log_plot(cycl_results, use = c("N_S", "N_A", "N_B", "N_AB", "R")) +
    annotate("text", x = 0, y = Inf, label = "B", hjust = 0.8, vjust = 1.5,
        size = 15, fontface = "bold") +
    theme(legend.position = "right",
          legend.box = "vertical") +
    guides(color = guide_legend(order = 1))  # Strains legend first

# print as a pdf
pdf("figs/dynamics.pdf", width = 20, height = 10)
comb_graph + cycl_graph +
    plot_layout(widths = c(1, 1.05), guides = "collect")
dev.off()
