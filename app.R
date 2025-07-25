library(shiny)
library(shinyMatrix)
library(shinyjs)
library(rsconnect)

# Load the stochastic.R file
source("stochastic.R")

general_text <- "This is a model describing the evolution of antibiotic resistance
in a bacterial population. It uses the Adaptive Tau package to stochastify
an ODE model. For details see <br> Delaney, Letten, and Engelstaedter,
<i>Drug mode of action and resource constraints modulate antimicrobial
resistance evolution</i> <a href='https://doi.org/10.1101/2023.08.29.555413'
target='_blank'>https://doi.org/10.1101/2023.08.29.555413</a>. <br>
For any suggestions or comments contact Oscar Delaney on
<a href='mailto:o.delaney@uq.net.au'>o.delaney@uq.net.au</a> <br>"

drugs_text <- "A and B are arbitrary antibiotics. The influx the concentration
\\(C_j(0)\\) is the amount of drug \\(j\\) added at bottleneck events, in units
of \\(z\\) (see 'Bacteria' tab). The mutation rate \\(m_j\\) is the proportion
of genome replications that result in resistance to drug \\(j\\). The elimination
rate \\(\\gamma_j\\) is the rate at which the drug degenerates in the body,
in units of \\(hours^{-1}\\). The maximum bactericidal and bacteriostatic
activity of drug \\(j\\) are \\(\\theta_j\\) and \\(\\phi_j\\) respectively.
The pharmacodynamic shape parameter of of drug \\(j\\) is \\(\\beta_j\\)."

bacteria_text <- "The rows represent four bacterial strains: Susceptible, 
A-resistant, B-resistant, and double resistant. \\(N_i(0)\\) is the starting
population size of strain \\(i\\). \\(\\mu_i\\) is the growth rate of strain
\\(i\\) with unlimited resources. \\(K_i\\) is the resource concentration that
produces half the maximal growth rate in strain \\(i\\). \\(\\delta_i\\) is the
intrinsic death rate of strain \\(i\\). \\(z_{i,j}\\) is the concentration of
drug \\(j\\) that produces half its maximal effect in strain \\(i\\)."

events_text <- "The two possible time-discontinuous features of the model
are dosing with antibiotics, and bottlenecks, where the populations of each
strain are scaled down by dilution into fresh media. Resources can also be
supplied continuously to the system, as in a chemostat. Drugs can be set to
disappear when the next dose arrives, which is less realistic but sometimes
convenient, or to persist through dosing."

# Default values
drugs_default <- matrix(
    c(
        6, 6, # drug influx concentrations, MIC units
        1e-9, 1e-9, # mutation rates
        0.35, 0.35, # drug elimination rates
        1, 1, # maximum bactericidal activity
        0, 0, # maximum bacteriostatic activity
        1, 1 # pharmacodynamics shape parameter
    ),
    nrow = 2, ncol = 6,
    dimnames = list(
        c("A", "B"),
        c("\\(C_j(0)\\)", "\\(m_j\\)", "\\(\\gamma_j\\)", "\\(\\theta_j\\)", "\\(\\phi_j\\)", "\\(\\beta_j\\)")
    )
)

bacteria_default <- matrix(
    c(
        "1e+10", 0, 0, 0, # init: initial populations
        1, 1, 1, 1, # mu: growth rates
        rep(1e8, 4), # k: resources at half-maximal growth rate
        rep(0.35, 4), # delta: intrinsic death rate
        1, 28, 1, 28, # zeta_A: concentration of A at half-maximal death rate
        1, 1, 28, 28 # zeta_B: concentration of B at half-maximal death rate
    ),
    nrow = 4, ncol = 6,
    dimnames = list(
        c("S", "A", "B", "AB"),
        c("\\(N_i(0)\\)", "\\(\\mu_i\\)", "\\(K_i\\)", "\\(\\delta_i\\)", "\\(z_{i,A}\\)", "\\(z_{i,B}\\)")
    )
)

general_content <- tagList(
    HTML(general_text),
    numericInput("rep", "Number of Runs", value = 10, min = 1, step = 1),
    numericInput("time", "Simulation Time (hours)", value = 60, step = 1),
    numericInput("dt", "Granularity (hours)", value = 0.1, step = 0.1),
    checkboxInput("deterministic", "Deterministic Model", FALSE),
    checkboxInput("cycl", "Cycle between drugs", FALSE),
    numericInput("seed", "Random Seed", value = NULL),
    numericInput("HGT", "recombination rate", value = 0, min = 0, step = 1e-15)
)

drugs_content <- tagList(
    p(drugs_text),
    withMathJax(matrixInput("drugs", value = drugs_default, class = "numeric"))
)

bacteria_content <- tagList(
    p(bacteria_text),
    withMathJax(matrixInput("bacteria", value = bacteria_default, class = "numeric"))
)

events_content <- tagList(
    p(events_text),
    numericInput("tau", "Bottleneck period (hours)", value = "1e4", step = 1),
    numericInput("R0", "Media resource concentration", value = "1e11", step = 1e9),
    numericInput("supply", "Resource supply rate (hours^-1)", value = "1e8", step = 1e9),
    numericInput("D", "Bottleneck Dilution fraction", value = 1e-1, step = 0.1),
    numericInput("dose_gap", "Dosing period (hours)", value = 10, step = 1),
    numericInput("dose_rep", "Drug cycling period (doses)", value = 1, step = 1),
    checkboxInput("keep_old_drugs", "Old drugs persist through dosing", TRUE),
)

display_content <- tagList(
    radioButtons("display", "Uncertainty Display",
        choices = c(
            "Median + 25th and 75th percentiles" = "median",
            "Mean + 95% confidence interval" = "mean",
            "Plot all the runs" = "all"
        ),
        selected = "all"
    )
)

graph_content <- tagList(
    plotOutput("plot", height = "600", width = "100%")
)

# UI
ui <- fluidPage(
    useShinyjs(),
    titlePanel("Stochastic Simulation of Antibiotic Resistance"),
    p(
        actionButton("run_simulation", "Run Simulation"),
        actionButton("reset_all", "Reset all parameters to defaults")
    ),
    div(
        id = "everything",
        sidebarPanel(
            tabsetPanel(
                tabPanel("General", wellPanel(general_content)),
                tabPanel("Drugs", wellPanel(drugs_content)),
                tabPanel("Bacteria", wellPanel(bacteria_content)),
                tabPanel("Events", wellPanel(events_content)),
                tabPanel("Display", wellPanel(display_content))
            )
        ),
        mainPanel(graph_content)
    )
)

server <- function(input, output, session) {
    # Create a reactive expression for the simulation results
    simulation_result <- eventReactive(input$run_simulation, {
        names <- c("N_S", "N_A", "N_B", "N_AB")
        simulate(
            rep = input$rep,
            deterministic = input$deterministic,
            seed = if (is.na(input$seed)) NULL else input$seed,
            cycl = input$cycl,
            dose_rep = input$dose_rep,
            dose_gap = input$dose_gap,
            keep_old_drugs = input$keep_old_drugs,
            time = input$time,
            tau = input$tau,
            dt = input$dt,
            R0 = input$R0,
            supply = input$supply,
            D = input$D,
            HGT = input$HGT,
            m_A = input$drugs["A", "\\(m_j\\)"],
            m_B = input$drugs["B", "\\(m_j\\)"],
            d_A = input$drugs["A", "\\(\\gamma_j\\)"],
            d_B = input$drugs["B", "\\(\\gamma_j\\)"],
            influx = setNames(input$drugs[, "\\(C_j(0)\\)"], c("C_A", "C_B")),
            init = setNames(input$bacteria[, c("\\(N_i(0)\\)")], names),
            bcidal_A = input$drugs["A", "\\(\\theta_j\\)"],
            bstatic_A = input$drugs["A", "\\(\\phi_j\\)"],
            zeta_A = setNames(input$bacteria[, "\\(z_{i,A}\\)"], names),
            kappa_A = input$drugs["A", "\\(\\beta_j\\)"],
            bcidal_B = input$drugs["B", "\\(\\theta_j\\)"],
            bstatic_B = input$drugs["B", "\\(\\phi_j\\)"],
            zeta_B = setNames(input$bacteria[, "\\(z_{i,B}\\)"], names),
            kappa_B = input$drugs["B", "\\(\\beta_j\\)"],
            delta = input$bacteria[, "\\(\\delta_i\\)"],
            mu = input$bacteria[, "\\(\\mu_i\\)"],
            k = input$bacteria[, "\\(K_i\\)"]
        )
    })
    output$plot <- renderPlot({
        # Check if the simulation_result has been executed
        if (!is.null(simulation_result())) {
            # Create a plot of the simulation results
            log_plot(simulation_result()[[1]], type = input$display)
        }
    })
    # add code to reset all parameters to defaults
    observeEvent(input$reset_all, {
        shinyjs::reset("everything")
        # update the matrices to use default values
        updateMatrixInput(session, inputId = "drugs", value = drugs_default)
        updateMatrixInput(session, inputId = "bacteria", value = bacteria_default)
    })
}

# Run the application
shinyApp(ui = ui, server = server)

# To deploy the app use deployApp(appName = "AMR-evolution")