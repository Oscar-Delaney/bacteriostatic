# Drug mode of action and resource constraints modulate antimicrobial resistance evolution
We present a theoretical model comparing the efficacy of bacteriostatic and bactericidal drugs, and drugs of intermediate type, at preventing the evolutionary rescue of an initially susceptible bacterial population. We find that, all else equal, in resource-abundant environments bacteriostatic drugs are best, as they constrain cell divisions and thus allow fewer resistance mutations to occur.

We also share an interactive online simulation of AMR evolution at https://oscar-delaney.shinyapps.io/AMR-evolution/

The basic structure is:
* stochastic.R codes the core simulation
* bstatic_code.R creates the figures for our paper, available in the 'figs' folder
* app.R creates the Shiny App
