library(ggplot2)
path_1d_lognormal_0_1 = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/uncertainties_by_model/0.1_step_size/1d lognormal_GO_results_plus_uncertainties_0.1_step.csv'

data <- read.csv(file = path_1d_lognormal_0_1, header = T, stringsAsFactors = F)

# fix direction of x axis labels
# check if need to take logs, or are godambe stats weird?
ggplot(data, aes(x=data$GO.Term.Name, y=as.numeric(data$log_mu))) + geom_point() + geom_errorbar(aes(ymin = data$log_mu - log(data$log_mu.uncertainty), ymax = data$log_mu + log(data$log_mu.uncertainty))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
