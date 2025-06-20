# plotting-control-variates.R (cluster)
#
# Plots for Control Variates manuscript. Make sure directory is set to root of
# the project
#-------------------------------------------------------------------------------
library(tidyverse)

# Load in results. Separate cluster runs to keep sizes under GitHub limits
path <- 'output/sim_results/sim_results_2025-04-22_01-42-50.485214.csv' #  'sim_results/sim_results_2024-08-29_22-11-29.251738.csv'

res <- read.csv(path)
#-------------------------------------------------------------------------------
# Some pre-processing

# # Create CI coverage variable
res <- res %>% mutate(cicovmime = 1-as.numeric( cilowmime>1 | cihimime<1  ),
                      cicovoracle = 1-as.numeric( ciloworacle>1 | cihioracle<1  ),
                      cicovcv = 1-as.numeric( cilowcv>1 | cihicv<1  ),
                      cicovdrcmd = 1-as.numeric( cilowdrcmd>1 | cihidrcmd<1  ),
                      cicovnaive = 1-as.numeric( cilownaive>1 | cihinaive<1  )
                      )

# Get long form of results for ATE estimates
res_long <- res %>% select(-matches('^var1')) %>%
  pivot_longer(cols = starts_with('tau1'),
               names_to = 'method',
               values_to = 'tau_est') %>%
  mutate(method=case_when(method == 'tau1mime' ~ 'MIME',
                          method == 'tau1hat' ~ 'Val. data only',
                          method == 'tau1tilde' ~ 'C.V.',
                          method == 'tau1oracle' ~ 'Oracle',
                          method == 'tau1naive' ~ 'Naive',
                          method == 'tau1mimeunc' ~ 'MIME Unc.',
                          method == 'tau1cvgen' ~ 'C.V. Gen.',
                          method == 'tau1gen' ~ 'Generalization',
                          method== 'tau1cv' ~ 'Control variates',
                          method== 'tau1drcmd' ~ 'DRCMD')) %>%
  filter(!(method %in% c('MIME Unc.', 'C.V.')) )

# Same, but with CIs
res_long_cov <- res %>% select(-matches('^var1')) %>%
  select(-matches('^tau1')) %>%
  pivot_longer(cols = starts_with('cicov'),
               names_to = 'method',
               values_to = 'ci_coverage') %>%
  mutate(method=case_when(method == 'cicovmime' ~ 'MIME',
                          method == 'cicovoracle' ~ 'Oracle',
                          method == 'cicovtilde' ~ 'C.V.',
                          method == 'cicovcvgen' ~ 'C.V. Gen.',
                          method== 'cicovcv' ~ 'Control variates',
                          method== 'cicovdrcmd' ~ 'DRCMD',
                        method == 'cicovnaive' ~ 'Naive')) %>%
  filter(!(method %in% c('MIME Unc.', 'C.V.')) )
#-------------------------------------------------------------------------------
# Get summary levels across different varying parameters

grp_vec <- c('method','sens','n','rho','etas1')

# First, for tmt effect estimates
res_summ <- res_long %>% group_by(across(all_of(grp_vec))) %>%
  summarize(per_bias = mean(100*(tau_est-1)),
            rmse = sqrt(mean( (tau_est-1)^2 ) ) ) %>%
  mutate(val_nonr = ifelse(is.na(etas1),'S completely random', 'S depends on X'),
         rho_desc = paste('P(S=1) =',rho)) %>%
  pivot_longer(cols=c(per_bias,rmse),names_to='outcome',values_to='value')

# Next, for ci coverage
res_summ_cov <- res_long_cov %>% group_by(across(all_of(grp_vec))) %>%
  summarize(cov_rate = mean(100*ci_coverage,na.rm=T)) %>%
  mutate(val_nonr = ifelse(is.na(etas1),'S completely random', 'S depends on X'),
         rho_desc = paste('P(S=1) =',rho)) %>%
  pivot_longer(cols=c(cov_rate),names_to='outcome',values_to='value')

# Bind together
res_sum <- rbind(res_summ, res_summ_cov) %>%
  mutate(outcome=recode(outcome,'rmse'='RMSE','per_bias'='% Bias','cov_rate'='95% C.I. cov.')) %>%
  filter(!( (outcome=='RMSE') & (method %in% c('Val. data only')) & !(is.na(etas1)) ) ) %>%
  filter(!( (outcome=='RMSE') & (method %in% c('Naive'))  ) ) %>%
  mutate(method=replace(method,method=='C.V. Gen.', 'Control variates'))
#-------------------------------------------------------------------------------
# Create plots

nval=5000

# Set up axis limits
y_limits <- res_sum %>%
  filter(n == nval, is.na(etas1), rho == 0.3) %>%
  group_by(outcome) %>%
  summarize(min_value = min(value, na.rm = TRUE),
            max_value = max(value, na.rm = TRUE))

# Expand limits slightly so stuff isn't borderline cut off
y_limits <- y_limits %>%
  mutate(min_value = min_value - (0.05 * (max_value - min_value)),
         max_value = max_value + (0.05 * (max_value - min_value)))


# Make grid plot for outcomes under random S
scondrandom_plot <- res_sum %>% filter(n==nval,etas1==0.1) %>%
  filter(method != 'DRCMD') %>%
  filter(!(method=='Naive' & outcome=='95% C.I. cov.')) %>%
  ggplot(aes(x=(1-sens),y=value,group=as.factor(method),
                                                                          color=as.factor(method))) +
  geom_point() + geom_line() + facet_grid(outcome ~ as.factor(rho_desc),scales='free') +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(size=10,angle=45,hjust=1)) +
  labs(title='Validation data obtained conditionally on X, A* and Y',
       x='1-Sensitivity',
       y='Operating characteristics',
       color='') ; scondrandom_plot
ggsave(paste0("figures/cv-complexval-n",nval,".pdf"), plot = scondrandom_plot, width = 8, height = 5, units = "in")



