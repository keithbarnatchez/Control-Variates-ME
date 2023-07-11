 
path <- 'sim_results/sim_results_2023-07-06_18-00-51.csv'
res <- read.csv(path)

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
                          method == 'tau1cvgen' ~ 'C.V. Gen.'))

res_long_cov <- res %>% select(-matches('^var1')) %>% 
  select(-matches('^tau1')) %>%
  pivot_longer(cols = starts_with('cicov'),
               names_to = 'method',
               values_to = 'ci_coverage') %>%
  mutate(method=case_when(method == 'cicovmime' ~ 'MIME',
                          method == 'cicovoracle' ~ 'Oracle',
                          method == 'cicovtilde' ~ 'C.V.',
                          method == 'cicovcvgen' ~ 'C.V. Gen.'))


res_long %>% filter(n==5000,rho==0.2) %>% ggplot(aes(x=as.factor(method),y=tau_est,
                         fill=as.factor(method))) + geom_violin() + 
  geom_hline(yintercept = 1,color='red') +
  facet_wrap(~as.factor(sens),ncol=2) + theme_bw() +
  theme(legend.position = 'bottom') + theme_bw() +
  labs(title='Varying sensitivity',
       subtitle = 'Specificity fixed at 0.95',
       y='ATE estimate',
       x='Method') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')

res_long %>%  ggplot(aes(x=as.factor(method),y=tau_est,
                   fill=as.factor(method))) + geom_violin() + 
                    geom_hline(yintercept = 1,color='red') +
  facet_wrap(~as.factor(sig_u),ncol=2) + theme_bw() +
  theme(legend.position = 'bottom') +
  labs(title='Varying sensitivity',
       subtitle = 'Specificity fixed at 0.95',
       y='ATE estimate',
       x='Method') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')


#--------------
# LINE PLOTS

grp_vec <- c('method','sens','n','rho','etas1')

# First, for tmt effect estimates
res_summ <- res_long %>% group_by(across(all_of(grp_vec))) %>%
  summarize(per_bias = mean(100*(tau_est-1)),
            rmse = sqrt(mean( (tau_est-1)^2 ) ) ) %>% 
  mutate(val_nonr = ifelse(is.na(etas1),'S completely random', 'S depends on X'),
         rho_desc = paste('rho =',rho))

# Next, for ci coverage
res_summ_cov <- res_long_cov %>% group_by(across(all_of(grp_vec))) %>%
  summarize(cov_rate = mean(100*ci_coverage,na.rm=T)) %>% 
  mutate(val_nonr = ifelse(is.na(etas1),'S completely random', 'S depends on X'),
         rho_desc = paste('rho =',rho))

# Make grid plot for bias
bias_plot <- res_summ %>% filter(n==5000) %>% ggplot(aes(x=(1-sens),y=per_bias,group=as.factor(method),
                                                                              color=as.factor(method))) +
  geom_point() + geom_line() + facet_grid(val_nonr ~ as.factor(rho_desc),scales='free') + 
theme_bw() +
  theme(legend.position = 'bottom',
       axis.text.x = element_text(size=10,angle=45,hjust=1)) +
  labs(title='Percent bias, varying measurement error sensitivity',
       x='1-Sensitivity',
       y='Percent bias',
       color='') ; bias_plot
ggsave("figures/cv-proj-perbias-grid.pdf", plot = bias_plot, width = 8, height = 5, units = "in")

# RMSE with all methods
rmse_plot <- res_summ %>% filter(n==5000) %>%
  ggplot(aes(x=(1-sens),y=rmse,group=as.factor(method),
             color=as.factor(method))) +
  geom_point() + geom_line() + facet_grid(val_nonr ~ as.factor(rho_desc),scales='free') + 
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(size=10,angle=45,hjust=1)) +
  labs(title='RMSE, varying measurement error sensitivity',
       x='1-Sensitivity',
       y='RMSE',
       color='') ; rmse_plot
ggsave("figures/cv-proj-rmse-grid.pdf", plot = rmse_plot, width = 8, height = 5, units = "in")

# RMSE, excluding heavily-biased methods
rmse_plot <- res_summ %>% filter(n==5000,method!='Naive',method!='MIME Unc.') %>%
 ggplot(aes(x=(1-sens),y=rmse,group=as.factor(method),
                                                         color=as.factor(method))) +
  geom_point() + geom_line() + facet_grid(val_nonr ~ as.factor(rho_desc),scales='free') + 
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(size=10,angle=45,hjust=1)) +
  labs(title='RMSE, varying measurement error sensitivity',
       x='1-Sensitivity',
       y='RMSE',
       color='') ; rmse_plot
ggsave("figures/cv-proj-rmse-grid-nonaive.pdf", plot = rmse_plot, width = 8, height = 5, units = "in")

# CI coverage
cov_plot <- res_summ_cov %>% filter(n==5000,method!='Naive',method!='MIME Unc.') %>%
  ggplot(aes(x=(1-sens),y=cov_rate,group=as.factor(method),
             color=as.factor(method))) +
  geom_point() + geom_line() + geom_hline(yintercept = 95, color='red',size=.8) + facet_grid(val_nonr ~ as.factor(rho_desc),scales='free') + 
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(size=10,angle=45,hjust=1)) +
  labs(title='95% CI coverage, varying measurement error sensitivity',
       x='1-Sensitivity',
       y='95% CI Coverage Rate',
       color='') ; cov_plot




