 
path <- 'sim_results/sim_results_2023-03-13_21-22-21.csv'
res <- read.csv(path)


convert_res <- function(resdf) {
  nn <- nrow(resdf)
  newvec <- c(resdf$tau1hat,resdf$tau1tilde,resdf$tau1mime,
              resdf$tau1rcal,resdf$tau1naive,resdf$tau1oracle)
  sig_u_vec <- rep(res$n112,6)
  labs <- c(rep('Val. data only',nn),rep('C.V.',nn),rep('MIME',nn),
            rep('RC',nn),rep('Naive',nn),rep('Oracle',nn))
  
  df <- data.frame(tau_est=newvec,method=labs,sig_u=sig_u_vec)
  return(df)
}

res_long <- convert_res(res)

res_long %>%  ggplot(aes(x=as.factor(method),y=tau_est,
                   fill=as.factor(method))) + geom_violin() + 
  geom_hline(yintercept=mean(tau_est)) +
  facet_wrap(~as.factor(sig_u),ncol=2) + theme_bw() +
  theme(legend.position = 'bottom') 

res_long %>% filter(sig_u=0.25) %>% 
             
                                                                         
