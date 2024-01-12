# inla results 
library('tidyverse')
library('INLA')
library('sf')

# easy reporting results 
library('brinla')


df <- read.csv("local_ecology_data/local_ecology_df.csv")[-1]

adj  <- INLA::inla.read.graph(filename = "local_ecology_data/cuy_map.adj")

n_t <- length(unique(df$tract))


## priors 
# check and set priors 
inla.models()$latent$besag$hyper
inla.models()$latent$iid$hyper
inla.models()$latent$iid2d$hyper
inla.models()$latent$ar1$hyper
prec_prior <- list(prec = list(prior = "loggamma", param = c(.1, .1)))
theta_prior <- list(theta = list(prior = "loggamma", param = c(.1, .1)))
besag_prior <-  list(prec = list(prior = "loggamma", param = c(.1, .1)))
prioriid2d <-  list(theta = list(prior = "wishart2d", 
                                 param =c(4,1, 1)))
####
## fit models 
###

# base model
mod.inla0 <- inla(total_founded ~  1, family= "poisson", 
                  data=df,
                  control.family=list(link='log'),
                  control.fixed=list(mean=0, prec = 1),
                  control.predictor=list(link=1, compute=TRUE),
                  control.compute=list( cpo=TRUE, waic=TRUE))

summary(mod.inla0)

# covariates 
mod.inla1 <- inla(total_founded ~   year +   scale(log(total_population)) + 
                    vacancy_avg_std + perc_mob_std + total_disorder_std + 
                    org_std + div_std + perc_black_std + perc_latinx_std + disadvantage_std, 
                  family= "poisson", 
                  data=df,
                  control.family=list(link='log'),
                  control.fixed=list(mean=0, prec = 1),
                  control.predictor=list(link=1, compute=TRUE),
                  control.compute=list( cpo=TRUE, waic=TRUE))

summary(mod.inla1)

## correlated RE
mod.inla2 <- inla(total_founded ~ year +   scale(log(total_population)) + 
                    vacancy_avg_std + perc_mob_std + total_disorder_std + 
                    org_std + div_std + perc_black_std + perc_latinx_std + disadvantage_std +
                    f(id1, year, copy = "id") + 
                    f(id, n = 2*n_t, model = "iid2d", hyper = prioriid2d), family= "poisson", 
                  data=df,
                  control.family=list(link='log'),
                  control.fixed=list(mean=0, prec = 1),
                  control.predictor=list(link=1, compute=TRUE),
                  control.compute=list( cpo=TRUE, waic=TRUE))

summary(mod.inla2)

# time constant ICAR 
mod.inla3 <- inla(total_founded ~ year +   scale(log(total_population)) + 
                    vacancy_avg_std + perc_mob_std + total_disorder_std + 
                    org_std + div_std + perc_black_std + perc_latinx_std + disadvantage_std +
                    f(id, n = 2*n_t, model = "iid2d", hyper = prioriid2d) +
                    f(id1, year, copy = "id") +
                    f(id2, model = "besag", graph = adj, hyper =  besag_prior),
                  family= "poisson",
                  data=df,
                  control.family=list(link='log'),
                  control.fixed=list(mean=0, prec = 1),
                  control.predictor=list(link=1, compute=TRUE),
                  control.compute=list( cpo=TRUE, waic=TRUE))

summary(mod.inla3)

# AR1 ICAR
mod.inla4 <- inla(total_founded ~ year +   scale(log(total_population)) + 
                    vacancy_avg_std + perc_mob_std + total_disorder_std + 
                    org_std + div_std + perc_black_std + perc_latinx_std + disadvantage_std +
                    f(id, n = 2*n_t, model = "iid2d", hyper = prioriid2d) +  f(id1, year, copy = "id") +
                    f(id2 ,group = year_id, model = "besag", graph = adj,control.group = list(model = "ar1", hyper = theta_prior), hyper = besag_prior ),
                  family= "poisson", 
                  data=df,
                  control.family=list(link='log'),
                  control.fixed=list(mean=0, prec = 1),
                  control.predictor=list(link=1, compute=TRUE),
                  control.compute=list( cpo=TRUE, waic=TRUE))
summary(mod.inla4)
# spline 
nsx <-  splines::ns(df$total_orgs, knots = quantile(df$total_orgs, prob = c(.25,.75)), Boundary.knots = c(min(df$total_orgs), max(df$total_orgs))) 

df$nx1 <- nsx[,1]
df$nx2 <- nsx[,2]
df$nx3 <- nsx[,3]

mod.inlaf.spl <- inla(total_founded ~ year +  scale(log(total_population)) + 
                        vacancy_avg_std + perc_mob_std + total_disorder_std + 
                        div_std + perc_black_std + perc_latinx_std + disadvantage_std + 
                        nx1 + nx2 + nx3 +
                        f(id, n = 2*n_t, model = "iid2d", hyper = prioriid2d) +  f(id1, year, copy = "id") +
                        f(id2 ,group = year_id, model = "besag", graph = adj,control.group = list(model = "ar1", hyper = theta_prior), hyper = besag_prior),
                      family= "poisson",
                      data=df,
                      control.family=list(link='log'),
                      control.fixed=list(mean=0, prec = 1),
                      control.predictor=list(link=1, compute=TRUE),
                      control.compute=list( cpo=TRUE, waic=TRUE))
summary(mod.inlaf.spl)

(waic_list <- c("0" = mod.inla0$waic$waic,
                "1" = mod.inla1$waic$waic,
                "2" = mod.inla2$waic$waic,
                "3" = mod.inla3$waic$waic,
                "4" = mod.inla4$waic$waic,
                "5" = mod.inlaf.spl$waic$waic))



(ic_dat <- waic_list  |> as.data.frame() |>  
    mutate(delta_waic = waic_list - min(waic_list),
                             exp_delta_waic = exp(-.5*delta_waic),
                             total_waic = sum(exp_delta_waic),
                             w_waic = exp_delta_waic / total_waic))

# plot curve 
df$meanp <- mod.inlaf.spl$summary.fitted.values$mean
df$uu <- mod.inlaf.spl$summary.fitted.values$`0.975quant`
df$ll <- mod.inlaf.spl$summary.fitted.values$`0.025quant`

modpred <- df |> select(total_orgs, meanp, uu, ll) |>  
  pivot_longer(cols = !total_orgs, names_to = "prediction")
(fig2 <- ggplot(data = modpred, aes(x = total_orgs, y = value, linetype = prediction)) + 
    geom_smooth( color = "black", formula = y ~ splines::ns(x,
                                                            knots = quantile(x, prob = c(.25,.75)), Boundary.knots = c(min(x), max(x))), se = FALSE) +
    scale_linetype_manual(values=c("dotted", "solid", "dotted"), labels=c('2.5 %', 'Mean', '97.5 %'))+
    coord_cartesian(ylim=c(0,15), xlim = c(0,500)) + 
    labs(y = "Nonprofit founding events", x = "Nonprofit density",  linetype = "Prediction") + 
    theme_bw() +
    theme(text=element_text(family="Times New Roman",
                            size=12,),
          legend.position = "none") )


