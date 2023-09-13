library('INLA')
# sensitivity analysis 

####
# remove most dense tracts 
####
df <- read.csv("local_ecology_data/local_ecology_df.csv")[-1]
# bring in spatial information with tracts removed 
adj_sen  <- INLA::inla.read.graph(filename = "local_ecology_data/cuy_map_queensen.adj")
dfsen <- df[df$tract != 156101 & df$tract != 107701, ]
n_t <- length(unique(df$tract))
dfsen$id <- as.numeric(as.factor(dfsen$tract))
dfsen$id1 <- dfsen$id
dfsen$id2 <- dfsen$id
dfsen$year_id <- as.numeric(as.factor(dfsen$year))
nsx <-  splines::ns(dfsen$total_orgs, knots = quantile(dfsen$total_orgs, prob = c(.25,.75)), Boundary.knots = c(min(dfsen$total_orgs), max(dfsen$total_orgs))) 

dfsen$nx1 <- nsx[,1]
dfsen$nx2 <- nsx[,2]
dfsen$nx3 <- nsx[,3]
# priors 
prec_prior <- list(prec = list(prior = "loggamma", param = c(.1, .1)))
theta_prior <- list(theta = list(prior = "loggamma", param = c(.1, .1)))
besag_prior <-  list(prec = list(prior = "loggamma", param = c(.1, .1)))
prioriid2d <-  list(theta = list(prior = "wishart2d", 
                                 param =c(4,1, 1)))

mod.inlaf.spls <- inla(total_founded ~ year +  scale(log(total_population)) + 
                        vacancy_avg_std + perc_mob_std + total_disorder_std + 
                        div_std + perc_black_std + perc_latinx_std + disadvantage_std + 
                        nx1 + nx2 + nx3 +
                        f(id, n = 2*n_t, model = "iid2d", hyper = prioriid2d) +  f(id1, year, copy = "id") +
                        f(id2 ,group = year_id, model = "besag", graph = adj_sen,control.group = list(model = "ar1", hyper = theta_prior), hyper = besag_prior),
                      family= "poisson",
                      data=dfsen,
                      control.family=list(link='log'),
                      control.fixed=list(mean=0, prec = 1),
                      control.predictor=list(link=1, compute=TRUE),
                      control.compute=list( cpo=TRUE, waic=TRUE))
summary(mod.inlaf.spls)
#####
# check against diffuse priors 
#####
rm(list = ls())
df <- read.csv("local_ecology_data/local_ecology_df.csv")[-1]
nsx <-  splines::ns(df$total_orgs, knots = quantile(df$total_orgs, prob = c(.25,.75)), Boundary.knots = c(min(df$total_orgs), max(df$total_orgs))) 

df$nx1 <- nsx[,1]
df$nx2 <- nsx[,2]
df$nx3 <- nsx[,3]

adj  <- INLA::inla.read.graph(filename = "local_ecology_data/cuy_map.adj")

n_t <- length(unique(df$tract))
# set diffuse priors 
prec_priors <- list(prec = list(prior = "loggamma", param = c(1, 1)))
theta_priors <- list(theta = list(prior = "loggamma", param = c(1, 1)))
besag_priors <-  list(prec = list(prior = "loggamma", param = c(1, 1)))
prioriid2ds <-  list(theta = list(prior = "wishart2d", 
                                 param =c(4,5, 5)))

mod.inlaf.splsen <- inla(total_founded ~ year +  scale(log(total_population)) + 
                        vacancy_avg_std + perc_mob_std + total_disorder_std + 
                        div_std + perc_black_std + perc_latinx_std + disadvantage_std + 
                        nx1 + nx2 + nx3 +
                        f(id, n = 2*n_t, model = "iid2d", hyper = prioriid2ds) +  f(id1, year, copy = "id") +
                        f(id2, group = year_id, model = "besag", graph = adj,control.group = list(model = "ar1", hyper = theta_priors), hyper = besag_priors),
                      family= "poisson",
                      data=df,
                      control.family=list(link='log'),
                      control.fixed=list(mean=0, prec = 2),
                      control.predictor=list(link=1, compute=TRUE),
                      control.compute=list( cpo=TRUE, waic=TRUE))
summary(mod.inlaf.splsen)
