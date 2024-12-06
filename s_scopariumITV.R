#######################################################################################
###     "Intraspecific trait variation and coordination in Schizachyrium scoparium" ###
###                                                                                 ###
###      Data preparation and model specification                                   ###
###      J. Zeldin, 2024                                                            ###
###                                                                                 ###
#######################################################################################

#### Dependencies ####
library(tidyverse);library(brms)
library(brms);library(tidybayes)
library(posterior);library(bayesplot)
library(ggokabeito);library(ggdist)
library(ggrepel);library(ggthemes)
library(ggpubr)

#### Read source data ####

plant_dat <- read.csv("s_scoparium_ITV_whole_plant_data.csv")

#### Scale and center response traits #### 

plant_dat$rdmc_scale <- scale(plant_dat$rdmc)
plant_dat$srl_scale <- scale(plant_dat$srl)
plant_dat$cn_scale<- scale(plant_dat$C_N_ratio)
plant_dat$diam_scale <- scale(plant_dat$root_diam)
plant_dat$ldmc_scale <- scale(plant_dat$ldmc_av)
plant_dat$sla_scale <- scale(plant_dat$sla_av)


#### Model specification for univariate trait models with bayesian hierarchical models.
   # Model construction detailed for each response trait, example validations and
   # posterior predictive checks given for root dry matter content,
   # similar modeling approach and evaluation procedure for remaining response traits 

  # prior specification - weakly informative, half-cauchy priors on scale coefficients
    gen_priors <- c(prior(cauchy(0,3),class = "Intercept", dpar = "sigma"), 
                    prior(cauchy(0,3),class = "b", dpar = "sigma"), 
                    prior(cauchy(0,3), class = "sd"))
    
  ### Root dry matter content

  # model structure 
  rdmc_bf <- bf(rdmc_scale ~ 0 + (1|population) + 
                  # among-population variation term
                  (1|gr(uniq_gen, by = population)), 
                  # within population, among genotype variation term
                  # modeled by population (distributional approach)
                  sigma ~ population) 
                  # within genotype scale parameter - resid var.- 
                  # modeled by location effect -population-
  
  # model building
  rdmc_mod <- brm(rdmc_bf, data = plant_dat, prior = gen_priors,
                  iter = 10000, warmup = 2500, chains = 4, cores = 4, 
                  # 10,000 iterations with 2,500 warm up across 4 chains
                  control = list(adapt_delta = 0.9999))

  # Model summaries, validation and evaluation
    # model summaries
    summary(rdmc_mod)
    prior_summary(rdmc_mod)
    summary(rdmc_mod,prob = 0.89) %>% 
      .$fixed %>% 
      mutate(across(Estimate:`u-89% CI`,exp)) %>% 
      rbind(sum$random$population,sum$random$uniq_gen)
    
    # check the chains, traceplots
    plot(rdmc_mod, regex = TRUE)
    
    # posterior pred Checks
    y <- plant_dat$rdmc_scale[,1]; yrep<- posterior_predict(rdmc_mod, draws = 500)
    ppc_stat_grouped(y, yrep, group = plant_dat$population,stat = "mean" )
    ppc_stat_grouped(y, yrep, group = plant_dat$population,stat = "sd" )
    
    
  ### Specific Root Length (SRL) ###
    
    # Model structure
    srl_bf <-  bf(srl_scale ~ 0 + (1|population) + (1|gr(uniq_gen, by = population)),
                  sigma ~ population)
    # Run model
    srl_mod <- brm(srl_bf,data = plant_dat,warmup = 2500,
                   prior = gen_priors,
                   iter = 10000,chains = 4,cores  = 4,
                   control = list(adapt_delta=0.9999))
    
    ### Root diameter ###
    
    # Model structure
    diam_bf <-  bf(diam_scale ~0 + (1|population) + (1|gr(uniq_gen, by = population)),
                  sigma ~ population)
    # Run model
    diam_mod <- brm(diam_bf,data = plant_dat,warmup = 2500,
                    prior = gen_priors,
                    iter = 10000, chains = 4,cores  = 4,
                    control = list(adapt_delta=0.9999))
    
    ### CN Ratio ###
    
    # Model structure
    cn_bf <-  bf(cn_scale ~ 0 + (1|population) + (1|gr(uniq_gen, by = population)), 
                 sigma ~ population)
    # Run model
    cn_mod <- brm(cn_bf,data = plant_dat,warmup = 2500,
                  prior = gen_priors,
                  iter = 10000,chains = 4,cores  = 4,
                  control = list(adapt_delta=0.9999))
    
    ### Specific leaf area (SLA)
    
    # Model structure
    sla_bf <-  bf(sla_scale ~ 0 + (1|population) + (1|gr(uniq_gen, by = population)),
                  sigma ~ population)
    # Run model
    sla_mod <- brm(sla_bf,data = plant_dat,warmup = 2500,
                   prior = gen_priors,
                   iter = 10000,chains = 4,cores  = 4,
                   control = list(adapt_delta=0.9999))
    
    ### Leaf dry matter content (LDMC) ###
    
    # Model structure
    ldmc_bf <-  bf(ldmc_scale ~ 0 + (1|population) + (1|gr(uniq_gen, by = population)),
                   sigma ~ population)
    # Run model
    ldmc_mod <- brm(ldmc_bf,data = plant_dat,warmup = 2500,
                    prior = gen_priors,
                    iter = 10000,chains = 4,cores  = 4,
                    control = list(adapt_delta=0.9999))
    
    
#### Multivariate Analysis  - principal component analysis and ordination ####
    
    # Generate PCA prep table
    pca_prep <- na.omit(dplyr::select(plant_dat,population,uniq_gen,cn_scale,ldmc_scale,sla_scale,
                                      srl_scale,diam_scale,rdmc_scale))
    # PCA
    pca <- prcomp((dplyr::select(pca_prep,-population,-uniq_gen)))
    
    # component importance summary
    summary(pca)$importance
    
    # generate trait loadings
    PCloadings <- data.frame(variables = rownames(pca$rotation), pca$rotation) 
    PCloadings$variables <- fct_recode(PCloadings$variables, 
                                        CN_ratio = "cn_scale",
                                        LDMC = "ldmc_scale",
                                        SLA = "sla_scale",
                                        SRL = "srl_scale",
                                        Root_diam = "diam_scale",
                                        RDMC = "rdmc_scale")
    # Calculate population centroids
    pca_cents <- pca_prep %>% 
      mutate(pc1 = pca$x[,1], pc2 = pca$x[,2],
             pc3 = pca$x[,3], pc4 = pca$x[,4]) %>% 
      group_by(population) %>% 
      summarize(cent1 = mean(pc1), cent2 = mean(pc2), se1 = sd(pc1)/sqrt(length(pc2)),
                       se2 = sd(pc2)/sqrt(length(pc2)),
                       cent3 = mean(pc3), cent4 = mean(pc4), se3 = sd(pc3)/sqrt(length(pc3)),
                       se4 = sd(pc4)/sqrt(length(pc4))) %>% 
      mutate(population = fct_recode(population, D2 = "NACHPP", D1 = "NACHIK", C = "HOSAH",
                                     B = "ALBANY", A = "DRUMLIN")) %>% 
      mutate(population = fct_relevel(population, "A","B","C","D1","D2"))
    

    # build ordination
    pca_prep %>% 
      # attach PCA scores
      mutate(pc1 = pca$x[,1], pc2 = pca$x[,2], pc3 = pca$x[,3], pc4 = pca$x[,4]) %>% 
      mutate(population = fct_recode(population, D2 = "NACHPP", D1 = "NACHIK", C = "HOSAH",
                                     B = "ALBANY", A = "DRUMLIN")) %>% 
      # relabel populations to match manuscript
      mutate(population = fct_relevel(population, "A","B","C","D1","D2")) %>% 
      # plot
      ggplot(aes(pc1,pc2,color = population))+
      # individual observations
      geom_point(size = 2, alpha = 0.15)+
      # trait loadings
      geom_segment(data = PCloadings, 
                   aes(x = 0, y = 0, xend = (PC1*4),yend = (PC2*4)), arrow = arrow(length = unit(1/2, "picas")),
                   color = "black")+
      geom_label_repel(data = PCloadings,aes(x = PC1*4, y = PC2*4,label = variables),
                       alpha = 0.7,inherit.aes = F)+
      stat_stars(alpha = 0.25)+
      # centroid errorbars
      geom_errorbar(data = pca_cents,aes(x = cent1,ymin = cent2 - se2, ymax = cent2 + se2),width = 0.1, inherit.aes = F)+
      geom_errorbarh(data =pca_cents,aes(y = cent2,xmin = cent1 - se1, xmax = cent1 + se1),height =0.1, inherit.aes = F)+
      # centroid points
      geom_point(data = pca_cents,aes(cent1,cent2, fill = population),size = 4,
                 inherit.aes = F, color = "black",shape = 21)+
      # aesthetics
      scale_color_okabe_ito(name = "Population")+ scale_fill_okabe_ito(name = "Population")+
      xlab("PC1 (37%)") + ylab("PC2 (26%)")+ theme_few()
    
    
#### Multivariate correlation modeling ####
    
    # prior specification
    mult_cor_priors <- c(prior(lkj(1), class = rescor))
    
    # model structure
    cor_bf <- bf(mvbind(cn_scale,ldmc_scale,sla_scale,
                        rdmc_scale,srl_scale,diam_scale) ~ 0 + (1|p|population))
    
    # model building
    cor_mod <- brm(cor_bf+set_rescor(TRUE),prior = mult_cor_priors, 
                   data = plant_dat,chains = 4,
                   cores = 4, iter = 10000,warmup = 2500,
                   control = list(adapt_delta = 0.99999))
    
    # model summary
    summary(cor_mod)
    
    # traceplots
    plot(cor_mod)

    
    # Posterior plot
    cor_mod %>% 
      # gather draws from posterior
      gather_draws(`rescor.*`,regex = T) %>% 
      # calculate 89 and 95 HDIs
      mean_hdci(.width = c(.89,.95)) %>%
      pivot_wider(names_from = .width,values_from = c(.lower,.upper)) %>% 
      # assign pos/neg/ns tags
      mutate(sign = c("x","x","p","n","n",
                      "x","x","n","x","n","n",
                      "x","n","p","n")) %>% 
      # plot
      ggplot(aes(.variable,.value, color = sign))+
      # dashed line at 0
      geom_hline(yintercept = 0,linetype = "dashed", color = "grey60")+
      geom_point(size = 3)+
      geom_errorbar(aes(ymin = .lower_0.89, ymax = .upper_0.89), linewidth = 0.9,width = 0)+
      geom_errorbar(aes(ymin = .lower_0.95, ymax = .upper_0.95), width = 0)+
      scale_y_continuous(name = "Correlation coefficent", limits = c(-.75,.75))+
      # rescale to name trait pairs
      scale_x_discrete(name = NULL, 
                       labels = c('CN ratio - Diameter','CN ratio - LDMC',
                                  'CN ratio - RDMC', 'CN ratio - SLA',
                                  'CN ratio - SRL','LDMC - Diameter', 'LDMC - RDMC',
                                  'LDMC - SLA','LDMC - SRL','RDMC - Diameter',
                                  'RDMC - SRL', 'SLA - Diameter','SLA - RDMC','SLA - SRL', 
                                  'SRL - Diameter'))+
      scale_color_manual(guide = "none",values = c("tomato2","skyblue2","black"))+
      coord_flip()+
      theme_few()+
      theme(axis.text.y = element_text(hjust = 0))
    