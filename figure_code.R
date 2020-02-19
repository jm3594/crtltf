library(ggpubr)
library(tidyverse)
library(zeallot)
library(latex2exp)
library(directlabels)
library(gridExtra)

source("R/cpa.did.normal.ltf.R")
source("R/get_ltf_var.R")

# function for figure 1, to get variance components
get_vars <- function(vart,icc, rho.c, rho.s){
  varc <- rho.c*icc*vart
  varct <- (1 - rho.c)*icc*vart
  vars <- rho.s*(1 - icc)*vart
  varst <- (1 - rho.s)*(1 - icc)*vart
  list(varc = varc, varct = varct,
       vars = vars, varst = varst)
}

# Figure 1 ------------------------------------------------------------------------------
vart <- 1
icc <- 0.05
rho.c <- 0.3
rho.s <- seq(0.1,0.9,by=0.1)
J <- 30
K <- 100

ltf_0 <-
  c(varc,varct,vars,varst) %<-% get_vars(vart=vart, icc=icc, rho.c=rho.c, rho.s=rho.s)
L1 <- L2 <- seq(from = 0, to = 90, by = 1)
G1 <- G2 <- 0
vardf <- get_var(vart,icc,rho.c,rho.s,J,K,L1,L2,G1,G2) %>%
  mutate_at(vars(ltfnr:cohort), list(~./cohort)) %>%
  mutate(ltfnr_by_mix = ltfnr/mix,xsec_by_mix = xsec/mix)
vardf2 <- vardf %>%  filter(L1 == L2) %>%
  select(rho.c,rho.s,L1,ltfnr:xsec) %>%
  gather(design,reff,ltfnr:xsec) %>%
  mutate(design = factor(design),
         design = fct_relevel(design,"mix","xsec","ltfnr","reduced"),
         design = fct_recode(design,
                             Replacement="mix","Cross-sectional"="xsec",
                             "Non-replacement"="ltfnr","Reduced"="reduced"))

# for paper
text_size <- 10

caption1 <- str_wrap("Variance of cross-sectional, replacement, non-replacement, and
                     reduced designs relative to the complete design as a function of percent
                     loss to follow-up.", width = 80)

fig1 <- ggplot(vardf2,aes(x=L1,y=reff,color=design)) +
  geom_line() +
  scale_y_continuous(limits = c(1,1.75)) +
  facet_wrap(vars(rho.s), nrow = 3,
             labeller=label_bquote(cols = rho[S] == .(rho.s))) +
  # labs(color = "Design",
  #      x = "Percent Loss to Follow-up",
  #      y = "Relative Variance",
  #      title = "Figure 1",
  #      caption = caption1) +
  # theme(legend.position = "bottom",
  #       plot.caption = element_text(hjust=0, size = 10))
  labs(color = "Design",
       x = "Percent Loss to Follow-up",
       y = "Relative Variance",
       title = "Moyer Figure 1 TOP") +
  theme(legend.position = "bottom")

#ggsave("figure1.pdf",fig1,width=6,height=7)
ggsave("figures/figure1.pdf",fig1,width=6,height=5,dpi=800)
#fig1

# Figure 2 ------------------------------------------------------------------------------

which.median <- function(x) {
  floor(median(x))
}

label_equal <- function(labels){
  label_both(labels,sep = " = ")
}

brks <- c(0,10,30,50,70,90)
#brks <- seq(0,90,by=10)

# for paper
text_size <- 10
key_width_unit <- 0.25

# for posterFfig
# text_size <- 24
# key_width_unit <- 1

vardf <- vardf %>% filter(rho.s %in% c(0.2,0.4,0.6,0.8))

make_contour_plot <- function(design,title,color,vbrks,brks,
                              text_size,key_width_unit,dat=vardf){
  ggplot(dat, aes_string(x = "L1", y = "L2", z = design)) +
    geom_raster(aes_string(fill = design)) + #, colour = "gray90") +
    facet_grid(.~rho.s, labeller = label_bquote(cols = rho[S] == .(rho.s))) +
    scale_x_continuous(breaks = brks) +
    scale_y_continuous(breaks = brks) +
    scale_fill_gradient(low = "white", high = color) +
    scale_color_manual(values = c("white", color)) +
    labs(title = title,
         x = "Arm 1 Percent LTF",
         y = "Arm 2 Percent LTF",
         fill = "Variance Ratio") +
    geom_contour(color=color,breaks=vbrks) +
    geom_dl(aes(label=..level..),method=list("bottom.pieces",cex=0.3),
            stat="contour", breaks=vbrks) +
    theme(legend.position = "right",
          plot.title = element_text(size = text_size),
          text = element_text(size = text_size),
          #legend.text = element_text(angle = -45, hjust=0.25),
          legend.key.width = unit(key_width_unit,"cm"),
          strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")))
}

v1brks <- seq(1.02,1.2,by=0.04)
v1 <- make_contour_plot("mix","Replacement relative to Complete","red",
                        v1brks,brks,text_size,key_width_unit)

v2brks <- c(1,1.25,1.5,1.75,2)#seq(1,1.8,by=.2)
v2 <- make_contour_plot("ltfnr","Non-replacement relative to Complete","blue",
                        v2brks,brks,text_size,key_width_unit)

v3brks <- seq(1,2,by=.2)
v3 <- make_contour_plot("ltfnr_by_mix","Non-replacement relative to Replacement","green",
                        v3brks,brks,text_size,key_width_unit)

v4brks <- seq(1,1.2,by = .04)
v4 <- make_contour_plot("xsec_by_mix","Cross-sectional relative to Mixture","purple",
                        v4brks,brks,text_size,key_width_unit)

fig2title <- "Moyer Figure 2 TOP"

fig2 <- arrangeGrob(v1,v2,v3,v4,nrow=4,ncol=1,
                    top=text_grob(fig2title, hjust=0, x =0, size = 14))

# fig2 <- ggarrange(v1, v2, v3, labels = c("a", "b","c"), ncol = 3, nrow = 1)
# fig2 <- ggarrange(v1, v2, v3, v4,
#                   labels = c("a", "b","c","d"), ncol = 1, nrow = 4)
#ggsave("fig2.png",width=14,height=9.33)
# annotate_figure(fig2, left = "Arm 2 Percent LTF", bottom= "Arm 1 Percent LTF",
#                 top = "Figure 2")
# annotate_figure(fig2,
#                 top = "Figure 2")
caption2 <- "Relative variances as functions of loss to follow-up in treatment arms."

# fig2ann <- annotate_figure(fig2,bottom = text_grob(label=caption2,hjust=0,x=0,size=10))
fig2ann <- fig2

#ggsave("figure2.pdf",fig2,width=6.5,height=9)
ggsave("figures/figure2.pdf",fig2ann,width=6,height=7,dpi=800)


# Figure 3 ------------------------------------------------------------------------------
alpha <- 0.05
power <- NA
nclusters <- c(5,10,30)
nsubjects <- 100
d <- 0.20
icc <- 0.05
rho_c <- c(0.1,0.3,0.5)
rho_s <- seq(0.05,0.9,by=0.05)
vart <- 1
ltf_0 <- ltf_1 <- c(.50)
gtf_0 <- gtf_1 <- 0

fig3_text <- 12
fig3_xtext <- 10

caption3 <- "Power as a function of subject autocorrelation for the no replacement and
reduced cohort designs."

fig3dfa <- expand.grid(alpha=alpha,power=power,nclusters=nclusters,nsubjects=nsubjects,
                       d=d,icc=icc,rho_c=rho_c,rho_s=rho_s,vart=vart,
                       ltf_0=ltf_0,ltf_1=ltf_1,gtf_0=gtf_0,gtf_1=gtf_1) %>%
  filter(ltf_0==ltf_1) %>%
  mutate(power = pmap_dbl(.,cpa.did.normal.ltf),
         type = "Non-replacement")
fig3dfb <- expand.grid(alpha=alpha,power=power,nclusters=nclusters,nsubjects=nsubjects,
                       d=d,icc=icc,rho_c=rho_c,rho_s=rho_s,vart=vart,
                       ltf_0=ltf_0,ltf_1=ltf_1,gtf_0=gtf_0,gtf_1=gtf_1) %>%
  filter(ltf_0==ltf_1) %>%
  mutate(nsubjects=(1-ltf_0)*nsubjects,
         ltf_0=0,ltf_1=0) %>%
  mutate(power = pmap_dbl(.,cpa.did.normal.ltf),
         type = "Reduced",
         ltf_0=fig3dfa$ltf_0)
fig3df <- rbind(fig3dfa,fig3dfb)
fig3<-ggplot(fig3df, aes(x=rho_s,y=power,color = type)) +
  #geom_point() +
  geom_line() +
  facet_grid(rho_c~nclusters,
             labeller = label_bquote(rows = rho[C] == .(rho_c),
                                     cols = J == .(nclusters))) +
  # labs(color = "Design",
  #      y = "Power",
  #      x=TeX("Subject Autocorrelation $\\rho_S$"),
  #      title = "Figure 3",
  #      caption = str_wrap(caption3,width=90)) +
  # theme(plot.title = element_text(size = fig3_text),
  #       text = element_text(size=fig3_text),
  #       plot.caption = element_text(hjust= 0, size=10),
  #       legend.position = "bottom",
  #       axis.text.x = element_text(size = fig3_xtext))
  labs(color = "Design",
       y = "Power",
       x=TeX("Subject Autocorrelation $\\rho_S$"),
       title = "Moyer Figure 3 TOP") +
  theme(plot.title = element_text(size = fig3_text),
        text = element_text(size=fig3_text),
        legend.position = "bottom",
        axis.text.x = element_text(size = fig3_xtext))
#legend.text = element_text(angle = -45, hjust=0.25),
# legend.key.width = unit(1,"cm"))
# ggsave("fig3.png",width=9,height=6)
ggsave("figures/figure3.pdf",fig3,width=6.5,height=6.5,dpi=800)
# fig3

# Figure 4 ------------------------------------------------------------------------------
# function for figure 4
get_K1 <- function(K,K2,rhos) { K*K2*(3 - 4*rhos)/(K2*(2 - 4*rhos) + K)}

fig4_text <- 12
fig4_xtext <- 10
rhos <- c(0.60, 0.65, 0.70, 0.74, 0.749, 0.7499)
K2 <- seq(from=0,to=100,by=1)
K <- 100
caption4 <- "Blue and red regions indicate where the no replacement and reduced cohort designs
are more powerful, respectively. As subject autocorrelation increases, the reduced cohort exhibits
more power as represented by the increased red area."
d <- expand.grid(K=K,K2=K2,rhos=rhos) %>%
  mutate(K1 = pmap_dbl(select(., K,K2,rhos),get_K1)) %>%
  mutate(K1r = K1/K, K2r = K2/K) %>%
  mutate(L1r = 1 - K1r, L2r = 1 - K2r) %>%
  filter(K1 > 0, K1 <= K2)
label_both_equal <- function(labels) label_both(labels, sep = " = ")
fig4 <- ggplot(d, aes(x=K2r,y=K1r)) +
  geom_line() +
  geom_point(aes(x=0.95,y=0.84)) +
  geom_segment(aes(x=0,y=0,xend=1,yend=1), lty = "dotted") +
  geom_ribbon(aes(ymin=K1r,ymax=K2r, fill = "Reduced Cohort more powerful"), alpha = 0.3) +
  geom_area(aes(y = K1r, fill = "No Replacement more powerful"), alpha = 0.3) +
  scale_fill_manual(values=c("blue","red")) +
  # labs(x = TeX("Follow-up Rate in Arm 2 ($1-\\lambda_2$)"),
  #      y = TeX("Follow-up Rate in Arm 1 ($1-\\lambda_1$)"),
  #      title = "Figure 4",
  #      caption = str_wrap(caption4,width=80),
  #      fill="") +
  # facet_wrap(.~rhos,nrow=2, labeller = label_bquote(rho[S] == .(rhos))) +
  # theme(plot.title = element_text(size = fig4_text),
  #       text = element_text(size=fig4_text),
  #       axis.text.x = element_text(size = fig4_xtext),
  #       legend.position = "bottom",
  #       plot.caption = element_text(hjust=0, size=10),
  #       #legend.text = element_text(angle = -45, hjust=0.25),
  #       legend.key.width = unit(0.75,"cm"))
  labs(x = TeX("Follow-up Rate in Arm 2 ($1-\\lambda_2$)"),
       y = TeX("Follow-up Rate in Arm 1 ($1-\\lambda_1$)"),
       title = "Moyer Figure 4 TOP",
       fill="") +
  facet_wrap(.~rhos,nrow=2, labeller = label_bquote(rho[S] == .(rhos))) +
  theme(plot.title = element_text(size = fig4_text),
        text = element_text(size=fig4_text),
        axis.text.x = element_text(size = fig4_xtext),
        legend.position = "bottom",
        plot.caption = element_text(hjust=0, size=10),
        #legend.text = element_text(angle = -45, hjust=0.25),
        legend.key.width = unit(0.75,"cm"))
#ggsave("fig4.png",width=12,height=8)
#fig4
ggsave("figures/figure4.pdf",fig4,width=6,height=5,dpi=800)

# Figure 5 ------------------------------------------------------------------------------
alpha <- 0.05
power <- NA
nclusters <- 15
nsubjects <- 40:300
d <- 0.12
varc <- 0.0218
varct <- 0.0047
vars <- 0.3342
varst <- 0.2567
vart <- varc+varct+vars+varst
icc <- (varc+varct)/(vart)
rho_c <- varc/(varc+varct)
rho_s <- vars/(vars+varst)
ltf_0 <- 0.05
ltf_1 <- 0.16
gtf_0 <- 0
gtf_1 <- 0

fig5adfa <- expand.grid(alpha=alpha,power=power,nclusters=nclusters,nsubjects=nsubjects,
                        d=d,icc=icc,rho_c=rho_c,rho_s=rho_s,vart=vart,
                        ltf_0=ltf_0,ltf_1=ltf_1,gtf_0=gtf_0,gtf_1=gtf_1) %>%
  mutate(power = pmap_dbl(.,cpa.did.normal.ltf)) %>%rename(Power = power) %>%
  mutate(Design = "Non-replacement")

gtf_0 <- 0.05
gtf_1 <- 0.16
fig5adfb <- expand.grid(alpha=alpha,power=power,nclusters=nclusters,nsubjects=nsubjects,
                        d=d,icc=icc,rho_c=rho_c,rho_s=rho_s,vart=vart,
                        ltf_0=ltf_0,ltf_1=ltf_1,gtf_0=gtf_0,gtf_1=gtf_1) %>%
  mutate(power = pmap_dbl(.,cpa.did.normal.ltf)) %>%
  rename(Power = power) %>%
  mutate(Design = "Replacement")

fig5adf <- rbind(fig5adfa,fig5adfb) %>%
  mutate(type = 'a')

fig5a <- ggplot(fig5adf,aes(x=nsubjects,y=Power,color=Design)) +
  geom_line() + geom_hline(yintercept = 0.8, lty = "dashed") +
  labs(title = "Figure 5") + scale_y_continuous(limits=c(0.5,.9)) +
  theme(panel.grid.minor.x = element_blank())

ltf_0 <- 0.50
ltf_1 <- 0.50
gtf_0 <- 0
gtf_1 <- 0
fig5bdfa <- expand.grid(alpha=alpha,power=power,nclusters=nclusters,nsubjects=nsubjects,
                        d=d,icc=icc,rho_c=rho_c,rho_s=rho_s,vart=vart,
                        ltf_0=ltf_0,ltf_1=ltf_1,gtf_0=gtf_0,gtf_1=gtf_1) %>%
  mutate(power = pmap_dbl(.,cpa.did.normal.ltf)) %>%rename(Power = power) %>%
  mutate(Design = "Non-replacement")

gtf_0 <- 0.50
gtf_1 <- 0.50
fig5bdfb <- expand.grid(alpha=alpha,power=power,nclusters=nclusters,nsubjects=nsubjects,
                        d=d,icc=icc,rho_c=rho_c,rho_s=rho_s,vart=vart,
                        ltf_0=ltf_0,ltf_1=ltf_1,gtf_0=gtf_0,gtf_1=gtf_1) %>%
  mutate(power = pmap_dbl(.,cpa.did.normal.ltf)) %>%
  rename(Power = power) %>%
  mutate(Design = "Replacement")

fig5bdf <- rbind(fig5bdfa,fig5bdfb) %>%
  mutate(type = "b")

fig5b <- ggplot(fig5bdf,aes(x=nsubjects,y=Power,color=Design)) +
  geom_line() + geom_hline(yintercept = 0.8, lty = "dashed") +
  labs(title = "Figure 5") + scale_y_continuous(limits=c(0.5,.9)) +
  theme(panel.grid.minor.x = element_blank())

fig5df <- rbind(fig5adf,fig5bdf)

to_label <- as_labeller(c(a = 'list(lambda[1] == 0.05, lambda[2] == 0.16)',
                          b = 'lambda[1] == ~ lambda[2] == 0.50'),
                        label_parsed)

caption5 <- "Plot of power as a function of cluster size using parameters given by
Miller et al (2012). The left plot uses observed pre/post loss to follow-up rates
while the right plot uses 50% loss to follow-up."

fig5 <- ggplot(fig5df, aes(x = nsubjects, y = Power, color = Design)) +
  geom_line() + geom_hline(yintercept = 0.8, lty = "dashed") +
  facet_grid(.~type, labeller = to_label) +
  # labs(title = "Figure 5",
  #      x = "Number of Subjects per Cluster",
  #      caption = str_wrap(caption5, width = 90)) +
  # scale_y_continuous(limits=c(0.5,.9)) +
  # theme(legend.position = "bottom",
  #       plot.title = element_text(size = 12),
  #       plot.caption = element_text(size = 10, hjust = 0)) +
  # theme(panel.grid.minor.x = element_blank())
  labs(title = "Moyer Figure 5 TOP",
       x = "Number of Subjects per Cluster") +
  scale_y_continuous(limits=c(0.5,.9)) +
  theme(legend.position = "bottom") +
  theme(panel.grid.minor.x = element_blank())

# fig5 <- ggarrange(fig5a, fig5b,
#                   labels = c("a", "b"), ncol = 2, nrow = 1)
# #ggsave("fig5.png",width=14,height=9.33)
# annotate_figure(fig5, top = "Figure 5")
ggsave("figures/figure5.pdf",fig5,width=6,height=4,dpi=800)
#fig5
