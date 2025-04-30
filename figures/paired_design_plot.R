library(tidyverse)
library(patchwork)
library(kimma)

#Model
model <- kmFit(example.voom, 
               model="~virus + (1|ptID)",
               run_lm = TRUE, run_lme = TRUE)

#Genes only signif when paired
lm_ns <- model$lm %>%
  filter(variable=="virusHRV" & pval > 0.2)
lme_s <- model$lme %>%
  filter(variable=="virus" & pval < 0.05)

gene.OI <- intersect(lm_ns$gene, lme_s$gene)[1]

dat <- as.data.frame(example.voom$E) %>% 
  rownames_to_column() %>% 
  filter(rowname == gene.OI) %>% 
  pivot_longer(-rowname, names_to = "libID") %>% 
  left_join(example.voom$targets)

p1 <- dat %>% 
  ggplot(aes(x=virus, y=value)) +
  geom_jitter(width=0.2, height=0) +
  geom_point(stat="summary", fun="mean", color="red", shape="square") + 
  geom_errorbar(stat="summary", color="red",
                fun.min=function(x) {mean(x)-sd(x)/sqrt(length(x))}, 
                fun.max=function(x) {mean(x)+sd(x)/sqrt(length(x))}) +
  theme_classic() +
  labs(y="Normalized log2 expression", 
       title = paste0("P = ", 
                      round(lm_ns[lm_ns$gene==gene.OI, "pval"], digits=2)))

p2 <- dat %>% 
  ggplot(aes(x=virus, y=value)) +
  geom_point() +
  geom_line(aes(group = ptID)) +
  theme_classic() +
  labs(y="Normalized log2 expression", 
       title = paste0("P = ", 
                      round(lme_s[lme_s$gene==gene.OI, "pval"], digits=3)))

#Randomize ptID
set.seed(42)
ptID2 <- sample(unique(dat$ptID))
dat2 <- dat %>% 
  arrange(virus, ptID) %>% 
  mutate(ptID2 = c(.$ptID[1:6], ptID2))

lme2 <- lme4::lmer(value~virus+(1|ptID2), dat2)
p.lme <- broom::tidy(car::Anova(lme2))

p3 <- dat2 %>% 
  ggplot(aes(x=virus, y=value)) +
  geom_point() +
  geom_line(aes(group = ptID2)) +
  theme_classic() +
  labs(y="Normalized log2 expression", 
       title = paste0("P = ", 
                      round(p.lme$p.value, digits=3)))

p_all <- p1+p2+p3
p_all

ggsave("~/Documents/GitHub/BIGslu/R_packages/kimma_vignette/paired_design.png", p_all, width=6, height=3)
