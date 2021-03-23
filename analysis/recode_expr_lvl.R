setwd("~/Documents/nathansonlab/tecac-manuscript/")

#pluta 12/10/20

library(ggplot2)
library(dplyr)

# file containing genes and expression levels
dat <- read.table("gene_expr_lvl.txt", header = T)
dat$Fetal.testis <- as.factor(dat$Fetal.testis)



# make the tertiles and attach
dat <- dat %>% mutate(tertiles = ntile(TPM.counts, 3)) %>% 
              mutate(tertiles = if_else(tertiles == 1, 
                            'No', if_else(tertiles == 2, "Low", "High"))) %>% arrange(TPM.counts)

dat <- dat[ order(dat$Scored.Gene..closest.gene.bolded.),]

p1 <- ggplot(data = dat, aes(x = 1:dim(dat)[1], y = TPM.counts, color = as.factor(tertiles))) + 
  geom_point() +
  xlab("Genes") + 
  guides(fill = guide_legend("Expression Level")) +
  geom_hline(yintercept = 698) +
  geom_hline(yintercept = 2348) +
  theme_minimal()

png("expression_tertiles.png")
print(p1)
dev.off()