# set options
options(stringsAsFactors = F)
Sys.setlocale("LC_ALL", locale="de_DE")
setwd("/Users/stefanhartmann/Dropbox/Privat/SiGS/MultivariateModelGesamtkorpus/model_final")

# install / load packages
sapply(c("dplyr", "tidyr", "lme4", "lmerTest", "gvlma", "Hmisc", "lattice", "knitr", "ggplot2", "effects"), 
       function(x) if(!is.element(x, installed.packages())) install.packages(x, dependencies = T))
lapply(list("dplyr", "tidyr", "lme4", "lmerTest", "gvlma", "Hmisc", "lattice", "knitr", "ggplot2", "effects"), 
       require, character.only=T)


# read data
kwic <- read.table("sigsN.csv", sep="\t", head=T, encoding="UTF-8")


# HELPER FUNCTIONS --------------------------------------------------
# function for getting variance inflation factors
# adopted from Stackoverflow
# Link: http://stackoverflow.com/questions/26633483/collinearity-after-accounting-for-random-mixed-effects
# User: colin (http://stackoverflow.com/users/2777850/colin)
vif.lme <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v 
}

# function for getting model evaluation
get.c <- function(mymodel, mybelege=kwic) {
  return(somers2(fitted(mymodel), mybelege$upperCase))
}


# DATA MANIPULATION ----------------------------------------------------

# Animacy as factor
kwic$Animacy3 <- as.factor(kwic$Animacy3)
kwic$Animacy <- as.factor(kwic$Animacy)



# MODEL FITTING -----------------------------------------------------

# fit model
m1 <- glmer(upperCase~Animacy+LogFreq+Animacy:LogFreq+
              (1+Animacy|Place)+(1|graphtoken_lemma),
             data=kwic,family="binomial", glmerControl(optimizer = "bobyqa"))
summary(m1)


# MODEL CHECKING ----------------------------------------------------

get.c(m1)

# do random slopes lead to a better model fit?
m1_no_slo <- glmer(upperCase~factor(Animacy)+log10(LemmaFreq)+factor(Animacy):log10(LemmaFreq)+
                     (1|Place)+(1|graphtoken_lemma),
                   data=kwic,family="binomial", glmerControl(optimizer = "bobyqa"))

anova(m1, m1_no_slo) # yes


# do we need the random effects?
m1_no_rand1 <- glmer(upperCase~factor(Animacy)+log10(LemmaFreq)+factor(Animacy):log10(LemmaFreq)+
                     (1|graphtoken_lemma),
                   data=kwic,family="binomial", glmerControl(optimizer = "bobyqa"))

m1_no_rand2 <- glmer(upperCase~factor(Animacy)+log10(LemmaFreq)+factor(Animacy):log10(LemmaFreq)+
                     (1+Animacy|Place),
                   data=kwic,family="binomial", glmerControl(optimizer = "bobyqa"))

anova(m1, m1_no_rand2) # yes
anova(m1, m1_no_rand1) # yes


# do we need the interaction Frequency - Animacy?
m1_no_int <- glmer(upperCase~Animacy+log10(LemmaFreq)+
                     (1+Animacy|Place)+(1|graphtoken_lemma),
                   data=kwic,family="binomial", glmerControl(optimizer = "bobyqa"))

anova(m1, m1_no_int) #yes




# null model without animacy and frequency
m0 <- glmer(upperCase~
              (1|Place)+(1|graphtoken_lemma),
            data=kwic,family="binomial", glmerControl(optimizer = "bobyqa"))

anova(m1, m0)

# do we need animacy?
m0_no_anim <- glmer(upperCase~LogFreq+
                      (1|Place)+(1|graphtoken_lemma),
                    data=kwic,family="binomial", glmerControl(optimizer = "bobyqa"))

anova(m1, m0_no_anim)

# explorative visualization
dotplot(ranef(m1))
dotplot(ranef(m1))[1]

plot(effect("Animacy", m1), main = "Animacy", ylab = "uppercase", xlab="Animacy")
plot(effect("Animacy:LogFreq", m1))

table(kwic$Animacy, kwic$upperCase) %>% t() %>% 
  prop.table(mar=2) %>% barplot(names.arg=c("abstract", "concrete", "animal", "human", "superhuman"), cex.names=0.5,
                                main="Animacy and capitalization")

# predicted vs. fitted
kwic$predicted <- round(fitted(m1, 0))

comp_kwic <- kwic %>% group_by(graphtoken_lemma) %>% summarise(
  Lemma=unique(graphtoken_lemma),
  Observed_upper=length(which(upperCase==1)),
  Pred_upper = sum(predicted),
  # Pred_upper=unique(predicted),
  Freq=n(),
  Diff = Observed_upper - Pred_upper,
  Diff2 = (Observed_upper - Pred_upper) / n()
)


# highest deviations between predicted and observed forms
head(comp_kwic[order(comp_kwic$Diff),], 20)
tail(comp_kwic[order(comp_kwic$Diff),], 20)

comp_kwic[order(comp_kwic$Diff2),] %>% filter(Freq>5) %>% head(20)
comp_kwic[order(comp_kwic$Diff2),] %>% filter(Freq>5) %>%  tail(20)

