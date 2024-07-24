library(tidyverse)
library(readxl)
library(logKDE)

chk_pcb <- read_xlsx("data/CB_Stilly_Puyallup_monitoring_data.xlsx") |>
  separate(Species, c("lifestage", "species")) |>
  mutate(river = str_to_title(RiverSystem),
         unit = "ng/g (wet)") |>
  select(river, species, lifestage,
         pcbs = `Max of CONC_FOUND`,
         unit) |>
  mutate(gt_lla = pcbs >= 100)


p_mort <- function(p_exp, mort_reg) {

}


mort_qreg <- function(pcb) {
  pmax(0.1702 + 0.221 * log10(pcb), 0)
}

growth_qreg <- function(pcb) {
  pmax(0.15 + 0.0938 * log10(pcb), 0)
}

## f(d) 5 [exp(21.69 1. The0.0329d)] / [1 1 exp(21.69 1 0.0329d)
## growth_surv <- function()


puy_pcb <- chk_pcb |>
  filter(river == "Puyallup",
         lifestage == "juvenile") |>
  mutate(pcbs = 1e-3 * pcbs,
         unit = "μg/g (wet)")

puy_pcb_dens <- density(puy_pcb$pcbs)

plot(puy_pcb_dens, ylim = c(0, 0.05))
hist(puy_pcb$pcbs, probability = TRUE, add = TRUE, color = NA)

ld <- logdensity(puy_pcb$pcbs, from = 1e-7, to = 500e-3)

hist(puy_pcb$pcbs, breaks = seq(0, 200e-3, 5e-3), probability = TRUE)
lines(ld)
## lines(puy_pcb_dens)
##

plot(ld, xlim = c(0, 0.2))
hist(puy_pcb$pcbs, breaks = seq(0, 200e-3, 5e-3), probability = TRUE, add = TRUE)
abline(v = 0.1, lty = "dashed")

ld_fun <- approxfun(ld)

integrate(\(x) ld_fun(x) * mort_qreg(x), 1e-7, 0.2)

puy_pcb |>
  arrange(pcbs) |>
  mutate(qtl = quantile(pcbs))

boot_mort <- function(n, pcbs) {
  s <- sample(pcbs, length(pcbs), replace = TRUE)
  mort_q
}

boot_ex <- replicate(10000, sample(puy_pcb$pcbs, replace = TRUE), simplify = FALSE)
boot_mort <- map_dbl(boot_ex, \(x) mean(mort_qreg(x)))

## boot_ld <- map(boot_ex, logdensity, from = 1e-10, to = 500e-3)
## boot_ldf <- map(boot_ld, approxfun)
## boot_int <- map(boot_ldf, \(ldf) integrate(\(x) ldf(x) * mort_qreg(x), 1e-10, 200e-3))
## boot_surv <- map_dbl(boot_int, pluck, "value")


hist(boot_surv)
hist(boot_sum)

plot(density(boot_surv))
lines(density(boot_sum), col = "red")
curve(dnorm(x, mean(boot_sum), sd(boot_sum)), add = TRUE)

plot(boot_ld[[1]], col = rgb(0, 0, 0, 0.1), ylim = c(0, 100))
walk(boot_ld[2:200], lines, col = rgb(0, 0, 0, 0.1))
abline(v= 0.1, lty = "dashed")

chk_pcb |>
  filter(lifestage == "juvenile") |>
  mutate(river = ifelse(river == "Puyallup", "Puyallup/White", river)) |>
ggplot(aes(x = pcbs)) +
  geom_histogram(breaks = seq(0, 175, 5)) +
  ## geom_density(bounds = c(0, Inf)) +
  geom_vline(xintercept = 100, linetype = "dashed") +
  facet_wrap(~ river, nrow = 2, scales = "free_y") +
  scale_y_continuous(expand = expansion(c(0, 0.05), 0)) +
  scale_x_continuous(expand = expansion(0, 0)) +
  labs(x = "Total PCBs (ng/g ww)",
       y = "Number of observations") +
  theme_bw()
ggsave("PCB_obs.png")
