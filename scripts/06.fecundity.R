# ---- FECUNDITY ANALYSIS ----

# Extract relevant columns - same structure as recombination
fecund <- by_vialday[, c("Treatment", "PaternalStock", "MaternalVial", "vial_num", "vial_letter", "num_offspring", "num_moms")]

# Calculate per vial per letter fecundity
fecund$fecundity <- fecund$num_offspring / fecund$num_moms

# Summarize by MaternalVial + vial_letter (mirrors recomb_rate summaryBy)
fecund_by_vial <- summaryBy(fecundity + num_offspring + num_moms ~ 
                              MaternalVial + vial_letter + PaternalStock + Treatment,
                            data = fecund, FUN = sum, na.rm = T)

# Recalculate fecundity after summing (total offspring / total moms)
fecund_by_vial$fecundity <- fecund_by_vial$num_offspring.sum / fecund_by_vial$num_moms.sum

# Save full vial + brood level fecundity
write.csv(fecund_by_vial, "output/fecund_by_vial.csv", row.names = FALSE, quote = FALSE)

# Summary by Treatment + PaternalStock
fecund_summary <- aggregate(fecundity ~ Treatment + PaternalStock, 
                            data = fecund_by_vial, mean)

# Make factors
fecund_by_vial$vial_letter  <- as.factor(fecund_by_vial$vial_letter)
fecund_by_vial$MaternalVial <- as.factor(fecund_by_vial$MaternalVial)

# ---- STATS ----
# Negative Binomial model - mirrors recombination GLMER structure
fit_fec <- glmer.nb(round(fecundity) ~ (1|MaternalVial) + Treatment * PaternalStock * vial_letter, 
                    data = fecund_by_vial)
fecundity_stats <- Anova(fit_fec)
fecundity_stats
write.csv(fecundity_stats, "output/fecundity-stats_anova.csv")

# ---- FIGURES ----
# Overall fecundity by treatment
fecund_figure <- ggplot(fecund_by_vial, aes(x = Treatment, y = fecundity, col = factor(PaternalStock))) +
  xlab("Caloric Density") + ylab("# progeny per mom") + ggtitle("Fecundity") +
  ylim(0, 80) +
  theme(axis.text.x = element_text(angle = 0)) +
  geom_boxplot() + scale_color_discrete(name = "Paternal Stock")
fecund_figure
ggsave("images/fecundity-figure.png")

# Fecundity by brood/day
fecund_figure_day <- ggplot(fecund_by_vial, aes(x = Treatment, y = fecundity, col = factor(PaternalStock))) + 
  facet_wrap(~vial_letter) +
  xlab("Caloric Density") + ylab("# progeny per mom") + ggtitle("Fecundity by Brood") +
  ylim(0, 80) +
  theme(axis.text.x = element_text(angle = 0)) +
  geom_boxplot() + scale_color_discrete(name = "Paternal Stock")
fecund_figure_day
ggsave("images/fecundity-figure-byDay.png")

# ---- SUMMARY STATS ----
mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock == "42"])
mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock == "217"])
mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock == "42"  & fecund_by_vial$Treatment == "0.5x"])
mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock == "217" & fecund_by_vial$Treatment == "0.5x"])
mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock == "42"  & fecund_by_vial$Treatment == "1x"])
mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock == "217" & fecund_by_vial$Treatment == "1x"])
mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock == "42"  & fecund_by_vial$Treatment == "2x"])
mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock == "217" & fecund_by_vial$Treatment == "2x"])

# Zero progeny vials
length(fecund_by_vial$MaternalVial[fecund_by_vial$fecundity == 0 & fecund_by_vial$PaternalStock == "217"])
length(fecund_by_vial$MaternalVial[fecund_by_vial$fecundity == 0 & fecund_by_vial$PaternalStock == "42"])