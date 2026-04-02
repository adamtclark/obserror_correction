# The following code replicates Fig. S1 in the manuscript

#setwd("~/Dropbox/Projects/117_ObservationError/src")
rm(list = ls())

d = read.csv("data/NutNet_ReSurvey_processed_250821.csv")

cover = c(d$cover_1, d$cover_2, d$cover_3, d$cover_4)/100
cover = cover[cover!=0]
cover = cover[!is.na(cover)]
#plot(density(cover, from = 0, to = 1, bw = 0.05), xlab = "Cover Value", main = "")

pdf("figures/cover_hist.pdf", width = 8, height = 4)
#png("figures/cover_hist.png", width = 8, height = 4, units = "in", res = 200)
par(mar=c(4,4,2,2), mfrow = c(1,2))
hist(cover, xlab = "Cover Value, Fraction", main = "")
box()
title("A. All Species", adj = 0)

tb = table(cover[cover <= 0.01])
par(mar=c(4.8,4,6,2))
tmp = barplot(tb, log = "y", names.arg = FALSE, axes = FALSE)
box()
axis(2, las=2)
axis(1, at = tmp, names(tb), las = 2)
axis(3, at = tmp, as.numeric(names(tb))*100, las = 2)
mtext(side = 1, line = 3.8, "Cover Value, Fraction")
mtext(side = 3, line = 3, "Cover, Percent")
mtext(side = 2, line = 2.9, "Frequency, Log Scale")
title("B. All Species with Cover ≤ 1%", adj = 0, line = 4.7)
dev.off()