## Requires a PCC object from sensitivity analysis run

# Define color gradient
pal <- colorRampPalette(c("brown1", "wheat", "steelblue"))
x <- PCC[,1]
normalized_x <- (x - min(x)) / diff(range(x))
colors <- pal(length(x))[as.numeric(cut(normalized_x, breaks = length(pal(length(x))), include.lowest = TRUE))]


barplot(PCC[,1], horiz = TRUE, xlim = c(-0.61,0.61), xaxt='n', ann=FALSE, yaxt='n',
        xlab = "Parameter Value", ylab = "PRCC", col = colors)

abline(v=0, lwd = 2)
yes <- c(0.7, 1.9, 3.1, 4.3, 5.5, 6.7, 7.9, 9.1, 10.3, 11.5, 12.7, 13.9,
         15.1, 16.3, 17.5, 18.7, 19.9, 21.1)  # hand-roll plot spots?
segments(y0=yes, y1=yes, x0 = PCC[,4], x1=PCC[,5])        # horizonta whiskers
segments(y0=yes-.1, y1=yes+.1, x0 = PCC[,4], x1=PCC[,4])  # bottom horizontal
segments(y0=yes-.1, y1=yes+.1, x0 = PCC[,5], x1=PCC[,5])  # top horizontal

## NOTE: HARD-CODED!!!
labels = expression(nu[H], omega, alpha[H], nu[M], mu[A], alpha[M], tau,
                    rho[RM], beta[M], rho[RH], rho[VM], mu[E], beta[AM],
                    beta[MA], beta[HA], rho[VH], zeta, beta[AH])

text(y = yes, x = c(rep(0.06, 7), rep(0.15, 3), rep(-0.06, 8)), labels = labels)
axis(1, at = seq(-0.6, 0.6, by = 0.2), label = seq(-0.6, 0.6, by = 0.2),
     las = 1)
