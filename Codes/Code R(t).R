###Méthode de Cori et al.(2013) pour estimer le NOMBRE DE REPRODUCTION EFFECTIF(Rt):

#Installer et charger les packages nécessaires:
if (!require(EpiEstim)) install.packages("EpiEstim")
if (!require(lubridate)) install.packages("lubridate")
library(EpiEstim)
library(lubridate)

#Lire le fichier CSV:
tab <- read.csv("C:/Users/milic/OneDrive/ULB 2025-2026/Q2/BINGF432/TP/Test genomic CSV.csv",
                stringsAsFactors = FALSE)

#Convertir la colonne de dates au format Date:
tab$collection_date <- as.Date(tab$collection_date)

# Trier les données par ordre chronologique:
tab <- tab[order(tab$collection_date), ]

# Transformer les dates en jours numérotés depuis la première observation:
days <- as.integer(tab$collection_date - min(tab$collection_date)) + 1
total_number_of_days <- as.integer(max(tab$collection_date) - min(tab$collection_date)) + 1

# Compter le nombre de cas observés chaque jour:
daily_cases <- rep(0, total_number_of_days)
for (i in seq_along(daily_cases)) {
  daily_cases[i] <- sum(days == i)
}

# Définir les paramètres de simulation pour l'estimation de Rt:
n <- 1000
mean_range <- c(9, 23)
sd_range <- c(4, 8)
t_start <- seq(2, length(daily_cases) - 6)
t_end <- seq(8, length(daily_cases))
all_Rt <- matrix(NA, nrow = length(t_start), ncol = n)

# Répéter l’estimation de Rt avec différentes valeurs aléatoires de l’intervalle sériel:
for (i in 1:n) {
  mean_si_i <- runif(1, mean_range[1], mean_range[2])
  sd_si_i <- runif(1, sd_range[1], sd_range[2])
  
  res_i <- estimate_R(
    incid = daily_cases,
    method = "parametric_si",
    config = make_config(list(
      mean_si = mean_si_i,
      std_si = sd_si_i,
      t_start = t_start,
      t_end = t_end
    ))
  )
  
  all_Rt[, i] <- res_i$R$`Mean(R)`
}

# Calculer la médiane et les bornes de l’intervalle de confiance de Rt:
R_median <- apply(all_Rt, 1, median, na.rm = TRUE)
R_lower <- apply(all_Rt, 1, quantile, probs = 0.025, na.rm = TRUE)
R_upper <- apply(all_Rt, 1, quantile, probs = 0.975, na.rm = TRUE)

# Associer chaque estimation de Rt à une date:
Rt_days <- floor((t_start + t_end) / 2)
R_dates <- min(tab$collection_date) + Rt_days

# Ouvrir une fenêtre graphique et définir la mise en forme du graphique:
windows(width = 18, height = 6)
par(oma=c(0,0,0,0), mar=c(5.0,5.0,1.0,1.0), lwd=0.3, bty="o",
    col="gray30", col.axis="gray30", fg="gray30")

# Initialiser le graphique de Rt dans le temps:
plot(R_dates, R_median, lwd=0.7, type="n", axes=FALSE, xlab=NA, ylab=NA,
     main="Evolution de Rt avec propagation de l'incertitude du SI",
     xlim=c(min(R_dates), max(R_dates)),
     ylim=c(min(R_lower, na.rm = TRUE) - 0.2, max(R_upper, na.rm = TRUE) + 0.3),
     xaxs="i")

# Ajouter la zone d’incertitude autour de Rt:
xx_l <- c(R_dates, rev(R_dates))
yy_l <- c(R_lower, rev(R_upper))
polygon(xx_l, yy_l, col=rgb(187/255, 187/255, 187/255, 0.35), border=0)

# Tracer la courbe médiane de Rt:
lines(R_dates, R_median, lwd=0.8, col="gray30")

# Ajouter une ligne de référence à Rt = 1:
abline(h=1, lty=2, lwd=0.5, col="red")

# Ajouter les axes et le label de l’axe Y:
axis.Date(1, at=seq.Date(from=min(R_dates), to=max(R_dates), by="2 weeks"),
          format="%Y-%m", las=2, cex.axis=0.7)
axis(side = 2, lwd = 0.5, cex.axis = 0.7, mgp = c(0, 0.6, 0), lwd.tick = 0.5,
     col.lab = "gray30", col = "gray30", tck = -0.03, las = 1, padj = 0.4,
     at = seq(0, ceiling(max(R_upper, na.rm = TRUE)), by = 1))
mtext("Effective reproduction number (Rt)", side=2, col="gray30", cex=0.9, line=1.7)





