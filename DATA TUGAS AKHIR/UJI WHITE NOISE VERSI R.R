# ============================================================
# Uji White Noise (Ljung-Box) - Semua Model ARIMA
# Model  : ARIMA(3,1,0), ARIMA(3,1,2), ARIMA(1,1,1)
# Data   : M2SL
# Lag    : 6, 12, 18, 24, 36
# ============================================================

# install.packages("forecast")
library(forecast)

# ---- Data M2SL ----
m2sl <- c(
  13900.3, 13949.3, 14012.7, 14035.9, 14092.4, 14137.2, 14174.6, 14223.2,
  14245.5, 14256.6, 14282.8, 14386.8, 14459.7, 14506.3, 14557.5, 14594.7,
  14704.4, 14802.8, 14886.8, 14974.4, 15047.0, 15184.1, 15296.0, 15347.5,
  15427.1, 15493.4, 16033.0, 17065.3, 17932.9, 18184.3, 18335.4, 18411.2,
  18602.5, 18757.6, 18997.0, 19115.7, 19379.0, 19641.8, 19873.2, 20176.9,
  20460.4, 20477.8, 20633.3, 20840.9, 20981.0, 21160.6, 21335.5, 21498.5,
  21650.1, 21725.2, 21787.2, 21769.4, 21719.9, 21651.5, 21652.9, 21638.4,
  21536.7, 21461.7, 21406.7, 21290.9, 21279.9, 21243.7, 20947.4, 20761.4,
  20829.6, 20799.4, 20791.4, 20777.5, 20747.8, 20732.9, 20746.4, 20777.8,
  20843.7, 20926.9, 20973.2, 20958.8, 21019.9, 21067.7, 21094.4, 21185.8
)

# ---- Definisi model ----
models <- list(
  "ARIMA(3,1,0)" = list(order = c(3, 1, 0), fitdf = 3),
  "ARIMA(3,1,2)" = list(order = c(3, 1, 2), fitdf = 5),
  "ARIMA(1,1,1)" = list(order = c(1, 1, 1), fitdf = 2)
)

lags <- c(6, 12, 18, 24, 36)

# ---- Fungsi uji Ljung-Box ----
uji_ljungbox <- function(resid, lags, fitdf, nama_model) {
  cat("\n------------------------------------------\n")
  cat(sprintf(" Model: %s\n", nama_model))
  cat("------------------------------------------\n")
  cat(sprintf("%-6s %-14s %-12s %-20s\n", "Lag", "Q-statistik", "p-value", "Kesimpulan"))
  cat(rep("-", 56), "\n", sep = "")
  
  rows <- list()
  for (lag in lags) {
    uji <- Box.test(resid, lag = lag, type = "Ljung-Box", fitdf = fitdf)
    p   <- uji$p.value
    Q   <- uji$statistic
    ket <- ifelse(p > 0.05, "White noise (lulus)", "TIDAK white noise")
    cat(sprintf("%-6d %-14.4f %-12.4f %-20s\n", lag, Q, p, ket))
    rows[[length(rows) + 1]] <- data.frame(
      Model = nama_model, Lag = lag,
      Q_stat = round(Q, 4), p_value = round(p, 4),
      Kesimpulan = ket, stringsAsFactors = FALSE
    )
  }
  cat(rep("-", 56), "\n", sep = "")
  do.call(rbind, rows)
}

# ============================================================
# FITTING & UJI SEMUA MODEL
# ============================================================
cat("==========================================\n")
cat("   UJI WHITE NOISE LJUNG-BOX - M2SL\n")
cat("   Lag: 6, 12, 18, 24, 36 | alpha=0.05\n")
cat("==========================================\n")

semua_hasil <- list()
fitted_models <- list()

for (nama in names(models)) {
  cfg   <- models[[nama]]
  model <- Arima(m2sl, order = cfg$order)
  fitted_models[[nama]] <- model
  
  cat(sprintf("\n[%s] AIC=%.2f | BIC=%.2f\n", nama, AIC(model), BIC(model)))
  
  hasil <- uji_ljungbox(residuals(model), lags, cfg$fitdf, nama)
  semua_hasil[[nama]] <- hasil
}

# ============================================================
# REKAPITULASI
# ============================================================
rekap <- do.call(rbind, semua_hasil)
rownames(rekap) <- NULL

cat("\n==========================================\n")
cat("         REKAPITULASI SEMUA MODEL\n")
cat("==========================================\n")
print(rekap, row.names = FALSE)

# Hitung jumlah lag yang lulus per model
cat("\n--- Ringkasan Kelulusan White Noise ---\n")
for (nama in names(models)) {
  sub   <- rekap[rekap$Model == nama, ]
  lulus <- sum(sub$p_value > 0.05)
  cat(sprintf("%-15s : %d/%d lag lulus white noise\n", nama, lulus, length(lags)))
}

# ============================================================
# PLOT DIAGNOSTIK - ACF & PACF RESIDUAL
# ============================================================
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

warna <- c("ARIMA(3,1,0)" = "#185FA5",
           "ARIMA(3,1,2)" = "#1D9E75",
           "ARIMA(1,1,1)" = "#D85A30")

for (nama in names(models)) {
  resid <- residuals(fitted_models[[nama]])
  acf(resid,  lag.max = 36, main = paste("ACF Residual -",  nama), col = warna[nama])
  pacf(resid, lag.max = 36, main = paste("PACF Residual -", nama), col = warna[nama])
}

par(mfrow = c(1, 1))

# ============================================================
# PLOT P-VALUE LJUNG-BOX SEMUA MODEL (satu grafik)
# ============================================================
plot(NULL, xlim = c(4, 38), ylim = c(0, 1.05),
     xlab = "Lag", ylab = "p-value",
     main = "Ljung-Box p-value - Semua Model ARIMA",
     xaxt = "n")
axis(1, at = lags)
abline(h = 0.05, col = "red", lty = 2, lwd = 1.5)

pch_list <- c(19, 17, 15)
lty_list <- c(1, 2, 3)

for (i in seq_along(names(models))) {
  nama <- names(models)[i]
  sub  <- rekap[rekap$Model == nama, ]
  lines(sub$Lag, sub$p_value, col = warna[nama], lty = lty_list[i], lwd = 1.5)
  points(sub$Lag, sub$p_value, col = warna[nama], pch = pch_list[i], cex = 1.2)
}

legend("bottomright",
       legend = c(names(models), "alpha = 0.05"),
       col    = c(unname(warna), "red"),
       lty    = c(lty_list, 2),
       pch    = c(pch_list, NA),
       cex    = 0.85, bg = "white")