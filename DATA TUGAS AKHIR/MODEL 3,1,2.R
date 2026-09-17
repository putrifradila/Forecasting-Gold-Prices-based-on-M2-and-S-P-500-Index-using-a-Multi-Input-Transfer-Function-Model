# ============================================================
# Uji White Noise (Ljung-Box) untuk ARIMA(3,1,2)
# Data: M2SL | Lag: 6, 12, 18, 24, 36
# ============================================================

# Install package jika belum ada
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

# ---- Fitting ARIMA(3,1,2) ----
model <- Arima(m2sl, order = c(3, 1, 2))

cat("==========================================\n")
cat("       ARIMA(3,1,2) - Model Summary\n")
cat("==========================================\n")
print(summary(model))

# ---- Ekstrak Residual ----
resid <- residuals(model)

# ---- Uji Ljung-Box pada lag 6, 12, 18, 24, 36 ----
lags <- c(6, 12, 18, 24, 36)

cat("\n==========================================\n")
cat("   Uji Ljung-Box - White Noise Test\n")
cat("==========================================\n")
cat(sprintf("%-6s %-12s %-12s %-20s\n", "Lag", "Statistik Q", "p-value", "Kesimpulan"))
cat(rep("-", 55), "\n", sep = "")

hasil <- data.frame(
  Lag         = integer(),
  Q_statistik = numeric(),
  p_value     = numeric(),
  Kesimpulan  = character(),
  stringsAsFactors = FALSE
)

for (lag in lags) {
  uji <- Box.test(resid, lag = lag, type = "Ljung-Box", fitdf = 5) # fitdf = p+q = 3+2
  p   <- uji$p.value
  Q   <- uji$statistic
  ket <- ifelse(p > 0.05, "White noise (lulus)", "TIDAK white noise")
  
  cat(sprintf("%-6d %-12.4f %-12.4f %-20s\n", lag, Q, p, ket))
  
  hasil <- rbind(hasil, data.frame(
    Lag         = lag,
    Q_statistik = round(Q, 4),
    p_value     = round(p, 4),
    Kesimpulan  = ket,
    stringsAsFactors = FALSE
  ))
}

cat(rep("-", 55), "\n", sep = "")
cat("Catatan: H0 = residual white noise | alpha = 0.05\n")
cat("         fitdf = 5 (p=3, q=2) sesuai orde ARIMA(3,1,2)\n\n")

# ---- Plot Diagnostik ----
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# 1. Residual plot
plot(resid, type = "l", main = "Residual ARIMA(3,1,2)",
     ylab = "Residual", xlab = "Waktu", col = "#185FA5")
abline(h = 0, col = "red", lty = 2)

# 2. ACF residual
acf(resid, lag.max = 36, main = "ACF Residual", col = "#185FA5")

# 3. PACF residual
pacf(resid, lag.max = 36, main = "PACF Residual", col = "#185FA5")

# 4. Histogram residual
hist(resid, breaks = 15, main = "Distribusi Residual",
     xlab = "Residual", col = "#B5D4F4", border = "white")
curve(dnorm(x, mean = mean(resid), sd = sd(resid)) * length(resid) * diff(range(resid)) / 15,
      add = TRUE, col = "#185FA5", lwd = 2)

par(mfrow = c(1, 1))

# ---- p-value plot ----
plot(hasil$Lag, hasil$p_value,
     type = "b", pch = 19, col = ifelse(hasil$p_value > 0.05, "#1D9E75", "#D85A30"),
     ylim = c(0, 1), xlab = "Lag", ylab = "p-value",
     main = "Ljung-Box p-value per Lag (ARIMA 3,1,2)",
     xaxt = "n")
axis(1, at = lags)
abline(h = 0.05, col = "red", lty = 2, lwd = 1.5)
text(lags, hasil$p_value + 0.04, labels = round(hasil$p_value, 3), cex = 0.8)
legend("bottomright", legend = c("p-value", "alpha = 0.05"),
       col = c("#1D9E75", "red"), lty = c(1, 2), pch = c(19, NA), cex = 0.8)

# ---- Tampilkan tabel hasil akhir ----
cat("==========================================\n")
cat("          Tabel Hasil Lengkap\n")
cat("==========================================\n")
print(hasil, row.names = FALSE)