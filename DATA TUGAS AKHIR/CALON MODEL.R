# ============================================================
# Uji White Noise (Ljung-Box) + Signifikansi Parameter
# Model  : ARIMA(3,1,0), (3,1,1), (3,1,2),
#          ARIMA(1,1,1), (2,1,1), (2,1,2)
# Data   : M2SL_1
# Lag    : 6, 12, 18, 24, 36
# ============================================================

# install.packages("forecast")
library(forecast)

# ---- Data M2SL_1 ----
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

# ---- Definisi model (nama, order, fitdf=p+q) ----
models <- list(
  "ARIMA(3,1,0)" = list(order = c(3, 1, 0), fitdf = 3),
  "ARIMA(3,1,1)" = list(order = c(3, 1, 1), fitdf = 4),
  "ARIMA(3,1,2)" = list(order = c(3, 1, 2), fitdf = 5),
  "ARIMA(1,1,1)" = list(order = c(1, 1, 1), fitdf = 2),
  "ARIMA(2,1,1)" = list(order = c(2, 1, 1), fitdf = 3),
  "ARIMA(2,1,2)" = list(order = c(2, 1, 2), fitdf = 4)
)

lags <- c(6, 12, 18, 24, 36)

# ============================================================
# FUNGSI SIGNIFIKANSI PARAMETER
# ============================================================
cek_parameter <- function(model, nama_model) {
  coef_val <- coef(model)
  se_val   <- sqrt(diag(vcov(model)))
  t_val    <- coef_val / se_val
  p_val    <- 2 * (1 - pnorm(abs(t_val)))
  
  cat(sprintf("\n  Parameter ARIMA %s:\n", nama_model))
  cat(sprintf("  %-12s %10s %10s %10s %10s %s\n",
              "Parameter", "Estimasi", "Std.Error", "t-hitung", "p-value", "Signifikan"))
  cat("  ", rep("-", 68), "\n", sep = "")
  
  params_df <- data.frame(
    Model     = nama_model,
    Parameter = names(coef_val),
    Estimasi  = round(coef_val, 4),
    StdError  = round(se_val, 4),
    t_hitung  = round(t_val, 4),
    p_value   = round(p_val, 4),
    Signifikan = ifelse(p_val < 0.001, "***",
                        ifelse(p_val < 0.01,  "**",
                               ifelse(p_val < 0.05,  "*",
                                      ifelse(p_val < 0.1,   ".",  "Tidak")))),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(params_df))) {
    cat(sprintf("  %-12s %10.4f %10.4f %10.4f %10.4f %s\n",
                params_df$Parameter[i],
                params_df$Estimasi[i],
                params_df$StdError[i],
                params_df$t_hitung[i],
                params_df$p_value[i],
                params_df$Signifikan[i]))
  }
  cat("  Signif: *** p<0.001  ** p<0.01  * p<0.05  . p<0.1\n")
  
  invisible(params_df)
}

# ============================================================
# FUNGSI UJI LJUNG-BOX
# ============================================================
uji_ljungbox <- function(resid, lags, fitdf, nama_model) {
  cat(sprintf("\n  Ljung-Box ARIMA %s (fitdf=%d):\n", nama_model, fitdf))
  cat(sprintf("  %-6s %-14s %-12s %-20s\n", "Lag", "Q-statistik", "p-value", "Kesimpulan"))
  cat("  ", rep("-", 54), "\n", sep = "")
  
  rows <- list()
  for (lag in lags) {
    uji <- Box.test(resid, lag = lag, type = "Ljung-Box", fitdf = fitdf)
    p   <- uji$p.value
    Q   <- uji$statistic
    ket <- ifelse(p > 0.05, "White noise (lulus)", "TIDAK white noise")
    cat(sprintf("  %-6d %-14.4f %-12.4f %-20s\n", lag, Q, p, ket))
    rows[[length(rows) + 1]] <- data.frame(
      Model      = nama_model,
      Lag        = lag,
      Q_stat     = round(Q, 4),
      p_value    = round(p, 4),
      Kesimpulan = ket,
      stringsAsFactors = FALSE
    )
  }
  cat("  ", rep("-", 54), "\n", sep = "")
  do.call(rbind, rows)
}

# ============================================================
# LOOP UTAMA: FITTING, PARAMETER, WHITE NOISE
# ============================================================
cat("============================================================\n")
cat("  CEK PARAMETER & WHITE NOISE - M2SL_1\n")
cat("  Model: ARIMA(3,1,0) (3,1,1) (3,1,2) (1,1,1) (2,1,1) (2,1,2)\n")
cat("  Lag  : 6, 12, 18, 24, 36 | alpha = 0.05\n")
cat("============================================================\n")

rekap_lb     <- list()
rekap_params <- list()
fitted_list  <- list()
ringkasan    <- data.frame()

for (nama in names(models)) {
  cfg   <- models[[nama]]
  model <- Arima(m2sl, order = cfg$order)
  fitted_list[[nama]] <- model
  
  rmse <- sqrt(mean(residuals(model)^2))
  
  cat("\n============================================================\n")
  cat(sprintf(" %s\n", nama))
  cat(sprintf(" AIC=%.4f | BIC=%.4f | RMSE=%.4f\n",
              AIC(model), BIC(model), rmse))
  cat("============================================================\n")
  
  # Parameter
  p_df <- cek_parameter(model, nama)
  rekap_params[[nama]] <- p_df
  
  # White noise
  lb_df <- uji_ljungbox(residuals(model), lags, cfg$fitdf, nama)
  rekap_lb[[nama]] <- lb_df
  
  # Ringkasan
  lulus <- sum(lb_df$p_value > 0.05)
  sig_param <- sum(p_df$p_value < 0.05 & p_df$Parameter != "intercept")
  tot_param  <- sum(p_df$Parameter != "intercept")
  
  ringkasan <- rbind(ringkasan, data.frame(
    Model          = nama,
    AIC            = round(AIC(model), 4),
    BIC            = round(BIC(model), 4),
    RMSE           = round(rmse, 4),
    WN_lulus       = paste0(lulus, "/", length(lags)),
    Param_signif   = paste0(sig_param, "/", tot_param),
    stringsAsFactors = FALSE
  ))
}

# ============================================================
# REKAPITULASI AKHIR
# ============================================================
cat("\n============================================================\n")
cat("                  REKAPITULASI AKHIR\n")
cat("============================================================\n")
cat(sprintf("%-15s %10s %10s %10s %12s %14s\n",
            "Model", "AIC", "BIC", "RMSE", "WN Lulus", "Param Signif"))
cat(rep("-", 75), "\n", sep = "")
for (i in seq_len(nrow(ringkasan))) {
  r <- ringkasan[i, ]
  cat(sprintf("%-15s %10.4f %10.4f %10.4f %12s %14s\n",
              r$Model, r$AIC, r$BIC, r$RMSE, r$WN_lulus, r$Param_signif))
}
cat(rep("-", 75), "\n", sep = "")
cat("WN Lulus     : jumlah lag yang lulus uji white noise (dari 5 lag)\n")
cat("Param Signif : jumlah parameter signifikan p<0.05 (diluar intercept)\n")

# ============================================================
# PLOT P-VALUE LJUNG-BOX SEMUA MODEL
# ============================================================
rekap_lb_all <- do.call(rbind, rekap_lb)

warna <- c(
  "ARIMA(3,1,0)" = "#185FA5",
  "ARIMA(3,1,1)" = "#1D9E75",
  "ARIMA(3,1,2)" = "#D85A30",
  "ARIMA(1,1,1)" = "#993556",
  "ARIMA(2,1,1)" = "#BA7517",
  "ARIMA(2,1,2)" = "#534AB7"
)
pch_list <- c(19, 17, 15, 18, 8, 10)
lty_list <- c(1, 2, 3, 4, 5, 6)

plot(NULL, xlim = c(4, 38), ylim = c(0, 1.05),
     xlab = "Lag", ylab = "p-value",
     main = "Ljung-Box p-value - Semua Model ARIMA",
     xaxt = "n")
axis(1, at = lags)
abline(h = 0.05, col = "red", lty = 2, lwd = 1.5)

for (i in seq_along(names(models))) {
  nama <- names(models)[i]
  sub  <- rekap_lb_all[rekap_lb_all$Model == nama, ]
  lines(sub$Lag, sub$p_value, col = warna[nama], lty = lty_list[i], lwd = 1.5)
  points(sub$Lag, sub$p_value, col = warna[nama], pch = pch_list[i], cex = 1.2)
}

legend("bottomright",
       legend = c(names(models), "alpha = 0.05"),
       col    = c(unname(warna), "red"),
       lty    = c(lty_list, 2),
       pch    = c(pch_list, NA),
       cex    = 0.75, bg = "white")

# ============================================================
# PLOT ACF & PACF RESIDUAL SEMUA MODEL (2 kolom x 6 baris)
# ============================================================
par(mfrow = c(6, 2), mar = c(3, 4, 2.5, 1))
for (nama in names(models)) {
  resid <- residuals(fitted_list[[nama]])
  acf(resid,  lag.max = 36, main = paste("ACF -",  nama), col = warna[nama])
  pacf(resid, lag.max = 36, main = paste("PACF -", nama), col = warna[nama])
}
par(mfrow = c(1, 1))
# ============================================================
# CEK STASIONERITAS & INVERTIBILITAS AR & MA
# ============================================================

# Fungsi untuk mengecek akar-akar polinomial AR dan MA
cek_stasioner_invertibel <- function(model, nama_model) {
  # Ekstrak koefisien AR dan MA
  coef_ar <- NULL
  coef_ma <- NULL
  
  if (length(model$coef) > 0) {
    # Ambil koefisien AR (dimulai dengan "ar")
    ar_names <- grep("^ar", names(model$coef), value = TRUE)
    if (length(ar_names) > 0) {
      coef_ar <- model$coef[ar_names]
    }
    
    # Ambil koefisien MA (dimulai dengan "ma")
    ma_names <- grep("^ma", names(model$coef), value = TRUE)
    if (length(ma_names) > 0) {
      coef_ma <- model$coef[ma_names]
    }
  }
  
  cat(sprintf("\n  --- STASIONERITAS & INVERTIBILITAS: %s ---\n", nama_model))
  
  # Cek stasioneritas (akar-akar polinomial AR)
  if (!is.null(coef_ar) && length(coef_ar) > 0) {
    # Polinomial AR: 1 - phi1*B - phi2*B^2 - ... - phip*B^p
    # Akar-akar dari persamaan karakteristik: 1 - phi1*z - phi2*z^2 - ... = 0
    ar_poly <- c(1, -coef_ar)
    ar_roots <- tryCatch(polyroot(ar_poly), error = function(e) NULL)
    
    if (!is.null(ar_roots)) {
      ar_mod <- Mod(ar_roots)
      cat("  Koefisien AR:", paste(round(coef_ar, 4), collapse=", "), "\n")
      cat("  Akar-akar AR (|akar|):", paste(round(ar_mod, 4), collapse=", "), "\n")
      
      if (all(ar_mod > 1)) {
        cat("  ✓ STASIONER: Semua |akar| > 1 (akar di luar lingkaran unit)\n")
        ar_stationary <- TRUE
      } else {
        cat("  ✗ TIDAK STASIONER: Terdapat |akar| <= 1\n")
        ar_stationary <- FALSE
      }
    } else {
      cat("  ✗ Gagal menghitung akar AR\n")
      ar_stationary <- NA
    }
  } else {
    cat("  Model tanpa komponen AR → stasioner (differencing)\n")
    ar_stationary <- TRUE
  }
  
  # Cek invertibilitas (akar-akar polinomial MA)
  if (!is.null(coef_ma) && length(coef_ma) > 0) {
    # Polinomial MA: 1 + theta1*B + theta2*B^2 + ... + thetaq*B^q
    ma_poly <- c(1, coef_ma)
    ma_roots <- tryCatch(polyroot(ma_poly), error = function(e) NULL)
    
    if (!is.null(ma_roots)) {
      ma_mod <- Mod(ma_roots)
      cat("  Koefisien MA:", paste(round(coef_ma, 4), collapse=", "), "\n")
      cat("  Akar-akar MA (|akar|):", paste(round(ma_mod, 4), collapse=", "), "\n")
      
      if (all(ma_mod > 1)) {
        cat("  ✓ INVERTIBEL: Semua |akar| > 1 (akar di luar lingkaran unit)\n")
        ma_invertible <- TRUE
      } else {
        cat("  ✗ TIDAK INVERTIBEL: Terdapat |akar| <= 1\n")
        ma_invertible <- FALSE
      }
    } else {
      cat("  ✗ Gagal menghitung akar MA\n")
      ma_invertible <- NA
    }
  } else {
    cat("  Model tanpa komponen MA → invertibel\n")
    ma_invertible <- TRUE
  }
  
  return(list(
    Model = nama_model,
    AR_Stationary = ar_stationary,
    MA_Invertible = ma_invertible,
    AR_Roots = if (exists("ar_mod")) round(ar_mod, 4) else NULL,
    MA_Roots = if (exists("ma_mod")) round(ma_mod, 4) else NULL
  ))
}

# ============================================================
# FUNGSI PLOT AKAR-AKAR (UNIT ROOT PLOT)
# ============================================================
plot_unit_roots <- function(model, nama_model) {
  coef_ar <- grep("^ar", names(model$coef), value = TRUE)
  coef_ma <- grep("^ma", names(model$coef), value = TRUE)
  
  # Inisialisasi plot
  plot(NA, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), 
       xlab = "Real", ylab = "Imaginary",
       main = paste("Unit Root Plot -", nama_model),
       asp = 1)
  abline(h = 0, v = 0, col = "gray", lty = 2)
  symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, 
          fg = "black", lty = 2, lwd = 1.5)
  
  # Plot akar AR (merah)
  if (length(coef_ar) > 0) {
    ar_poly <- c(1, -model$coef[coef_ar])
    ar_roots <- polyroot(ar_poly)
    points(Re(ar_roots), Im(ar_roots), col = "red", pch = 19, cex = 1.2)
  }
  
  # Plot akar MA (biru)
  if (length(coef_ma) > 0) {
    ma_poly <- c(1, model$coef[coef_ma])
    ma_roots <- polyroot(ma_poly)
    points(Re(ma_roots), Im(ma_roots), col = "blue", pch = 17, cex = 1.2)
  }
  
  legend("bottomright", 
         legend = c("AR Roots (stasioneritas)", "MA Roots (invertibilitas)", "Unit Circle"),
         col = c("red", "blue", "black"),
         pch = c(19, 17, NA),
         lty = c(NA, NA, 2),
         cex = 0.8, bg = "white")
}

# ============================================================
# LOOP CEK STASIONERITAS & INVERTIBILITAS
# ============================================================
cat("\n============================================================\n")
cat("         CEK STASIONERITAS & INVERTIBILITAS\n")
cat("============================================================\n")

ringkasan_akar <- data.frame()

for (nama in names(models)) {
  model <- fitted_list[[nama]]
  hasil <- cek_stasioner_invertibel(model, nama)
  
  ringkasan_akar <- rbind(ringkasan_akar, data.frame(
    Model = nama,
    Stasioner = ifelse(hasil$AR_Stationary, "Ya", 
                       ifelse(is.na(hasil$AR_Stationary), "NA", "Tidak")),
    Invertibel = ifelse(hasil$MA_Invertible, "Ya",
                        ifelse(is.na(hasil$MA_Invertible), "NA", "Tidak"))
  ))
  
  # Plot unit roots untuk setiap model (optional, bisa dinonaktifkan)
  # plot_unit_roots(model, nama)
}

# ============================================================
# REKAP AKHIR SEMUA KRITERIA
# ============================================================
cat("\n============================================================\n")
cat("              REKAPITULASI LENGKAP\n")
cat("============================================================\n")

# Gabungkan ringkasan sebelumnya dengan ringkasan akar
ringkasan_lengkap <- merge(ringkasan, ringkasan_akar, by = "Model")
ringkasan_lengkap <- ringkasan_lengkap[order(ringkasan_lengkap$AIC), ]

cat(sprintf("%-15s %10s %10s %10s %12s %14s %10s %10s\n",
            "Model", "AIC", "BIC", "RMSE", "WN Lulus", "Param Signif", 
            "Stasioner", "Invertibel"))
cat(rep("-", 95), "\n", sep = "")
for (i in seq_len(nrow(ringkasan_lengkap))) {
  r <- ringkasan_lengkap[i, ]
  cat(sprintf("%-15s %10.4f %10.4f %10.4f %12s %14s %10s %10s\n",
              r$Model, r$AIC, r$BIC, r$RMSE, r$WN_lulus, r$Param_signif,
              r$Stasioner, r$Invertibel))
}
cat(rep("-", 95), "\n", sep = "")

# ============================================================
# TABEL KOMPARASI AKAR-akar
# ============================================================
cat("\n============================================================\n")
cat("         TABEL AKAR-AKAR KARAKTERISTIK\n")
cat("============================================================\n")

for (nama in names(models)) {
  model <- fitted_list[[nama]]
  
  coef_ar <- grep("^ar", names(model$coef), value = TRUE)
  coef_ma <- grep("^ma", names(model$coef), value = TRUE)
  
  cat(sprintf("\n--- %s ---\n", nama))
  
  if (length(coef_ar) > 0) {
    ar_poly <- c(1, -model$coef[coef_ar])
    ar_roots <- polyroot(ar_poly)
    ar_mod <- Mod(ar_roots)
    cat("  Akar AR invers (1/|akar|):", paste(round(1/ar_mod, 4), collapse=", "), "\n")
    if (all(ar_mod > 1)) {
      cat("  → Stasioner (semua |1/akar| < 1)\n")
    } else {
      cat("  → Tidak stasioner (terdapat |1/akar| >= 1)\n")
    }
  } else {
    cat("  Tidak ada akar AR\n")
  }
  
  if (length(coef_ma) > 0) {
    ma_poly <- c(1, model$coef[coef_ma])
    ma_roots <- polyroot(ma_poly)
    ma_mod <- Mod(ma_roots)
    cat("  Akar MA invers (1/|akar|):", paste(round(1/ma_mod, 4), collapse=", "), "\n")
    if (all(ma_mod > 1)) {
      cat("  → Invertibel (semua |1/akar| < 1)\n")
    } else {
      cat("  → Tidak invertibel (terdapat |1/akar| >= 1)\n")
    }
  } else {
    cat("  Tidak ada akar MA\n")
  }
}