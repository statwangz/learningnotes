library(ggplot2)
library(latex2exp)

df_11_10 <- data.frame(time = c(res_als_11[[1]]$time,
                                res_si_11[[1]]$time,
                                res_si_als_11[[1]]$time),
                       rel_obj = c(res_als_11[[1]]$`relative objevtive`,
                                   res_si_11[[1]]$`relative objevtive`,
                                   res_si_als_11[[1]]$`relative objevtive`),
                       method = factor(c(rep("ALS", res_als_11[[1]]$iteration),
                                         rep("softImpute", res_si_11[[1]]$iteration),
                                         rep("softImpute-ALS", res_si_als_11[[1]]$iteration))))
plot_11_10 <- ggplot(df_11_10, aes(x = time, y = rel_obj, color = method)) +
  geom_point() +
  ylab("Relative Objective (log scale)") +
  xlab("Time in Seconds") +
  ggtitle(TeX("(1000, 1000), 80% NAs, true rank = 30, r = 50, $\\lambda$ = 10, SNR = 1, iterations: ALS = 14, softImpute = 98, softImpute-ALS = 54")) +
  scale_y_log10() +
  theme_set(theme_bw()) +
  theme(legend.position = c(.92, .9),
        legend.title = element_blank())

df_11_30 <- data.frame(time = c(res_als_11[[2]]$time,
                                res_si_11[[2]]$time,
                                res_si_als_11[[2]]$time),
                       rel_obj = c(res_als_11[[2]]$`relative objevtive`,
                                   res_si_11[[2]]$`relative objevtive`,
                                   res_si_als_11[[2]]$`relative objevtive`),
                       method = factor(c(rep("ALS", res_als_11[[2]]$iteration),
                                         rep("softImpute", res_si_11[[2]]$iteration),
                                         rep("softImpute-ALS", res_si_als_11[[2]]$iteration))))
plot_11_30 <- ggplot(df_11_30, aes(x = time, y = rel_obj, color = method)) +
  geom_point() +
  ylab("Relative Objective (log scale)") +
  xlab("Time in Seconds") +
  ggtitle(TeX("(1000, 1000), 80% NAs, true rank = 30, r = 50, $\\lambda$ = 30, SNR = 1, iterations: ALS = 6, softImpute = 45, softImpute-ALS = 33")) +
  scale_y_log10() +
  theme_set(theme_bw()) +
  theme(legend.position = c(.92, .9),
        legend.title = element_blank())

df_11_50 <- data.frame(time = c(res_als_11[[3]]$time,
                                res_si_11[[3]]$time,
                                res_si_als_11[[3]]$time),
                       rel_obj = c(res_als_11[[3]]$`relative objevtive`,
                                   res_si_11[[3]]$`relative objevtive`,
                                   res_si_als_11[[3]]$`relative objevtive`),
                       method = factor(c(rep("ALS", res_als_11[[3]]$iteration),
                                         rep("softImpute", res_si_11[[3]]$iteration),
                                         rep("softImpute-ALS", res_si_als_11[[3]]$iteration))))
plot_11_50 <- ggplot(df_11_50, aes(x = time, y = rel_obj, color = method)) +
  geom_point() +
  ylab("Relative Objective (log scale)") +
  xlab("Time in Seconds") +
  ggtitle(TeX("(1000, 1000), 80% NAs, true rank = 30, r = 50, $\\lambda$ = 50, SNR = 1, iterations: ALS = 6, softImpute = 27, softImpute-ALS = 24")) +
  scale_y_log10() +
  theme_set(theme_bw()) +
  theme(legend.position = c(.92, .9),
        legend.title = element_blank())

plot_11_10
plot_11_30
plot_11_50