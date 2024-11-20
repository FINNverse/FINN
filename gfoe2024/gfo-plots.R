library(ggplot2)
library(data.table)
# read rds file
runtimes <- data.table(readRDS("gfoe2024/results_runtime_21-8-24.RDS"))

runtimes[,method := factor(NGPUs, levels = c(0,1,4), labels = c("CPU", "1 GPU", "4 GPUs")),]

ggplot(runtimes, aes(x=sites, color=method)) +
  geom_line(aes(y=mean_time)) +
  labs(x="sites",
       y="runtime (s)") +
  ggthemes::theme_base() +
  theme(
    legend.position = c(0, 1),
    legend.justification = c(-.3, 1.1)
    )+
  scale_color_discrete(name = "Hardware Configuration")+
  coord_cartesian(ylim = c(0, NA), xlim = c(min(runtimes$sites), NA))+
  scale_x_continuous(breaks = c(min(runtimes$sites),seq(500, max(runtimes$sites), 500)))+
  scale_y_continuous(breaks = c(0, seq(0, max(runtimes$mean_time), 60)))


