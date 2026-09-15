# Concordance metrics schematic

source(here::here("code/00_config.R"))
theme_set(theme_classic(base_size = 10))

data <- read_rds(here(paste0("output/time_varying_estimates/SC2_run1_NDs-imputed_Rt.rds")))

draws_a <- extract_Rt_draws(get_samples(data$Boston_MA$clinical)) %>% mutate(source="a")
draws_b <- extract_Rt_draws(get_samples(data$Boston_MA$untargeted_ww)) %>% mutate(source="b")
draws_b_perm <- permute_draws(draws_b) %>% mutate(source="b_perm")
all_draws <- rbind(draws_a, draws_b, draws_b_perm)
all_draws_summary <- all_draws %>% 
  group_by(source,date) %>% 
  summarise(ymed = median(value),
             ymin = quantile(value, 0.25),
             ymax = quantile(value,0.75))

ggplot(all_draws_summary %>% filter(date >= as.Date("2024-10-01")), 
       aes(x = date,y = ymed,ymin = ymin, ymax = ymax)) +
  geom_hline(yintercept = 1,linetype="dashed",linewidth=0.5) +
  geom_line(aes(col = source, linetype = source), linewidth = 1) +
  # geom_ribbon(aes(fill = source), alpha = 0.25) +
  labs(x="",y="Time-varying estimate of transmission dynamics",col="Data source",fill="Data source") +
  scale_color_manual(values = c("#e41a1c","#377eb8","#377eb8")) +
  scale_fill_manual(values = c("#e41a1c","#377eb8","#377eb8")) +
  scale_linetype_manual(values = c("solid","solid","dashed")) +
  ylim(c(0.88,1.11)) +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right")

ggplot(all_draws_summary %>% filter(date >= as.Date("2024-10-01"), source != "b_perm"), 
       aes(x = date,y = ymed,ymin = ymin, ymax = ymax)) +
  geom_hline(yintercept = 1,linetype="dashed",linewidth=0.5) +
  geom_line(aes(col = source), linewidth = 1, alpha = 0.5) +
  # geom_ribbon(aes(fill = source), alpha = 0.25) +
  labs(x="",y="Time-varying estimate of transmission dynamics",col="Data source",fill="Data source") +
  scale_color_manual(values = c("#e41a1c","#377eb8","#377eb8")) +
  scale_fill_manual(values = c("#e41a1c","#377eb8","#377eb8")) +
  # scale_linetype_manual(values = c("solid","solid","dashed")) +
  ylim(c(0.88,1.11)) +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right")
ggsave(here("output/figures/a_vs_b.png"),width = 6, height = 4, units = "in", dpi = 400)

ggplot(all_draws_summary %>% filter(date >= as.Date("2024-10-01"), source != "b"), 
       aes(x = date,y = ymed,ymin = ymin, ymax = ymax)) +
  geom_hline(yintercept = 1,linetype="dashed",linewidth=0.5) +
  geom_line(aes(col = source), linewidth = 1, alpha = 0.5) +
  # geom_ribbon(aes(fill = source), alpha = 0.25) +
  labs(x="",y="Time-varying estimate of transmission dynamics",col="Data source",fill="Data source") +
  scale_color_manual(values = c("#e41a1c","#377eb8","#377eb8")) +
  scale_fill_manual(values = c("#e41a1c","#377eb8","#377eb8")) +
  # scale_linetype_manual(values = c("solid","solid","dashed")) +
  ylim(c(0.88,1.11)) +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right")
ggsave(here("output/figures/a_vs_b_perm.png"),width = 6, height = 4, units = "in", dpi = 400)
