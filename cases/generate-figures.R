width <- 6.5

ggsave(paste0(fig_folder, "/le_both.pdf"), le_both, width = width, height = 6.475/5.55*width)
ggsave(paste0(fig_folder, "/le_optimal.pdf"), le2.plot, width = width, height=3.4/5.55*width)
ggsave(paste0(fig_folder, "/delay_both.pdf"), delay_both, width = width, height = 6.475/5.55*width)

ggsave(paste0(fig_folder, "/delay_switch.pdf"), delay_switch + theme(text = element_text(size=9)), width = width, height = 3.7/5.55*width)
ggsave(paste0(fig_folder, "/g2.pdf"), fig_g2, width = width, height=width)
ggsave(paste0(fig_folder, "/g2_reductions_only.pdf"), g2b, width = width, height=0.6*width)
ggsave(paste0(fig_folder, "/sgg_age.pdf"), sgg_age + theme(text = element_text(size=9)), width = width, height=1.85/5.55*width)
ggsave(paste0(fig_folder, "/g1_joint.pdf"),g1_joint, width = width, height=1.85/5.55*width)


# FDF figures
ggsave(paste0(fig_folder, "/fdf2.pdf"), fig_fdf2, width = width, height=3.4/5.55*width)
ggsave(paste0(fig_folder, "/fdf1.pdf"),fig_fdf1 + theme(text = element_text(size=9)), width = width, height=2.3/5.55*width)
ggsave(paste0(fig_folder, "/sfdf.pdf"),fig_sfdf, width = width, height=7.4/5.55*width)
# ggsave(paste0(fig_folder, "/sfdf2.pdf"),fig_sfdf2, width = 5.55, height=7.4)

ggsave(paste0(fig_folder, "/fig_kappa.pdf"),fig_kappa + theme(text = element_text(size=8)), width = width, height=3/5.55*width)

ggsave(paste0(fig_folder, "/pre-epidemic.pdf"),benefits_gg, width = width, height=3.24/5.55*width)

ggsave(paste0(fig_folder, "/vax_rate.pdf"), fig.vax_rate, width = width, height=3.4/5.55*width)
ggsave(paste0(fig_folder, "/delay_impact.pdf"), fig_delay, width = width, height=3.4/5.55*width)
ggsave(paste0(fig_folder, "/supply_impact.pdf"), fig_extra_supply, width = width, height=3.4/5.55*width)

if (all_k)
  ggsave("figures/fdf_best_policy_allk.pdf", fig2.all_k, width = width, height=0.6*width)
