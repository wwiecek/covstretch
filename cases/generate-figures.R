ggsave(paste0(fig_folder, "/le_both.pdf"), le_both, width = 5.55, height = 6.475)
ggsave(paste0(fig_folder, "/delay_both.pdf"), delay_both, width = 5.55, height = 6.475)

ggsave(paste0(fig_folder, "/delay_switch.pdf"), delay_switch + theme(text = element_text(size=9)), width = 5.55, height = 3.7)
ggsave(paste0(fig_folder, "/g2.pdf"), fig_g2, width = 5.55, height=5.55)
ggsave(paste0(fig_folder, "/sgg_age.pdf"), sgg_age + theme(text = element_text(size=7)), width = 5.55, height=1.85)
ggsave(paste0(fig_folder, "/g1_joint.pdf"),g1_joint, width = 5.55, height=1.85)


# FDF figures
ggsave(paste0(fig_folder, "/fdf2.pdf"), fig_fdf2, width = 5.55, height=3.4)
ggsave(paste0(fig_folder, "/fdf1.pdf"),fig_fdf1 + theme(text = element_text(size=9)), width = 5.55, height=2.3)
ggsave(paste0(fig_folder, "/sfdf.pdf"),fig_sfdf, width = 5.55, height=7.4)
# ggsave(paste0(fig_folder, "/sfdf2.pdf"),fig_sfdf2, width = 5.55, height=7.4)

ggsave(paste0(fig_folder, "/fig_kappa.pdf"),fig_kappa + theme(text = element_text(size=8)), width = 5.55, height=3)

ggsave(paste0(fig_folder, "/pre-epidemic.pdf"),benefits_gg, width = 5.55, height=3.24)

ggsave(paste0(fig_folder, "/vax_rate.pdf"), fig.vax_rate, width = 5.55, height=3.4)
