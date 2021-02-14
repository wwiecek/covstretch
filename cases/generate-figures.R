ggsave(paste0(fig_folder, "/le_both.pdf"), le_both, width = 6, height = 7)
ggsave(paste0(fig_folder, "/delay_both.pdf"), delay_both, width = 6, height = 7)

ggsave(paste0(fig_folder, "/delay_switch.pdf"), delay_switch, width = 6, height = 4)
ggsave(paste0(fig_folder, "/g2.pdf"), fig_g2, width = 6, height=6)
ggsave(paste0(fig_folder, "/sgg_age.pdf"), sgg_age, width = 7.5, height=2.5)
ggsave(paste0(fig_folder, "/g1_joint.pdf"),g1_joint, width = 7.5, height=2.5)
ggsave(paste0(fig_folder, "/fdf2.pdf"), fig_fdf2 , width = 6.5, height=4)
ggsave(paste0(fig_folder, "/fdf1.pdf"),fig_fdf1, width = 6, height=2.5)

ggsave(paste0(fig_folder, "/fig_kappa.pdf"),fig_kappa, width = 7.5, height=4)

ggsave(paste0(fig_folder, "/pre-epidemic.pdf"),benefits_gg, width = 6, height=3.5)
