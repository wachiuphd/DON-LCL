library("ggpubr")
LCL <- ggarrange(pfit, pec10,
                 labels = c("A", "B"),
                 ncol = 2, nrow = 1)
ggsave("LCL_figure.pdf", plot=LCL, scale=0.65, width=24, height=15)

LCLposter <- ggarrange(pfit, pec10, TDVF05,
                 labels = c("A", "B", "C"),
                 ncol = 3, nrow = 1)
ggsave("LCL_figure for poster.pdf", plot=LCLposter, scale=0.7, width=36, height=12)
