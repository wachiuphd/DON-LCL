library("ggpubr")
TDVF <- ggarrange(TDVF05, TDVF01, 
                  labels = c("A", "B"),
                  ncol = 2, nrow = 1)
ggsave("TDVF_figure.pdf", plot=TDVF, scale=0.5, width=20, height=8)
