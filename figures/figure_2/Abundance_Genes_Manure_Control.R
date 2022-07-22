library(readxl)
library(tidyverse)

processed_data <- read_excel("Armstrong_qPCR_processed.xlsx")

manure_control <- processed_data %>%
  filter(Treatment == "Manure" | Treatment == "Strip + No Manure")

ARG_bottom <- manure_control %>% filter(Primer %in% c("tnpA-02", "IS6100", "sul1", "intl1(clinical)", "IS1247")) %>%
  pivot_longer(cols = Primer)

levels1 = ordered(c("tet(M)", "tet(Q)", "tet(44)", "tnpA-06", "tet(T)", "tet(W)", "ermF", "tet(36)", "aphA3",
                    "ermQ", "aadE", "ermB", "tet(O)", "sat4", "erm(35)", "aadD", "tet(L)", "tet(X)",
                    "sul2", "tnpA-05", "intl2", "ermT", "tnpA-02", "IS6100", "sul1",
                    "intl1(clinical)", "IS1247"))

levels2 = ordered(c("tnpA-02", "IS6100", "sul1", "intl1(clinical)", "IS1247"))

levels3 = ordered(c("Manure", "Water", "Sediment"))

all_ARG <- ggplot(manure_control, aes(ordered(Primer, levels = levels1), Relative_Abundance,
                                      color = ordered(Matrix, levels = levels3))) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  labs(y = "Relative Abundance") +
  scale_color_manual(labels = c("Manure (n = 1)", "Control Runoff Water (n = 18)", "Control Runoff Sediment (n = 18)"), values = c("#999999", "#56B4E9", "#E69F00")) +
  scale_x_discrete(labels = c("tet(M)" = "tet(M)", "tet(Q)" = "tet(Q)", "tet(44)" = "tet(44)", "tnpA-06" = "tnpA-06", "tet(T)" = "tet(T)",
                              "tet(W)" = "tet(W)", "ermF" = "ermF", "tet(36)" = "tet(36)", "aphA3" = "aphA3", "ermQ" = "ermQ", "aadE" = "aadE",
                              "ermB" = "ermB", "tet(O)" = "tet(O)", "sat4" = "sat4", "erm(35)" = "erm(35)", "aadD" = "aadD", "tet(L)" = "tet(L)",
                              "tet(X)" = "tet(X)", "sul2" = "sul2", "tnpA-05" = "tnpA-05", "intl2" = "intl2", "ermT" = "ermT", "tnpA-02" = "*tnpA-02",
                              "IS6100" = "*IS6100", "sul1" = "sul1", "intl1(clinical)" = "*intl1(clinical)", "IS1247" = "*IS1247")) +
  theme_light() +
  theme(axis.text.x = element_text(family = "Times New Roman", face = "italic", size = 16, color = "black",
                                   angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(family = "Times New Roman", size = 16, color = "black"),
        axis.title.y = element_text(family = "Times New Roman", face = "bold", size = 18),
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman", size = 18, color = "black"),
        legend.position = "top",
        panel.background = element_rect(colour = "black", linetype = "solid", size = 0.75))

all_ARG

bottom_ARG <- ggplot(ARG_bottom, aes(ordered(value, levels = levels2), Relative_Abundance,
                                     color = ordered(Matrix, levels = levels3))) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  scale_color_manual(values = c("#999999", "#56B4E9", "#E69F00")) +
  scale_x_discrete(labels = c("tnpA-02" = "*tnpA-02", "IS6100" = "*IS6100", "sul1" = "sul1", "intl1(clinical)" = "*intl1(clinical)", "IS1247" = "*IS1247")) +
  theme_light() +
  theme(axis.text.x = element_text(family = "Times New Roman", face = "italic", size = 14, color = "black",
                                   angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(family = "Times New Roman", size = 14, color = "black"),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.background = element_rect(colour = "black", linetype = "solid", size = 0.75),
        plot.background = element_rect(colour = "black", size = 0.75))

bottom_ARG

combined_plot <- all_ARG + annotation_custom(ggplotGrob(bottom_ARG), xmin = "sat4", xmax = "IS1247", ymin = 0.035, ymax = 0.22)

combined_plot

ggsave("figure_2.tiff", plot = combined_plot,
       width = 12.5, height = 7.5, dpi = 300, units = c("in"))
