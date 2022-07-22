library(readxl)
library(ggplot2)
library(scales)
theme_set(theme_bw())

data <- read_xlsx("Flow_Data.xlsx")

figure <- ggplot(data, aes(x = as.factor(Time), y = Flow, color = as.factor(Plot))) +
  geom_line(aes(group = Plot)) +
  geom_point(aes(shape = Treatment), size = 2.5) +
  xlab("Time after runoff (min)") +
  ylab(bquote("Flow Rate"~(m^3~s^-1))) +
  theme(legend.title = element_text(family = "Times New Roman", size = 18),
        legend.text = element_text(family = "Times New Roman", size = 16),
        axis.text.x = element_text(family = "Times New Roman", size = 14, color = "black"),
        axis.text.y = element_text(family = "Times New Roman", size = 14, color = "black"),
        axis.title.x = element_text(family = "Times New Roman", size = 18),
        axis.title.y = element_text(family = "Times New Roman", size = 18),
        panel.background = element_rect(colour = "black", linetype = "solid", size = 0.75)) +
  scale_color_manual(name = "Plot", values=c('#9ECAE1', '#2171B5', '#08306B',
                              '#A1D99B', '#238B45', '#00441B',
                              '#FC9272', '#CB181D', '#67000D')) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )

figure

ggsave("figure_S1.tiff", plot = figure,
       width = 10, height = 7, units = c("in"))
