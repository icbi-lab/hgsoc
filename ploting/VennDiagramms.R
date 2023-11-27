install.packages("ggVennDiagram")
library(ggvenn) 

up = read.csv('Up_untreated.csv')
down = read.csv('Down_untreated.csv')
up_ola = read.csv('VennUp.csv')
down_ola = read.csv('Venn_down.csv')

a <- list(`OVC3 vs. UWB` = up$OVC3_UWB,
          `UWB iso vs. UWB` = up$UWB_UWBiso)
ggvenn(a,show_percentage = FALSE,  fill_color = c('red','orangered'),text_size = 10) 
ggsave("VennDiagramm_UP_CellLines.svg")

b <- list(`OVC3 vs. UWB` = down$OVC3_UWB,
          `UWB iso vs. UWB` = down$UWB_UWBiso)
ggvenn(b,show_percentage = FALSE,  fill_color = c('blue','lightblue'),text_size = 10) 
ggsave("VennDiagramm_Down_CellLines.svg")

c <- list(`OVC3 olaparib vs. DMSO` = up_ola$OVC3_ola,
          `UWB olaparib vs. DMSO` = up_ola$UWB_ola)
ggvenn(c,show_percentage = FALSE,  fill_color = c('red','orangered'),text_size = 10) 
ggsave("VennDiagramm_UP_Treatement.svg")

d <- list(`OVC3 olaparib vs. DMSO` = down_ola$OVC_ola,
      `UWB olaparib vs. DMSO` = down_ola$UWB_ola)
ggvenn(d,show_percentage = FALSE,  fill_color = c('blue','lightblue'),text_size = 10) 
ggsave("VennDiagramm_Down_Treatement.svg")
