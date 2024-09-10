#R3.5
library(ggplot2)


############## 最终版————无字版本
df <- read.csv('/home/wus/2023-4-25_012_gmk_Mammary-Gland/GSE108097/SRR6428483/fgseaRes-5-7-mac-m2.csv')


p<-ggplot(
     data=df, # 定义源数据
     mapping=aes( # 定义美学元素
         x=FDR, # x轴代表wt取值
         y=NES,
         colour=factor(color)
 # y轴代表mpg取值
     )
 )

p2 <- p+geom_point(size=9)+


  scale_color_manual(values = c('#Ec221f','#bec1c4'),
                     name = 'signature',
                     labels = df$name)+theme_classic()+

scale_y_continuous(limits = c(1, 2.6),breaks = c(0, 0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6 ))
p2

 p3 <- p2+ theme(
    legend.position =  'none',##删除图例
    axis.title.x =element_text(size=30), 
    axis.title.y=element_text(size=30),
    axis.text=element_text(size=24),
    axis.ticks=element_line(size=2), #刻度线粗细
  	axis.ticks.length=unit(0.5, "cm"), ##刻度线长度
  	axis.line = element_line(size=2) #边框粗细
    )

#ggsave(p2,filename = "mac1-notext.png",width = 12,height = 9)
ggsave(p3,filename = "mac_1-2.6_nolegend.png",width = 10,height = 9,dpi=600)

p4 <- p+geom_point(size=9)+


  scale_color_manual(values = c('#Ec221f','#bec1c4'),
                     name = 'signature',
                     labels = df$name)+theme_classic()+

scale_y_continuous(limits = c(0, 4.2),breaks = c(0,0.8,1.6,2.4,3.2,4.0,4.2 ))


 p5 <- p4+ theme(
    legend.position =  'none',##删除图例
    axis.title.x =element_text(size=30), 
    axis.title.y=element_text(size=30),
    axis.text=element_text(size=24), 
    #panel.border = element_rect(size=2),
    axis.ticks=element_line(size=2), #刻度线粗细
  	axis.ticks.length=unit(0.5, "cm"), ##刻度线长度
  	axis.line = element_line(size=2) #边框粗细
    )

#ggsave(p2,filename = "mac1-notext.png",width = 12,height = 9)
setwd('/home/wus/2023-4-25_012_gmk_Mammary-Gland/GSE108097/SRR6428483/')
ggsave(p5,filename = "mac_0-2.6_nolegend-size9.png",width = 10,height = 9,dpi=600)



########
