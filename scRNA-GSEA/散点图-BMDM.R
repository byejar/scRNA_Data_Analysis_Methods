#R3.5
library(ggplot2)
df <- read.csv('/home/wus/12-15_004_zc_TPL2022112807_RNA/TPL2022111843/analysis/8_GSEA/2-14/cp10000-v3.csv')
p<-ggplot(
     data=df, # 定义源数据
     mapping=aes( # 定义美学元素
         x=FDR, # x轴代表wt取值
         y=NES,colour=factor(color)
 # y轴代表mpg取值
     )
 )
#p+geom_point(aes(colour=factor(colour),shape=factor(shape))) # 以散点图形式呈现
p+geom_point()+stat_smooth(method = 'loess', span = 0.3, se = FALSE)
p+geom_point()+geom_line(size=2)
#等价形式

x=df$FDR
y=df$NES
plot(x,y)
sp<-spline(x,y)
lines(sp)



############
df <- read.csv('/home/wus/12-15_004_zc_TPL2022112807_RNA/TPL2022111843/analysis/8_GSEA/2-14/cp10000-v3.csv')
p<-ggplot(
     
 )
#p+geom_point(aes(colour=factor(colour),shape=factor(shape))) # 以散点图形式呈现
p+geom_point(data=df[1:13,],aes( x=df$FDR[1:13],y=df$NES[1:13]))+stat_smooth(method = 'lm', span = 0.3, se = FALSE)+
geom_point(data=df[14:252,],aes( x=df$FDR[14:252],y=df$NES[14:252]))+stat_smooth()


##########
df <- read.csv('/home/wus/12-15_004_zc_TPL2022112807_RNA/TPL2022111843/analysis/8_GSEA/2-14/cp10000-v3.csv')
df2 <- read.csv('/home/wus/12-15_004_zc_TPL2022112807_RNA/TPL2022111843/analysis/8_GSEA/2-14/cp10000-v4-mac.csv')

p<-ggplot(
     data=df, # 定义源数据
     mapping=aes( # 定义美学元素
         x=FDR, # x轴代表wt取值
         y=NES
 # y轴代表mpg取值
     )
 )

p2 <- p+geom_point(size=0.1)+geom_line(size=2)+geom_point(data=df2,aes(x=FDR,y=NES,colour=color),size=10)+
##标注排序修改breaks,breaks——labels,
 scale_colour_discrete(values=c('#Ec221f','#3582c8','#41ba42','#835988','#Ff7420','#Ac5021','#Fa78bd','#979696','#65bea0','#Fd8757','#8f99c4','#31475e','#f0a500'),
    #name="Experimental",
    labels = df2$NAME)+theme_classic()

 p3 <- p2+ theme(
    legend.position = c(.95, .95),#plot内位置
    legend.justification = c("right", "top"),#固定右上角
    legend.background = element_blank(),#图例背景色
    legend.key = element_blank(),#图标背景色element_rect(colour="black")
    legend.box.background = element_rect(fill=NA,color = "black",linetype = 1), #图例外框和背景色默认填充白色（删除fill=NA）
    #legend.box.just = "right",
    #legend.margin = margin(6, 6, 6, 6)#边框大小调整
    legend.text=element_text(size=20),#图例字体大小
    axis.title.x =element_text(size=30), axis.title.y=element_text(size=30),
    axis.text=element_text(size=24)
    )

ggsave(p2,filename = "mac1.png",width = 12,height = 9)
ggsave(p3,filename = "mac2.png",width = 12,height = 9,dpi=600)




########## 再测试
df <- read.csv('/home/wus/12-15_004_zc_TPL2022112807_RNA/TPL2022111843/analysis/8_GSEA/2-14/cp10000-v3.csv')
df2 <- read.csv('/home/wus/12-15_004_zc_TPL2022112807_RNA/TPL2022111843/analysis/8_GSEA/2-14/cp10000-v4-mac.csv')

p<-ggplot(
     data=df, # 定义源数据
     mapping=aes( # 定义美学元素
         x=FDR, # x轴代表wt取值
         y=NES
 # y轴代表mpg取值
     )
 )

p2 <- p+geom_point(size=0.1)+geom_line(size=1.5)+geom_point(data=df2,aes(x=FDR,y=NES,colour=factor(color)),size=10)+


  scale_color_manual(values = c('#Ec221f','#3582c8','#41ba42','#835988','#Ff7420','#Ac5021','#Fa78bd','#979696','#65bea0','#Fd8757','#8f99c4','#31475e','#f0a500'),
                     name = 'macrophage signature',
                     labels = df2$NAME)+theme_classic()


 p3 <- p2+ theme(
    legend.position = c(.95, .95),#plot内位置
    legend.justification = c("right", "top"),#固定右上角
    legend.background = element_blank(),#图例背景色
    legend.key = element_blank(),#图标背景色element_rect(colour="black")
    legend.box.background = element_rect(fill=NA,color = "black",linetype = 1), #图例外框和背景色默认填充白色（删除fill=NA）
    #legend.box.just = "right",
    #legend.margin = margin(6, 6, 6, 6)#边框大小调整
    legend.text=element_text(size=20),#图例字体大小
    axis.title.x =element_text(size=30), axis.title.y=element_text(size=30),
    axis.text=element_text(size=24)
    )

ggsave(p2,filename = "mac1.png",width = 12,height = 9)
ggsave(p3,filename = "mac2-v2.png",width = 12,height = 9,dpi=600)

############ 最终版————有字版本

df <- read.csv('/home/wus/12-15_004_zc_TPL2022112807_RNA/TPL2022111843/analysis/8_GSEA/2-14/2-15-cp.csv')


p<-ggplot(
     data=df, # 定义源数据
     mapping=aes( # 定义美学元素
         x=FDR, # x轴代表wt取值
         y=NES,
         colour=factor(color)
 # y轴代表mpg取值
     )
 )

p2 <- p+geom_point(size=4)+


  scale_color_manual(values = c('#Ec221f','#3582c8','#41ba42','#bec1c4'),
                     name = 'signature',
                     labels = df$name)+theme_classic()+

scale_y_continuous(limits = c(1, 2.6),breaks = c(0, 0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6 ))
p2

 p3 <- p2+ theme(
    legend.position = c(.95, .95),#plot内位置
    legend.justification = c("right", "top"),#固定右上角
    legend.background = element_blank(),#图例背景色
    legend.key = element_blank(),#图标背景色element_rect(colour="black")
    legend.box.background = element_rect(fill=NA,color = "black",linetype = 1), #图例外框和背景色默认填充白色（删除fill=NA）
    #legend.box.just = "right",
    #legend.margin = margin(6, 6, 6, 6)#边框大小调整
    legend.text=element_text(size=20),#图例字体大小
    axis.title.x =element_text(size=30), axis.title.y=element_text(size=30),
    axis.text=element_text(size=24)
    )

ggsave(p2,filename = "mac1.png",width = 12,height = 9)
ggsave(p3,filename = "mac2-v3.png",width = 10,height = 9,dpi=600)

############## 最终版————无字版本
df <- read.csv('/home/wus/12-15_004_zc_TPL2022112807_RNA/TPL2022111843/analysis/8_GSEA/2-14/2-15-cp.csv')


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


  scale_color_manual(values = c('#Ec221f','#bec1c4','#bec1c4','#bec1c4'),
                     name = 'signature',
                     labels = df$name)+theme_classic()+

scale_y_continuous(limits = c(1, 2.6),breaks = c(0, 0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6 ))


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
ggsave(p5,filename = "BMDM_0-2.6_nolegend-5-10-size9.png",width = 10,height = 9,dpi=600)



########
