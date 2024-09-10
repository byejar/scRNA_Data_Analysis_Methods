


###macs 8-27
library(dplyr)
df <- read.csv('/home/wus/for_seurat/8_2/fig-8-26-macs.csv', sep=",", header=T)

df %>%
  group_by(group) %>%
  mutate(percent = round(Num/sum(Num),4)) ->df2
df2
ggplot(df2,aes(x=group,y=percent,fill=type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=percent),position=position_stack(vjust=0.5))

df2$output<-paste(as.character(df2$percent),'(',as.character(df2$Num),')')
df2
ggplot(df2,aes(x=group,y=percent,fill=type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=output),position=position_stack(vjust=0.5))

### tcell 8-27
library(dplyr)
rm(list=ls())
df <- read.csv('/home/wus/for_seurat/8_2/fig-8-26-tcell.csv', sep=",", header=T)

df %>%
  group_by(group) %>%
  mutate(percent = round(Num/sum(Num),4)) ->df2
df2
ggplot(df2,aes(x=group,y=percent,fill=type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=percent),position=position_stack(vjust=0.5))

df2$output<-paste(as.character(df2$percent),'(',as.character(df2$Num),')')
df2
ggplot(df2,aes(x=group,y=percent,fill=type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=output),position=position_stack(vjust=0.5))



### neu
library(dplyr)
rm(list=ls())
df <- read.csv('/home/wus/for_seurat/8_2/fig-8-26-neu.csv', sep=",", header=T)

df %>%
  group_by(group) %>%
  mutate(percent = round(Num/sum(Num),4)) ->df2
df2
ggplot(df2,aes(x=group,y=percent,fill=type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=percent),position=position_stack(vjust=0.5))

df2$output<-paste(as.character(df2$percent),'(',as.character(df2$Num),')')
df2
ggplot(df2,aes(x=group,y=percent,fill=type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=output),position=position_stack(vjust=0.5))


  ### immune
library(dplyr)
rm(list=ls())
df <- read.csv('/home/wus/for_seurat/8_2/fig-8-26-immune.csv', sep=",", header=T)

df %>%
  group_by(group) %>%
  mutate(percent = round(Num/sum(Num),4)) ->df2
df2
ggplot(df2,aes(x=group,y=percent,fill=type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=percent),position=position_stack(vjust=0.5))

df2$output<-paste(as.character(df2$percent),'(',as.character(df2$Num),')')
df2
ggplot(df2,aes(x=group,y=percent,fill=type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=output),position=position_stack(vjust=0.5))

