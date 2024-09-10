# scRNA
##################### 在win上处理之后传回服务器，在服务器转换为长数据
library(tidyverse)
library(ggplot2)
library(ggpubr)

a<-read.csv('/home/wus/016-2023-6-15_PTEN-mutation-TCGA-CIBERSORT/file/GBM_box-7-17.csv')
head(a)
#a[, ] <- sapply(a[, ], as.character)


library(tidyr)
chji_long <- gather(a, key = "CELL", value = "SCORE",
                    -'SAMPLE', -'TYPE')

chji_long$CELL<- as.factor(chji_long$CELL)

bxp <- ggboxplot(
    chji_long, x="CELL", y="SCORE", color = "TYPE",bxp.errorbar =TRUE,outlier.shape = NA)

bxp<-bxp+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

bxp
#### test1 
# bxp<-bxp+stat_compare_means(aes(group=CELL),label = "p.format",method = "t.test",angle=45,label.y.npc='center')+ 
#  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
# bxp
#### test2
 bxp <- ggboxplot(
    chji_long, x="CELL", y="SCORE", color = "TYPE",bxp.errorbar =TRUE,outlier.shape = NA)+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

bxp + stat_compare_means(aes(group = TYPE),hide.ns=T,position='jitter',label = "p.format")



# #### 加p值
# library(rstatix)
# stat.test <- chji_long %>% 
#   group_by(CELL) %>%
#   t_test(SCORE ~ TYPE) %>%
#   adjust_pvalue(method = "bonferroni") %>% 
#   add_significance("p.adj")


# stat.test <- stat.test %>% 
#   add_xy_position(x='CELL',dodge = 1)
#   #add_y_position() 
# as.data.frame(stat.test)

# bxp + stat_pvalue_manual(
#   stat.test,  label = "p", tip.length = 0
#   ) + 
# scale_y_continuous(expand = expansion(mult = c(0, .5)))+ 
#  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))


##### 多组 
chji_long$TYPE <- factor(chji_long$TYPE, levels=c('MUT', "HOMDEL ", "no alteration"))
comparisons <- list(c("HOMDEL ", "no alteration"), c('MUT', 'no alteration'))

stat.test <- chji_long %>% 
  group_by(CELL) %>%
  t_test(SCORE ~ TYPE, comparisons = comparisons) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p")

stat.test <- stat.test %>%
  add_xy_position(x = "CELL", dodge = 0.8)
print(stat.test, width=Inf)





#############
library(ggpubr)
p <- ggboxplot(chji_long, x="CELL", y="SCORE", color = "TYPE",
    palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
    trim = TRUE) + 
    labs(x=NULL, fill=NULL)
p8 <- p + stat_pvalue_manual(stat.test, label = "{p} {p.signif}", tip.length = 0, step.increase = 0)+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
p8

########################## UCEC

a<-read.csv('/home/wus/016-2023-6-15_PTEN-mutation-TCGA-CIBERSORT/file/UCEC_box-7-17.csv')
head(a)
#a[, ] <- sapply(a[, ], as.character)


library(tidyr)
chji_long <- gather(a, key = "CELL", value = "SCORE",
                    -'SAMPLE', -'TYPE')

chji_long$CELL<- as.factor(chji_long$CELL)
 bxp <- ggboxplot(
    chji_long, x="CELL", y="SCORE", color = "TYPE",bxp.errorbar =TRUE,outlier.shape = NA)

bxp<-bxp+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

bxp
#### test1 
# bxp<-bxp+stat_compare_means(aes(group=CELL),label = "p.format",method = "t.test",angle=45,label.y.npc='center')+ 
#  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
# bxp
#### test2
 # bxp <- ggboxplot(
 #    chji_long, x="CELL", y="SCORE", color = "TYPE",bxp.errorbar =TRUE,outlier.shape = NA)+ 
 # theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

# bxp + stat_compare_means(aes(group = TYPE),hide.ns=T,position='jitter',label = "p.format")



# #### 加p值
# library(rstatix)
# stat.test <- chji_long %>% 
#   group_by(CELL) %>%
#   t_test(SCORE ~ TYPE) %>%
#   adjust_pvalue(method = "bonferroni") %>% 
#   add_significance("p.adj")


# stat.test <- stat.test %>% 
#   add_xy_position(x='CELL',dodge = 1)
#   #add_y_position() 
# as.data.frame(stat.test)

# bxp + stat_pvalue_manual(
#   stat.test,  label = "p", tip.length = 0
#   ) + 
# scale_y_continuous(expand = expansion(mult = c(0, .5)))+ 
#  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))



##### 多组 
chji_long$TYPE <- factor(chji_long$TYPE, levels=c('MUT', "HOMDEL ", "no alteration"))
comparisons <- list(c("HOMDEL ", "no alteration"), c('MUT', 'no alteration'))

stat.test <- chji_long %>% 
  group_by(CELL) %>%
  t_test(SCORE ~ TYPE, comparisons = comparisons) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p")

stat.test <- stat.test %>%
  add_xy_position(x = "CELL", dodge = 0.8)
print(stat.test, width=Inf)





#############
library(ggpubr)
p <- ggboxplot(chji_long, x="CELL", y="SCORE", color = "TYPE",
    palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
    trim = TRUE) + 
    labs(x=NULL, fill=NULL)
p8 <- p + stat_pvalue_manual(stat.test, label = "{p} {p.signif}", tip.length = 0, step.increase = 0)+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
p8