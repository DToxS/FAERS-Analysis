###calaculate davis selectivity score using kinobeads 
library(data.table)
library(magrittr)
library(ggplot2)
library(iheatmapr)
library(scales)
library(tidyr)

##read data
kd = fread("data/aan4368_Table_S3-kd.csv", sep = ",", header = T)
ki.smiles = fread("data/KI-smiles.2.csv")
ki = read.csv("data/KI-ecfp2.csv", sep = ",", header = F)
#kd = fread("data/41587_2011_BFnbt1990_MOESM5_ESM-davis.csv", sep = ",", header = T)
kd.names = names(kd)
coen.scores = fread("data/KI-risk-score-training.csv", sep = ",", header = T)
###figure 4d 
library(ggrepel)
coen.scores.testing = coen.scores[type=="testing"]
ggplot(coen.scores.testing,aes( observed, predicted, label= threelettercode)) + 
  geom_point(color = "blue", size = 3) + 
  geom_text_repel(nudge_x = .05, size = 10) + 
  expand_limits(x = c(0,1.5), y = c(0,1.5)) + 
  geom_abline( intercept= 0 , linetype = 2) + theme_bw() +
  theme(axis.text = element_text(size = 20))

plot(coen.scores[type == "testing"]$observed, coen.scores[type == "testing"]$predicted)
##functions
ggplotRegression <- function (fit) {
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
kd_selectivity = function(df, threshold = 100, startcol = 3 )
{
  return.df = data.frame()
  for ( i in startcol:ncol(df))
  {
    col = df[,..i]
    drug = names(col)
    selectivity = length(which(col < threshold)) / nrow(col)
    temp.df = cbind(drug, as.numeric(selectivity), as.numeric(threshold)) %>% as.data.frame
    names(temp.df) = c("Drugnames", "selectivity", "threshold")
    return.df = rbind(return.df, temp.df)
  }
  return(return.df)
}
##parse data 
kd = kd[Kinase == "yes"]
selected.drugs = which(kd.names %in% coen.scores$Drugnames)
kd.kinase = kd[,c(1,..selected.drugs)]
selectivities = kd_selectivity(kd.kinase, 500, 2)
select.m = merge(selectivities, coen.scores, by = "Drugnames")
select.m$selectivity = select.m$selectivity %>% as.character() %>% as.numeric()
select.m$rescaled = rescale(select.m$observed, c(-1,1))
select.m$rescaledpred = rescale(select.m$predicted, c(-1,1))
select.m$drug =NULL
### correlation 
ggplotRegression(glm(selectivity ~ 0 + observed , data=select.m))
### heatmap 
kd.kinase.df = kd.kinase %>%  as.data.frame()
row.names(kd.kinase.df) = kd.kinase.df$`Gene name`
kd.kinase.df$`Gene name` = NULL
####remove no information kinases rows 
noinfo = 2600000
rowsums = apply(kd.kinase.df , 1 , function (x) {sum(x)}) %>% as.numeric  
kd.kinase.df = kd.kinase.df[which(rowsums != noinfo),] 

### KI chemical similarity 
ki.wide = spread(ki, key =  V1, value = V3)
row.names(ki.wide) = ki.wide$V2 
ki.wide$V2 = NULL
clust = hclust(dist(data.matrix(ki.wide)))
select.m$drug = NULL
coen.scores$Drugnames = as.factor(as.character(coen.scores$Drugnames))
coen.scores.2 = coen.scores[,c(6,2,3,4)] %>% as.data.frame()
rescale(c(coen.scores.2$observed , coen.scores.2$predicted), c(-1,1))[1:26]

z = merge(selectivities, coen.scores.2 , by = "Drugnames")
z$observed_scaled = rescale(c(z$observed,z$predicted),c(-1,1))[1:26]
z$predicted_scaled = rescale(c(z$observed, z$predicted),c(-1,1))[27:52]
kd.kinase.df.2 = -1*log((t(data.matrix(kd.kinase.df))/100000))%>% as.data.frame

main_heatmap(data.matrix(kd.kinase.df.2) ,
               layout = list(font = list(size = 8))) %>%
  add_row_clusters(factor(select.m$PrimaryType), colors = 
                     c("f0134d", "ff6f5e", "f8b195", "c06c84", "40bfc1", "6c5b7b")) %>%
  add_row_dendro(clust) %>%
  add_col_clustering() %>%
  add_row_labels(font = list(size = 20)) %>% 
  add_row_barplot(y = z$Drugnames, z$selectivity, 
                  layout = 
                    list(title = "Kinase inhibitor<br>selectivity<br>S(500nM)")) %>%
  add_row_barplot(y = z$Drugnames, z$observed_scaled, color = "green",  
                  layout = 
                    list(title = "Observed<br>risk score")) %>%
  add_row_barplot(y = z$Drugnames, z$predicted_scaled, color = "blue",  
                  layout = 
                    list(title = "Predicted<br>risk score")) %>%
  add_row_barplot(y = z$Drugnames, z$abserror, color = "red", 
                  layout = 
                    list(title = "Error"))

ggplot(select.m[select.m$Drugnames != "Ibrutinib",], aes(type, abserror )) + geom_boxplot()
####elastic net regression to identify kinase targets of cardiotoxicity 
library(caret)
library(glmnet)

kd.kinase.df.2 = kd.kinase.df
kd.kinase.df.2[ kd.kinase.df.2 < 500 ] = 1 
kd.kinase.df.2[ kd.kinase.df.2 >= 500  ] = 0 

kd.kinase.df.3 = kd.kinase.df.2[which(apply(kd.kinase.df.2, 1 , function(x) {sum(x)}) > 3),]  

kd.kinase.scaled = t(data.matrix(kd.kinase.df.3)) %>% as.data.frame()
kd.kinase.scaled$Drugnames = row.names(kd.kinase.scaled)
kd.kinase.scaled.m = merge(kd.kinase.scaled, z, by = "Drugnames")
model = glm(formula = observed_scaled ~ . ,data =  kd.kinase.scaled.m[1:25,c(2:19, 25)] )
summary(model)
apply(kd.kinase.scaled.m[,c(2:19)], 2, function(x) {sum(x)})
ggplotRegression(ll)
ll = lm(selectivity ~ 0 + observed , data=select.m)
ggplot(select.m, aes( select.m$observed,
                      select.m$selectivity, 
                      label = select.m$Drugnames)) + 
  geom_point() + 
  geom_abline(intercept=0, slope=ll$coefficients[1], color='#2C3E50', size=1.1) + 
  coord_cartesian(xlim = c(0,2)) + 
  geom_text_repel(size = 10)+ 
  theme_bw() + 
  theme(axis.text = element_text(size = 20))
lm_fit <- lm()

heatmap(log(data.matrix(t(kd.kinase[,3:ncol(kd.kinase)]))))

### Structure active similarity map

ki = read.csv("data/KI-tanimoto.csv", sep = ",", header = T)
ki = merge(ki, coen.scores.2[,1:2], by.x = c("D1"), by.y = c("Drugnames"))
ki = merge(ki, coen.scores.2[,1:2], by.x = c("D2"), by.y = c("Drugnames"))
ki$absdiff = abs(ki$observed.x - ki$observed.y)
ki = ki[ki$ECFP4 != 1,] 

remove_duplicates= function(long_data, key1, key2 ) { 
    seen = c()
    df.final = data.frame()
    for ( i in 1:nrow(long_data) )
    {
      seen = c(seen, as.character(long_data[i,key1]))
      test = long_data[i,key2] %>% as.character()
      if ( test %in% seen )
      {
        next 
      }
      else 
      {
        df.final = rbind(df.final, long_data[i,])
      }
    }
    return(df.final)
}
ki.single  = remove_duplicates(ki, 'D2', 'D1')
quantile(unique(ki.single$Weighted),probs = c(0,.1,.25,.5,.75,1))
max(unique(ki.single$absdiff)) / 2 
library(ggrepel)
xlimits = c(.36, NA)
ggplot(ki.single, aes( 
    x = Weighted , 
    y = absdiff, 
    label = paste(D1, ",",  D2, sep = "")
  )) + 
  geom_point() + 
  geom_hline(yintercept = 0.8198644, color = 'red', size = .5) + 
  geom_vline(xintercept =  0.35, color = 'red', size = .5) + 
  geom_text_repel(data = ki.single[ki.single$Weighted >  0.35, ],nudge_x = .01, size = 5, 
                  xlim = xlimits) + 
  theme_bw() +
  theme(axis.text = element_text(size = 15))



