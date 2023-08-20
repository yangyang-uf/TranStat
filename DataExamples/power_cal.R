library(lme4)
library(geepack)
n_sim <- 100
ORs_TRUE <- ORi_TRUE <- '0.5'
ORs <- as.numeric(ORs_TRUE)
ORi <- as.numeric(ORi_TRUE)

#setwd(paste0('/Users/yy70876/teach/TranStat/Examples/power/individual_randomization/', n_sim, '/ORs_', ORs_TRUE, '_ORi_', ORi_TRUE)) 
setwd(paste0('/Users/yy70876/teach/TranStat/Examples/power/cluster_randomization/', n_sim, '/ORs_', ORs_TRUE, '_ORi_', ORi_TRUE)) 
out <- read.table('output.txt')
col.names<-c('sn', 'inc', 'inf', 'b', 'b.se', 'b.cil', 'b.ciu', 'p', 'p.se', 'p.cil', 'p.ciu', 'or1', 'or1.se', 'or1.cil', 'or1.ciu', 'or2', 'or2.se', 'or2.cil', 'or2.ciu', 'or3', 'or3.se', 'or3.cil', 'or3.ciu')
col.names <- c(col.names, paste('var', 1:15, sep=''))
colnames(out) <- col.names
#c(median(out$b), median(out$p), median(out$or1), median(out$or3))
#c(sd(out$b), sd(out$p), sd(out$or1), sd(out$or3))
c((ORs-mean(out$or1))/ORs, mean((out$or1-ORs)^2), mean(out$or1.ciu < 1 | out$or1.cil > 1), mean(out$or1.ciu >= ORs & out$or1.cil <= ORs), 
  (ORi-mean(out$or3))/ORi, mean((out$or3-ORi)^2), mean(out$or3.ciu < 1 | out$or3.cil > 1), mean(out$or3.ciu >= ORi & out$or3.cil <= ORi))



store.gee <- store.lmm <- c()
for(i in 0:(n_sim-1))
{
   print(i)
   dat <- read.table(paste0('sim_pop_', i, '.txt'))
   colnames(dat) <- c('id', 'hh', 'idx', 'infection', 'day_ill', 'vacc')
   index_cases <- subset(dat, idx==1, select=c(hh, vacc))
   colnames(index_cases) <- c('hh', 'idx_vacc')
   contacts <- subset(dat, idx==0)
   contacts <- merge(contacts, index_cases, by='hh', all.x=TRUE)
      
   mf <- formula(infection ~ idx_vacc + vacc)
   gee <- geeglm(mf, data=contacts, id=hh, family=binomial(link = "logit"), corstr="exchangeable")
   if(summary(gee)$error == 0)
   {
      out<-as.data.frame(coef(summary(gee)))[,1:2]
      out$cil<-out[,1] - 1.96 * out[,2]
      out$ciu<-out[,1] + 1.96 * out[,2]
      store.gee <- rbind(store.gee, exp(cbind(out[2,-2], out[3,-2])))
   }
   lmm <- glmer( infection ~ idx_vacc + vacc + (1 | hh), data=contacts, family = binomial)
   if(summary(lmm)$optinfo$conv$op == 0)
   {
      out<-as.data.frame(summary(lmm)$coefficients)[,1:2]
      out$cil<-out[,1] - 1.96 * out[,2]
      out$ciu<-out[,1] + 1.96 * out[,2]
      store.lmm <- rbind(store.lmm, exp(cbind(out[2,-2], out[3,-2])))
   }  
}
save(store.gee, store.lmm, file='gee_lmm.RData')


load('gee_lmm.RData')
colnames(store.gee) <- c('ORi', 'ORi.cil', 'ORi.ciu', 'ORs', 'ORs.cil', 'ORs.ciu')
rownames(store.gee) <- NULL
c((ORs-mean(store.gee$ORs))/ORs, mean((store.gee$ORs-ORs)^2), mean(store.gee$ORs.ciu < 1 | store.gee$ORs.cil > 1), 
  mean(store.gee$ORs.ciu >= ORs & store.gee$ORs.cil <= ORs), 
  (ORi-mean(store.gee$ORi))/ORi, mean((store.gee$ORi-ORi)^2), mean(store.gee$ORi.ciu < 1 | store.gee$ORi.cil > 1),
  mean(store.gee$ORi.ciu >= ORi & store.gee$ORi.cil <= ORi))
colnames(store.lmm) <- c('ORi', 'ORi.cil', 'ORi.ciu', 'ORs', 'ORs.cil', 'ORs.ciu')
rownames(store.lmm) <- NULL
c((ORs-mean(store.lmm$ORs))/ORs, mean((store.lmm$ORs-ORs)^2), mean(store.lmm$ORs.ciu < 1 | store.lmm$ORs.cil > 1), 
  mean(store.lmm$ORs.ciu >= ORs & store.lmm$ORs.cil <= ORs), 
  (ORi-mean(store.lmm$ORi))/ORi, mean((store.lmm$ORi-ORi)^2), mean(store.lmm$ORi.ciu < 1 | store.lmm$ORi.cil > 1),
  mean(store.lmm$ORi.ciu >= ORi & store.lmm$ORi.cil <= ORi))

   
   