
source("0-pkgs-funcs-data.R")

### c-score of all sites ----
(cscores = cs_sim(veg, nullmethod = "tswap", n = 5000)) 

# siginificant overall co-oc pattern. z > 2, suggesting competition?

### pairwise co-oc ----
## who co-oc with who?

pairs.all = pair.ff(comm = veg, nullmethod = "tswap", formu = "(A-J)*(B-J)", nsim = 5000, burn = 5000)
pairs.all = p.adj(pairs.all)
head(pairs.all)

pairs.pos = filter(pairs.all, SES > 0 & p.value < 0.05)
pairs.neg = filter(pairs.all, SES < 0 & p.value < 0.05)

pairs.all.summary = data.frame(sp.rich = length(unique(c(pairs.all$sp1_name, pairs.all$sp2_name))),
            all.pos = sum(pairs.all$SES>0), all.neg = sum(pairs.all$SES<0), 
            sig.pos = sum(pairs.all$SES>0 & pairs.all$p.value < 0.05),
            sig.neg = sum(pairs.all$SES<0 & pairs.all$p.value < 0.05),
            all.pos.prop = sum(pairs.all$SES>0)/choose(length(unique(c(pairs.all$sp1_name, pairs.all$sp2_name))), 2), 
            all.neg.prop = sum(pairs.all$SES<0)/choose(length(unique(c(pairs.all$sp1_name, pairs.all$sp2_name))), 2), 
            sig.pos.prop = sum(pairs.all$SES>0 & pairs.all$p.value < 0.05)/choose(length(unique(c(pairs.all$sp1_name, pairs.all$sp2_name))), 2),
            sig.neg.prop = sum(pairs.all$SES<0 & pairs.all$p.value < 0.05)/choose(length(unique(c(pairs.all$sp1_name, pairs.all$sp2_name))), 2))

#### envi diff of pairs ----
explain.pos = get_together(veg = veg, cooc = pairs.pos, envi = envi, pos.neg = "pos", m = 5)
explain.neg = get_together(veg = veg, cooc = pairs.neg, envi = envi, pos.neg = "neg", m = 5)


save(cscores, pairs.all, pairs.all.summary, explain.pos, explain.neg, file = "testresult.RData")

load("testresult.RData")

explain.all.summary = rbind(explain.pos[[1]], explain.neg[[1]])
names(explain.all.summary)[2] = "Explanation"
names(explain.all.summary)[3] = "Proportion"
ggplot(explain.all.summary, aes(x=variable, y=Proportion, fill = Explanation)) + 
  #   scale_fill_grey()+
  geom_bar(stat="identity") +
  facet_wrap(~posneg) +
  theme(legend.position="top", axis.title.x=element_blank()) +
  guides(fill = guide_legend(title = "Causal explanation"))
