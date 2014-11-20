# Blois et al 2014 test
D Li  
11/20/2014  

This project tests the effectiveness of Blois *et al* [2014](http://onlinelibrary.wiley.com/doi/10.1111/ecog.00779/abstract) Ecography paper's framework to disentangle different mechanisms for co-occurrent species. It has the following steps:

1. Test the randomness of overall co-occurrence patterns. Results were stored in `cscores`.
2. Test the significance of co-occurrence for all possible species pairs. In the test dataset, we have 100 speciese, so we have `choose(100, 2) = 4950` pairs. We will collect all significant positive (`pairs.pos`) and negative (`pairs.neg`) co-occurrent species pairs.
3. Test the underlying mechanisms of each of these significant pairs. Particularly, we only tested *environmental filtering* and *species interactions* here, without testing *dispersal limitation*. For positive pairs, we compared sites that they both present and sites that they both absent. If these two groups of sites are significantly different in environmental conditions, then we infer the possible mechanism as environmental filtering. If not, we infer as species interaction.
  - However, they are some issues:
      + we may not include all important environmental varibles.
      + for the same pair, we may conclude that it was shaped by environmental filtering by one environmental variable and species interaction by another environmental variable.
4. Get the proportion of species pairs that can be explained by either mechanism for each environmental variable.

See R codes for details.

### Data

```r
## data----
load("xy.Rdata")
envi = data.frame(x)
row.names(envi) = paste0("site", 1:2500)
colnames(envi) = paste0("envi", 1:5)
envi$site = row.names(envi) 

veg = data.frame(y)
row.names(veg) = paste0("site", 1:2500)
colnames(veg) = paste0("sp", 1:100)
rm(x, y)

# sort(colSums(veg))
# 
# ## plot
par(mfrow=c(3,2))
hist(envi$envi1)
hist(envi$envi2)
hist(envi$envi3)
hist(envi$envi4)
hist(envi$envi5)
hist(rowSums(veg))
```

<img src="http://i.imgur.com/lonbhcV.png" title="plot of chunk data" alt="plot of chunk data" style="display: block; margin: auto;" />

```r
# hist(colSums(veg))
# par(mfrow=c(1,1))
# plot(envi$envi1, rowSums(veg))
```



```r
# load the results so can avoid re-run
load("testresult.RData")
```

### Overall co-occurrence pattern (site-level)
Then let's look at the overall randomness of co-occurrence patterns. C-score is the average of all checkboad units. Thus, higher cscores suggest stronger competition among species. But we cannot naively say that species competition is the main mechanisms. Because things like habitat preferences, dispersal ability, etc. can also get the same pattern.

```r
cscores = cs_sim_fixsp_fixsite(veg)
```

```r
cscores
```

```
## $z
## [1] 2.143
## 
## $means
## [1] 78693
## 
## $pval
## [1] 0.000999
## 
## $statistic
## statistic 
##     78801
```
It is highly significant, which suggests that the overall co-occurrence pattern is non-random. Potentially, species competition is very important.

### Pairwise co-occurrence patterns (species-level)
Then we will test the significance of co-occurrence for each species pair.

```r
pairs.all = pair.ff(comm = veg, nsim = 5000, burn = 5000)
pairs.all = p.adj(pairs.all)

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
```

```r
head(pairs.all)
```

```
##   pair_sp obs sim.mean      SES   p.value sp1_name sp2_name
## 1 sp2 sp1  27    20.79  3.42874 0.0003124      sp2      sp1
## 2 sp3 sp1  33    28.56  2.21702 0.1768612      sp3      sp1
## 3 sp4 sp1  31    31.11 -0.05172 0.9948682      sp4      sp1
## 4 sp5 sp1  14    15.38 -1.51110 0.2964070      sp5      sp1
## 5 sp6 sp1  21    22.03 -0.35432 1.0000000      sp6      sp1
## 6 sp7 sp1  26    31.64 -1.46131 0.4097048      sp7      sp1
```

```r
head(pairs.pos)
```

```
##    pair_sp obs sim.mean   SES   p.value sp1_name sp2_name
## 1  sp2 sp1  27    20.79 3.429 0.0003124      sp2      sp1
## 2  sp8 sp1  33    28.23 2.251 0.0003124      sp8      sp1
## 3  sp9 sp1  37    33.37 1.539 0.0003124      sp9      sp1
## 4 sp14 sp1  29    18.04 9.377 0.0003124     sp14      sp1
## 5 sp18 sp1  25    18.29 2.200 0.0003124     sp18      sp1
## 6 sp23 sp1  41    33.68 3.670 0.0003124     sp23      sp1
```

```r
head(pairs.neg)
```

```
##    pair_sp obs sim.mean    SES   p.value sp1_name sp2_name
## 1 sp16 sp1  18    27.78 -6.551 0.0003124     sp16      sp1
## 2 sp20 sp1  16    22.21 -2.502 0.0003124     sp20      sp1
## 3 sp21 sp1  17    22.59 -2.990 0.0003124     sp21      sp1
## 4 sp24 sp1  21    27.52 -3.058 0.0003124     sp24      sp1
## 5 sp30 sp1  13    20.15 -5.056 0.0003124     sp30      sp1
## 6 sp34 sp1  27    35.53 -2.009 0.0003124     sp34      sp1
```

```r
pairs.all.summary
```

```
##   sp.rich all.pos all.neg sig.pos sig.neg all.pos.prop all.neg.prop sig.pos.prop
## 1     100    2275    2675    1518    1771       0.4596       0.5404       0.3067
##   sig.neg.prop
## 1       0.3578
```
I found 65% of species pairs to be significant, which is very high compare with real dataset.

### Get explanation for co-occurrence pairs

```r
#### envi diff of pairs ----
explain.pos = get_together(veg = veg, cooc = pairs.pos, envi = envi, pos.neg = "pos", m = 5)
explain.neg = get_together(veg = veg, cooc = pairs.neg, envi = envi, pos.neg = "neg", m = 5)
```

```r
explain.all.summary = rbind(explain.pos[[1]], explain.neg[[1]])
names(explain.all.summary)[2] = "Explanation"
names(explain.all.summary)[3] = "Proportion"
explain.all.summary 
```

```
##    variable             Explanation Proportion posneg
## 1   p.envi1 Environmental filtering     0.6779    pos
## 2   p.envi1     Species interaction     0.3221    pos
## 3   p.envi2 Environmental filtering     0.6831    pos
## 4   p.envi2     Species interaction     0.3169    pos
## 5   p.envi3 Environmental filtering     0.7220    pos
## 6   p.envi3     Species interaction     0.2780    pos
## 7   p.envi4 Environmental filtering     0.6963    pos
## 8   p.envi4     Species interaction     0.3037    pos
## 9   p.envi5 Environmental filtering     0.7075    pos
## 10  p.envi5     Species interaction     0.2925    pos
## 11  p.envi1 Environmental filtering     0.1976    neg
## 12  p.envi1     Species interaction     0.8024    neg
## 13  p.envi2 Environmental filtering     0.3128    neg
## 14  p.envi2     Species interaction     0.6872    neg
## 15  p.envi3 Environmental filtering     0.3456    neg
## 16  p.envi3     Species interaction     0.6544    neg
## 17  p.envi4 Environmental filtering     0.3840    neg
## 18  p.envi4     Species interaction     0.6160    neg
## 19  p.envi5 Environmental filtering     0.3196    neg
## 20  p.envi5     Species interaction     0.6804    neg
```

```r
library(ggplot2)
ggplot(explain.all.summary, aes(x=variable, y=Proportion, fill = Explanation)) + 
  geom_bar(stat="identity") +
  facet_wrap(~posneg) +
  theme(legend.position="top", axis.title.x=element_blank()) +
  guides(fill = guide_legend(title = "Causal explanation"))
```

<img src="http://i.imgur.com/PO3lrOr.png" title="plot of chunk envi2" alt="plot of chunk envi2" style="display: block; margin: auto;" />
For positive co-occurrent pairs, about 70% of them were driven by environmental filtering. Other ~30% were shaped by potential interactions, e.g. facilitation, among species. For negative co-occurrent pairs, about 75% of them were driven by biotic competition. 

### How good is this framework?


