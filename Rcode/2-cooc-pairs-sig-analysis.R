# will run on Darwin00
source("Rcode/0-pkgs-funcs-data.R")

csp.1950.long = subset(all, date == "1950s" & type2 == "CSP")
csp.2000.long = subset(all, date == "2000s" & type2 == "CSP")
nuf.1950.long = subset(all, date == "1950s" & type2 == "NUF")
nuf.2000.long = subset(all, date == "2000s" & type2 == "NUF")
suf.1950.long = subset(all, date == "1950s" & type2 == "SUF")
suf.2000.long = subset(all, date == "2000s" & type2 == "SUF")

### ---- CSP 1950 ----
## change to s_q by sp matrix
# csp.1950.wide.sq = to.wide.sq(csp.1950.long)

## calculate who co-oc with who and their significance
# cooc.csp.1950.wide.sq=pair.ff(csp.1950.wide.sq, nsim = 5000, burn = 5000)
# cooc.csp.1950.wide.sq = p.adj(cooc.csp.1950.wide.sq) # adjust for multiple comparison
cooc.csp.1950.wide.sq = read.csv("cooc_subset/csp1950coocSig.all2.csv", stringsAsFactors = F)
  
pos.cooc.csp.1950.wide.sq = subset(cooc.csp.1950.wide.sq, SES>0 & p.value < 0.05)

## Phylo/func distance
  # phylo and trait data from Rcode/0-pkgs-and-data.R 

splist.csp.1950 = unique(csp.1950.long$sp) # all sp list of 30 sites 20 quad

csp.1950.pos.dist = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                   sp.list = splist.csp.1950, 
                                   pairs = pos.cooc.csp.1950.wide.sq, N = 5000)

###--- negative co-oc sp
neg.cooc.csp.1950.wide.sq = subset(cooc.csp.1950.wide.sq, SES < 0 & p.value < 0.05)
csp.1950.neg.dist = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                           sp.list = splist.csp.1950, 
                                           pairs = neg.cooc.csp.1950.wide.sq, N = 5000)

### ---- CSP 2000 all data ----
## change to s_q by sp matrix
# csp.2000.wide.sq = to.wide.sq(csp.2000.long)
# 
# ## calculate who co-oc with who and their significance
# cooc.csp.2000.wide.sq=pair.ff(csp.2000.wide.sq, nsim = 5000, burn = 5000)
# cooc.csp.2000.wide.sq = p.adj(cooc.csp.2000.wide.sq) # adjust for multiple comparison
cooc.csp.2000.wide.sq = read.csv("cooc_subset/csp2000coocSig.all2.csv", stringsAsFactors = F)
pos.cooc.csp.2000.wide.sq = subset(cooc.csp.2000.wide.sq, SES>0 & p.value < 0.05)

## Phylo/func distance
# phylo and trait data from Rcode/0-pkgs-and-data.R 

splist.csp.2000 = unique(csp.2000.long$sp) # all sp list of 30 sites 20 quad

csp.2000.pos.dist = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                           sp.list = splist.csp.2000, 
                                           pairs = pos.cooc.csp.2000.wide.sq, N = 5000)

###--- negative co-oc sp
neg.cooc.csp.2000.wide.sq = subset(cooc.csp.2000.wide.sq, SES < 0 & p.value < 0.05)
csp.2000.neg.dist = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                           sp.list = splist.csp.2000, 
                                           pairs = neg.cooc.csp.2000.wide.sq, N = 5000)

### ---- CSP 2000 subset data ----
## change to s_q by sp matrix
# csp.2000.long.subset = long_subset_by_1950s(csp.2000.long)
# csp.2000.wide.sq.subset = to.wide.sq(csp.2000.long.subset)
# 
# ## calculate who co-oc with who and their significance
# cooc.csp.2000.wide.sq.subset=pair.ff(csp.2000.wide.sq.subset, nsim = 5000, burn = 5000)
# cooc.csp.2000.wide.sq.subset = p.adj(cooc.csp.2000.wide.sq.subset) # adjust for multiple comparison
cooc.csp.2000.wide.sq.subset = read.csv("cooc_subset/csp2000coocSig.all2.subset.csv", stringsAsFactors = F)
pos.cooc.csp.2000.wide.sq.subset = subset(cooc.csp.2000.wide.sq.subset, SES>0 & p.value < 0.05)

## Phylo/func distance
# phylo and trait data from Rcode/0-pkgs-and-data.R 

splist.csp.2000.subset = unique(csp.2000.long.subset$sp) # all sp list of 30 sites 20 quad

csp.2000.pos.dist.subset = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                           sp.list = splist.csp.2000.subset, 
                                           pairs = pos.cooc.csp.2000.wide.sq.subset,
                                           N = 5000)

###--- negative co-oc sp
neg.cooc.csp.2000.wide.sq.subset = subset(cooc.csp.2000.wide.sq.subset, SES < 0 & p.value < 0.05)
csp.2000.neg.dist.subset = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                           sp.list = splist.csp.2000.subset, 
                                           pairs = neg.cooc.csp.2000.wide.sq.subset,
                                           N = 5000)

### ---- nuf 1950 ----
## change to s_q by sp matrix
# nuf.1950.wide.sq = to.wide.sq(nuf.1950.long)
# 
# ## calculate who co-oc with who and their significance
# cooc.nuf.1950.wide.sq=pair.ff(nuf.1950.wide.sq, nsim = 5000, burn = 5000)
# cooc.nuf.1950.wide.sq = p.adj(cooc.nuf.1950.wide.sq) # adjust for multiple comparison
cooc.nuf.1950.wide.sq = read.csv("cooc_subset/nuf1950coocSig.all2.csv", stringsAsFactors = F)
pos.cooc.nuf.1950.wide.sq = subset(cooc.nuf.1950.wide.sq, SES>0 & p.value < 0.05)

## Phylo/func distance
# phylo and trait data from Rcode/0-pkgs-and-data.R 

splist.nuf.1950 = unique(nuf.1950.long$sp) # all sp list of 30 sites 20 quad

nuf.1950.pos.dist = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                           sp.list = splist.nuf.1950, 
                                           pairs = pos.cooc.nuf.1950.wide.sq, N = 5000)

###--- negative co-oc sp
neg.cooc.nuf.1950.wide.sq = subset(cooc.nuf.1950.wide.sq, SES < 0 & p.value < 0.05)
nuf.1950.neg.dist = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                           sp.list = splist.nuf.1950, 
                                           pairs = neg.cooc.nuf.1950.wide.sq, N = 5000)

### ---- nuf 2000 all data ----
## change to s_q by sp matrix
# nuf.2000.wide.sq = to.wide.sq(nuf.2000.long)
# 
# ## calculate who co-oc with who and their significance
# cooc.nuf.2000.wide.sq=pair.ff(nuf.2000.wide.sq, nsim = 5000, burn = 5000)
# cooc.nuf.2000.wide.sq = p.adj(cooc.nuf.2000.wide.sq) # adjust for multiple comparison
cooc.nuf.2000.wide.sq = read.csv("cooc_subset/nuf2000coocSig.all2.csv", stringsAsFactors = F)

pos.cooc.nuf.2000.wide.sq = subset(cooc.nuf.2000.wide.sq, SES>0 & p.value < 0.05)

## Phylo/func distance
# phylo and trait data from Rcode/0-pkgs-and-data.R 

splist.nuf.2000 = unique(nuf.2000.long$sp) # all sp list of 30 sites 20 quad

nuf.2000.pos.dist = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                           sp.list = splist.nuf.2000, 
                                           pairs = pos.cooc.nuf.2000.wide.sq, N = 5000)

###--- negative co-oc sp
neg.cooc.nuf.2000.wide.sq = subset(cooc.nuf.2000.wide.sq, SES < 0 & p.value < 0.05)
nuf.2000.neg.dist = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                           sp.list = splist.nuf.2000, 
                                           pairs = neg.cooc.nuf.2000.wide.sq, N = 5000)

### ---- nuf 2000 subset data ----
## change to s_q by sp matrix
# nuf.2000.long.subset = long_subset_by_1950s(nuf.2000.long)
# nuf.2000.wide.sq.subset = to.wide.sq(nuf.2000.long.subset)
# 
# ## calculate who co-oc with who and their significance
# cooc.nuf.2000.wide.sq.subset=pair.ff(nuf.2000.wide.sq.subset, nsim = 5000, burn = 5000)
# cooc.nuf.2000.wide.sq.subset = p.adj(cooc.nuf.2000.wide.sq.subset) # adjust for multiple comparison
cooc.nuf.2000.wide.sq.subset = read.csv("cooc_subset/nuf2000coocSig.all2.subset.csv", stringsAsFactors = F)

pos.cooc.nuf.2000.wide.sq.subset = subset(cooc.nuf.2000.wide.sq.subset, SES>0 & p.value < 0.05)

## Phylo/func distance
# phylo and trait data from Rcode/0-pkgs-and-data.R 

splist.nuf.2000.subset = unique(nuf.2000.long.subset$sp) # all sp list of 30 sites 20 quad

nuf.2000.pos.dist.subset = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                                  sp.list = splist.nuf.2000.subset, 
                                                  pairs = pos.cooc.nuf.2000.wide.sq.subset,
                                                  N = 5000)

###--- negative co-oc sp
neg.cooc.nuf.2000.wide.sq.subset = subset(cooc.nuf.2000.wide.sq.subset, SES < 0 & p.value < 0.05)
nuf.2000.neg.dist.subset = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                                  sp.list = splist.nuf.2000.subset, 
                                                  pairs = neg.cooc.nuf.2000.wide.sq.subset,
                                                  N = 5000)

### ---- suf 1950 ----
## change to s_q by sp matrix
# suf.1950.wide.sq = to.wide.sq(suf.1950.long)
# 
# ## calculate who co-oc with who and their significance
# cooc.suf.1950.wide.sq=pair.ff(suf.1950.wide.sq, nsim = 5000, burn = 5000)
# cooc.suf.1950.wide.sq = p.adj(cooc.suf.1950.wide.sq) # adjust for multiple comparison
cooc.suf.1950.wide.sq = read.csv("cooc_subset/suf1950coocSig.all2.csv", stringsAsFactors = F)

pos.cooc.suf.1950.wide.sq = subset(cooc.suf.1950.wide.sq, SES>0 & p.value < 0.05)

## Phylo/func distance
# phylo and trait data from Rcode/0-pkgs-and-data.R 

splist.suf.1950 = unique(suf.1950.long$sp) # all sp list of 30 sites 20 quad

suf.1950.pos.dist = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                           sp.list = splist.suf.1950, 
                                           pairs = pos.cooc.suf.1950.wide.sq, N = 5000)

###--- negative co-oc sp
neg.cooc.suf.1950.wide.sq = subset(cooc.suf.1950.wide.sq, SES < 0 & p.value < 0.05)
suf.1950.neg.dist = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                           sp.list = splist.suf.1950, 
                                           pairs = neg.cooc.suf.1950.wide.sq, N = 5000)

### ---- suf 2000 all data ----
## change to s_q by sp matrix
# suf.2000.wide.sq = to.wide.sq(suf.2000.long)
# 
# ## calculate who co-oc with who and their significance
# cooc.suf.2000.wide.sq=pair.ff(suf.2000.wide.sq, nsim = 5000, burn = 5000)
# cooc.suf.2000.wide.sq = p.adj(cooc.suf.2000.wide.sq) # adjust for multiple comparison
cooc.suf.2000.wide.sq = read.csv("cooc_subset/suf2000coocSig.all2.csv", stringsAsFactors = F)

pos.cooc.suf.2000.wide.sq = subset(cooc.suf.2000.wide.sq, SES>0 & p.value < 0.05)

## Phylo/func distance
# phylo and trait data from Rcode/0-pkgs-and-data.R 

splist.suf.2000 = unique(suf.2000.long$sp) # all sp list of 30 sites 20 quad

suf.2000.pos.dist = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                           sp.list = splist.suf.2000, 
                                           pairs = pos.cooc.suf.2000.wide.sq, N = 5000)

###--- negative co-oc sp
neg.cooc.suf.2000.wide.sq = subset(cooc.suf.2000.wide.sq, SES < 0 & p.value < 0.05)
suf.2000.neg.dist = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                           sp.list = splist.suf.2000, 
                                           pairs = neg.cooc.suf.2000.wide.sq, N = 5000)

### ---- suf 2000 subset data ----
## change to s_q by sp matrix
# suf.2000.long.subset = long_subset_by_1950s(suf.2000.long)
# suf.2000.wide.sq.subset = to.wide.sq(suf.2000.long.subset)
# 
# ## calculate who co-oc with who and their significance
# cooc.suf.2000.wide.sq.subset=pair.ff(suf.2000.wide.sq.subset, nsim = 5000, burn = 5000)
# cooc.suf.2000.wide.sq.subset = p.adj(cooc.suf.2000.wide.sq.subset) # adjust for multiple comparison
cooc.suf.2000.wide.sq.subset = read.csv("cooc_subset/suf2000coocSig.all2.subset.csv", stringsAsFactors = F)

pos.cooc.suf.2000.wide.sq.subset = subset(cooc.suf.2000.wide.sq.subset, SES>0 & p.value < 0.05)

## Phylo/func distance
# phylo and trait data from Rcode/0-pkgs-and-data.R 

splist.suf.2000.subset = unique(suf.2000.long.subset$sp) # all sp list of 30 sites 20 quad

suf.2000.pos.dist.subset = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                                  sp.list = splist.suf.2000.subset, 
                                                  pairs = pos.cooc.suf.2000.wide.sq.subset,
                                                  N = 5000)

###--- negative co-oc sp
neg.cooc.suf.2000.wide.sq.subset = subset(cooc.suf.2000.wide.sq.subset, SES < 0 & p.value < 0.05)
suf.2000.neg.dist.subset = pairs.trait.null.model(trait.dist.matrix=trait.dist.all.matrix, 
                                                  sp.list = splist.suf.2000.subset, 
                                                  pairs = neg.cooc.suf.2000.wide.sq.subset,
                                                  N = 5000)
# save(cooc.csp.1950.wide.sq,
# pos.cooc.csp.1950.wide.sq,
# neg.cooc.csp.1950.wide.sq,
# csp.1950.pos.dist,
# csp.1950.neg.dist,
# cooc.csp.2000.wide.sq,
# pos.cooc.csp.2000.wide.sq,
# neg.cooc.csp.2000.wide.sq,
# csp.2000.pos.dist,
# csp.2000.neg.dist,
# cooc.csp.2000.wide.sq.subset,
# pos.cooc.csp.2000.wide.sq.subset,
# neg.cooc.csp.2000.wide.sq.subset,
# csp.2000.pos.dist.subset,
# csp.2000.neg.dist.subset,
# cooc.nuf.1950.wide.sq,
# pos.cooc.nuf.1950.wide.sq,
# neg.cooc.nuf.1950.wide.sq,
# nuf.1950.pos.dist,
# nuf.1950.neg.dist,
# cooc.nuf.2000.wide.sq,
# pos.cooc.nuf.2000.wide.sq,
# neg.cooc.nuf.2000.wide.sq,
# nuf.2000.pos.dist,
# nuf.2000.neg.dist,
# cooc.nuf.2000.wide.sq.subset,
# pos.cooc.nuf.2000.wide.sq.subset,
# neg.cooc.nuf.2000.wide.sq.subset,
# nuf.2000.pos.dist.subset,
# nuf.2000.neg.dist.subset,
# cooc.suf.1950.wide.sq,
# pos.cooc.suf.1950.wide.sq,
# neg.cooc.suf.1950.wide.sq,
# suf.1950.pos.dist,
# suf.1950.neg.dist,
# cooc.suf.2000.wide.sq,
# pos.cooc.suf.2000.wide.sq,
# neg.cooc.suf.2000.wide.sq,
# suf.2000.pos.dist,
# suf.2000.neg.dist,
# cooc.suf.2000.wide.sq.subset,
# pos.cooc.suf.2000.wide.sq.subset,
# neg.cooc.suf.2000.wide.sq.subset,
# suf.2000.pos.dist.subset,
# suf.2000.neg.dist.subset, file = "coocSig.all.RData")