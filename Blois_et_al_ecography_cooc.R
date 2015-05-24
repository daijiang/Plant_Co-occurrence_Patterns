library(plyr)
library(dplyr)
library(reshape2)
library(stringr)
library(vegan)

## classify sites for each pairs
# site classification, x is a vector of three, will be used in extract_site function
class_site = function(x){
  if (x[1] == 0 & x[2] == 0) {
    x[3] = "none"
  } else {
    if (x[1] == 0 & x[2] > 0) {
      x[3] = "b"
    } else {
      if (x[1] > 0 & x[2] == 0) {
        x[3] = "a"
      } else {
        x[3] = "both"
      }
    }
  }
}

# this function extract envi var, co-oc sites, and classify site types for a pair of species
# this function will be used in the get_p_envi_pairs() function below.
# veg is site by sp matrix
# cooc.pair has at least sp1_name and sp2_name columns
# envi is a dataframe with site, long, lat, and envi var columns
extract_site = function(veg, cooc.pair, envi){
  x = veg[c(cooc.pair$sp1_name, cooc.pair$sp2_name)]
  x$site.class = NA # three columns now
  x$site.class = apply(x, 1, class_site)
  # combine with site envi data
  x$site = if(any(str_detect(string = row.names(x), pattern = "_"))) { # some of my site names have _.
    gsub(x = row.names(x), "^([0-9.]*)_.*", "\\1")} else row.names(x)
  join(x, envi, by="site", type = "left")
}


# this function works on dataframe. It will compare envi var of groups of sites for all pairs.
# veg.df is a site by sp matrix
# cooc.df is a data frame for all sig pairs, including at least three columns:
# pair_sp, which is a combination of two sp, e.g. Acer rubrum Pinus strobus
# sp1_name
# sp2_name
# envi.df is a site by envi matrix, first 3 columns must be site lat and long, plus at least one other column.
# n how many sites at least one pair co-oc and not co-oc to be included in analysis? (to be statistically meaningful)

# output data frame will have these columns: sp_pair, a, b, both, none, p.geod, and p.values for other env variables.
get_p_envi_pairs = function(veg.df, cooc.df, envi.df, posneg, n = 5){
  data.pairs = dlply(cooc.df, .(pair_sp), function(x) extract_site(veg = veg.df, cooc.pair = x, envi = envi.df))
  #data.pairs$site.class = as.character(data.pairs$site.class)
  analysis.pairs = ldply(data.pairs, function(x){
    # print(names(x)[1:2])
    lx = data.frame(a = sum(x$site.class == "a"), b = sum(x$site.class == "b"), # how many sites in each type?
                    both = sum(x$site.class == "both"), none = sum(x$site.class == "none"))
    envi.var = names(x)[-c(1:6)] # no sp1 sp2 site.class site lat long
    if(posneg == "pos"){ # positive cooc: compare both and none
      data2 =  subset(x, site.class %in% c("both", "none"))
      if (lx$both >= n & lx$none >= n) { # at least n sites for both and none
        # geo distance compare
        p.geod = adonis(data2[c("long", "lat")] ~ data2$site.class, method = "euclidean")$aov.tab[6][1,1]
        p.envi = vector(length = length(envi.var))
        for (i in seq_along(envi.var)){
          if(min(ddply(select_(data2, "site.class", envi.var[i]), .(site.class), function(x) sum(!is.na(x[,2])))$V1) < 3){
            p.envi[i] = NA
          } else {
            p.envi[i] = t.test(as.formula(paste0(envi.var[i], "~ site.class")), data = data2)$p.value
          }
        }
        names(p.envi) = paste0("p.", envi.var)
        mean.envi.both.a = vector(length = length(envi.var))
        mean.envi.none.b = vector(length = length(envi.var))
        for (i in seq_along(envi.var)){
          # print(i)
          if(min(ddply(select_(data2, "site.class", envi.var[i]), .(site.class), function(x) sum(!is.na(x[,2])))$V1) < 3){
            mean.envi.both.a[i] = mean.envi.none.b[i] = NA
          } else {
            mean.envi.both.a[i] = t.test(as.formula(paste0(envi.var[i], "~ site.class")), data = data2, paired = F, var.equal = T)$estimate[1]
            mean.envi.none.b[i] = t.test(as.formula(paste0(envi.var[i], "~ site.class")), data = data2, paired = F, var.equal = T)$estimate[2]
          }
        }
        names(mean.envi.both.a) = paste0("mean.envi.both.a.", envi.var)
        names(mean.envi.none.b) = paste0("mean.envi.none.b.", envi.var)
        data.frame(lx, p.geod = p.geod, t(p.envi), t(mean.envi.both.a), t(mean.envi.none.b))
      } else { # if less than n sites for both and none, set as NA.
        mean.envi.both.a= mean.envi.none.b=p.envi = rep(NA, length(envi.var))
        names(p.envi) = paste0("p.", envi.var)
        names(mean.envi.both.a) = paste0("mean.envi.both.a.", envi.var)
        names(mean.envi.none.b) = paste0("mean.envi.none.b.", envi.var)
        data.frame(lx, p.geod = NA, t(p.envi), t(mean.envi.both.a), t(mean.envi.none.b))
      }} else { # negative cooc: compare a and b
        data2 = filter(x, site.class %in% c("a", "b"))
        if (lx$a >= n & lx$b >= n) {
          p.geod = adonis(data2[c("long", "lat")] ~ data2$site.class, method = "euclidean")$aov.tab[6][1,1]
          p.envi = vector(length = length(envi.var))
          for (i in seq_along(envi.var)){
            if(min(ddply(select_(data2, "site.class", envi.var[i]), .(site.class), function(x) sum(!is.na(x[,2])))$V1) < 3){
              p.envi[i] = NA
            } else {
              p.envi[i] = t.test(as.formula(paste0(envi.var[i], "~ site.class")), data = data2)$p.value
            }
          }
          names(p.envi) = paste0("p.", envi.var)
          mean.envi.both.a = vector(length = length(envi.var))
          mean.envi.none.b = vector(length = length(envi.var))
          for (i in seq_along(envi.var)){
            # print(i)
            if(min(ddply(select_(data2, "site.class", envi.var[i]), .(site.class), function(x) sum(!is.na(x[,2])))$V1) < 3){
              mean.envi.both.a[i] = mean.envi.none.b[i] = NA
            } else {
              mean.envi.both.a[i] = t.test(as.formula(paste0(envi.var[i], "~ site.class")), data = data2, paired = F, var.equal = T)$estimate[1]
              mean.envi.none.b[i] = t.test(as.formula(paste0(envi.var[i], "~ site.class")), data = data2, paired = F, var.equal = T)$estimate[2]
            }
          }
          names(mean.envi.both.a) = paste0("mean.envi.both.a.", envi.var)
          names(mean.envi.none.b) = paste0("mean.envi.none.b.", envi.var)
          data.frame(lx, p.geod = p.geod, t(p.envi), t(mean.envi.both.a), t(mean.envi.none.b))
          
        } else {
          mean.envi.both.a= mean.envi.none.b=p.envi = rep(NA, length(envi.var))
          names(p.envi) = paste0("p.", envi.var)
          names(mean.envi.both.a) = paste0("mean.envi.both.a.", envi.var)
          names(mean.envi.none.b) = paste0("mean.envi.none.b.", envi.var)
          data.frame(lx, p.geod = NA, t(p.envi), t(mean.envi.both.a), t(mean.envi.none.b))
        }
      }
  })
  analysis.pairs # a data frame, summarize for all pairs.
}

# this function is used to get explaination for each pair of species for each env variable
# p.wide is the result data frame from the above function: get_p_envi_pairs()
get_explain_envi_pairs = function(p.wide){
  n.var = (dim(p.wide)[2] - 6)/3
  p.long = melt(p.wide, id.vars = 1:6, measure.vars = 7:(7+n.var-1), value.name = "p.value") # 1:6 columns: sp-pair a b both none p.geod
  p.long2 = ddply(p.long, .(variable), function(x){
    x2 = adply(x, 1, function(x1){
      if(!is.na(x1$p.geod)){
        if(is.na(x1$p.value)){x1$explain = NA} else{
          if(x1$p.geod <= 0.05){
            if(x1$p.value <= 0.05) x1$explain = "Dispersal limitation and/or Environmental filtering" else x1$explain = "Dispersal limitation"
          }else{ # p.geod > 0.05
            if(x1$p.value <= 0.05) x1$explain = "Environmental filtering" else x1$explain = "Species interaction"
          }
        }}else{
          x1$explain = NA
        }
    })
    x2$explain = x2$V1; x2$V1=NULL
    x2
  })
  mean.envi.both.a = melt(p.wide, id.vars = 1:6, measure.vars = (7+n.var):(7+2*n.var-1), 
                          variable.name = "both.a.variable", value.name = "mean.envi.both.a")
  mean.envi.none.b = melt(p.wide, id.vars = 1:6, measure.vars = (7+2*n.var):(7+3*n.var-1), 
                          variable.name = "none.b.variable", value.name = "mean.envi.none.b")
  p.long2$mean.envi.both.a = mean.envi.both.a$mean.envi.both.a
  p.long2$mean.envi.none.b = mean.envi.none.b$mean.envi.none.b
  p.long2
}

# this function is a wrap of all above ones.
# veg site by species matrix, site as rows.
# cooc cooc pairs dataframe (aggregated OR seggregated), must have three columns: pair_sp, sp1_name, sp2_name
# envi environmental dataframe, the first three columns must be: site, long, lat, other columns for other envi variables.
# climate.period which subset of envi data to use? I have 1950s and 2000s envi data in one dataframe.
# pos.neg "pos" for aggregated pairs, "neg" for seggregated pairs.
# m how many sites at least to be included in data analysis? Some species only co-oc in one or two sites and not meaningful to be included.
# veg.type: I have three veg types in my analysis...
get_together = function(veg, cooc, envi, pos.neg, m = 5){
  sitediff.p = get_p_envi_pairs(veg.df = veg, cooc.df = cooc,
                                envi.df = envi,
                                posneg = pos.neg, n = m)
  p.explain0 = get_explain_envi_pairs(sitediff.p)
  p.explain = na.omit(p.explain0)
  p.explain.prop = ddply(p.explain[c("variable", "explain")], 
                         .(variable), 
                         function(x) data.frame(prop.table(table(x$explain))))
  #   p.explain.prop$Var1 = factor(p.explain.prop$Var1, 
  #                                levels = c("Dispersal limitation", "Environmental filtering", 
  #                                           "Dispersal limitation and/or Environmental filtering", "Species interaction"))
  p.explain.prop$posneg = pos.neg
  list(p.prop = p.explain.prop, pairs.explain = dcast(p.explain0[, -c(8,10,11)], ...~variable, value.var = "explain"), 
       p.mean.detail = p.explain0)
}

## examples ---------
get_together(veg = veg.csp.1950.wide.site, cooc = coocsig.csp1950.pos, veg.type = "csp",
             envi = envi.csp.all, climate.period = "55_59", pos.neg = "pos", m = 5)
# examples of my data frames. veg.type can be any name.
# climate.period: since I have climate data for two time periods.. you may want to change the function above for your own work.
# pos.neg: I am using this function on positive co-oc pairs now.
# m: for positive pairs, how many sites do they need to both present or absent to be included ? 
# If the # sites they both occur is to low, say, 1, it is meaningless.

# head(veg.csp.1950.wide.site)
#          Acer rubrum Achillea millefolium Alnus incana Amelanchier spp Amphicarpaea bracteata Andropogon gerardii Anemone quinquefolia Antennaria spp
# 3293           0                    0            0               0                      0                   0                    0              0
# 3294           3                    0            0               1                      0                   0                    0              0
# 3295           8                    0            1               0                      0                   0                    0              0
# 3296           0                    1            0               0                      0                   5                    0              5
# 3297           0                    1            0               0                      0                   0                    3              3
# 3299          15                    0            0               2                      0                   0                    0              0

# head(coocsig.csp1950.pos)
#                           pair_sp obs sim.mean       SES      p.value           sp1_name    sp2_name
# 7   Aralia nudicaulis Acer rubrum  30  17.3288  6.341080 0.0008650752  Aralia nudicaulis Acer rubrum
# 11  Betula papyrifera Acer rubrum   9   2.8084  6.749735 0.0008650752  Betula papyrifera Acer rubrum
# 14 Clintonia borealis Acer rubrum  15   6.6922  5.687225 0.0008650752 Clintonia borealis Acer rubrum
# 16    Coptis trifolia Acer rubrum  14   3.1592  5.814496 0.0008650752    Coptis trifolia Acer rubrum
# 17  Cornus canadensis Acer rubrum  38  16.1272 10.517854 0.0008650752  Cornus canadensis Acer rubrum
# 21     Epigaea repens Acer rubrum   7   1.7630  4.444030 0.0008650752     Epigaea repens Acer rubrum

# head(envi.csp.all)
#           site      long      lat   soil_pc1   soil_pc2 shade   precip     tmax       tmin  date
# 3293_1958 3293 -90.22775 44.18436  1.5294029 -0.3692148  20.9 1.902397 13.13436  0.1288019 55_59
# 3294_1958 3294 -90.27711 44.18011  0.9239977 -2.2331857  45.1 1.902397 13.13436  0.1288019 55_59
# 3295_1958 3295 -90.30947 44.21922 -3.1552649  0.7124214  58.9 1.922669 12.95795 -0.1061756 55_59
# 3296_1958 3296 -90.31997 44.21583 -3.5212085  0.9327750  19.7 1.922669 12.95795 -0.1061756 55_59
# 3297_1958 3297 -90.30694 44.21197 -3.2778733 -1.8698072  56.6 1.922669 12.95795 -0.1061756 55_59
# 3299_1958 3299 -90.29964 44.14242  1.4882514  2.4954042  87.0 1.899665 13.19993  0.1736423 55_59
