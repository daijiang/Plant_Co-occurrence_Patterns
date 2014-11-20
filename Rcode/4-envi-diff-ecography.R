# will run on Darwin00
source("Rcode/0-pkgs-funcs-data.R")

## results from darwin00 ----
# get from codes in 3-func-dist.R
# 
# load("coocSig.csp2.RData")
# load("coocSig.nuf2.RData")
# load("coocSig.suf2.RData")

## ---------
coocsig.csp1950.pos = filter(cooc.csp.1950.wide.sq, SES > 0 & p.value < 0.05)
coocsig.csp1950.neg = filter(cooc.csp.1950.wide.sq, SES < 0 & p.value < 0.05)
coocsig.csp2000.pos = filter(cooc.csp.2000.wide.sq, SES > 0 & p.value < 0.05)
coocsig.csp2000.neg = filter(cooc.csp.2000.wide.sq, SES < 0 & p.value < 0.05)
coocsig.nuf1950.pos = filter(cooc.nuf.1950.wide.sq, SES > 0 & p.value < 0.05)
coocsig.nuf1950.neg = filter(cooc.nuf.1950.wide.sq, SES < 0 & p.value < 0.05)
coocsig.nuf2000.pos = filter(cooc.nuf.2000.wide.sq, SES > 0 & p.value < 0.05)
coocsig.nuf2000.neg = filter(cooc.nuf.2000.wide.sq, SES < 0 & p.value < 0.05)
coocsig.suf1950.pos = filter(cooc.suf.1950.wide.sq, SES > 0 & p.value < 0.05)
coocsig.suf1950.neg = filter(cooc.suf.1950.wide.sq, SES < 0 & p.value < 0.05)
coocsig.suf2000.pos = filter(cooc.suf.2000.wide.sq, SES > 0 & p.value < 0.05)
coocsig.suf2000.neg = filter(cooc.suf.2000.wide.sq, SES < 0 & p.value < 0.05)

coocsig.csp2000.pos.subset = filter(cooc.csp.2000.wide.sq.subset, SES > 0 & p.value < 0.05)
coocsig.csp2000.neg.subset = filter(cooc.csp.2000.wide.sq.subset, SES < 0 & p.value < 0.05)
coocsig.nuf2000.pos.subset = filter(cooc.nuf.2000.wide.sq.subset, SES > 0 & p.value < 0.05)
coocsig.nuf2000.neg.subset = filter(cooc.nuf.2000.wide.sq.subset, SES < 0 & p.value < 0.05)
coocsig.suf2000.pos.subset = filter(cooc.suf.2000.wide.sq.subset, SES > 0 & p.value < 0.05)
coocsig.suf2000.neg.subset = filter(cooc.suf.2000.wide.sq.subset, SES < 0 & p.value < 0.05)



## actual test for veg types and date -----
# csp  -------
csp.1950.pos.explain.prop = get_together(veg = veg.csp.1950.wide.site, cooc = coocsig.csp1950.pos, veg.type = "csp",
                                         envi = envi.csp.all, climate.period = "55_59", pos.neg = "pos", m = 5)
csp.1950.neg.explain.prop = get_together(veg = veg.csp.1950.wide.site, cooc = coocsig.csp1950.neg, veg.type = "csp",
                                         envi = envi.csp.all, climate.period = "55_59", pos.neg = "neg", m = 5)
csp.2000.pos.explain.prop = get_together(veg = veg.csp.2000.wide.site, cooc = coocsig.csp2000.pos, veg.type = "csp",
                                         envi = envi.csp.all, climate.period = "02_06", pos.neg = "pos", m = 5)
csp.2000.neg.explain.prop = get_together(veg = veg.csp.2000.wide.site, cooc = coocsig.csp2000.neg, veg.type = "csp",
                                         envi = envi.csp.all, climate.period = "02_06", pos.neg = "neg", m = 5)

csp.2000.pos.explain.prop.subset = get_together(veg = veg.csp.2000.wide.site.subset, cooc = coocsig.csp2000.pos.subset, veg.type = "csp",
                                                envi = envi.csp.all, climate.period = "02_06", pos.neg = "pos", m = 5)
csp.2000.neg.explain.prop.subset = get_together(veg = veg.csp.2000.wide.site.subset, cooc = coocsig.csp2000.neg.subset, veg.type = "csp",
                                                envi = envi.csp.all, climate.period = "02_06", pos.neg = "neg", m = 5)

# suf  -------
suf.1950.pos.explain.prop = get_together(veg = veg.suf.1950.wide.site, cooc = coocsig.suf1950.pos, veg.type = "suf",
                                         envi = envi.ns.all, climate.period = "50_56", pos.neg = "pos", m = 5)
suf.1950.neg.explain.prop = get_together(veg = veg.suf.1950.wide.site, cooc = coocsig.suf1950.neg, veg.type = "suf",
                                         envi = envi.ns.all, climate.period = "50_56", pos.neg = "neg", m = 5)
suf.2000.pos.explain.prop = get_together(veg = veg.suf.2000.wide.site, cooc = coocsig.suf2000.pos, veg.type = "suf",
                                         envi = envi.ns.all, climate.period = "00_06", pos.neg = "pos", m = 5)
suf.2000.neg.explain.prop = get_together(veg = veg.suf.2000.wide.site, cooc = coocsig.suf2000.neg, veg.type = "suf",
                                         envi = envi.ns.all, climate.period = "00_06", pos.neg = "neg", m = 5)

suf.2000.pos.explain.prop.subset = get_together(veg = veg.suf.2000.wide.site.subset, cooc = coocsig.suf2000.pos.subset, veg.type = "suf",
                                                envi = envi.ns.all, climate.period = "00_06", pos.neg = "pos", m = 5)
suf.2000.neg.explain.prop.subset = get_together(veg = veg.suf.2000.wide.site.subset, cooc = coocsig.suf2000.neg.subset, veg.type = "suf",
                                                envi = envi.ns.all, climate.period = "00_06", pos.neg = "neg", m = 5)


# nuf  -------
nuf.1950.pos.explain.prop = get_together(veg = veg.nuf.1950.wide.site, cooc = coocsig.nuf1950.pos, veg.type = "nuf",
                                         envi = envi.ns.all, climate.period = "50_56", pos.neg = "pos", m = 5)
nuf.1950.neg.explain.prop = get_together(veg = veg.nuf.1950.wide.site, cooc = coocsig.nuf1950.neg, veg.type = "nuf",
                                         envi = envi.ns.all, climate.period = "50_56", pos.neg = "neg", m = 5)
nuf.2000.pos.explain.prop = get_together(veg = veg.nuf.2000.wide.site, cooc = coocsig.nuf2000.pos, veg.type = "nuf",
                                         envi = envi.ns.all, climate.period = "00_06", pos.neg = "pos", m = 5)
nuf.2000.neg.explain.prop = get_together(veg = veg.nuf.2000.wide.site, cooc = coocsig.nuf2000.neg, veg.type = "nuf",
                                         envi = envi.ns.all, climate.period = "00_06", pos.neg = "neg", m = 5)

nuf.2000.pos.explain.prop.subset = get_together(veg = veg.nuf.2000.wide.site.subset, cooc = coocsig.nuf2000.pos.subset, veg.type = "nuf",
                                                envi = envi.ns.all, climate.period = "00_06", pos.neg = "pos", m = 5)
nuf.2000.neg.explain.prop.subset = get_together(veg = veg.nuf.2000.wide.site.subset, cooc = coocsig.nuf2000.neg.subset, veg.type = "nuf",
                                                envi = envi.ns.all, climate.period = "00_06", pos.neg = "neg", m = 5)

# 
# save(csp.1950.pos.explain.prop, csp.1950.neg.explain.prop,
#      csp.2000.pos.explain.prop, csp.2000.neg.explain.prop,
#      suf.1950.pos.explain.prop, suf.1950.neg.explain.prop,
#      suf.2000.pos.explain.prop, suf.2000.neg.explain.prop,
#      nuf.1950.pos.explain.prop, nuf.1950.neg.explain.prop,
#      nuf.2000.pos.explain.prop, nuf.2000.neg.explain.prop,
#      csp.2000.pos.explain.prop.subset, csp.2000.neg.explain.prop.subset,
#      suf.2000.pos.explain.prop.subset, suf.2000.neg.explain.prop.subset,
#      nuf.2000.pos.explain.prop.subset, nuf.2000.neg.explain.prop.subset, file = "detailEcography2.RData")

## combine all together into one data frame -------------------
csp.1950.pos.explain.prop$p.prop$data.source = "alldata"
csp.1950.neg.explain.prop$p.prop$data.source = "alldata"
csp.2000.pos.explain.prop$p.prop$data.source = "alldata"
csp.2000.neg.explain.prop$p.prop$data.source = "alldata"
csp.2000.pos.explain.prop.subset$p.prop$data.source = "subset"
csp.2000.neg.explain.prop.subset$p.prop$data.source = "subset"
suf.1950.pos.explain.prop$p.prop$data.source = "alldata"
suf.1950.neg.explain.prop$p.prop$data.source = "alldata"
suf.2000.pos.explain.prop$p.prop$data.source = "alldata"
suf.2000.neg.explain.prop$p.prop$data.source = "alldata"
suf.2000.pos.explain.prop.subset$p.prop$data.source = "subset"
suf.2000.neg.explain.prop.subset$p.prop$data.source = "subset"
nuf.1950.pos.explain.prop$p.prop$data.source = "alldata"
nuf.1950.neg.explain.prop$p.prop$data.source = "alldata"
nuf.2000.pos.explain.prop$p.prop$data.source = "alldata"
nuf.2000.neg.explain.prop$p.prop$data.source = "alldata"
nuf.2000.pos.explain.prop.subset$p.prop$data.source = "subset"
nuf.2000.neg.explain.prop.subset$p.prop$data.source = "subset"

all.explain.prop = rbind(csp.1950.pos.explain.prop$p.prop,
                         csp.1950.neg.explain.prop$p.prop,
                         csp.2000.pos.explain.prop$p.prop,
                         csp.2000.neg.explain.prop$p.prop,
                         csp.2000.pos.explain.prop.subset$p.prop,
                         csp.2000.neg.explain.prop.subset$p.prop,
                         suf.1950.pos.explain.prop$p.prop,
                         suf.1950.neg.explain.prop$p.prop,
                         suf.2000.pos.explain.prop$p.prop,
                         suf.2000.neg.explain.prop$p.prop,
                         suf.2000.pos.explain.prop.subset$p.prop,
                         suf.2000.neg.explain.prop.subset$p.prop,
                         nuf.1950.pos.explain.prop$p.prop,
                         nuf.1950.neg.explain.prop$p.prop,
                         nuf.2000.pos.explain.prop$p.prop,
                         nuf.2000.neg.explain.prop$p.prop,
                         nuf.2000.pos.explain.prop.subset$p.prop,
                         nuf.2000.neg.explain.prop.subset$p.prop)
all.explain.prop$posneg[all.explain.prop$posneg == "pos"] = "Aggregated"
all.explain.prop$posneg[all.explain.prop$posneg == "neg"] = "Segregated"
all.explain.prop$date[all.explain.prop$date %in% c("50_56", "55_59")] = "1950s"
all.explain.prop$date[all.explain.prop$date %in% c("00_06", "02_06")] = "2000s"
all.explain.prop$variable = str_sub(all.explain.prop$variable, 3)

## -------------------

all.explain.prop2 = filter(all.explain.prop, variable != "soil_pc2" & variable != "landsc_pc2")
names(all.explain.prop2)[2] = "Explanation"
names(all.explain.prop2)[3] = "Proportion"
all.explain.prop2$variable[all.explain.prop2$variable == "soil_pc1"] = "Soil"
all.explain.prop2$variable[all.explain.prop2$variable == "landsc_pc1"] = "Landscape"
all.explain.prop2$variable[all.explain.prop2$variable == "climate_pc1"] = "Climate"
all.explain.prop2$variable[all.explain.prop2$variable == "shade"] = "Shade"
all.explain.prop2$variable = factor(all.explain.prop2$variable, levels = c("Landscape", "Shade", "Soil", "Climate"))
all.explain.prop2$vegtype = toupper(all.explain.prop2$vegtype)

q1 = ggplot(subset(all.explain.prop2, vegtype == "csp"), aes(x=variable, y=Proportion, fill = Explaination)) + 
  #   scale_fill_grey()+
  geom_bar(stat="identity") + facet_grid(vegtype~date+posneg) + 
  theme(legend.position="top", axis.title.x=element_blank()) +
  scale_fill_manual(values=c("grey80", "grey50", "grey20", "black"),
                    guide=guide_legend(title="Inferring"))
q2 = ggplot(subset(all.explain.prop2, vegtype != "csp"), aes(x=variable, y=Proportion, fill = Explaination)) + 
  geom_bar(stat="identity") + facet_grid(vegtype~date+posneg)+ 
  theme(legend.position="none", axis.title.x=element_blank()) +
  scale_fill_manual(values=c("grey80", "grey50", "grey20", "black"),
                    guide=guide_legend(title="Inferring"))

library(dli55)
pdf("output/sitediff_geodistance_bw.pdf", width=10, height=10)
multiplot(q1,q2,layout = matrix(c(1,1,2,2,2), nrow = 5))
dev.off()

# all.explain.prop2$Explaination = as.character(all.explain.prop2$Explaination)

## plot using all data ----

q1.color = ggplot(subset(all.explain.prop2, vegtype == "CSP" & data.source == "alldata"), aes(x=variable, y=Proportion, fill = Explanation)) + 
  #   scale_fill_grey()+
  geom_bar(stat="identity") + facet_grid(vegtype~date+posneg) + 
  theme(legend.position="top", axis.title.x=element_blank()) +
  guides(fill = guide_legend(title = "Causal explanation"))
q2.color = ggplot(filter(all.explain.prop2, vegtype != "CSP", data.source == "alldata"), aes(x=variable, y=Proportion, fill = Explanation)) + 
  geom_bar(stat="identity") + facet_grid(vegtype~date+posneg)+ 
  theme(legend.position="none", axis.title.x=element_blank()) 
pdf("../Figs/sitediff_geodistance_alldata.pdf", width=10, height=10)
multiplot(q1.color,q2.color,layout = matrix(c(1,1,2,2,2), nrow = 5))
dev.off()

## plot using subset data ----
all.explain.prop2.subset = all.explain.prop2[-which(all.explain.prop2$date == "2000s" & 
                                                      all.explain.prop2$data.source == "alldata"), ]
q1.color.subset = ggplot(subset(all.explain.prop2, vegtype == "CSP" & data.source == "alldata"),
                         aes(x=variable, y=Proportion, fill = Explanation)) + 
  #   scale_fill_grey()+
  geom_bar(stat="identity") + facet_grid(vegtype~date+posneg) + 
  theme(legend.position="top", axis.title.x=element_blank()) +
  guides(fill = guide_legend(title = "Causal explanation"))
q2.color.subset = ggplot(filter(all.explain.prop2.subset, vegtype != "CSP"),
                  aes(x=variable, y=Proportion, fill = Explanation)) + 
  geom_bar(stat="identity") + facet_grid(vegtype~date+posneg)+ 
  theme(legend.position="none", axis.title.x=element_blank()) 
pdf("../Figs/sitediff_geodistance_subset.pdf", width=10, height=10)
multiplot(q1.color.subset,q2.color.subset,layout = matrix(c(1,1,2,2,2), nrow = 5))
dev.off()

## all details for pairs, envi filtering by diff variables?
## thus overall envi filtering is strong?
csp.1950.pos.explain.prop$pairs.explain$type = "csp.1950.pos"
csp.1950.neg.explain.prop$pairs.explain$type = "csp.1950.neg"
csp.2000.pos.explain.prop$pairs.explain$type = "csp.2000.pos"
csp.2000.neg.explain.prop$pairs.explain$type = "csp.2000.neg"
suf.1950.pos.explain.prop$pairs.explain$type = "suf.1950.pos"
suf.1950.neg.explain.prop$pairs.explain$type = "suf.1950.neg"
suf.2000.pos.explain.prop$pairs.explain$type = "suf.2000.pos"
suf.2000.neg.explain.prop$pairs.explain$type = "suf.2000.neg"
nuf.1950.pos.explain.prop$pairs.explain$type = "nuf.1950.pos"
nuf.1950.neg.explain.prop$pairs.explain$type = "nuf.1950.neg"
nuf.2000.pos.explain.prop$pairs.explain$type = "nuf.2000.pos"
nuf.2000.neg.explain.prop$pairs.explain$type = "nuf.2000.neg"

all.pairs.explain.csp = rbind(
  csp.1950.pos.explain.prop$pairs.explain,
  csp.1950.neg.explain.prop$pairs.explain,
  csp.2000.pos.explain.prop$pairs.explain,
  csp.2000.neg.explain.prop$pairs.explain)

all.pairs.explain.nsuf = rbind(
  suf.1950.pos.explain.prop$pairs.explain,
  suf.1950.neg.explain.prop$pairs.explain,
  suf.2000.pos.explain.prop$pairs.explain,
  suf.2000.neg.explain.prop$pairs.explain,
  nuf.1950.pos.explain.prop$pairs.explain,
  nuf.1950.neg.explain.prop$pairs.explain,
  nuf.2000.pos.explain.prop$pairs.explain,
  nuf.2000.neg.explain.prop$pairs.explain)

library(tidyr)
all.pairs.explain.csp =all.pairs.explain.csp %>%
  separate(type, into = c("veg", "date", "posneg"), sep = "\\.")
all.pairs.explain.nsuf = all.pairs.explain.nsuf %>%
  separate(type, into = c("veg", "date", "posneg"), sep = "\\.")

names(all.pairs.explain.csp)
all.pairs.explain.csp.naomit = ddply(na.omit(all.pairs.explain.csp), 
                                     .(date, posneg), mutate, 
                                     pair.num = seq_len(length(veg)))
all.pairs.explain.csp.naomit.long = melt(all.pairs.explain.csp.naomit, measure.vars = 7:10)
ggplot(filter(all.pairs.explain.csp.naomit.long, variable != "p.soil_pc2"), 
       aes(x = variable, y = pair.num)) + 
  geom_tile(aes(fill = value)) + 
  facet_grid(date~posneg, scales = "free") +
  theme(legend.position="top")

all.pairs.explain.nsuf.naomit = ddply(na.omit(all.pairs.explain.nsuf), 
                                      .(veg, date, posneg), mutate, 
                                      pair.num = seq_len(length(veg)))
all.pairs.explain.nsuf.naomit.long = melt(all.pairs.explain.nsuf.naomit,
                                          measure.vars = 7:11)
ggplot(na.omit(filter(all.pairs.explain.nsuf.naomit.long, 
                      variable != "p.soil_pc2" & variable != "p.landsc_pc2")), 
       aes(x = variable, y = pair.num)) + 
  geom_tile(aes(fill = value)) + 
  facet_grid(date+ veg~posneg , scales = "free") +
  theme(legend.position="top")

## relationship between envi diff and trait diff?
