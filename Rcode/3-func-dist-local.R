library(tidyr)

## results from darwin00 ----
# load("coocSig.all.RData")
# load("coocSig.csp2.RData")
# load("coocSig.nuf.RData")
# load("coocSig.suf.RData")

# read cooc pairs from cooc_subset file folder, which is came from 
# Darwin00:private/cooc/coocff.all2.R


## combine all functional distance summary for sig pairs into one ----
csp.1950.pos.dist$dist.obs.summary$type = "CSP"
csp.1950.neg.dist$dist.obs.summary$type = "CSP"
csp.2000.pos.dist$dist.obs.summary$type = "CSP"
csp.2000.neg.dist$dist.obs.summary$type = "CSP"
csp.2000.pos.dist.subset$dist.obs.summary$type = "CSP"
csp.2000.neg.dist.subset$dist.obs.summary$type = "CSP"
nuf.1950.pos.dist$dist.obs.summary$type = "NUF"
nuf.1950.neg.dist$dist.obs.summary$type = "NUF"
nuf.2000.pos.dist$dist.obs.summary$type = "NUF"
nuf.2000.neg.dist$dist.obs.summary$type = "NUF"
nuf.2000.pos.dist.subset$dist.obs.summary$type = "NUF"
nuf.2000.neg.dist.subset$dist.obs.summary$type = "NUF"
suf.1950.pos.dist$dist.obs.summary$type = "SUF"
suf.1950.neg.dist$dist.obs.summary$type = "SUF"
suf.2000.pos.dist$dist.obs.summary$type = "SUF"
suf.2000.neg.dist$dist.obs.summary$type = "SUF"
suf.2000.pos.dist.subset$dist.obs.summary$type = "SUF"
suf.2000.neg.dist.subset$dist.obs.summary$type = "SUF"

csp.1950.pos.dist$dist.obs.summary$date = "1950s"
csp.1950.neg.dist$dist.obs.summary$date = "1950s"
csp.2000.pos.dist$dist.obs.summary$date = "2000s"
csp.2000.neg.dist$dist.obs.summary$date = "2000s"
csp.2000.pos.dist.subset$dist.obs.summary$date = "2000s"
csp.2000.neg.dist.subset$dist.obs.summary$date = "2000s"
nuf.1950.pos.dist$dist.obs.summary$date = "1950s"
nuf.1950.neg.dist$dist.obs.summary$date = "1950s"
nuf.2000.pos.dist$dist.obs.summary$date = "2000s"
nuf.2000.neg.dist$dist.obs.summary$date = "2000s"
nuf.2000.pos.dist.subset$dist.obs.summary$date = "2000s"
nuf.2000.neg.dist.subset$dist.obs.summary$date = "2000s"
suf.1950.pos.dist$dist.obs.summary$date = "1950s"
suf.1950.neg.dist$dist.obs.summary$date = "1950s"
suf.2000.pos.dist$dist.obs.summary$date = "2000s"
suf.2000.neg.dist$dist.obs.summary$date = "2000s"
suf.2000.pos.dist.subset$dist.obs.summary$date = "2000s"
suf.2000.neg.dist.subset$dist.obs.summary$date = "2000s"

csp.1950.pos.dist$dist.obs.summary$data.source = "alldata"
csp.1950.neg.dist$dist.obs.summary$data.source = "alldata"
csp.2000.pos.dist$dist.obs.summary$data.source = "alldata"
csp.2000.neg.dist$dist.obs.summary$data.source = "alldata"
csp.2000.pos.dist.subset$dist.obs.summary$data.source = "subset"
csp.2000.neg.dist.subset$dist.obs.summary$data.source = "subset"
nuf.1950.pos.dist$dist.obs.summary$data.source = "alldata"
nuf.1950.neg.dist$dist.obs.summary$data.source = "alldata"
nuf.2000.pos.dist$dist.obs.summary$data.source = "alldata"
nuf.2000.neg.dist$dist.obs.summary$data.source = "alldata"
nuf.2000.pos.dist.subset$dist.obs.summary$data.source = "subset"
nuf.2000.neg.dist.subset$dist.obs.summary$data.source = "subset"
suf.1950.pos.dist$dist.obs.summary$data.source = "alldata"
suf.1950.neg.dist$dist.obs.summary$data.source = "alldata"
suf.2000.pos.dist$dist.obs.summary$data.source = "alldata"
suf.2000.neg.dist$dist.obs.summary$data.source = "alldata"
suf.2000.pos.dist.subset$dist.obs.summary$data.source = "subset"
suf.2000.neg.dist.subset$dist.obs.summary$data.source = "subset"

csp.1950.pos.dist$dist.obs.summary$posneg = "pos"
csp.1950.neg.dist$dist.obs.summary$posneg = "neg"
csp.2000.pos.dist$dist.obs.summary$posneg = "pos"
csp.2000.neg.dist$dist.obs.summary$posneg = "neg"
csp.2000.pos.dist.subset$dist.obs.summary$posneg = "pos"
csp.2000.neg.dist.subset$dist.obs.summary$posneg = "neg"
nuf.1950.pos.dist$dist.obs.summary$posneg = "pos"
nuf.1950.neg.dist$dist.obs.summary$posneg = "neg"
nuf.2000.pos.dist$dist.obs.summary$posneg = "pos"
nuf.2000.neg.dist$dist.obs.summary$posneg = "neg"
nuf.2000.pos.dist.subset$dist.obs.summary$posneg = "pos"
nuf.2000.neg.dist.subset$dist.obs.summary$posneg = "neg"
suf.1950.pos.dist$dist.obs.summary$posneg = "pos"
suf.1950.neg.dist$dist.obs.summary$posneg = "neg"
suf.2000.pos.dist$dist.obs.summary$posneg = "pos"
suf.2000.neg.dist$dist.obs.summary$posneg = "neg"
suf.2000.pos.dist.subset$dist.obs.summary$posneg = "pos"
suf.2000.neg.dist.subset$dist.obs.summary$posneg = "neg"

func.dist.summary = rbind(csp.1950.pos.dist$dist.obs.summary,
                          csp.1950.neg.dist$dist.obs.summary,
                          csp.2000.pos.dist$dist.obs.summary,
                          csp.2000.neg.dist$dist.obs.summary,
                          csp.2000.pos.dist.subset$dist.obs.summary,
                          csp.2000.neg.dist.subset$dist.obs.summary,
                          nuf.1950.pos.dist$dist.obs.summary,
                          nuf.1950.neg.dist$dist.obs.summary,
                          nuf.2000.pos.dist$dist.obs.summary,
                          nuf.2000.neg.dist$dist.obs.summary,
                          nuf.2000.pos.dist.subset$dist.obs.summary,
                          nuf.2000.neg.dist.subset$dist.obs.summary,
                          suf.1950.pos.dist$dist.obs.summary,
                          suf.1950.neg.dist$dist.obs.summary,
                          suf.2000.pos.dist$dist.obs.summary,
                          suf.2000.neg.dist$dist.obs.summary,
                          suf.2000.pos.dist.subset$dist.obs.summary,
                          suf.2000.neg.dist.subset$dist.obs.summary)
func.dist.summary.subset = func.dist.summary[-which(func.dist.summary$date == "2000s" &
                                                      func.dist.summary$data.source == "alldata"),]

## combine all mean func distance of random pair together ----
func.dist.detail = data.frame(distance = c(csp.1950.pos.dist$mean.dist.rand,
                                           csp.2000.pos.dist$mean.dist.rand,
                                           csp.2000.pos.dist.subset$mean.dist.rand,
                                           nuf.1950.pos.dist$mean.dist.rand,
                                           nuf.2000.pos.dist$mean.dist.rand,
                                           nuf.2000.pos.dist.subset$mean.dist.rand,
                                           suf.1950.pos.dist$mean.dist.rand,
                                           suf.2000.pos.dist$mean.dist.rand,
                                           suf.2000.pos.dist.subset$mean.dist.rand),
                              ty = rep(c("CSP.1950s.NA",
                                         "CSP.2000s.NA",
                                         "CSP.2000s.NA",
                                         "NUF.1950s.NA",
                                         "NUF.2000s.NA",
                                         "NUF.2000s.NA",
                                         "SUF.1950s.NA",
                                         "SUF.2000s.NA",
                                         "SUF.2000s.NA"), each = 5000), 
                              data.source = rep(rep(rep( c("alldata", "subset"), times = c(2,1)), 3), each = 5000)
)

func.dist.detail = separate(func.dist.detail, col = ty, into = c("type", "date", "Pairs"))
func.dist.detail$rand.obs = "rand"
func2 = select(func.dist.summary, dist.obs, type, date, posneg, data.source)
names(func2)[1] = "distance"
names(func2)[4] = "Pairs"
func2$rand.obs = "obs"
func2$Pairs[func2$Pairs == "pos"] = "Aggregated"
func2$Pairs[func2$Pairs == "neg"] = "Segregated"
func.dist.detail = rbind(func.dist.detail, func2)

## plot functional distance with all data ----
func.plot.alldata = filter(func.dist.detail, data.source == "alldata")
func.plot.alldata$type = factor(func.plot.alldata$type, levels = c("NUF", "CSP", "SUF"))

sp.func.diff.alldata = ggplot(data = func.plot.alldata) + geom_histogram(data= subset(func.plot.alldata, rand.obs == "rand"),
                                                             aes(x = distance, y = ..density..), fill = "lightskyblue") + 
  #   geom_density(data= subset(all.rand.dist.detail.long2, obs.rand == "rand"),aes(x=value), color = "blue") +
  geom_vline(data = subset(func.plot.alldata, rand.obs == "obs"), aes(xintercept = distance, linetype = Pairs),
             show_guide = TRUE)+xlab("Functional distance")+ylab("Count") +theme(legend.position = "top")+
  facet_grid(date~type)
ggsave(filename = "../Figs/dist.obs.rand.alldata.pdf", plot = sp.func.diff.alldata, width = 10, height = 6, units="in")

## plot functional distance with subset data ----
func.plot.subset = func.dist.detail[-which(func.dist.detail$date == "2000s" & func.dist.detail$data.source == "alldata"), ]
func.plot.subset$type = factor(func.plot.subset$type, levels = c("NUF", "CSP", "SUF"))

sp.func.diff.subset = ggplot(data = func.plot.subset) + geom_histogram(data= subset(func.plot.subset, rand.obs == "rand"),
                                                                         aes(x = distance, y = ..density..), fill = "lightskyblue") + 
  #   geom_density(data= subset(all.rand.dist.detail.long2, obs.rand == "rand"),aes(x=value), color = "blue") +
  geom_vline(data = subset(func.plot.subset, rand.obs == "obs"), aes(xintercept = distance, linetype = Pairs),
             show_guide = TRUE)+xlab("Functional distance")+ylab("Count") +theme(legend.position = "top")+
  facet_grid(date~type)
ggsave(filename = "../Figs/dist.obs.rand.subset.pdf", plot = sp.func.diff.subset, width = 10, height = 6, units="in")
