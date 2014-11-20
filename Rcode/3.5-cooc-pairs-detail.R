## combine all cooc pairs together
cooc.csp.1950.wide.sq$type = "CSP"
cooc.csp.2000.wide.sq$type = "CSP"
cooc.csp.2000.wide.sq.subset$type = "CSP"
cooc.nuf.1950.wide.sq$type = "NUF"
cooc.nuf.2000.wide.sq$type = "NUF"
cooc.nuf.2000.wide.sq.subset$type = "NUF"
cooc.suf.1950.wide.sq$type = "SUF"
cooc.suf.2000.wide.sq$type = "SUF"
cooc.suf.2000.wide.sq.subset$type = "SUF"

cooc.csp.1950.wide.sq$date = "1950s"
cooc.csp.2000.wide.sq$date = "2000s"
cooc.csp.2000.wide.sq.subset$date = "2000s"
cooc.nuf.1950.wide.sq$date = "1950s"
cooc.nuf.2000.wide.sq$date = "2000s"
cooc.nuf.2000.wide.sq.subset$date = "2000s"
cooc.suf.1950.wide.sq$date = "1950s"
cooc.suf.2000.wide.sq$date = "2000s"
cooc.suf.2000.wide.sq.subset$date = "2000s"

cooc.csp.1950.wide.sq$data.source = "alldata"
cooc.csp.2000.wide.sq$data.source = "alldata"
cooc.csp.2000.wide.sq.subset$data.source = "subset"
cooc.nuf.1950.wide.sq$data.source = "alldata"
cooc.nuf.2000.wide.sq$data.source = "alldata"
cooc.nuf.2000.wide.sq.subset$data.source = "subset"
cooc.suf.1950.wide.sq$data.source = "alldata"
cooc.suf.2000.wide.sq$data.source = "alldata"
cooc.suf.2000.wide.sq.subset$data.source = "subset"

cooc.all = rbind(cooc.csp.1950.wide.sq,
                 cooc.csp.2000.wide.sq,
                 cooc.csp.2000.wide.sq.subset,
                 cooc.nuf.1950.wide.sq,
                 cooc.nuf.2000.wide.sq,
                 cooc.nuf.2000.wide.sq.subset,
                 cooc.suf.1950.wide.sq,
                 cooc.suf.2000.wide.sq,
                 cooc.suf.2000.wide.sq.subset)
cooc.all = tbl_df(cooc.all)
str(cooc.all)


## summarize cooc pairs subset data------
cooc.all.subset = cooc.all[-which(cooc.all$date=="2000s" & cooc.all$data.source == "alldata"),]
str(cooc.all.subset)
cooc.all.subset.groupby = group_by(cooc.all.subset, type, date)
cooc.all.subset.summary = cooc.all.subset.groupby %>%
  summarize(sp.more.5q = length(unique(c(sp1_name, sp2_name))),
            all.agg = sum(SES>0), all.seg = sum(SES<0), 
            sig.agg = sum(SES>0 & p.value < 0.05),
            sig.seg = sum(SES<0 & p.value < 0.05),
            all.agg.prop = sum(SES>0)/choose(length(unique(c(sp1_name, sp2_name))), 2), 
            all.seg.prop = sum(SES<0)/choose(length(unique(c(sp1_name, sp2_name))), 2), 
            sig.agg.prop = sum(SES>0 & p.value < 0.05)/choose(length(unique(c(sp1_name, sp2_name))), 2),
            sig.seg.prop = sum(SES<0 & p.value < 0.05)/choose(length(unique(c(sp1_name, sp2_name))), 2))

cooc.all.subset.summary$sp.rich = c(sum(colSums(veg.csp.1950.wide.site) > 0),
#                 sum(colSums(veg.csp.2000.wide.site) > 0),
                sum(colSums(veg.csp.2000.wide.site.subset) > 0),
                sum(colSums(veg.nuf.1950.wide.site) > 0),
#                 sum(colSums(veg.nuf.2000.wide.site) > 0),
                sum(colSums(veg.nuf.2000.wide.site.subset) > 0),
                sum(colSums(veg.suf.1950.wide.site) > 0),
#                 sum(colSums(veg.suf.2000.wide.site) > 0),
                sum(colSums(veg.suf.2000.wide.site.subset) > 0))
cooc.all.subset.summary = select(cooc.all.subset.summary, type, date, sp.rich, sp.more.5q:sig.seg.prop)
cooc.all.subset.summary
