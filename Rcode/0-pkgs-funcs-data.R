## pkgs and functions ----
source("funs_include_geodist_envi.R")

## data----
load("xy.Rdata")
envi = data.frame(x)
row.names(envi) = paste0("site", 1:2500)
colnames(envi) = paste0("envi", 1:5)
envi$lat = runif(n = 2500, min = 40.00000, max = 45.999999)
envi$long = runif(n = 2500, min = -90.00000, max = -80.999999)
envi$site = row.names(envi) 
envi = select(envi, site, long, lat, envi1:envi5)

veg = data.frame(y)
row.names(veg) = paste0("site", 1:2500)
colnames(veg) = paste0("sp", 1:100)
rm(x, y)

## get 100 sites and 30 sp as code-test data
envi.test = envi[1:100, ]
veg.test = veg[1:100, 1:50]

# sort(colSums(veg))
# 
# ## plot
# par(mfrow=c(3,2))
# hist(envi$envi1)
# hist(envi$envi2)
# hist(envi$envi3)
# hist(envi$envi4)
# hist(envi$envi5)
# hist(rowSums(veg))
# # hist(colSums(veg))
# par(mfrow=c(1,1))
# plot(envi$envi1, rowSums(veg))
