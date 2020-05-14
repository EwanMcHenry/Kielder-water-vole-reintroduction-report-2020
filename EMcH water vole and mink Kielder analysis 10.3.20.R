
# load libraries and set WD ----
require(lubridate) # for date formatting
require(ggplot2) # for graphics
require(viridis) # colour blidn pallet
require(mgcv) # for gams
library(cowplot) # for pannel plots in ggplot
library (foreign) # for reading and writting .dbf files


setwd("D:\\Kielder water voles\\Kielder water voles")

# load data: mink raft deteciton, WV raft detections, survey locations ----

raft.monitor = read.csv("Ewan analysis\\mink.monitoring.csv")

raft.alts = read.csv("Ewan analysis\\raft.altitudes.csv")
v.rafts.dh = read.csv("Ewan analysis\\water.vole.rafts.csv") # vole raft detection history 
wv.surveys = read.csv("Ewan analysis\\water.vole.surveys.csv")

wv.colext = read.csv("Ewan analysis\\colvsext.csv") # estiamted extinctions and 

rv.points = read.csv ("Shape files\\split_points.csv")
rv.lines = read.csv ("Shape files\\split_lines.csv")

soils = read.dbf("Shape files\\soils_clip\\soils_clip2.dbf", as.is = T)

wv.raft.detcts = read.dbf ("Shape files\\water.detections.by.year_v1.01.1.dbf", as.is = T)

wv.d.hists = read.csv ("Ewan analysis\\patch survey detection histories.csv")

# data curation ----

raft.monitor = raft.monitor[, names(raft.monitor) %in% c("Raft", "Eastings", "Northings", "Tracks.present.", "Droppings.present.", "Tracks", 
                                                         "Droppings", "Mink.Present.", "Water.Vole.Present." , "Does.the.raft.require.maintenance", 
                                                         "Maintenance.required", "Date", "Grouped.Date")]

raft.monitor$Date = as.Date(raft.monitor$Date)
raft.monitor$Grouped.Date = as.Date(raft.monitor$Grouped.Date)
temp = as.numeric(raft.monitor$Does.the.raft.require.maintenance) *-1 +2
raft.monitor$Does.the.raft.require.maintenance = temp

# give missing location data ----
raft.monitor$Eastings[raft.monitor$Raft == "Raft 11" ] = 364900
raft.monitor$Northings[raft.monitor$Raft == "Raft 11" ] = 583600

raft.monitor$Eastings[raft.monitor$Raft == "Raft 12" ] = 368640
raft.monitor$Northings[raft.monitor$Raft == "Raft 12" ] = 585490

raft.monitor$Eastings[raft.monitor$Raft == "Raft 26" ] = 374080
raft.monitor$Northings[raft.monitor$Raft == "Raft 26" ] = 594780

raft.monitor$Eastings[raft.monitor$Raft == "Raft 22" ] = 363106
raft.monitor$Northings[raft.monitor$Raft == "Raft 22" ] = 592522

# raft.monitor$Eastings[raft.monitor$Raft == "Raft 6" ] = 
# raft.monitor$Northings[raft.monitor$Raft == "Raft 6" ] = 

# remove the incorect water vole raft detection
raft.monitor$Water.Vole.Present.[18] = 0


# check internal consitency of raft locations
raft.loc.table = table (paste(raft.monitor$Eastings, raft.monitor$Northings), raft.monitor$Raft)
raft.loc.table.01= which(raft.loc.table >0, arr.ind = T)
# no rafts with different locations, some with same location = all location 0, remove them from spatial analysis
as.numeric(raft.loc.table.01[,1])[which(duplicated(as.numeric(raft.loc.table.01[,1])) == T )]
raft.loc.table[,which(duplicated(raft.loc.table.01[,1]) == T)]

 ## # # # ## # # # # ## # # # #
# water vole surveys

wv.surveys = wv.surveys[, names(wv.surveys) %in% c("Site.Name", "Eastings", "Northings" , "Date" , "Length" , "Sightings" , "Latrines",
                                                   "Watervole.present.", "Feeding.Remains") ]
wv.surveys$Date = as.Date(wv.surveys$Date)


# add more date and variable options ----
raft.monitor$year = year(raft.monitor$Date) 
raft.monitor$month = month(raft.monitor$Date)
raft.monitor$julian = yday(raft.monitor$Date)
raft.monitor$month.day = format (raft.monitor$Date, format="%B %d")

raft.monitor$grouped.julian = yday(raft.monitor$Grouped.Date) 

raft.monitor$days.since.start =  as.numeric(raft.monitor$Date -  min(raft.monitor$Date))

wv.surveys$year = year(wv.surveys$Date) 
wv.surveys$month = month(wv.surveys$Date)
wv.surveys$julian = yday(wv.surveys$Date)
wv.surveys$month.day = format (wv.surveys$Date, format="%B %d")

wv.surveys$days.since.start =  as.numeric(wv.surveys$Date -  min(wv.surveys$Date))



dum.raft = (substring(raft.monitor$Raft, 6))
dum.raft [dum.raft == "27W"] = 27
dum.raft [dum.raft == "28W"] = 28
dum.raft [dum.raft == "28W"] = 28

raft.monitor$Raft.num = (substring(raft.monitor$Raft, 6))
raft.monitor$Raft.num[ raft.monitor$Raft.num =="27SW "] = "27"
raft.monitor$Raft.num[ raft.monitor$Raft.num =="28SW"] = "28"


# ploting variables ----
julian.breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 365)
mid.month.breaks = julian.breaks[1:(length(julian.breaks)-1)]+14
dist.to.nearest.wv.breaks = seq(0,20, 2)
n.checked.breaks = seq(0,35, 5)
long.term.breaks = c(0, 366, 731, 1096)
long.term.lables = 2016:2019

# create dispersal season variable ----
raft.monitor$dispersal.season = (raft.monitor$julian >= 213 & raft.monitor$julian <= max(raft.monitor$julian) )*1 # between 1st august and end of year
# 213 = 1st august 


#- - - ## - - - - ---- ---- ---- ----#######
#   distance calculations ----


# version of data removing observations without coordinates and rename raft 11
space.rafts = raft.monitor[raft.monitor$Eastings >0,]
space.rafts$Raft.num [space.rafts$Raft.num == "11a"] = "11"
space.rafts$Raft.num = as.numeric(space.rafts$Raft.num)
# 
raft.survey.locs = rbind (space.rafts[,  names(space.rafts) %in% c("Eastings", "Northings" )] ,
                          wv.surveys[,  names(wv.surveys) %in% c("Eastings", "Northings" )] 
                          )
raft.survey.dists =as.matrix (dist(raft.survey.locs))

space.rafts$dists = list(NA)
space.rafts$dist.dates = list(NA)
space.rafts$dist.wv.present = list(NA)

for (i in 1:dim(space.rafts)[1]){
  space.rafts$dists[[i]] = raft.survey.dists[i,] 
  space.rafts$dist.dates [[i]] = c(space.rafts$Date, wv.surveys$Date )
  space.rafts$dist.wv.present [[i]] = c(space.rafts$Water.Vole.Present., wv.surveys$Watervole.present.)
}


# ---- ---- ---- ---- ---- ###
# how often is raft maintainence required ----
table(1-(as.numeric (raft.monitor$Does.the.raft.require.maintenance)-1))


maintanence.glm = glm(raft.monitor$Does.the.raft.require.maintenance ~1, family = "binomial")
estimates = coef(summary(maintanence.glm))
plogis( estimates[, "Estimate"])
plogis( estimates[, "Estimate"] + 1.96*(estimates[, "Std. Error"] ))
plogis( estimates[, "Estimate"] - 1.96*(estimates[, "Std. Error"] ))


# where do malfunctions happen?
hist(table(1-(as.numeric (raft.monitor$Does.the.raft.require.maintenance)-1), raft.monitor$Raft)[2,]/
  table(1-(as.numeric (raft.monitor$Does.the.raft.require.maintenance)-1), raft.monitor$Raft)[1,])

# how often are mink detected ----
table(as.numeric (raft.monitor$Mink.Present.))

detect.glm.00 = glm(raft.monitor$Mink.Present. ~1, family = "binomial")
estimates = coef(summary(detect.glm.00))
plogis( estimates[, "Estimate"])
plogis( estimates[, "Estimate"] + 1.96*(estimates[, "Std. Error"] ))
plogis( estimates[, "Estimate"] - 1.96*(estimates[, "Std. Error"] ))

# how many rafts are checked each 2 week window ----
sort(raft.monitor$Grouped.Date)[1]
sort(raft.monitor$Grouped.Date)[length(raft.monitor$Grouped.Date)]
n.nochecks = 86- length(unique(raft.monitor$Grouped.Date))

length(unique(raft.monitor$Grouped.Date))/86



n.checked = c(as.numeric(table(raft.monitor$Grouped.Date[!duplicated(paste(raft.monitor$Grouped.Date, raft.monitor$Raft))]))
                , rep (0,n.nochecks))

n.checked.when.checked = as.numeric(table(raft.monitor$Grouped.Date[!duplicated(paste(raft.monitor$Grouped.Date, raft.monitor$Raft))]))
mean(n.checked.when.checked)
sd(n.checked.when.checked)

rafts.checked = ggplot()+
  geom_histogram(aes(x = n.checked.when.checked, ..density..),
                   breaks = n.checked.breaks, 
                 fill = "grey", color = "grey", alpha = 1)+  
  geom_density(aes(x = n.checked.when.checked, ..density..),
               fill = "grey", color = "black", alpha = 0.4)  +
    scale_x_continuous( name = "Number of rafts checked") +
  scale_y_continuous(lim = c(0, 0.07), name = "Freq. density") +
  theme_minimal()

# look at raft 30 detections ----
raft.monitor[raft.monitor$Raft == "Raft 30",]

# seasonality of detections ----
grouped.raft.monitor = raft.monitor[!duplicated(paste(raft.monitor$Grouped.Date, raft.monitor$Raft)),]

julian.densplot.m0 = 
  ggplot(raft.monitor)+
  geom_histogram(data = subset(grouped.raft.monitor, Mink.Present. ==0),
                 aes(x = julian, y =  ..density..),
                 breaks = julian.breaks, 
                 fill = "grey", color = "grey", alpha = 1)+  
  # geom_density(data=subset(grouped.raft.monitor,Mink.Present. == 0),
  #              aes(x = julian, ..density..),
  #              fill = "grey", color = "black", alpha = 0.4)  +

  scale_x_continuous(lim = c(0, max(julian.breaks) ),
                     labels = month.abb[seq(1,12,2)], breaks = mid.month.breaks[seq(1,12,2)], name = "") +
  scale_y_continuous(lim = c(0, 0.008), name = "Freq. density") +
  theme_minimal()

julian.densplot.m1 = ggplot(raft.monitor)+
    # geom_density(data=subset(raft.monitor,Mink.Present. == 1),
    #            aes(x = julian, ..density..),
    #            fill = "red", color = "red", alpha = 0.2)  +  
  geom_histogram(data = subset(raft.monitor, Mink.Present. ==1), 
                 aes(x = julian, y =  ..density..),
                 breaks = julian.breaks, 
                 fill = "red", color = "black", alpha = 0.6)+  
  scale_x_continuous(lim = c(0, max(julian.breaks) ),
                     labels = month.abb[seq(1,12,2)], breaks = mid.month.breaks[seq(1,12,2)], name = "") +
  scale_y_continuous(lim = c(0, 0.008), name = "") +
  theme_minimal()

plot_grid(julian.densplot.m0, julian.densplot.m1, labels = "AUTO")

# time since start for monitoring frequency ----

  ggplot(grouped.raft.monitor)+
  geom_histogram(data =grouped.raft.monitor,
                 aes(x = days.since.start, y =  ..density..),
                 fill = "grey", color = "grey", alpha = 1)+  
  geom_density(data=  grouped.raft.monitor,
               aes(x = days.since.start, ..density..),
               fill = "grey", color = "black", alpha = 0.4)  +
  scale_x_continuous(name = "Year", breaks = long.term.breaks, labels = long.term.lables) +
  scale_y_continuous( name = "Freq. density") +
  theme_minimal()






# by year + season tables ----
table(raft.monitor$dispersal.season, raft.monitor$year )
table(raft.monitor$dispersal.season[raft.monitor$Mink.Present.==1], raft.monitor$year[raft.monitor$Mink.Present.==1] )

table(raft.monitor$dispersal.season[raft.monitor$Mink.Present.==1], raft.monitor$year[raft.monitor$Mink.Present.==1] )[,1:2]/table(raft.monitor$dispersal.season, raft.monitor$year )[,1:2]

# distance from each mink detection to closest water vole deteciton in last year ----

for (i in 1:dim(space.rafts)[1]){
  which.wv.in.1yr = which (space.rafts$dist.wv.present[[i]] ==1 & 
                                   space.rafts$dist.dates[[i]]- space.rafts$Date[i] <= 0 &
                                   space.rafts$dist.dates[[i]]- space.rafts$Date[i] >= -365 &
                                   space.rafts$dists[[i]] > 0)
  space.rafts$dist.to.1yr.nearest.wv[i] = min( space.rafts$dists[[i]][which.wv.in.1yr])
}
space.rafts$dist.to.1yr.nearest.wv[space.rafts$dist.to.1yr.nearest.wv == "Inf"] = NA

interested.rafts = !is.na(space.rafts$dist.to.1yr.nearest.wv) & space.rafts$dist.to.1yr.nearest.wv < 20000
space.rafts$dist.to.1yr.nearest.wv = space.rafts$dist.to.1yr.nearest.wv/1000

dist.to.1yr.nearest.wv.densplot.m0 = 
  ggplot(raft.monitor)+
  geom_histogram(data = subset(space.rafts, Mink.Present. ==0),
                 aes(x = dist.to.1yr.nearest.wv, y =  ..density..),
                 breaks = dist.to.nearest.wv.breaks, 
                 fill = "grey", color = "grey", alpha = 1)+  
  geom_density(data=subset(space.rafts,Mink.Present. == 0),
               aes(x = dist.to.1yr.nearest.wv, ..density..),
               fill = "grey", color = "black", alpha = 0.4)  +
    scale_x_continuous(name = "Distance to nearest water vole sign (km)") +
  scale_y_continuous(lim = c(0,0.3) , name = "Freq. density") +
  theme_minimal()

dist.to.1yr.nearest.wv.densplot.m1 = ggplot(space.rafts)+
  geom_density(data=subset(space.rafts,Mink.Present. == 1),
               aes(x = dist.to.1yr.nearest.wv, ..density..),
               fill = "red", color = "red", alpha = 0.2)  +  
  geom_histogram(data = subset(space.rafts, Mink.Present. ==1), 
                 aes(x = dist.to.1yr.nearest.wv, y =  ..density..),
                 breaks = dist.to.nearest.wv.breaks, 
                 fill = "red", color = "black", alpha = 0.6)+  
  scale_x_continuous( name = "") +
  scale_y_continuous(lim = c(0, 0.3), name = "") +
  theme_minimal()

plot_grid(dist.to.1yr.nearest.wv.densplot.m0, dist.to.1yr.nearest.wv.densplot.m1, labels = "AUTO")

glm.dist.to.1yr.nearest.wv = glm(Mink.Present. ~ dist.to.1yr.nearest.wv, family = "binomial", data = space.rafts)
summary(glm.dist.to.1yr.nearest.wv)

# distance from each mink detection to closest water vole deteciton in past ----

for (i in 1:dim(space.rafts)[1]){
  which.wv.in.past = which (space.rafts$dist.wv.present[[i]] ==1 & 
                               space.rafts$dist.dates[[i]]- space.rafts$Date[i] <= 0 &
                               space.rafts$dists[[i]] > 0)
  space.rafts$dist.to.past.nearest.wv[i] = min( space.rafts$dists[[i]][which.wv.in.past])
}
space.rafts$dist.to.past.nearest.wv[space.rafts$dist.to.past.nearest.wv == "Inf"] = NA

interested.rafts = !is.na(space.rafts$dist.to.past.nearest.wv) & space.rafts$dist.to.past.nearest.wv < 20000
space.rafts$dist.to.past.nearest.wv = space.rafts$dist.to.past.nearest.wv/1000

dist.to.past.nearest.wv.densplot.m0 = 
  ggplot(raft.monitor)+
  geom_histogram(data = subset(space.rafts, Mink.Present. ==0),
                 aes(x = dist.to.past.nearest.wv, y =  ..density..),
                 breaks = dist.to.nearest.wv.breaks, 
                 fill = "grey", color = "grey", alpha = 1)+  
  geom_density(data=subset(space.rafts,Mink.Present. == 0),
               aes(x = dist.to.past.nearest.wv, ..density..),
               fill = "grey", color = "black", alpha = 0.4)  +
  scale_x_continuous(name = "Distance to nearest water vole sign") +
  scale_y_continuous(lim = c(0,0.3) , name = "Freq. density") +
  theme_minimal()

dist.to.past.nearest.wv.densplot.m1 = ggplot(space.rafts)+
  geom_density(data=subset(space.rafts,Mink.Present. == 1),
               aes(x = dist.to.past.nearest.wv, ..density..),
               fill = "red", color = "red", alpha = 0.2)  +  
  geom_histogram(data = subset(space.rafts, Mink.Present. ==1), 
                 aes(x = dist.to.past.nearest.wv, y =  ..density..),
                 breaks = dist.to.nearest.wv.breaks, 
                 fill = "red", color = "black", alpha = 0.6)+  
  scale_x_continuous() +
  scale_y_continuous(lim = c(0, 0.3), name = "") +
  theme_minimal()

plot_grid(dist.to.past.nearest.wv.densplot.m0, dist.to.past.nearest.wv.densplot.m1, labels = "AUTO")

glm.dist.to.past.nearest.wv = glm(Mink.Present. ~ dist.to.past.nearest.wv, family = "binomial", data = space.rafts)
summary(glm.dist.to.past.nearest.wv)

# distance from each mink detection to closest water vole deteciton in last 6 months ----

for (i in 1:dim(space.rafts)[1]){
  which.wv.in.6month = which (space.rafts$dist.wv.present[[i]] ==1 & 
                             space.rafts$dist.dates[[i]]- space.rafts$Date[i] <= 0 &
                             space.rafts$dist.dates[[i]]- space.rafts$Date[i] >= -183 &
                             space.rafts$dists[[i]] > 0)
  space.rafts$dist.to.6month.nearest.wv[i] = min( space.rafts$dists[[i]][which.wv.in.6month])
}
space.rafts$dist.to.6month.nearest.wv[space.rafts$dist.to.6month.nearest.wv == "Inf"] = NA

interested.rafts = !is.na(space.rafts$dist.to.6month.nearest.wv) & space.rafts$dist.to.6month.nearest.wv < 20000
space.rafts$dist.to.6month.nearest.wv = space.rafts$dist.to.6month.nearest.wv/1000

dist.to.6month.nearest.wv.densplot.m0 = 
  ggplot(raft.monitor)+
  geom_histogram(data = subset(space.rafts, Mink.Present. ==0),
                 aes(x = dist.to.6month.nearest.wv, y =  ..density..),
                 breaks = dist.to.nearest.wv.breaks, 
                 fill = "grey", color = "grey", alpha = 1)+  
  geom_density(data=subset(space.rafts,Mink.Present. == 0),
               aes(x = dist.to.6month.nearest.wv, ..density..),
               fill = "grey", color = "black", alpha = 0.4)  +
  scale_x_continuous(name = "Distance to nearest water vole sign") +
  scale_y_continuous(lim = c(0,0.3) , name = "Freq. density") +
  theme_minimal()

dist.to.6month.nearest.wv.densplot.m1 = ggplot(space.rafts)+
  geom_density(data=subset(space.rafts,Mink.Present. == 1),
               aes(x = dist.to.6month.nearest.wv, ..density..),
               fill = "red", color = "red", alpha = 0.2)  +  
  geom_histogram(data = subset(space.rafts, Mink.Present. ==1), 
                 aes(x = dist.to.6month.nearest.wv, y =  ..density..),
                 breaks = dist.to.nearest.wv.breaks, 
                 fill = "red", color = "black", alpha = 0.6)+  
  scale_x_continuous() +
  scale_y_continuous(lim = c(0, 0.3), name = "") +
  theme_minimal()

plot_grid(dist.to.6month.nearest.wv.densplot.m0, dist.to.6month.nearest.wv.densplot.m1, labels = "AUTO")

glm.dist.to.6month.nearest.wv = glm(Mink.Present. ~ dist.to.6month.nearest.wv, family = "binomial", data = space.rafts)
summary(glm.dist.to.6month.nearest.wv)

# distance from each mink detection to closest water vole deteciton in last 2 months ----

for (i in 1:dim(space.rafts)[1]){
  which.wv.in.2month = which (space.rafts$dist.wv.present[[i]] ==1 & 
                                space.rafts$dist.dates[[i]]- space.rafts$Date[i] <= 0 &
                                space.rafts$dist.dates[[i]]- space.rafts$Date[i] >= -61 &
                                space.rafts$dists[[i]] > 0)
  space.rafts$dist.to.2month.nearest.wv[i] = min( space.rafts$dists[[i]][which.wv.in.2month])
}
space.rafts$dist.to.2month.nearest.wv[space.rafts$dist.to.2month.nearest.wv == "Inf"] = NA

interested.rafts = !is.na(space.rafts$dist.to.2month.nearest.wv) & space.rafts$dist.to.2month.nearest.wv < 20000
space.rafts$dist.to.2month.nearest.wv = space.rafts$dist.to.2month.nearest.wv/1000

dist.to.2month.nearest.wv.densplot.m0 = 
  ggplot(raft.monitor)+
  geom_histogram(data = subset(space.rafts, Mink.Present. ==0),
                 aes(x = dist.to.2month.nearest.wv, y =  ..density..),
                 breaks = dist.to.nearest.wv.breaks, 
                 fill = "grey", color = "grey", alpha = 1)+  
  geom_density(data=subset(space.rafts,Mink.Present. == 0),
               aes(x = dist.to.2month.nearest.wv, ..density..),
               fill = "grey", color = "black", alpha = 0.4)  +
  scale_x_continuous(name = "") +
  scale_y_continuous(lim = c(0,0.3) , name = "Freq. density") +
  theme_minimal()

dist.to.2month.nearest.wv.densplot.m1 = ggplot(space.rafts)+
  geom_density(data=subset(space.rafts,Mink.Present. == 1),
               aes(x = dist.to.2month.nearest.wv, ..density..),
               fill = "red", color = "red", alpha = 0.2)  +  
  geom_histogram(data = subset(space.rafts, Mink.Present. ==1), 
                 aes(x = dist.to.2month.nearest.wv, y =  ..density..),
                 breaks = dist.to.nearest.wv.breaks, 
                 fill = "red", color = "black", alpha = 0.6)+  
  scale_x_continuous(name = "") +
  scale_y_continuous(lim = c(0, 0.3), name = "") +
  theme_minimal()

plot_grid(dist.to.2month.nearest.wv.densplot.m0, dist.to.2month.nearest.wv.densplot.m1, labels = "AUTO")

glm.dist.to.2month.nearest.wv = glm(Mink.Present. ~ dist.to.2month.nearest.wv, family = "binomial", data = space.rafts)
summary(glm.dist.to.2month.nearest.wv)

# grid plot of dist to nearest ----
plot_grid(dist.to.2month.nearest.wv.densplot.m0, dist.to.2month.nearest.wv.densplot.m1, 
          dist.to.1yr.nearest.wv.densplot.m0, dist.to.1yr.nearest.wv.densplot.m1, labels = "AUTO")
summary(glm.dist.to.2month.nearest.wv)
summary(glm.dist.to.1yr.nearest.wv)


# Ccolonisation extinction frequencies.




#====#==== # ====# ====#### #### ####
# shapefile work ----

# editign data file for soils polygon
# note that I copied the .dbf file and it is this new one tht I load, manipulate and overwrite the original with
# so that the same operations arent applied to a file multiple times etc

#ANYTHING THT CONTAINS "Peat" ----
soils$peat = F
soils$peat[grep("Peat", soils$full_soil1, value = F)] = T



# ==== # ==== #==== # ==== #==== # ==== #==== # ==== #==== # ==== #==== # ==== #==== # ==== #
# ==== # ==== #==== # ==== #==== # ==== #==== # ==== #==== # ==== #==== # ==== #==== # ==== #
# river sections attributes ----


# parametarising conectivity model ----
juv.per.patch = 3.25
pop.level.effective.dispersal = 0.015
alpha = 0.33
log(2)/alpha

# adding location, striaght line distance and mean altitude drop ----
for(i in 1:length(rv.lines$ID)){
  i.interest.locations = list(NA) # create new list to fill interesting locatiosn in for each year for each
  
  alts2m = rv.points$X2M.dtm[rv.points$ID == rv.lines$ID[i]]
  # taking the negative of alt frop as high numebrs are "bad" for voles
  rv.lines$mean.Alt.drop[i] = - abs(alts2m[1]- alts2m[2])/rv.lines$length[i]
  
  locs = cbind(rv.points$Eastings[rv.points$ID == rv.lines$ID[i]], rv.points$Northings[rv.points$ID == rv.lines$ID[i]])
  rv.lines$straight.line[i] = dist(locs)
  rv.lines$mid.x[i] = mean(rv.points$Eastings[rv.points$ID == rv.lines$ID[i]])
  rv.lines$mid.y[i] = mean(rv.points$Northings[rv.points$ID == rv.lines$ID[i]])
}

rv.lines$drop.per100m = (rv.lines$mean.Alt.drop/rv.lines$length)*100
# make rugosity - ratio betwwen river segment length and straight line length
rv.lines$rugosity = (rv.lines$length)/(rv.lines$straight.line) # taken as a log ratio to better tese out the medium quality habitat
# standardising rugosity and altitude drop
rv.lines$stnd_rugosity =  (rv.lines$rugosity- mean(rv.lines$rugosity))/sd(rv.lines$rugosity)
rv.lines$stnd_altdrop =  (rv.lines$drop.per100m- mean(rv.lines$drop.per100m))/sd(rv.lines$drop.per100m)
# not correlated
cor(rv.lines$stnd_rugosity, rv.lines$stnd_altdrop)
cor(rv.lines$rugosity, rv.lines$drop.per100m)
# create combined metric to show 
rv.lines$rug_alt = ((rv.lines$stnd_rugosity + rv.lines$stnd_altdrop)- mean(rv.lines$stnd_rugosity + rv.lines$stnd_altdrop))/
  sd(rv.lines$stnd_rugosity + rv.lines$stnd_altdrop)
#  plot: all high altitude drop is associated with low rigosity ----
ggplot(rv.lines)+
  geom_point(aes(x =drop.per100m, y = rugosity )) +
  scale_x_continuous(name = "Altitude change (m)") +
  scale_y_continuous(name = "Rugosity") +
  theme_minimal()
# flat.rugosity = rugosity, but with a 0 when altitude drop is >0.5
rv.lines$flat.rug = rv.lines$rugosity
rv.lines$flat.rug[rv.lines$drop.per100m < (-0.5)] = 0

rv.lines$flat.rug [ rv.lines$flat.rug < median(rv.lines$rugosity[rv.lines$drop.per100m >-0.5])]= 0

rv.lines$rug.order =0
rv.lines$rug.order[rv.lines$flat.rug>0] = rank(rv.lines$flat.rug[rv.lines$flat.rug>0])
    
# percentile where rug starts to get good
min(rv.lines$rug.order[rv.lines$flat.rug>2])/max(rv.lines$rug.order)
#id number of top percentile
max(rv.lines$rug.order)*99


                    
# adding vole connectivity in differnent years ----

rv.conect = rv.lines # create this to add connectivity values to, just done because I want to save the origianl rv.lines to csv, the lists in the conectiv df might make that

rv.conect$vole.dists = list(NA)
rv.conect$conectivity = list(NA)
for (ii in 1: length(unique(wv.surveys$year))){
  rv.conect$vole.dists[[ii]] = list(NA)
  rv.conect$conectivity[[ii]] =rep(NA, times = length(rv.lines$ID))
  rv.conect$col.prob[[ii]]= rep(NA, times = length(rv.lines$ID))
  for(i in 1:length(rv.lines$ID)){
    # for each year of interest 
    # calculated all the distances to all occupied water vole surveys in each year
    # store in a vector of distances within list of years within the element of the river section in the river sections dataframe 
    i.interest.locations[[ii]] = cbind(c(rv.lines$mid.x[i], wv.surveys$Eastings[wv.surveys$Watervole.present.==1 & wv.surveys$year==unique(wv.surveys$year)[ii] ]),
                                       c(rv.lines$mid.y[i], wv.surveys$Northings[wv.surveys$Watervole.present.==1 & wv.surveys$year==unique(wv.surveys$year)[ii] ]))
    rv.conect$vole.dists[[ii]][[i]] = as.matrix(dist(i.interest.locations[[ii]]))[2:dim(i.interest.locations[[ii]])[1],1]/1000
    
    rv.conect$conectivity[[ii]][i] = pop.level.effective.dispersal * sum(juv.per.patch* exp((-alpha) *rv.conect$vole.dists[[ii]][[i]]))
    rv.conect$col.prob[[ii]][i] = (1 - exp(- rv.conect$conectivity[[ii]][i]  ) )
    
  }
  print(ii)
}

rv.lines$cnct.2017 = rv.conect$conectivity[[1]]
rv.lines$cnct.2018 = rv.conect$conectivity[[2]]
rv.lines$cnct.2019 = rv.conect$conectivity[[3]]

rv.lines$col.p.2017 = rv.conect$col.prob[[1]]
rv.lines$col.p.2018 = rv.conect$col.prob[[2]]
rv.lines$col.p.2019 = rv.conect$col.prob[[3]]

ggplot()+
  geom_point(data = rv.lines, aes(x = mid.x, y = mid.y, col = col.p.2019)) +
  geom_point(aes(x = wv.surveys$Eastings[wv.surveys$Watervole.present.==1 & wv.surveys$year==2019 ], 
                 y = wv.surveys$Northings[wv.surveys$Watervole.present.==1 & wv.surveys$year==2019 ])) +
  geom_point(aes(x = rv.lines$mid.x[1842],
                 y = rv.lines$mid.y[1842]), col = "red")

    

# mink raft coverage ----
for(i in 1:length(rv.lines$ID)){
    # mink raft coverage
  #locaitons - this river section and all mink rafts
    i.interest.locations = cbind(c(rv.lines$mid.x[i], space.rafts$Eastings),
                                       c(rv.lines$mid.y[i], space.rafts$Northings))
  # distances to all mink rafts  
    raft.dists = c(as.matrix(dist(i.interest.locations))[2:dim(i.interest.locations)[1],1]/1000)
   # number of 2 week periods of raft checking within 2 km and....
      #in settled season of 2019
     rv.lines$settled2019.ncov[i] = length(unique(space.rafts$Grouped.Date[raft.dists<2 &
                                                                            space.rafts$dispersal.season ==0 &
                                                                            space.rafts$year == 2019
                                                                          ] ) )
      # in 2019
    rv.lines$y2019.ncov[i] = length(unique(space.rafts$Grouped.Date[raft.dists<2 &
                                                                      space.rafts$year == 2019
                                                                    ] ) )
    # in settled season of 2018
    rv.lines$settled2018.ncov[i] = length(unique(space.rafts$Grouped.Date[raft.dists<2 &
                                                                            space.rafts$dispersal.season ==0 &
                                                                            space.rafts$year == 2018
                                                                          ] ) )
    # in dispersal season of 2018
    rv.lines$dispers2018.ncov[i] = length(unique(space.rafts$Grouped.Date[raft.dists<2 &
                                                                            space.rafts$dispersal.season ==1 &
                                                                            space.rafts$year == 2018
                                                                          ] ) )
      # in 2018
    rv.lines$y2018.ncov[i] = length(unique(space.rafts$Grouped.Date[raft.dists<2 &
                                                                      space.rafts$year == 2018
                                                                    ] ) )
      # in 2017
    rv.lines$y2017.ncov[i] = length(unique(space.rafts$Grouped.Date[raft.dists<2 &
                                                                      space.rafts$year == 2017
                                                                    ] ) )
      # in 2016
    rv.lines$y2016.ncov[i] = length(unique(space.rafts$Grouped.Date[raft.dists<2 &
                                                                      space.rafts$year == 2016
                                                                    ] ) )
}

rv.lines$prp.19.set.cov = rv.lines$settled2019.ncov/12

rv.lines$prp.18.disp.cov = rv.lines$dispers2018.ncov/10
rv.lines$prp.18.cov = rv.lines$y2018.ncov/24
rv.lines$prp.18.set.cov = rv.lines$settled2018.ncov/14
rv.lines$prp.17.cov = rv.lines$y2017.ncov/24
rv.lines$prp.16.cov = rv.lines$y2016.ncov/24





#calculating areas outside river coverage ----
# rafts 30,20,36 and 18

outside.rafts = data.frame (year = space.rafts$year[space.rafts$Raft %in% c("Raft 30", "Raft 20", "Raft 36", "Raft 18")],
                           raft = droplevels(space.rafts$Raft[space.rafts$Raft %in% c("Raft 30", "Raft 20", "Raft 36", "Raft 18")]),
                           disp = space.rafts$dispersal.season[space.rafts$Raft %in% c("Raft 30", "Raft 20", "Raft 36", "Raft 18")],
                           grou = space.rafts$Grouped.Date[space.rafts$Raft %in% c("Raft 30", "Raft 20", "Raft 36", "Raft 18")]
)

outside.rafts= outside.rafts[!duplicated(paste(outside.rafts$grou, outside.rafts$raft)),]

table(outside.rafts$disp, outside.rafts$raft, outside.rafts$year )
table(outside.rafts$raft, outside.rafts$year )/24

# mink caught each trap, each season
caught.rafts = data.frame(year = space.rafts$year[space.rafts$Mink.Present.==1],
                          raft = droplevels(space.rafts$Raft[space.rafts$Mink.Present.==1]),
                          disp = space.rafts$dispersal.season[space.rafts$Mink.Present.==1]
                          )

table(caught.rafts$raft ,  caught.rafts$disp, caught.rafts$year)


# make data for water vole detections shapefile ----

for ( i in 1: dim (wv.raft.detcts)[1]){
  wv.raft.detcts$vole2017 [i] =  max(space.rafts$Water.Vole.Present.[space.rafts$Raft == wv.raft.detcts$Row.Labels[i] & 
                                                                       space.rafts$year == 2017])
  wv.raft.detcts$vole2018 [i] =  max(space.rafts$Water.Vole.Present.[space.rafts$Raft == wv.raft.detcts$Row.Labels[i] & 
                                                                       space.rafts$year == 2018])
  wv.raft.detcts$vole2019 [i] =  max(space.rafts$Water.Vole.Present.[space.rafts$Raft == wv.raft.detcts$Row.Labels[i] & 
                                                                       space.rafts$year == 2019])
}
wv.raft.detcts$vole2017[wv.raft.detcts$vole2017 == -Inf] = NA
wv.raft.detcts$vole2018[wv.raft.detcts$vole2018 == -Inf] = NA
wv.raft.detcts$vole2019[wv.raft.detcts$vole2019 == -Inf] = NA

# water vole survey stats ----

sum(wv.surveys$year == 2017,  na.rm = T) - sum(is.na(wv.surveys$Length[wv.surveys$year == 2017]))
sum(wv.surveys$Length[wv.surveys$year == 2017], na.rm = T)
median(wv.surveys$Length[wv.surveys$year == 2017], na.rm = T)
range(wv.surveys$Length[wv.surveys$year == 2017], na.rm = T)
sum(!is.na(wv.d.hists$X2017))
sum(wv.d.hists$X2017, na.rm = T)

length(unique(space.rafts$Raft[space.rafts$year == 2017 &
                                 space.rafts$Water.Vole.Present. == 1]))
length(unique(space.rafts$Raft[space.rafts$year == 2017]))


sum(wv.surveys$year == 2018,  na.rm = T) - sum(is.na(wv.surveys$Length[wv.surveys$year == 2018]))
wv.surveys[wv.surveys$year == 2018,]
median(wv.surveys$Length[wv.surveys$year == 2018], na.rm = T)
range(wv.surveys$Length[wv.surveys$year == 2018], na.rm = T)
sum(wv.surveys$Length[wv.surveys$year == 2018], na.rm = T)
sum(!is.na(wv.d.hists$X2018))
sum(wv.d.hists$X2018, na.rm = T)
wv.surveys$Date[wv.surveys$year == 2018]

length(unique(space.rafts$Raft[space.rafts$year == 2018 &
                                 space.rafts$Water.Vole.Present. == 1]))
length(unique(space.rafts$Raft[space.rafts$year == 2018]))

#wv.surveys[order(wv.surveys$Date),c(1,4:8)]

singl.surv.2019 = wv.surveys [wv.surveys$year == 2019 &
                                !duplicated(paste(wv.surveys$Date, wv.surveys$Site.Name)) ]
sum(wv.surveys$year == 2019,  na.rm = T) - sum(is.na(wv.surveys$Length[wv.surveys$year == 2019])) # survey number
wv.surveys[wv.surveys$year == 2019,]
median(wv.surveys$Length[wv.surveys$year == 2019], na.rm = T) # survey length median
range(wv.surveys$Length[wv.surveys$year == 2019], na.rm = T)
sum(wv.surveys$Length[wv.surveys$year == 2019], na.rm = T)

sum(!is.na(wv.d.hists$X2019)) # survey patch sum
sum(wv.d.hists$X2019, na.rm = T) # survey patch occupied
wv.surveys$Date[wv.surveys$year == 2019]

length(unique(space.rafts$Raft[space.rafts$year == 2019 &
                                 space.rafts$Water.Vole.Present. == 1]))
length((space.rafts$Raft[space.rafts$year == 2019 &
                                 space.rafts$Water.Vole.Present. == 1]))
length(unique(space.rafts$Raft[space.rafts$year == 2019]))
max(space.rafts$Date)


# writting files ----

write.dbf(soils, "Shape files\\soils_clip\\soils_clip.dbf", factor2char = TRUE, max_nchar = 254)


write.csv(rv.lines,"Shape files\\split_lines_added metrics.csv" )

write.dbf(rv.lines,"Shape files\\split_lines.dbf" )

write.dbf(wv.raft.detcts,"Shape files\\water.detections.by.year_v1.01.dbf" )



