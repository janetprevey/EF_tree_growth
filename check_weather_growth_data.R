

### Cheack and clean SSMT weather data - JUly 2019

#select working directory

setwd("C:/Users/jprevey/Desktop/NEFI_2019")

#load packages


# load data 


all_precip <- readRDS(file ="fin.cleaned.ppt.rds")

all_temps <- readRDS(file ="fin.cleaned.temps.rds")

### aggregate to daily temps

buckhorn_temps <- subset(all_temps, SITE_NAME=="BUCKHORN2")

buckhorn_temps <- separate(buckhorn_temps, TIME, into = c("date", "time"), sep = " ")

buckhorn_day_temps <- aggregate(clean.temps ~ date, data=buckhorn_temps, mean)

buckhorn_daily_temps <- separate(buckhorn_day_temps, date, into = c("year", "month", "day"), sep = "-")

write.csv(buckhorn_daily_temps, file="buckhorn_daily_temps.csv")

### aggregate to daily precip

buckhorn_precip <- subset(all_precip, SITE_NAME=="BUCKHORN2")

buckhorn_precip <- separate(buckhorn_precip, TIME, into = c("date", "time"), sep = " ")

buckhorn_day_precip <- aggregate(PPT.clean ~ date, data=buckhorn_precip, sum)

buckhorn_daily_precip <- separate(buckhorn_day_precip, date, into = c("year", "month", "day"), sep = "-")

write.csv(buckhorn_daily_precip, file="buckhorn_daily_precip.csv")

##check to see if the data looks logical 

plot(data=buckhorn_daily_temps,clean.temps ~ year)

plot(data=buckhorn_day_temps,clean.temps ~ date)

buckhorn_day_temps <- separate(buckhorn_temps, TIME, into = c("date", "time"), sep = " ")



#### Look at summaries of Buckhorn growth data


buck <- read.csv(file="buckhorn_growth_2007_2018.csv")

### check visually for differences between regions

ggplot(data = buck, aes(x = year, y = HT_cm)) + geom_point(aes(color = region)) + geom_smooth(aes(color = region), method = lm)



### what does the distribution of the y variable look like?

hist(buck$HT_cm)

buck2010 <-subset(buck, year=="2010")

hist(buck2010$HT_cm)

### Make a nice graph with the right colors 

theme_set(theme_bw(base_size=20))

bam <- ggplot(data = buck, aes(x = year, y = HT_cm)) + geom_point(aes(color = region)) + geom_smooth(aes(color = region), method = lm, size = 2)

blam <- bam + scale_color_manual(values=c("deepskyblue1", "purple","magenta","red","brown","darkblue","green","orange","yellow","darkslategray","black","bisque4"))

blam + labs(x="Year",y= "Height (cm)") 


