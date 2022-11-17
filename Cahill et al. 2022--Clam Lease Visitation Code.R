#Data and analysis for "Visitation patterns of two ray mesopredators at shellfish aquaculture leases in the Indian River Lagoon, Florida"
#Author: Brianna V. Cahill
#Purpose: Assess acoustic telemetry detections and visitation patterns for cownose and whitespotted eagle rays around Sebastian receivers (clam lease, inlet, river mouth)
#Input: clamleaseAcousticData.csv, CL_ReceiverPoints.csv, finaltides.csv, LOBOviz_EnvParameters.csv
#Output: publication figures and model results
#GitHub Repository: Clam-Lease-Visitation

#----------------------------------------- packages used -----------------------------------------#

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(lubridate)
library(suncalc)
library(FSA)
library(nortest)
library(viridis)
library(tidyverse)
library(stringr)
library(foreach)
library(VTrack)
library(mgcv)
library(ggpubr)
library(car)
library(GGally)
library(ellipse)
library(data.table)
library(stats)
library(MuMIn)
library(ggpubr)
library(lme4)
library(lmerTest)
library(nlme)
library(ggResidpanel)
library(MASS)
library(glmmTMB)
library(emmeans)
library(mgcViz)
library(scales)
library(gratia)

#----------------------------------------- preparing the data -----------------------------------------#

CompleteAcoData <- read.csv("Paper/Publication/clamleaseAcousticData.csv", header = TRUE)

#----- adjusting the structure -----#
CompleteAcoData$datecollected <- as.POSIXct(CompleteAcoData$datecollected, format = "%Y-%m-%d %H:%M:%S", tz = "UTC") 
CompleteAcoData$?..collectioncode <- as.factor(CompleteAcoData$?..collectioncode) 
CompleteAcoData$scientificname <- as.factor(CompleteAcoData$scientificname)
CompleteAcoData$commonname <- as.factor(CompleteAcoData$commonname)
CompleteAcoData$station <- as.factor(CompleteAcoData$station)
CompleteAcoData$receiver <- as.factor(CompleteAcoData$receiver)
CompleteAcoData$tagname <- as.factor(CompleteAcoData$tagname)
CompleteAcoData$timezone <- as.factor(CompleteAcoData$timezone)
CompleteAcoData$yearcollected <- as.factor(CompleteAcoData$yearcollected)
CompleteAcoData$monthcollected <- as.factor(CompleteAcoData$monthcollected)
CompleteAcoData$daycollected <- as.factor(CompleteAcoData$daycollected)
CompleteAcoData$julianday <- as.numeric(CompleteAcoData$julianday)

CompleteAcoData <- CompleteAcoData %>%
  separate(tagname, c("X", "Y", "ID"), "-")

#----- making additional column with EST time -----#
CompleteAcoData$EST <- with_tz(CompleteAcoData$datecollected, "America/New_York") 
CompleteAcoData <- subset(CompleteAcoData, CompleteAcoData$EST > as.POSIXct("2017-06-01 00:00:00", tz="EST") & CompleteAcoData$EST < as.POSIXct("2019-05-31 11:59:59", tz="EST"))

#----- prepping the acoustic data for VTrack format -----#
reformat_CompleteAcoData <- CompleteAcoData

names(reformat_CompleteAcoData)[names(reformat_CompleteAcoData) == "EST"] <- "Date.Time"
attr(reformat_CompleteAcoData, "tzone") <- "America/New_York"
names(reformat_CompleteAcoData)[names(reformat_CompleteAcoData) == "codespace"] <- "Code.Space"
names(reformat_CompleteAcoData)[names(reformat_CompleteAcoData) == "EST"] <- "Transmitter.S.N"
names(reformat_CompleteAcoData)[names(reformat_CompleteAcoData) == "receiver"] <- "Receiver.S.N"
names(reformat_CompleteAcoData)[names(reformat_CompleteAcoData) == "station"] <- "Station.Name"
names(reformat_CompleteAcoData)[names(reformat_CompleteAcoData) == "receiver_depth"] <- "Sensor.1"
names(reformat_CompleteAcoData)[names(reformat_CompleteAcoData) == "latitude"] <- "Station.Latitude"
names(reformat_CompleteAcoData)[names(reformat_CompleteAcoData) == "longitude"] <- "Station.Longitude"

reformat_CompleteAcoData$Units.1 <- "m"
reformat_CompleteAcoData$Sensor.2 <- NA
reformat_CompleteAcoData$Units.2 <- NA
reformat_CompleteAcoData$Transmitter.S.N <- NA

reformat_CompleteAcoData$Receiver.Name <- reformat_CompleteAcoData$Station.Name
reformat_CompleteAcoData$Transmitter.Name <- reformat_CompleteAcoData$ID

reformat_CompleteAcoData <- reformat_CompleteAcoData %>%
  select(Date.Time, Code.Space, ID, Sensor.1, Units.1, Sensor.2, Units.2, Transmitter.Name, Transmitter.S.N, Receiver.Name,
         Receiver.S.N, Station.Name, Station.Latitude, Station.Longitude) 

reformat_CompleteAcoData$Sensor.1 <- as.numeric(reformat_CompleteAcoData$Sensor.1)
reformat_CompleteAcoData$Units.1 <- as.character(reformat_CompleteAcoData$Units.1)
reformat_CompleteAcoData$Station.Name <- as.character(reformat_CompleteAcoData$Station.Name)

#----------------------------------------- VTrack, run residence extraction -----------------------------------------#

#----- receiver locations and detection radius -----#
receiverpoints <- read.csv ("CL_ReceiverPoints.csv", header = TRUE) #radius is in meters
receiverpoints_DM <- GenerateDirectDistance(receiverpoints)
colnames(receiverpoints_DM) <- c("DM", paste(receiverpoints$LOCATION))

#----- reading in the reformatted data -----#
Vclamlease <- ReadInputData(infile= reformat_CompleteAcoData, 
                            iHoursToAdd = 0, #added 0 hours because the date time has already been switched over to EST
                            fVemcoDualSensor = FALSE, 
                            dateformat =  NULL, 
                            sVemcoFormat = '1.0') #vemco format 1.0 is viewed as the "old" format, but since I chose to restructure the data based on the crocs formatting, I chose "1.0" instead of "default"

Vclamlease$TRANSMITTERID <- as.factor(Vclamlease$TRANSMITTERID)
Vclamlease$UNITS1 <- as.factor(Vclamlease$UNITS1)
Vclamlease$RECEIVERID <- as.factor(Vclamlease$RECEIVERID)
Vclamlease$STATIONNAME <- as.factor(Vclamlease$STATIONNAME)

#----- subsetting the data for the residence extraction -----#
V_SCLN <- subset(Vclamlease, Vclamlease$STATIONNAME == "SCLN")
V_SCLS <- subset(Vclamlease, Vclamlease$STATIONNAME == "SCLS")
V_IR16 <- subset(Vclamlease, Vclamlease$STATIONNAME == "IR16")
V_SINJ <- subset(Vclamlease, Vclamlease$STATIONNAME == "SINJ")
V_SISO <- subset(Vclamlease, Vclamlease$STATIONNAME == "SISO")
V_SIWP <- subset(Vclamlease, Vclamlease$STATIONNAME == "SIWP")

#----- determining the amount of time to pass through the 400 m detection radius -----#
#average ROM for SERs (DeGroot et al., 2020) was 1.4 km/hr, fastest ROM was 90 km/hr but was during the first 30 minutes of the track
#1400/(60*60) = 0.3888889 m/s
#400/0.3888889 = 1028.571 s to cross through this area, 17.14285 minutes

#----- residence extraction -----#
V_SCLN_extract <- RunResidenceExtraction(V_SCLN, "STATIONNAME", 2, 1028, sDistanceMatrix = receiverpoints_DM) #1028 seconds to pass through detection radius swimming at the average speed
V_SCLS_extract <- RunResidenceExtraction(V_SCLS, "STATIONNAME", 2, 1028, sDistanceMatrix = receiverpoints_DM)
V_IR16_extract <- RunResidenceExtraction(V_IR16, "STATIONNAME", 2, 1028, sDistanceMatrix = receiverpoints_DM)
V_SINJ_extract <- RunResidenceExtraction(V_SINJ, "STATIONNAME", 2, 1028, sDistanceMatrix = receiverpoints_DM)
V_SISO_extract <- RunResidenceExtraction(V_SISO, "STATIONNAME", 2, 1028, sDistanceMatrix = receiverpoints_DM)
V_SIWP_extract <- RunResidenceExtraction(V_SIWP, "STATIONNAME", 2, 1028, sDistanceMatrix = receiverpoints_DM)

all_extractions <- RunResidenceExtraction(Vclamlease, "STATIONNAME", 2, 1028, sDistanceMatrix = receiverpoints_DM)

#----- reading the residence extraction files -----#

SCLNlog <- V_SCLN_extract$residenceslog 
SCLSlog <- V_SCLS_extract$residenceslog 
IR16log <- V_IR16_extract$residenceslog 
SINJlog <- V_SINJ_extract$residenceslog 
SISOlog <- V_SISO_extract$residenceslog 
SIWPlog <- V_SIWP_extract$residenceslog 
all_extractions_log <- all_extractions$residenceslog

SCLNresid <- V_SCLN_extract$residences
SCLSresid <- V_SCLS_extract$residences
IR16resid <- V_IR16_extract$residences
SINJresid <- V_SINJ_extract$residences
SISOresid <- V_SISO_extract$residences
SIWPresid <- V_SIWP_extract$residences
all_extractions_resid <- all_extractions$residences
#residence event table that gives the start and end times of each event, the duration of the event, 
#the station and fish id, and the number of detections that are in each event

SCLNnonresid <- V_SCLN_extract$nonresidences
SCLSnonresid <- V_SCLS_extract$nonresidences 
IR16nonresid <- V_IR16_extract$nonresidences
SINJnonresid <- V_SINJ_extract$nonresidences
SISOnonresid <- V_SISO_extract$nonresidences
SIWPnonresid <- V_SIWP_extract$nonresidences
all_extractions_nonresid <- all_extractions$nonresidences
#nonresidences event table which gives the time a fish left a receiver, the time when it was detected on a different receiver, 
#and the time and distance (in meters) between those two detections

SCLNresid$minutes <- SCLNresid$DURATION /60
SCLSresid$minutes <- SCLSresid$DURATION /60
IR16resid$minutes <- IR16resid$DURATION /60
SINJresid$minutes <- SINJresid$DURATION /60
SISOresid$minutes <- SISOresid$DURATION /60
SIWPresid$minutes <- SIWPresid$DURATION /60

#----- merging the residency plots -----#
All_residency <- rbind.fill(SCLNresid, SCLSresid, IR16resid, SINJresid, SISOresid, SIWPresid)
attr(All_residency$STARTTIME, "tzone") <- "America/New_York"
attr(All_residency$ENDTIME, "tzone") <- "America/New_York"

#----- vlookup version in R to include individual animal info -----#
master <- read.csv("Master Acoustic Tag List_20210717.csv", header = TRUE)
master <- subset (master, master$Species == "Spotted Eagle Ray" | master$Species == "Cownose")
master <- subset (master, master$Coded.Continuous == "Coded")

tag_metadata <- read.csv("AcousticTagMetadata.csv", header = TRUE)
tag_metadata <- subset (tag_metadata, tag_metadata$Common.Name == "Spotted eagle ray" | 
                          tag_metadata$Common.Name == "Cownose ray" | 
                          tag_metadata$Common.Name == "Brazilian cownose ray")

names(tag_metadata)[names(tag_metadata) == "Tag.ID"] <- "TRANSMITTERID"
names(master)[names(master) == "OBJECTID"] <- "TRANSMITTERID"
tag_metadata$TRANSMITTERID <- as.factor(tag_metadata$TRANSMITTERID )

metadata <- merge(tag_metadata, unique(master)[, c("TRANSMITTERID", "Release.Date", "Expiration", "Region")], by="TRANSMITTERID", all.x=TRUE)
metadata$Release.Date<- as.Date(metadata$Release.Date, tz="America/New_York", format= "%m/%d/%y")
attr(metadata$Release.Date, "tzone") <- "America/New_York"
metadata$Expiration <- as.Date(metadata$Expiration, tz="America/New_York", format ="%m/%d/%y")
attr(metadata$Expiration, "tzone") <- "America/New_York"

All_residency_merge <- merge(All_residency, unique(metadata)[, c("TRANSMITTERID", "Sex", "Field.ID", "DW..cm.", "Tag.Serial.Number", "Common.Name", "Release.Date", "Expiration", "Region")], by="TRANSMITTERID", all.x=TRUE)

#----------------------------------------- sunrise/sunset info using suncalc -----------------------------------------#

#----- sun status for total detections -----#
suncalc_reformat <- reformat_CompleteAcoData

suncalc_reformat$Date <- as.Date(format(as.POSIXct(suncalc_reformat$Date.Time, tz="America/New_York"), "%Y-%m-%d"), tz="America/New_York")
attr(suncalc_reformat$Date, "tzone") <- "America/New_York" #the only way to seemingly establish what the time zone is for as.Date function

sebastian.sun <- getSunlightTimes(date= suncalc_reformat$Date, 
                                  lat = 27.86, #took the mean latitude
                                  lon = -80.46, #took the mean longitude
                                  keep = c("solarNoon", "nadir", "sunrise", "sunset", "sunriseEnd",
                                           "sunsetStart", "dawn", "dusk", "nauticalDawn", "nauticalDusk",
                                           "nightEnd", "night", "goldenHourEnd", "goldenHour"), 
                                  tz = "EST")
suncalc_reformat$Sunrise<-sebastian.sun$sunrise
suncalc_reformat$Dawn<-sebastian.sun$dawn
suncalc_reformat$Sunset<-sebastian.sun$sunset
suncalc_reformat$Dusk<-sebastian.sun$dusk

suncalc_reformat$DawnBuffer <- suncalc_reformat$Sunrise-suncalc_reformat$Dawn #determined time diff between sunrise and dawn, provided in minutes
suncalc_reformat$DuskBuffer <- suncalc_reformat$Dusk-suncalc_reformat$Sunset #determined time diff between sunset and dusk, provided in minutes

suncalc_reformat$PreDawn <- suncalc_reformat$Dawn-suncalc_reformat$DawnBuffer #used buffer time to determine pre-dawn timeframe
suncalc_reformat$PostDusk <- suncalc_reformat$Dusk + suncalc_reformat$DuskBuffer #used buffer time to determine post-dusk timeframe

suncalc_reformat$TOD<- with(suncalc_reformat, 
                            ifelse(suncalc_reformat$Date.Time >= suncalc_reformat$PreDawn & suncalc_reformat$Date.Time < suncalc_reformat$Sunrise, "Dawn",
                                   ifelse(suncalc_reformat$Date.Time >=suncalc_reformat$Sunrise & suncalc_reformat$Date.Time < suncalc_reformat$Sunset, "Day",
                                          ifelse(suncalc_reformat$Date.Time >= suncalc_reformat$Sunset & suncalc_reformat$Date.Time < suncalc_reformat$PostDusk, "Dusk", "Night"))))
#made a time of day column to effectively note the different periods of the day

detections <- subset(suncalc_reformat, select = -c(Sunrise, Dawn, Sunset, Dusk, DawnBuffer, DuskBuffer, PreDawn, PostDusk))
#made a new dataframe omitting the now unnecessary columns
attr(detections$Date, "tzone") <- "America/New_York" #the only way to seemingly establish what the time zone is for as.Date function

#----- sun status for residency detections -----#
All_residency_merge$Start_Date <- as.Date(format(as.POSIXct(All_residency_merge$STARTTIME)))
attr(All_residency_merge$Start_Date, "tzone") <- "America/New_York"

start.sun <- getSunlightTimes(date= All_residency_merge$Start_Date, 
                              lat = 27.86, #took the mean latitude
                              lon = -80.46, #took the mean longitude
                              keep = c("solarNoon", "nadir", "sunrise", "sunset", "sunriseEnd",
                                       "sunsetStart", "dawn", "dusk", "nauticalDawn", "nauticalDusk",
                                       "nightEnd", "night", "goldenHourEnd", "goldenHour"), 
                              tz = "EST")

All_residency_merge$Start.Sunrise<-start.sun$sunrise
All_residency_merge$Start.Sunset<-start.sun$sunset
All_residency_merge$Start.Dawn<-start.sun$dawn
All_residency_merge$Start.Dusk<-start.sun$dusk

All_residency_merge$Start.DawnBuffer <- All_residency_merge$Start.Sunrise-All_residency_merge$Start.Dawn #determined time diff between sunrise and dawn, provided in minutes
All_residency_merge$Start.DuskBuffer <- All_residency_merge$Start.Dusk-All_residency_merge$Start.Sunset #determined time diff between sunset and dusk, provided in minutes

All_residency_merge$Start.PreDawn <- All_residency_merge$Start.Dawn-All_residency_merge$Start.DawnBuffer #used buffer time to determine pre-dawn timeframe
All_residency_merge$Start.PostDusk <- All_residency_merge$Start.Dusk + All_residency_merge$Start.DuskBuffer #used buffer time to determine post-dusk timeframe

All_residency_merge$Start.TOD<- with(All_residency_merge, 
                                     ifelse(All_residency_merge$STARTTIME >= All_residency_merge$Start.PreDawn & All_residency_merge$STARTTIME < All_residency_merge$Start.Sunrise, "Dawn",
                                            ifelse(All_residency_merge$STARTTIME >=All_residency_merge$Start.Sunrise & All_residency_merge$STARTTIME < All_residency_merge$Start.Sunset, "Day",
                                                   ifelse(All_residency_merge$STARTTIME >= All_residency_merge$Start.Sunset & All_residency_merge$STARTTIME < All_residency_merge$Start.PostDusk, "Dusk", "Night"))))

visitations <- subset(All_residency_merge, select = -c(Start.Sunrise, Start.Dawn, Start.Sunset, Start.Dusk, 
                                                       Mid.Sunrise, Mid.Sunset, Mid.Dawn, Mid.Dusk,
                                                       End.Sunrise, End.Sunset, End.Dawn, End.Dusk,
                                                       Start.DawnBuffer, Start.DuskBuffer, Start.PreDawn, Start.PostDusk,
                                                       Mid.DawnBuffer, Mid.DuskBuffer, Mid.PreDawn, Mid.PostDusk,
                                                       End.DawnBuffer, End.DuskBuffer, End.PreDawn, End.PostDusk))

attr(visitations$Start_Date, "tzone") <- "America/New_York"

#----------------------------------------- restructuring dfs -----------------------------------------#

#----- visitations -----#
visitations$TRANSMITTERID <- as.factor(visitations$TRANSMITTERID)
visitations$RESIDENCEEVENT <- as.factor(visitations$RESIDENCEEVENT)
visitations$STATIONNAME <- as.factor(visitations$STATIONNAME)
visitations$ENDREASON <- as.factor(visitations$ENDREASON)
visitations$Sex <- as.factor(visitations$Sex)
visitations$Field.ID <- as.factor(visitations$Field.ID)
visitations$Tag.Serial.Number <- as.factor(visitations$Tag.Serial.Number)
visitations$Common.Name <- as.factor(visitations$Common.Name)
visitations$TRANSMITTERID <- as.factor(visitations$TRANSMITTERID)
visitations$Start.TOD<-as.factor(visitations$Start.TOD)
visitations$Mid.TOD<-as.factor(visitations$Mid.TOD)
visitations$End.TOD<-as.factor(visitations$End.TOD)

#----- detections -----#
detections$Code.Space <- as.factor(detections$Code.Space)
detections$ID <- as.factor(detections$ID)
detections$Transmitter.Name <- as.factor(detections$Transmitter.Name)
detections$Transmitter.S.N <- as.factor(detections$Transmitter.S.N)
detections$Station.Name <- as.factor(detections$Station.Name)
detections$TOD<-as.factor(detections$TOD)

#----------------------------------------- figuring out tide shit -----------------------------------------#
finaltides <- read.csv ("Tides/finaltides.csv", header = TRUE)
test <- finaltides

test$Date.Time <- with(test, as.POSIXct(paste(Date, Time), format= "%Y-%m-%d %H:%M"))
test$Date <- as.Date(as.character(test$Date), format= "%Y-%m-%d", tz="EST")
attr(test$Date, "tzone") <- "America/New_York"
names(test)[names(test) == "High.Low"] <- "Tide"

test$StartTide<- test$Date.Time - 1800
test$EndTide<- test$Date.Time + 1800

test <- test[ order(test$Date.Time , decreasing = FALSE ),]
rownames(test) <- NULL

empty2 <- test[rep(1:nrow(test), each = 2), ] #duplicates each row
empty2[1:nrow(empty2) %% 2 == 0, ] <- NA #makes each duplicate row NAs, edit tide shit in Excel

for (i in 1:dim(empty2)[1]) {
  #Grab data from dataframe
  #this section becomes the reference data
  if( (i%%2) != 0) {
    tide <- empty2[i,"Tide"]
    startTide <- empty2[i+2, "StartTide"]
    endTide <- empty2[i, "EndTide"]
  }
  
  if( (i%%2) == 0) {
    #check if tide is high on previous step, and assign following step as outgoing if TRUE
    #this section becomes the writing data
    ifelse(tide == "High", 
           empty2[i, "Tide"] <- "Outgoing",
           empty2[i, "Tide"] <- "Incoming")
    
    require(lubridate)
    empty2[i, "StartTide"] <- endTide + seconds(1)
    
    empty2[i, "EndTide"] <- startTide - seconds(1) ######### This is taking the previous starting tide, not the following starting tide
  }
}

rownames(empty2) <- NULL

alltides2 <- empty2
alltides2$Tide <- as.factor(alltides2$Tide)
alltides2$Pred <- as.numeric(alltides2$Pred)

#IR16 offset
#about a 2 hour and 20 minute offset to Sebastian downtown, time = (60*140) --> Time = 8400 seconds
#given distance of 3689 m, we'll say that the water is moving at about 0.439 m/s
#rate = distance/time, time = distance/rate --> time = 3689/8400, time = 8400 seconds
IR16_tides <- alltides2
IR16_tides$StartTide <- IR16_tides$StartTide - 8400
IR16_tides$EndTide <- IR16_tides$EndTide - 8400

#SCLN offset
#calculating time offset with water movement at 0.439 m/s and shortest distance of 6002 m
#rate = distance/time, time = distance/rate --> time = 6002/0.439, time = 13671.98 seconds
SCLN_tides <- alltides2
SCLN_tides$StartTide <- SCLN_tides$StartTide - 13672
SCLN_tides$EndTide <- SCLN_tides$EndTide - 13672

#SCLS offset
#calculating time offset with water movement at 0.439 m/s and shortest distance of 4344 m
#rate = distance/time, time = distance/rate --> time = 4344/0.428, time = 10149.53 seconds
SCLS_tides <- alltides2
SCLS_tides$StartTide <- SCLS_tides$StartTide - 10150
SCLS_tides$EndTide <- SCLS_tides$EndTide - 10150

#SINJ offset
#tide information is derived very close to this receiver, no offset was used
SINJ_tides <- alltides2
SINJ_tides$StartTide <- SINJ_tides$StartTide - 0
SINJ_tides$EndTide <- SINJ_tides$EndTide - 0

#SISO offset
#calculating time offset with water movement at 0.439 m/s and shortest distance of 431.5345 m
#rate = distance/time, time = distance/rate --> time = 432/0.428, time = 1009.346 seconds
SISO_tides <- alltides2
SISO_tides$StartTide <- SISO_tides$StartTide - 1009
SISO_tides$EndTide <- SISO_tides$EndTide - 1009

#SIWP offset  
#calculating time offset with water movement at 0.439 m/s and shortest distance of 84.9431 m
#rate = distance/time, time = distance/rate --> time = 85/0.428, time = 198.5981 seconds
SIWP_tides <- alltides2
SIWP_tides$StartTide <- SIWP_tides$StartTide - 199
SIWP_tides$EndTide <- SIWP_tides$EndTide - 199

#----------------------------------------------#
#----- Assigning Tides to Visitation Data -----#
#----------------------------------------------#
{
  #----- subsetting visitations -----#
  #visitations$Start.Tide <- NA
  visitations$Tide <- NA
  visitations$Tide <- as.character(visitations$Tide)
  test_visitations <- visitations
  IR16_visitations <- subset(test_visitations, test_visitations$STATIONNAME == "IR16")
  SCLN_visitations <- subset(test_visitations, test_visitations$STATIONNAME == "SCLN")
  SCLS_visitations <- subset(test_visitations, test_visitations$STATIONNAME == "SCLS")
  SINJ_visitations <- subset(test_visitations, test_visitations$STATIONNAME == "SINJ")
  SISO_visitations <- subset(test_visitations, test_visitations$STATIONNAME == "SISO")
  SIWP_visitations <- subset(test_visitations, test_visitations$STATIONNAME == "SIWP")
  
  #----- TEST Visitation -----#
  
  test_visit_tidecheck <- visitations[sample(nrow(visitations), 50), ]
  
  {
    for (i in 1:dim(test_visit_tidecheck)[1]) {
      #extract the date time
      dt <- test_visit_tidecheck[i,"Date.Time"]
      
      for(k in 1:dim(alltides2)[1]){
        #compare the date time data to the start and end times in the detections data
        start <- alltides2[k,"StartTide"]
        end <- alltides2[k,"EndTide"]
        t <- alltides2[k,"Tide"]
        if(dt > start & dt < end) {
          break
        }
      }
      
      test_visit_tidecheck[i,"Tide"] <- paste(t)
    }
  } #took about a minute for 50 data points
  check <- subset(test_visit_tidecheck, select = c(Date.Time, Tide))
  #----- IR16 Visitation -----#
  {
    for (i in 1:dim(IR16_visitations)[1]) {
      #extract the date time
      dt <- IR16_visitations$STARTTIME[i]
      
      for(k in 1:dim(IR16_tides)[1]){
        #compare the date time data to the start and end times in the detections data
        start <- IR16_tides$StartTide[k]
        end <- IR16_tides$EndTide[k]
        t <- IR16_tides$Tide[k]
        if(dt >= start & dt <= end) {
          break
        }
      }
      
      IR16_visitations$Tide[i] <- paste(t)
    }
  } 
  
  check <- subset(IR16_visitations, select = c(Date.Time, Tide))
  
  #----- SCLN Visitation TEST -----#
  {
    for (i in 1:dim(SCLN_visitations)[1]) {
      #extract the date time
      dt <- SCLN_visitations$STARTTIME[i]
      
      for(k in 1:dim(SCLN_tides)[1]){
        #compare the date time data to the start and end times in the detections data
        start <- SCLN_tides$StartTide[k]
        end <- SCLN_tides$EndTide[k]
        t <- SCLN_tides$Tide[k]
        if(dt >= start & dt <= end) {
          break
        }
      }
      
      SCLN_visitations$Tide[i] <- paste(t)
    }
  } 
  
  check <- subset(SCLN_visitations, select = c(Date.Time, Tide))
  
  #----- SCLS Visitation -----#
  {
    for (i in 1:dim(SCLS_visitations)[1]) {
      #extract the date time
      dt <- SCLS_visitations$STARTTIME[i]
      
      for(k in 1:dim(SCLS_tides)[1]){
        #compare the date time data to the start and end times in the detections data
        start <- SCLS_tides$StartTide[k]
        end <- SCLS_tides$EndTide[k]
        t <- SCLS_tides$Tide[k]
        if(dt >= start & dt <= end) {
          break
        }
      }
      
      SCLS_visitations$Tide[i] <- paste(t)
    }
  }
  check <- subset(SCLS_visitations, select = c(Date.Time, Tide))
  
  #----- SINJ Visitation -----#
  {
    for (i in 1:dim(SINJ_visitations)[1]) {
      #extract the date time
      dt <- SINJ_visitations$STARTTIME[i]
      
      for(k in 1:dim(SINJ_tides)[1]){
        #compare the date time data to the start and end times in the detections data
        start <- SINJ_tides$StartTide[k]
        end <- SINJ_tides$EndTide[k]
        t <- SINJ_tides$Tide[k]
        if(dt >= start & dt <= end) {
          break
        }
      }
      
      SINJ_visitations$Tide[i] <- paste(t)
    }
  } 
  check <- subset(SINJ_visitations, select = c(Date.Time, Tide))
  
  #----- SISO Visitation -----#
  {
    for (i in 1:dim(SISO_visitations)[1]) {
      #extract the date time
      dt <- SISO_visitations$STARTTIME[i]
      
      for(k in 1:dim(SISO_tides)[1]){
        #compare the date time data to the start and end times in the detections data
        start <- SISO_tides$StartTide[k]
        end <- SISO_tides$EndTide[k]
        t <- SISO_tides$Tide[k]
        if(dt >= start & dt <= end) {
          break
        }
      }
      
      SISO_visitations$Tide[i] <- paste(t)
    }
  }
  check <- subset(SISO_visitations, select = c(Date.Time, Tide))
  
  #----- SIWP Visitation -----#
  {
    for (i in 1:dim(SIWP_visitations)[1]) {
      #extract the date time
      dt <- SIWP_visitations$STARTTIME[i]
      
      for(k in 1:dim(SIWP_tides)[1]){
        #compare the date time data to the start and end times in the detections data
        start <- SIWP_tides$StartTide[k]
        end <- SIWP_tides$EndTide[k]
        t <- SIWP_tides$Tide[k]
        if(dt >= start & dt <= end) {
          break
        }
      }
      
      SIWP_visitations$Tide[i] <- paste(t)
    }
  }
  check <- subset(SIWP_visitations, select = c(Date.Time, Tide))
  
  visitations <- rbind(IR16_visitations, SCLN_visitations, SCLS_visitations, SINJ_visitations, SISO_visitations, SIWP_visitations)
  visitations$Tide <- as.factor(visitations$Tide)

}

#----------------------------------------- including additional env data -----------------------------------------#
LOBO <- read.csv("LOBOviz_EnvParameters.csv", header=TRUE)
names(LOBO)[names(LOBO) == "?..Date.Time"] <- "Date.Time"

#LOBO$Date.Time <- as.POSIXct(paste(LOBO$Date_EST, LOBO$Time_EST), format = "%Y-%m-%d %H:%M:%S")
LOBO$Date.Time <- as.POSIXct(LOBO$Date.Time, format = "%m/%d/%y %H:%M")
attr(LOBO$Date.Time, "tzone") <- "America/New_York"
LOBO$Date.Time.LOBO <-LOBO$Date.Time
attr(LOBO$Date.Time.LOBO, "tzone") <- "America/New_York"

LOBO$DO[LOBO$DO < 0] <- NA
LOBO$Nitrate[LOBO$Nitrate < 0] <- NA
LOBO$pH[LOBO$pH < 1 | LOBO$pH > 14] <- NA
LOBO$Phosphate[LOBO$Phosphate < 0] <- NA

LOBO_t <- as.data.table(LOBO)
setkeyv(LOBO_t, "Date.Time")

#----- visitations -----#
visitations$Date.Time <- visitations$STARTTIME
visitations_t <- as.data.table(visitations)
setkeyv(visitations_t, "STARTTIME")
visitations <- LOBO_t[visitations_t, roll="nearest", on="Date.Time"]

#----------------------------------------- adding in hours/decimal hours, seasons/julian day -----------------------------------------#

#----- detections seasons -----#
{
  detections$Season[detections$Date >= "2017-01-01" & detections$Date < "2017-03-20"] <- "Winter"
  detections$Season[detections$Date >= "2017-03-20" & detections$Date < "2017-06-20"] <- "Spring"
  detections$Season[detections$Date >= "2017-06-20" & detections$Date < "2017-09-22"] <- "Summer"
  detections$Season[detections$Date >= "2017-09-22" & detections$Date < "2017-12-21"] <- "Fall"
  detections$Season[detections$Date >= "2017-12-21" & detections$Date < "2018-03-20"] <- "Winter"
  detections$Season[detections$Date >= "2018-03-20" & detections$Date < "2018-06-21"] <- "Spring"
  detections$Season[detections$Date >= "2018-06-21" & detections$Date < "2018-09-22"] <- "Summer"
  detections$Season[detections$Date >= "2018-09-22" & detections$Date < "2018-12-21"] <- "Fall"
  detections$Season[detections$Date >= "2018-12-21" & detections$Date < "2019-03-20"] <- "Winter"
  detections$Season[detections$Date >= "2019-03-20" & detections$Date < "2019-06-21"] <- "Spring"
  detections$Season[detections$Date >= "2019-06-21" & detections$Date < "2019-09-23"] <- "Summer"
  detections$Season[detections$Date >= "2019-09-23" & detections$Date < "2019-12-21"] <- "Fall"
  detections$Season[detections$Date >= "2019-12-21" & detections$Date < "2020-03-19"] <- "Winter"
  detections$Season[detections$Date >= "2020-03-19" & detections$Date < "2020-06-20"] <- "Spring"
  detections$Season[detections$Date >= "2020-06-20" & detections$Date < "2020-09-22"] <- "Summer"
  detections$Season[detections$Date >= "2020-09-22" & detections$Date < "2020-12-21"] <- "Fall"
  detections$Season[detections$Date >= "2020-12-21" & detections$Date < "2021-03-20"] <- "Winter"
  detections$Season[detections$Date >= "2021-03-20" & detections$Date < "2021-06-20"] <- "Spring"
  detections$Season[detections$Date >= "2021-06-20" & detections$Date < "2021-09-22"] <- "Summer"
  detections$Season[detections$Date >= "2021-09-22" & detections$Date < "2021-12-21"] <- "Fall"
  
  detections$Season <- as.factor(detections$Season)
  detections$Season <- factor(detections$Season, levels = c("Winter", "Spring", "Summer", "Fall"))
  
  #----- visitation seasons -----#
  visitations$Season[visitations$Start_Date >= "2017-01-01" & visitations$Start_Date < "2017-03-20"] <- "Winter"
  visitations$Season[visitations$Start_Date >= "2017-03-20" & visitations$Start_Date < "2017-06-20"] <- "Spring"
  visitations$Season[visitations$Start_Date >= "2017-06-20" & visitations$Start_Date < "2017-09-22"] <- "Summer"
  visitations$Season[visitations$Start_Date >= "2017-09-22" & visitations$Start_Date < "2017-12-21"] <- "Fall"
  visitations$Season[visitations$Start_Date >= "2017-12-21" & visitations$Start_Date < "2018-03-20"] <- "Winter"
  visitations$Season[visitations$Start_Date >= "2018-03-20" & visitations$Start_Date < "2018-06-21"] <- "Spring"
  visitations$Season[visitations$Start_Date >= "2018-06-21" & visitations$Start_Date < "2018-09-22"] <- "Summer"
  visitations$Season[visitations$Start_Date >= "2018-09-22" & visitations$Start_Date < "2018-12-21"] <- "Fall"
  visitations$Season[visitations$Start_Date >= "2018-12-21" & visitations$Start_Date < "2019-03-20"] <- "Winter"
  visitations$Season[visitations$Start_Date >= "2019-03-20" & visitations$Start_Date < "2019-06-21"] <- "Spring"
  visitations$Season[visitations$Start_Date >= "2019-06-21" & visitations$Start_Date < "2019-09-23"] <- "Summer"
  visitations$Season[visitations$Start_Date >= "2019-09-23" & visitations$Start_Date < "2019-12-21"] <- "Fall"
  visitations$Season[visitations$Start_Date >= "2019-12-21" & visitations$Start_Date < "2020-03-19"] <- "Winter"
  visitations$Season[visitations$Start_Date >= "2020-03-19" & visitations$Start_Date < "2020-06-20"] <- "Spring"
  visitations$Season[visitations$Start_Date >= "2020-06-20" & visitations$Start_Date < "2020-09-22"] <- "Summer"
  visitations$Season[visitations$Start_Date >= "2020-09-22" & visitations$Start_Date < "2020-12-21"] <- "Fall"
  visitations$Season[visitations$Start_Date >= "2020-12-21" & visitations$Start_Date < "2021-03-20"] <- "Winter"
  visitations$Season[visitations$Start_Date >= "2021-03-20" & visitations$Start_Date < "2021-06-20"] <- "Spring"
  visitations$Season[visitations$Start_Date >= "2021-06-20" & visitations$Start_Date < "2021-09-22"] <- "Summer"
  visitations$Season[visitations$Start_Date >= "2021-09-22" & visitations$Start_Date < "2021-12-21"] <- "Fall"
  
  visitations$Season <- as.factor(visitations$Season)
  visitations$Season <- factor(visitations$Season, levels = c("Winter", "Spring", "Summer", "Fall"))
}
#----- hourly plots for detections -----#
detections$Hour <- as.factor(format(as.POSIXct(detections$Date.Time,format = "%Y-%m-%dT%H:%M:%S"),"%H"))
visitations$Start.Hour <- as.factor(format(as.POSIXct(visitations$STARTTIME,format = "%Y-%m-%dT%H:%M:%S"),"%H"))
#visitations$Mid.Hour <- as.factor(format(as.POSIXct(visitations$MIDTIME,format = "%Y-%m-%dT%H:%M:%S"),"%H"))
#visitations$End.Hour <- as.factor(format(as.POSIXct(visitations$ENDTIME,format = "%Y-%m-%dT%H:%M:%S"),"%H"))

#----- monthly and yearly plots for detections and visitations -----#
visitations$Month <- as.factor(format(as.POSIXct(visitations$STARTTIME,format = "%Y-%m-%dT%H:%M:%S"),"%m"))
visitations$Year <- as.factor(format(as.POSIXct(visitations$STARTTIME,format = "%Y-%m-%dT%H:%M:%S"),"%Y"))

detections$Month <- as.factor(format(as.POSIXct(detections$Date.Time,format = "%Y-%m-%dT%H:%M:%S"),"%m"))
detections$Year <- as.factor(format(as.POSIXct(detections$Date.Time,format = "%Y-%m-%dT%H:%M:%S"),"%Y"))

#----- year num -----#
visitations$YearNum[visitations$Start_Date >= "2017-06-01" & visitations$Start_Date < "2018-06-01"] <- 1
visitations$YearNum[visitations$Start_Date >= "2018-06-01" & visitations$Start_Date < "2019-06-01"] <- 2

detections$YearNum[detections$Date >= "2017-06-01" & detections$Date < "2018-06-01"] <- 1
detections$YearNum[detections$Date >= "2018-06-01" & detections$Date < "2019-06-01"] <- 2

#----- merging individual info with the detections data -----#
detections$TRANSMITTERID <- detections$Transmitter.Name

detections <- merge(detections, unique(metadata)[, c("TRANSMITTERID", "Sex", "Field.ID", "DW..cm.", "Tag.Serial.Number", "Common.Name", "Release.Date", "Expiration", "Region")], by="TRANSMITTERID", all.x=TRUE, allow.cartesian = TRUE)

#----- decimal hour -----#
visitations$Dec.Hour <- as.numeric(format(as.POSIXct(visitations$STARTTIME,format = "%Y-%m-%dT%H:%M:%S"),"%H"))
visitations$Dec.Min <- as.numeric(format(as.POSIXct(visitations$STARTTIME,format = "%Y-%m-%dT%H:%M:%S"),"%M"))
visitations$Dec.Min <- visitations$Dec.Min/60
visitations$Dec.Sec <- as.numeric(format(as.POSIXct(visitations$STARTTIME,format = "%Y-%m-%dT%H:%M:%S"),"%S"))
visitations$Dec.Sec <- visitations$Dec.Sec/(60*60)
visitations$Dec.Hour <- visitations$Dec.Hour + visitations$Dec.Min + visitations$Dec.Sec
visitations <- subset(visitations, select = -c(Dec.Min, Dec.Sec))

detections$Dec.Hour <- as.numeric(format(as.POSIXct(detections$Date.Time,format = "%Y-%m-%dT%H:%M:%S"),"%H"))
detections$Dec.Min <- as.numeric(format(as.POSIXct(detections$Date.Time,format = "%Y-%m-%dT%H:%M:%S"),"%M"))
detections$Dec.Min <- detections$Dec.Min/60
detections$Dec.Sec <- as.numeric(format(as.POSIXct(detections$Date.Time,format = "%Y-%m-%dT%H:%M:%S"),"%S"))
detections$Dec.Sec <- detections$Dec.Sec/(60*60)
detections$Dec.Hour <- detections$Dec.Hour + detections$Dec.Min + detections$Dec.Sec
detections <- subset(detections, select = -c(Dec.Min, Dec.Sec))

#----- animal code -----#
detections$AnimalCode[detections$Common.Name =="Spotted eagle ray"] <- "SER"
detections$AnimalCode[detections$Common.Name == "Cownose ray" | detections$Common.Name == "Brazilian cownose ray"] <- "CNR"
detections$AnimalCode <- as.factor(detections$AnimalCode)

visitations$AnimalCode[visitations$Common.Name == "Spotted eagle ray"] <- "SER"
visitations$AnimalCode[visitations$Common.Name == "Cownose ray" | visitations$Common.Name == "Brazilian cownose ray"] <- "CNR"
visitations$AnimalCode <- as.factor(visitations$AnimalCode)

#----- clam lease vs. not -----#
visitations$CL_Other <- with(visitations, ifelse(visitations$STATIONNAME == "SCLN" | 
                                                   visitations$STATIONNAME == "SCLS", "Clam Lease", "Other"))
visitations$CL_Other <- as.factor(visitations$CL_Other)

detections$CL_Other <- with(visitations, ifelse(detections$Station.Name == "SCLN" | 
                                                  detections$Station.Name == "SCLS", "Clam Lease", "Other"))
detections$CL_Other <- as.factor(detections$CL_Other)

#----------------------------------------- moonphase info, suncalc -----------------------------------------#
#----- detections -----#
sebastian.moon_detect <- getMoonIllumination(date= detections$Date, 
                                             keep = c("fraction", "phase", "angle"))

detections$Moon_Fraction<-sebastian.moon_detect$fraction
detections$Moon_Phase<-sebastian.moon_detect$phase
detections$Moon_Angle<-sebastian.moon_detect$angle

#----- visitations -----#
sebastian.moon_visit <- getMoonIllumination(date= visitations$Start_Date, 
                                            keep = c("fraction", "phase", "angle"))

visitations$Moon_Fraction<-sebastian.moon_visit$fraction
visitations$Moon_Phase<-sebastian.moon_visit$phase
visitations$Moon_Angle<-sebastian.moon_visit$angle

#----------------------------------------- matching transmitter ID to individual code -----------------------------------------#
transIDs <- c("16921", "16922", "16923", "16924", "46180", "14587", "14589", "14584", "14585", "14586",
              "13488", "13490", "84", "85", "86", "87", "7684", "7683", "7682", "7681", 
              "16927", "16929", "16930", "16942", "16943", "16944", "16945", "16946", "16947", "1706", 
              "1707", "1708", "17626", "15429", "15428", "17624", "14583", "15423")

animalNumber <- c("SER 1", "SER 2", "SER 3", "SER 4", "SER 5", "SER 6", "SER 7", "SER 8", "SER 9", "SER 10",
                  "SER 11", "SER 12", "SER 13", "SER 14", "SER 15", "SER 16", "SER 17", "SER 18", "SER 19", "SER 20", 
                  "CNR 1", "CNR 2", "CNR 3", "CNR 4", "CNR 5", "CNR 6", "CNR 7", "CNR 8", "CNR 9", "CNR 10", 
                  "CNR 11", "CNR 12", "CNR 13", "CNR 14", "CNR 15", "CNR 16", "CNR 17", "CNR 18")
tagAnimalNumber <- data.frame(transIDs, animalNumber)
names(tagAnimalNumber)[names(tagAnimalNumber) == "transIDs"] <- "TRANSMITTERID"
tagAnimalNumber$TRANSMITTERID <- as.factor(tagAnimalNumber$TRANSMITTERID)
tagAnimalNumber$animalNumber <- as.factor(tagAnimalNumber$animalNumber)

#----- for detections -----#
detections <- merge(detections, unique(tagAnimalNumber)[, c("TRANSMITTERID", "animalNumber")], by="TRANSMITTERID", all.x=TRUE, allow.cartesian = TRUE)
detections$animalNumber <- factor(detections$animalNumber, levels = c("SER 1", "SER 2", "SER 3", "SER 4", "SER 5", "SER 6", "SER 7", "SER 8", "SER 9", "SER 10",
                                                                      "SER 11", "SER 12", "SER 13", "SER 14", "SER 15", "SER 16", "SER 17", "SER 18", "SER 19", "SER 20", 
                                                                      "CNR 1", "CNR 2", "CNR 3", "CNR 4", "CNR 5", "CNR 6", "CNR 7", "CNR 8", "CNR 9", "CNR 10", 
                                                                      "CNR 11", "CNR 12", "CNR 13", "CNR 14", "CNR 15", "CNR 16", "CNR 17", "CNR 18"))

#----- for visitations -----#
visitations <- merge(visitations, unique(tagAnimalNumber)[, c("TRANSMITTERID", "animalNumber")], by="TRANSMITTERID", all.x=TRUE, allow.cartesian = TRUE)
visitations$animalNumber <- factor(visitations$animalNumber, levels = c("SER 1", "SER 2", "SER 3", "SER 4", "SER 5", "SER 6", "SER 7", "SER 8", "SER 9", "SER 10",
                                                                        "SER 11", "SER 12", "SER 13", "SER 14", "SER 15", "SER 16", "SER 17", "SER 18", "SER 19", "SER 20", 
                                                                        "CNR 1", "CNR 2", "CNR 3", "CNR 4", "CNR 5", "CNR 6", "CNR 7", "CNR 8", "CNR 9", "CNR 10", 
                                                                        "CNR 11", "CNR 12", "CNR 13", "CNR 14", "CNR 15", "CNR 16", "CNR 17", "CNR 18"))


#----------------------------------------- models -----------------------------------------#

#----- SER: time of day*station*sex GLMM -----#
test2 <- detections %>%
  filter(AnimalCode == "SER") %>%
  group_by(TRANSMITTERID, as.Date(Date.Time), Station.Name, TOD, Sex) %>%
  summarise(count = n())

test2 <- as.data.frame(test2)
test2$Sex <- as.factor(test2$Sex)
test2$TOD_cond<- with(test2, 
                      ifelse(test2$TOD == "Day", "Day",
                             ifelse(test2$TOD == "Night", "Night", "Crepuscular")))
test2$TOD_cond <- as.factor(test2$TOD_cond)

test2_glmm <- glmer(count ~ TOD_cond*Station.Name*Sex + (1 | TRANSMITTERID), data = test2, family = poisson())
summary(test2_glmm)
Anova(test2_glmm)
emmeans(test2_glmm, list(pairwise ~ TOD_cond), adjust = "tukey")
emmeans(test2_glmm, list(pairwise ~ Station.Name), adjust = "tukey")
emmeans(test2_glmm, list(pairwise ~ Sex), adjust = "tukey")
emmeans(test2_glmm, list(pairwise ~ TOD_cond*Station.Name), adjust = "tukey")
emmeans(test2_glmm, list(pairwise ~ TOD_cond*Sex), adjust = "tukey")
emmeans(test2_glmm, list(pairwise ~ Sex*Station.Name), adjust = "tukey")
emmeans(test2_glmm, list(pairwise ~ TOD_cond*Station.Name*Sex), adjust = "tukey")

#----- CNR: time of day*station*sex GLMM -----#
test3 <- detections %>%
  filter(AnimalCode == "CNR") %>%
  group_by(TRANSMITTERID, as.Date(Date.Time), Station.Name, TOD, Sex) %>%
  summarise(count = n())

test3 <- as.data.frame(test3)
test3$Sex <- as.factor(test3$Sex)
test3$TOD_cond<- with(test3, 
                      ifelse(test3$TOD == "Day", "Day",
                             ifelse(test3$TOD == "Night", "Night", "Crepuscular")))
test3$TOD_cond <- as.factor(test3$TOD_cond)

test3_glmm <- glmer(count ~ TOD_cond*Station.Name*Sex + (1 | TRANSMITTERID), data = test3, family = poisson())
summary(test3_glmm)
Anova(test3_glmm)
emmeans(test3_glmm, list(pairwise ~ TOD_cond), adjust = "tukey")
emmeans(test3_glmm, list(pairwise ~ Station.Name), adjust = "tukey")
emmeans(test3_glmm, list(pairwise ~ Sex), adjust = "tukey")
emmeans(test3_glmm, list(pairwise ~ TOD_cond*Station.Name), adjust = "tukey")
emmeans(test3_glmm, list(pairwise ~ TOD_cond*Sex), adjust = "tukey")
emmeans(test3_glmm, list(pairwise ~ Sex*Station.Name), adjust = "tukey")
emmeans(test3_glmm, list(pairwise ~ TOD_cond*Station.Name*Sex), adjust = "tukey")

#----- checking for collinearity -----#
nums <- subset(visitations, select = c("DO",  "Salinity",  "DO_sat", "Temp_C", "minutes", 
                                       "DW..cm.", "julian", "Dec.Hour", "Moon_Phase"))
cor(nums, use="pairwise.complete.obs") #'x' must be numeric error
model <- lm(minutes ~ ., nums)
summary(model)
vif(model) #really high VIF values with DO_sat, so removed from the model
mean(vif(model))
ggpairs(nums, title="correlogram with ggpairs()") 

#----- CNR: GAMM -----#
CNR <- subset(visitations, visitations$AnimalCode == "CNR")
CNRmodel <- uGamm(minutes ~ s(TRANSMITTERID, bs = 're', k=36) + factor(CL_Other) + s(Dec.Hour, bs="cc") + factor(Tide) + 
                    s(Moon_Fraction, bs = "cc") + s(Temp_C) + s(Salinity) + s(DO),
                  data = CNR, method = "REML", family = Gamma) #reran this at 8:15 PM on 9/28, finished at 8:37 PM
CNR1 <- dredge(CNRmodel)
b<- getViz(CNRmodel$gam) #residual plots for gamms
check(b)

CNR_optimal <- uGamm(minutes ~ s(TRANSMITTERID, bs = 're', k=36)+ s(Dec.Hour, bs="cc") + s(Moon_Fraction, bs = "cc"), 
                     data = SER, method = "REML", family = Gamma)

#----- SER model -----#
SER <- subset(visitations, visitations$AnimalCode == "SER")
SERmodel <- uGamm(minutes ~ s(TRANSMITTERID, bs = 're', k=36) + factor(CL_Other) + s(Dec.Hour, bs="cc") + factor(Tide) + 
                    s(Moon_Fraction, bs = "cc") + s(Temp_C) + s(Salinity) + s(DO),
                  data = SER, method = "REML", family = Gamma) 
SER1 <- dredge(SERmodel) 
c<- getViz(SERmodel$gam) #residual plots for gamms
check(c)

SER_optimal <- uGamm(minutes ~ s(TRANSMITTERID, bs = 're', k=36) + factor(CL_Other) + s(Dec.Hour, bs="cc") + factor(Tide) + 
                       s(Moon_Fraction, bs = "cc"), 
                     data = SER, method = "REML", family = Gamma)