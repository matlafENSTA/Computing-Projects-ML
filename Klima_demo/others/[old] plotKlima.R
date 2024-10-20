library(ggplot2)
library(readxl)
library(scales)
library(showtext)
library(thematic)
library(ggokabeito)
library(janitor)
#
# Sets up font, not necessary for running the model.
font_add_google("Fira Sans", "firasans")
showtext_auto()
colors <-  thematic::okabe_ito(8)

theme_add <- theme(text = element_text(family = 'firasans', size = 19),plot.title.position = 'plot',plot.title = element_text(face = 'bold', colour = thematic::okabe_ito(8)[6],))

theme_set(theme_bw()+theme_add)

#
# Read the climate file - contain baseline and the scenario for each climate model.
# Date, Baseline, CNRM_CCLM, .....
#
d <- read_excel("Runoff_long_term_Rinda.xlsx", sheet="Sheet1",na = "-9999")

#For saving files
#
model <- "CNRM_CCLM"
timeper <- "2071-00"

# Select the date column and the baseline
#
base <- d[,c(1,2)]

# Add the year and a filed containing a constant year with day and month (for wide formatting)
#
base <- transform(base,year=as.POSIXlt(Date)$year+1900)
base <- transform(base,mnd=as.POSIXlt(Date)$mon+1)
base <- transform(base,mday=as.POSIXlt(Date)$mday)
base <- transform(base,id=as.Date(paste("1900",mnd,mday,sep="-")))

#Reshape the data (column Baseline, id and year) so that we have one column for each year
#
base_w <- reshape(base[,c("Baseline","id","year")],timevar="year",idvar="id",direction="wide")

#Compute the mean across all years (columns 2 - 31)
#
base_w$middel <- apply(base_w[,2:31], 1, mean,na.rm=TRUE)

#Select the CNRM_CCLM climate model and the date
#
dcn <- d[,c(1,3)]

#Add id fields
#
dcn <- transform(dcn,year=as.POSIXlt(Date)$year+1900)
dcn <- transform(dcn,mnd=as.POSIXlt(Date)$mon+1)
dcn <- transform(dcn,mday=as.POSIXlt(Date)$mday)
dcn <- transform(dcn,id=as.Date(paste("1900",mnd,mday,sep="-")))

# Reshape into wide format and find mean
#
dcn_w <- reshape(dcn[,c("CNRM_CCLM","id","year")],timevar="year",idvar="id",direction="wide")
dcn_w$middel <- apply(dcn_w[,2:31], 1, mean,na.rm=TRUE)

#Plot individual years + baseline mean + CNRM_CCLM mean.
#
pall <- ggplot() + geom_line(data=dcn,aes(x=as.POSIXct(id),y=CNRM_CCLM,group=as.factor(year)),colour="grey80",alpha=0.5,size=0.7,show.legend = FALSE)+
  geom_line(data=dcn_w,aes(x=as.POSIXct(id),y=middel,colour="Mean future"))+
  geom_line(data=base_w,aes(x=as.POSIXct(id),y=middel,colour="Baseline"))+
  scale_colour_manual(values=c("Baseline"=colors[3],"Mean future"=colors[6], "Individual years - future"="grey80"))+
  labs(title=paste("Rinna",model,timeper),x=element_blank(), y=expression("Vassføring (m"^3%.%s^-1*")"))+
  scale_x_datetime(breaks=date_breaks("2 month"),labels = date_format("%b"))+
  theme(legend.title=element_blank(),legend.position="bottom",legend.margin=margin(t=-20))
print(pall)
ggsave(paste("All",model,"_",timeper,".png"),pall,device="png",unit="cm",width=30,height=20,dpi=200)
  
#Plot the realisations for all climate model.
#
dall <- d[,c(1,3:12)]

dall <- transform(dall,year=as.POSIXlt(Date)$year+1900)
dall <- transform(dall,mnd=as.POSIXlt(Date)$mon+1)
dall <- transform(dall,mday=as.POSIXlt(Date)$mday)
dall <- transform(dall,id=as.Date(paste("1900",mnd,mday,sep="-")))

#long format
#
dall_l <- pivot_longer(dall,cols=c(2:11))

#Plot a panel for each model
#
ppan <- ggplot() + geom_line(data=dall_l,aes(x=as.POSIXct(id),y=value,group=as.factor(year)),colour="grey80",alpha=0.5,size=0.7,show.legend = FALSE)+
  labs(title=element_blank(),x=element_blank(), y=expression("Vassføring (m"^3%.%s^-1*")"))+
  scale_x_datetime(breaks=date_breaks("2 month"),labels = date_format("%b"))+
  facet_wrap(.~name)
print(ppan)
#  theme(legend.title=element_blank(),legend.position="bottom",legend.margin=margin(t=-20))
print(pall)
ggsave(paste("All_",timeper,".png"),ppan,device="png",unit="cm",width=30,height=20,dpi=200)

#KALIBRERING HYDROLOGISK MODELL
#

k <- read_excel("RindaCalib.xlsx", sheet="RindaCalib_071",na = "-9999")
k <- clean_names(k)
k <- transform(k,RDate=as.Date(date,format="%d.%m.%Y"))

kq <- k[,c("x1obsrunoff","x1simrunoff","RDate")]

rk <- ggplot(kq) + geom_line(aes(x=RDate,y=x1obsrunoff,colour="Observert"))+
  geom_line(aes(x=RDate,y=x1simrunoff,colour="Simulert"))+
  scale_colour_manual(values=c("Observert"=colors[3],"Simulert"=colors[6]))+
  labs(title="Rinna kalibrert",x=element_blank(), y="Avrenning (mm)")+
  annotate("text", x=as.Date("2014-01-01"), y=42, label= "NSE = 0.71", size=4)+
  theme(legend.title=element_blank(),legend.position="bottom",legend.margin=margin(t=-20))
plot(rk)
ggsave("RinnaKalibrert.png",rk,device="png",unit="cm",width=30,height=20,dpi=200)
