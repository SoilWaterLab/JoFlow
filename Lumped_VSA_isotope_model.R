#IsoJoFlow hydrologic & stable water isotope model
#  - hydrology based on JoFlow (Archibald et al., 2014)
#  - snowmelt based on Walter (2005)
#  - PET from temperature model based on Priestly-Taylor approximation
#  - modified baseflow equation per Rossman et al (2009)
#  - stable water isotope fractionation of O (Gibson, 2002; Smith et al. 2016)
#  - plant growth model (Allen et al., 1998; 2005)
#  - dual crop coefficient partitioning of ET (Allen et al. 1998; 2005)
#
#Last Modified: JK 05/25/2016

#Model Inputs
#  - Parameters: daily precipitation, daily min temperature, daily max temperature, del O precip, del H precip, Relative humidity
#  - Units: precipitation (mm), temperatures (degrees C), del O and H (Vienna standard), Relative humidity (fraction)
#  - DataType: Date (integer), PRCP (double), TMAX (double), TMIN (double), del_O (double), del_H (double), RelHum (double)

precip= read.csv('C:/Isotope/CZO/Final/PrecipModelForcing_CKM.csv')

#R Libraries
setwd('C:/Isotope/CZO/')
library(EcoHydRology)
library(zoo)
source('IsoSnowMelt.R')

#1. Set Model Parameters
#Simulation Settings
outpath = 'C:/Isotope/CZO/Final/output/'
record = 1096 #Simulation length (days)

#Watershed model inputs
Depth = 750 # Average soil depth in watershed (mm)
SATper = 0.3 # Porosity (fraction)
AWCper = 0.13 # Available Water Capacity (AWC) (fraction)
no_wet_class = 1 # The number of wetness classes for saturated area designation (not currently used)
percentImpervious = 1 # Percent of the watershed that is impervious
Tp = 2 # Time to peak (hours)
latitudeDegrees = 40.67 #decimal degrees
albedo = 0.18 # Average albedo
StartCond = "avg" # Watershed conditions before first day of run ("wet", "dry", "avg")
SW1 = NULL # Soil water on the first day (depth, mm)
#BF1 = 10 #  mm/day can use nearby watershed baseflow, ConvertFlowUnits(cfs,WA=W_Area)
PETcap = 5 # Does not let PET get larger than this cut off (mm)
A1 = 0.01 #Baseflow recession coefficient (Rossman, 2009)
B = 2 #Baseflow recession exponent  (Rossman, 2009)
GW_min = 1000 #min GW before BF starts
Se_min = 78 #mm
C1 = 3.1
Ia_coef = 0.2 # range ~ 0.05 - 0.2
TEW = 25 #total evaporable water, based on Allen et al (2006) Table 19
REW = 15 #readily evaporable water, basde on Allen et al (2006) Table 19

#Isotope model inputs
standard = 2005/1000000 #Vienna standard
R_Soil_Init_O = 0.00198 #Initial R for unsaturated zone storage Oxygen
R_Soil_Init_H = 0.00196 #Initial R for unsaturated zone storage Hydrogen
R_GW_Init_O = 0.001985953 #Initial R for unsaturated zone storage Oxygen
R_GW_Init_H = 0.001890 #Initial R for unsaturated zone storage Hydrogen
E_ET_ratio = 2 #Ratio of evaporation to total ET (0,1) or 2 
#note: use E_ET_ratio = 2 to use dual crop coefficient model
mixing_model = 2 #unsaturated zone mixing model: 
#1 - no mixing (macro pore concept)
#2 - complete mixing
#3 - partial mixing (a weighting of the above models)
mew = 0.5 #weighting factor for partial mixing model (0,1)
SW_Volume = 0.2 #soil water contribution to baseflow (mm)

#Plant growth model inputs
Tb = 1 #base temperature for computing growing degree days (degrees C)
GDDmax = 2500 #maximum growing degree days in annual plant cycle
GDD_1 = 0*GDDmax #stage 1 growth
GDD_2 = 0.1*GDDmax #stage 2 growth
GDD_3 = 0.225*GDDmax #stage 3 growth
GDD_4 = 0.9*GDDmax #stage 4 growth
GDD_5 = 1*GDDmax #stage 5 growth
Kcmin = 0.25 #minimum crop coefficient (Kc) for deciduous forest
Kcmax = 1.05 #max crop coefficient (Kc)
Kcbini = 0.4 #initial basal crop coefficient
Kcbmid = 1 #mid basal crop coefficient
Kcbend = 0.3 #end basal crop coefficient
few = 0.9

#Constants to calculate
AWC = Depth*AWCper # AWC as a depth (mm)
SAT = Depth*SATper # porosity as a depth (mm)
runoff_breakdown = RunoffBreakdown(Tp, HrPrcDelay = (Tp/2-4)) 
latitude<-latitudeDegrees*pi/180 ## latitude in radians

#Clean up NCDC meteorological data
precip$PRCP[precip$PRCP == -9999] = 0

for (j in 1:length(precip$TMAX)) #Note: Some NCDC datasets have errors where TMIN > TMAX
{
  if (precip$TMIN[j] == -9999)
  {precip$TMIN[j] = precip$TMIN[j-1]}
  
  if (precip$TMIN[j] == -9999)
  {precip$TMAX[j] = precip$TMAX[j-1]}
  
  if (precip$TMIN[j]>precip$TMAX[j])
  {precip$TMAX[j] = precip$TMIN[j] + 1}
}

h = precip$RelHum

#Create date column for snowmelt model
precip$DATE=as.Date.numeric(precip$Date, origin="2009-01-01")

#3. Snowmelt Model (based on Walter et al (2005))
#Note: Accumulated snowpack del same as precip, snowpack melt determined from IsoSnowMelt
Ri_O = (precip$del_O/1000 + 1)*standard
Ri_O = Ri_O/(1 + Ri_O)
Ri_H = (precip$del_H/1000 + 1)*standard
Ri_H = Ri_H/(1 + Ri_H)
snowmelt=IsoSnowMelt(Date=precip$DATE, precip_mm=precip$PRCP, R_precip_O=Ri_O, R_precip_H=Ri_H, Tmax_C=precip$TMAX, Tmin_C=precip$TMIN, lat_deg=latitudeDegrees)

#Assign precipitation isotopic values
P <- (snowmelt$SnowMelt_mm[1:record]+snowmelt$Rain_mm[1:record])
Tmin <- precip$TMIN[1:record]
Tmax <- precip$TMAX[1:record]
Tav <- 0.5*(Tmin + Tmax)
P_heavy_O <- snowmelt$SnowMelt_mm*snowmelt$R_O + snowmelt$Rain_mm*Ri_O
P_light_O <- snowmelt$SnowMelt_mm*(1 - snowmelt$R_O) + snowmelt$Rain_mm*(1 - Ri_O)
del_p_O = 1000*((P_heavy_O/P_light_O)/standard - 1)

P_heavy_H <- snowmelt$SnowMelt_mm*snowmelt$R_H + snowmelt$Rain_mm*Ri_H
P_light_H <- snowmelt$SnowMelt_mm*(1 - snowmelt$R_H) + snowmelt$Rain_mm*(1 - Ri_H)
del_p_H = 1000*((P_heavy_H/P_light_H)/standard - 1)

#4. PET Model
# Create Julian Date format for use in PET model
x=as.POSIXlt(precip$DATE[1:record], format="Y%b%d")
day=x$yday+1

#Run PET from daily temp model (Priestly Taylor) 
PET<-PET_fromTemp(Jday=day,Tmax_C=Tmax,Tmin_C=Tmin,albedo=albedo,lat_radians=latitude)*10000## mm (Priestley-Taylor)
PET[which(PET>PETcap)]<-PETcap  #Sets a cap on PET estimates (usually ~ 5mm)

#5. Plant growth model to determine Kc
PlantModel = precip[1:record,]
PlantModel$Year = as.numeric(format(PlantModel$DATE,'%Y'))
PlantModel$Month = as.numeric(format(PlantModel$DATE,'%m'))
PlantModel$Avg_Temp = 0.5*(PlantModel$TMIN + PlantModel$TMAX)
PlantModel$DegreeDays = PlantModel$Avg_Temp - Tb
PlantModel[PlantModel$DegreeDays < 0,9] = 0
PlantModel$Accumulate_GDD = 0
PlantModel$DegreeDays_percent = 0

for (j in 2:nrow(PlantModel))
{
  #Calculate Growing Degree Days (GDD)
  PlantModel$Accumulate_GDD[j] = PlantModel$Accumulate_GDD[j-1] + PlantModel$DegreeDays[j]
  if (PlantModel$Year[j] > PlantModel$Year[j-1])
  {PlantModel$Accumulate_GDD[j] = PlantModel$DegreeDays[j]}
  #Calculate GDD as percent of GDDmax
  PlantModel$DegreeDays_percent[j] = 100*min(PlantModel$Accumulate_GDD[j] / GDDmax,1)
}

#Calculate alpha and daily Kcb
PlantModel$Kcb <- NA
for (j in 1:nrow(PlantModel))
{
  if (PlantModel$Accumulate_GDD[j] <= GDD_2)
     {PlantModel$Kcb[j] = Kcbini
    PlantModel$alphaK[j] = PlantModel$Accumulate_GDD[j] / GDD_4}
  if (PlantModel$Accumulate_GDD[j] > GDD_2 & PlantModel$Accumulate_GDD[j] <= GDD_3)
     {PlantModel$alphaK[j] = (GDD_2/GDD_4) + ((GDD_4 - GDD_2)/(GDD_3-GDD_2))*(PlantModel$Accumulate_GDD[j] - GDD_2)/GDD_4}
  if (PlantModel$Accumulate_GDD[j] > GDD_3 & PlantModel$Accumulate_GDD[j] <= GDD_4)
     {PlantModel$Kcb[j] = Kcbmid
    PlantModel$alphaK[j] = 1}
  if (PlantModel$Accumulate_GDD[j] > GDD_4 & PlantModel$Accumulate_GDD[j] < GDD_5)
     {PlantModel$alphaK[j] = 1 - 0.6*((PlantModel$Accumulate_GDD[j] - GDD_4)/(GDDmax - GDD_4))}
  if (PlantModel$Accumulate_GDD[j] > GDD_5)
     {PlantModel$Kcb[j] =  0.15 #set equal to minimum non-zero value per Allen et al (2005)
     PlantModel$alphaK[j] = 0}
}
PlantModel$Kcb <- na.approx(PlantModel$Kcb)

#Calculate Kc from alpha
PlantModel$Kc = Kcmin + PlantModel$alphaK*(Kcmax - Kcmin)

#Estimate AET from single crop (Kc) coefficient for impervious area calculation
ETo = PET*PlantModel$Kc
ET <- ETo    #  Initializing modeled actual ET, equal to ETo when deltaP > 0, but less than ETo on drier days

#6.Watershed Isotope Model
#Effective soil water storage coefficients, eswsc = Se*(2..... see Schneiderman et al 2007)
eswsc <- vector(mode="numeric", length=no_wet_class)
eqn16 <- function(x){(2*(sqrt(1-(1/no_wet_class)*(x-1))-sqrt(1-(1/no_wet_class)*x))/(1/no_wet_class))-1}
x <- seq(1,no_wet_class, by=1)
eswsc <- eqn16(x)

# Initializing vectors for Water Budget Loop, and Runoff estimation
Pe<-vector(length=length(P))# Effective Precipitation (non-impervious areas)

Ratio_O <- vector(length=length(P)) #mixing model ratio
Ratio_H <- vector(length=length(P)) #mixing model ratio
del_a_O <- vector(length=length(P)) #del atmospheric
del_a_H <- vector(length=length(P)) #del atmospheric
alpha_star <- vector(length=length(P))
epsilon <- vector(length=length(P))
del_e_O <- vector(length=length(P)) #del evaporation
del_e_H <- vector(length=length(P)) #del evaporation
R_evap_O<- vector(length=length(P)) #R of evaporation
R_evap_H<- vector(length=length(P)) #R of evaporation

PET_Evap <- vector(length=length(P))  
PET_Plant <- vector(length=length(P))
PET_Evap_heavy_O <- vector(length=length(P))
PET_Evap_light_O <- vector(length=length(P))  
PET_Plant_heavy_O <- vector(length=length(P))
PET_Plant_light_O <- vector(length=length(P))  
Evap_light_O <- vector(length=length(P))
Evap_heavy_O <- vector(length=length(P))
PET_Evap_heavy_H <- vector(length=length(P))
PET_Evap_light_H <- vector(length=length(P))  
PET_Plant_heavy_H <- vector(length=length(P))
PET_Plant_light_H <- vector(length=length(P))  
Evap_light_H <- vector(length=length(P))
Evap_heavy_H <- vector(length=length(P))
Revap <- vector(length=length(P))

SoilWater <- vector(length=length(P))##(mm)
SoilWater_heavy_O <-vector(length=length(P))##(mm)
SoilWater_light_O <-vector(length=length(P))##(mm)
SoilWater_heavy_H <-vector(length=length(P))##(mm)
SoilWater_light_H <-vector(length=length(P))##(mm)

del_SoilWater_O <- vector(length=length(P))
del_Transpiration_O <- vector(length=length(P))
del_TM_S_O<-vector(length=length(P))
del_baseflow_O<-vector(length=length(P))
del_totQ_O<-vector(length=length(P))
del_SoilWater_H <- vector(length=length(P))
del_Transpiration_H <- vector(length=length(P))
del_TM_S_H<-vector(length=length(P))
del_baseflow_H<-vector(length=length(P))
del_totQ_H<-vector(length=length(P))

baseflow<-vector(length=length(P))
baseflow_heavy_O <- vector(length=length(P))
baseflow_light_O <- vector(length=length(P))
baseflow_heavy_H <- vector(length=length(P))
baseflow_light_H <- vector(length=length(P))

excess<-vector(length=length(P))
excess_heavy_O<-vector(length=length(P))
excess_light_O<-vector(length=length(P))
excess_heavy_H<-vector(length=length(P))
excess_light_H<-vector(length=length(P))

TM_S<-vector(length=length(P))##This is the daily time-step T-M storage, used to calc baseflow
TM_S_heavy_O<-vector(length=length(P))
TM_S_light_O<-vector(length=length(P))
TM_S_heavy_H<-vector(length=length(P))
TM_S_light_H<-vector(length=length(P))
GW_Ratio_O<-vector(length=length(P))
GW_Ratio_H<-vector(length=length(P))
SW_Ratio_O<-vector(length=length(P))
SW_Ratio_H<-vector(length=length(P))

totQ<-vector(length=length(P))
totQ_heavy_H<-vector(length=length(P))
totQ_light_H<-vector(length=length(P))
totQ_heavy_O<-vector(length=length(P))
totQ_light_O<-vector(length=length(P))

Se<-vector(length=length(P))
Q2<-vector(length=length(P))

MaxWetClass<-vector(length=length(P))
sigmaS<-matrix(nrow=length(P), ncol=no_wet_class)  # Storage in each wetness class
runoff_by_class<-matrix(nrow=length(P), ncol=no_wet_class)

impervRunoff<-vector(length=length(P))
OverlandFlow<-vector(length=length(totQ))
OverlandFlow_heavy_O<-vector(length=length(totQ))
OverlandFlow_light_O<-vector(length=length(totQ))
OverlandFlow_heavy_H<-vector(length=length(totQ))
OverlandFlow_light_H<-vector(length=length(totQ))

ShallowInterflow<-vector(length=length(totQ))
ShallowInterflow_heavy_O<-vector(length=length(totQ))
ShallowInterflow_light_O<-vector(length=length(totQ))
ShallowInterflow_heavy_H<-vector(length=length(totQ))
ShallowInterflow_light_H<-vector(length=length(totQ))

#Imervious Area CN Runoff
deltaP <- P - ETo## (mm) neg values indicate net evap
impervIa<-Ia_coef * 5  ##Se = 5 mm for impervious areas (CN=98)
impervPe<-deltaP - impervIa##  Effective precipitation on impervious areas
impervPe[which(impervPe < 0)]<-0

#Model initial conditions
#these values generally do not affect the results, particularly if a spin-up year is used
del_p_O[1] = -10  #just setting this to something, doesn't affect results
del_p_H[1] = -10  #just setting this to something, doesn't affect results

TM_S[1]<-GW_min + 50
TM_S_heavy_O[1] = R_GW_Init_O*TM_S[1]
TM_S_light_O[1] = (1-R_GW_Init_O)*TM_S[1]
TM_S_heavy_H[1] = R_GW_Init_H*TM_S[1]
TM_S_light_H[1] = (1-R_GW_Init_H)*TM_S[1]

# Initialize soil wetness based on initial conditions ("wet", "dry", "avg") and initial isotopic ratio
SoilWater_heavy_O[1] <- ifelse(StartCond == "wet", AWC, ifelse(StartCond == "dry", R_Soil_Init_O*0.2*AWC, R_Soil_Init_O*0.75*AWC))
SoilWater_light_O[1] <- ifelse(StartCond == "wet", AWC, ifelse(StartCond == "dry", (1-R_Soil_Init_O)*0.2*AWC, (1-R_Soil_Init_O)*0.75*AWC))
SoilWater_heavy_H[1] <- ifelse(StartCond == "wet", AWC, ifelse(StartCond == "dry", R_Soil_Init_H*0.2*AWC, R_Soil_Init_H*0.75*AWC))
SoilWater_light_H[1] <- ifelse(StartCond == "wet", AWC, ifelse(StartCond == "dry", (1-R_Soil_Init_H)*0.2*AWC, (1-R_Soil_Init_H)*0.75*AWC))
SoilWater[1] = SoilWater_heavy_O[1] + SoilWater_light_O[1]

excess_heavy_O[1] <- 0# Assume no runoff from the day before
excess_light_O[1] <- 0# Assume no runoff from the day before
excess_heavy_H[1] <- 0# Assume no runoff from the day before
excess_light_H[1] <- 0# Assume no runoff from the day before

Se[1]<- Se_min+C1*(AWC-SoilWater[1])
sigmaS[1,] <- eswsc*Se[1]

totQ[1] <- Pe[1]*Pe[1]/(Pe[1]+Se[1])
totQ_heavy_O[1] = totQ[1]*Ri_O[1]
totQ_light_O[1] = totQ[1]*(1-Ri_O[1])
totQ_heavy_H[1] = totQ[1]*Ri_H[1]
totQ_light_H[1] = totQ[1]*(1-Ri_H[1])

Ia<-Ia_coef*Se  # Initial abstraction = depth of precip before runoff starts

PlantModel$De[1] = 0
PlantModel$Kr[1] = 0
PlantModel$Ke[1] = PlantModel$Kc[1]

#Thornthwaite-Mather Function and Runoff Generation (Pervious areas)
for(i in (2):length(deltaP))
{
  #Solve net P and daily AET from dual crop coefficient model for pervious areas
  deltaP[i] = P[i] - PET[i]*(PlantModel$Ke[i-1] + PlantModel$Kcb[i-1])
  
  #Calculate time step runoff with CN method
  if (deltaP[i]-Ia[i-1]>0)
  {Pe[i] <- deltaP[i]-Ia[i-1]} 
  else 
  {Pe[i]<-0}
  
  totQ[i] <- Pe[i]*Pe[i]/(Pe[i]+Se[i-1]) ## Effective storage is from previous day
  totQ_heavy_O[i] <- totQ[i]*Ri_O[i]
  totQ_light_O[i] <- totQ[i]*(1-Ri_O[i])
  totQ_heavy_H[i] <- totQ[i]*Ri_H[i]
  totQ_light_H[i] <- totQ[i]*(1-Ri_H[i])
  
  #Partition PET into transpiration and evaporation
  if (E_ET_ratio > 1)
  {
      PET_Evap[i] = PET[i]*PlantModel$Ke[i-1]
      PET_Plant[i] = PET[i]*PlantModel$Kcb[i-1]
  } else
  {
    PET_Evap[i] = E_ET_ratio*PET[i]
    PET_Plant[i] = (1 - E_ET_ratio)*PET[i]
  }
  #Determine how much transpiration pulls from each isotope pool
  #Assumes uniform vertical distribution of 18O stored in evaporable soil matrix
  PET_Plant_heavy_O[i] = max(0,PET_Plant[i]*SoilWater_heavy_O[i-1]/SoilWater[i-1])
  PET_Plant_light_O[i] = max(0,PET_Plant[i]*SoilWater_light_O[i-1]/SoilWater[i-1])
  if (is.nan(PET_Plant_light_O[i]))  {PET_Plant_light_O[i] = 0}
  if (is.nan(PET_Plant_heavy_O[i])) {PET_Plant_heavy_O[i] = 0}
  PET_Plant_heavy_H[i] = max(0,PET_Plant[i]*SoilWater_heavy_H[i-1]/SoilWater[i-1])
  PET_Plant_light_H[i] = max(0,PET_Plant[i]*SoilWater_light_H[i-1]/SoilWater[i-1])
  if (is.nan(PET_Plant_light_H[i]))  {PET_Plant_light_H[i] = 0}
  if (is.nan(PET_Plant_heavy_H[i])) {PET_Plant_heavy_H[i] = 0}
  
  #Determine evaporation fractionation (from Gibson, 2002)
  Ta = 0.5*(snowmelt$MaxT_C[i] + snowmelt$MinT_C[i])
  del_SoilWater_O[i-1] = ((SoilWater_heavy_O[i-1]/SoilWater_light_O[i-1])/standard - 1)*1000
  del_SoilWater_H[i-1] = ((SoilWater_heavy_H[i-1]/SoilWater_light_H[i-1])/standard - 1)*1000
  #If no rain today, assume precipitation del 18O same as day before
  if (is.nan(del_p_O[i]))  {del_p_O[i] = del_p_O[i-1]}
  if (is.nan(del_p_H[i]))  {del_p_H[i] = del_p_H[i-1]}
  #Determine liquid-vapor equilibrium (Gibson, 2002; Smith et al., 2016)
  term1 = exp(6.7123 *((10^3)/(Ta + 273.15))/1000)
  term2 = exp(0.35041*((10^9)/((Ta + 273.15)^3))/1000)
  term3 = exp(7.685/1000)
  term4 = exp(1.6664*((10^6)/((Ta + 273.15)^2)/1000))
  alpha_star[i] = (term1*term2)/(term3*term4)
  #Equilibrium and kinetic separation
  e_star = (alpha_star[i] - 1)*1000
  e_k = 1.047
  epsilon[i] = e_star + e_k
  #Estimate atmospheric del 18O
  del_a_O[i] = (del_p_O[i] - e_star)/alpha_star[i]
  del_a_H[i] = (del_p_H[i] - e_star)/alpha_star[i]
  #Estimate evaporation del 18O
  del_e_O[i] = (del_SoilWater_O[i-1] - h[i]*del_a_O[i] - epsilon[i])/(1-h[i]+e_k/1000)
  del_e_H[i] = (del_SoilWater_H[i-1] - h[i]*del_a_H[i] - epsilon[i])/(1-h[i]+e_k/1000)
  R_evap_O[i] = (del_e_O[i]/1000 + 1)*standard
  R_evap_H[i] = (del_e_H[i]/1000 + 1)*standard
  R_evap_total_O = R_evap_O[i]/(1 + R_evap_O[i])
  R_evap_total_H = R_evap_H[i]/(1 + R_evap_H[i])
  
  #PET partitioning based on previous time step, dual crop coefficient model (Allen et al., 1998)
  PET_Evap_heavy_O[i] = PET_Evap[i]*(R_evap_total_O)
  PET_Evap_light_O[i] = PET_Evap[i]*(1 - R_evap_total_O)
  PET_Evap_heavy_H[i] = PET_Evap[i]*(R_evap_total_H)
  PET_Evap_light_H[i] = PET_Evap[i]*(1 - R_evap_total_H)
  
  # Water Balance
  #a. DRYING CONDITION
  if(deltaP[i]<=0)
  {
      
      #Determine fractionation evaporation
      if (deltaP[i] != 0)
      {
            SoilWater[i] <- SoilWater[i-1]*exp(deltaP[i]/AWC)
            Losses = SoilWater[i-1] - SoilWater[i]
            Heavy_losses_O = Losses*(PET_Evap_heavy_O[i] + PET_Plant_heavy_O[i])/(PET_Evap_heavy_O[i] + PET_Plant_heavy_O[i] + PET_Evap_light_O[i] + PET_Plant_light_O[i])
            Light_losses_O = Losses - Heavy_losses_O
            Heavy_losses_H = Losses*(PET_Evap_heavy_H[i] + PET_Plant_heavy_H[i])/(PET_Evap_heavy_H[i] + PET_Plant_heavy_H[i] + PET_Evap_light_H[i] + PET_Plant_light_H[i])
            Light_losses_H = Losses - Heavy_losses_H
            
            SoilWater_light_O[i] = SoilWater_light_O[i-1] - Light_losses_O
            SoilWater_heavy_O[i] = SoilWater_heavy_O[i-1] - Heavy_losses_O
          
            SoilWater_light_H[i] = SoilWater_light_H[i-1] - Light_losses_H
            SoilWater_heavy_H[i] = SoilWater_heavy_H[i-1] - Heavy_losses_H

      }
      else
      {
        SoilWater[i] = SoilWater[i-1]
        SoilWater_light_O[i] = SoilWater_light_O[i-1]
        SoilWater_heavy_O[i] = SoilWater_heavy_O[i-1]
      
        SoilWater_light_H[i] = SoilWater_light_H[i-1]
        SoilWater_heavy_H[i] = SoilWater_heavy_H[i-1]
      }

      excess[i] <- 0
      excess_heavy_O[i] <- 0
      excess_light_O[i] <- 0
      excess_heavy_H[i] <- 0
      excess_light_H[i] <- 0
      
      #Cummulative evaporation losses (Allen et al., 2006)
      PlantModel$De[i] = PlantModel$De[i-1] - deltaP[i]
  }
  
  #b. WETTING BELOW FIELD CAPACITY
  else if (deltaP[i]+SoilWater[i-1] <= AWC)
  {  
      SoilWater_light_O[i] <- SoilWater_light_O[i-1] + P_light_O[i] - PET_Evap_light_O[i] - PET_Plant_light_O[i]
      SoilWater_heavy_O[i] <- SoilWater_heavy_O[i-1] + P_heavy_O[i] - PET_Evap_heavy_O[i] - PET_Plant_heavy_O[i]
      SoilWater_light_H[i] <- SoilWater_light_H[i-1] + P_light_H[i] - PET_Evap_light_H[i] - PET_Plant_light_H[i]
      SoilWater_heavy_H[i] <- SoilWater_heavy_H[i-1] + P_heavy_H[i] - PET_Evap_heavy_H[i] - PET_Plant_heavy_H[i]
      
      SoilWater[i]<- SoilWater_heavy_O[i] + SoilWater_light_O[i]
      excess[i] <- 0
      excess_heavy_O[i] <- 0
      excess_light_O[i] <- 0
      excess_heavy_H[i] <- 0
      excess_light_H[i] <- 0
      
      PlantModel$De[i] = min(0,PlantModel$De[i] - deltaP[i])
  }
  
  #c. WETTING ABOVE FIELD CAPACITY
  else { 
    
    SoilWater[i] <- AWC  ## So overall SW cannot exceed field capacity (AWC)

    #Determine total volume of percolated water
    excess[i] <- deltaP[i]+SoilWater[i-1]-AWC
    
    #Determine Isotopic Ratio of the Percolated Water
    #Model 1: No Mixing
    if (mixing_model == 1)
    {Ratio_O[i] = Ri_O[i]
     Ratio_H[i] = Ri_H[i]}
    
    #Model 2: Complete Mixing
    if (mixing_model == 2)
    {Ratio_O[i] = (SoilWater_heavy_O[i-1] +  P_heavy_O[i])/(SoilWater[i-1] + P_light_O[i] + P_heavy_O[i])
     Ratio_H[i] = (SoilWater_heavy_H[i-1] +  P_heavy_H[i])/(SoilWater[i-1] + P_light_H[i] + P_heavy_H[i])}  
    
    #Model 3: Partial Mixing       
    if (mixing_model == 3)
    {
      #Ratio of heavy isotope for no mixing
      Ratio_nomix_O = Ri_O[i]
      Ratio_nomix_H = Ri_H[i]
      #Ratio of heavy isotope for complete mixing
      Ratio_complete_O = (SoilWater_heavy_O[i-1] +  P_heavy_O[i])/(SoilWater[i-1] + P_light_O[i] + P_heavy_O[i])                
      Ratio_complete_H = (SoilWater_heavy_H[i-1] +  P_heavy_H[i])/(SoilWater[i-1] + P_light_H[i] + P_heavy_H[i])                
      #Ratio of heavy isotope for partial mixing    
      Ratio_O[i] = mew*Ratio_complete_O + (1-mew)*Ratio_nomix_O
      Ratio_H[i] = mew*Ratio_complete_H + (1-mew)*Ratio_nomix_H
    }
    
    #Compute new excess from the excess loss ratio
    excess_heavy_O[i] <- Ratio_O[i]*excess[i]
    excess_light_O[i] <- (1-Ratio_O[i])*excess[i]
    excess_heavy_H[i] <- Ratio_H[i]*excess[i]
    excess_light_H[i] <- (1-Ratio_H[i])*excess[i]
    
    #Compute new unsaturated zone storage from the excess loss ratio
    SoilWater_heavy_O[i] <- SoilWater_heavy_O[i-1] - excess_heavy_O[i] + deltaP[i]*Ri_O[i]
    SoilWater_light_O[i] <- SoilWater_light_O[i-1] - excess_light_O[i] + deltaP[i]*(1 - Ri_O[i])
    SoilWater_heavy_H[i] <- SoilWater_heavy_H[i-1] - excess_heavy_H[i] + deltaP[i]*Ri_H[i]
    SoilWater_light_H[i] <- SoilWater_light_H[i-1] - excess_light_H[i] + deltaP[i]*(1 - Ri_H[i])

    #reset cummulative evaporation losses under saturated conditions
    PlantModel$De[i] = 0
  }
  
  #Solve evaporation loss coefficient, Ke (Allen et al 1998)
  PlantModel$Kr[i] = 0.5*max(0,(TEW - PlantModel$De[i]))/(TEW - REW)
  PlantModel$Ke[i] = min(PlantModel$Kr[i]*(max(PlantModel$Kc) - PlantModel$Kcb[i]), few*max(PlantModel$Kc))
  
  #Compute baseflow
  #Baseflow modified to use an exponential, slightly better fit for summer months (JK 05/04/2016)
  if (TM_S[i-1] > GW_min)
  {baseflow[i] = min(A1*((TM_S[i-1] - GW_min)^B),TM_S[i-1]-1)}
  if (TM_S[i-1] <= GW_min)
  {baseflow[i] = 0}
  
  GW_Ratio_O[i] = TM_S_heavy_O[i-1] / TM_S[i-1]
  GW_Ratio_H[i] = TM_S_heavy_H[i-1] / TM_S[i-1]

  SW_Ratio_O[i] = SoilWater_heavy_O[i] / SoilWater[i]
  SW_Ratio_H[i] = SoilWater_heavy_H[i] / SoilWater[i]

  baseflow_heavy_O[i] <- baseflow[i]*(GW_Ratio_O[i]*baseflow[i] + SW_Ratio_O[i]*SW_Volume)/(baseflow[i] + SW_Volume)
  baseflow_light_O[i] <- baseflow[i]*((1-GW_Ratio_O[i])*baseflow[i] + (1 - SW_Ratio_O[i])*SW_Volume)/(baseflow[i] + SW_Volume)
  baseflow_heavy_H[i] <- baseflow[i]*(GW_Ratio_H[i]*baseflow[i] + SW_Ratio_H[i]*SW_Volume)/(baseflow[i] + SW_Volume)
  baseflow_light_H[i] <- baseflow[i]*((1-GW_Ratio_H[i])*baseflow[i] + (1 - SW_Ratio_H[i])*SW_Volume)/(baseflow[i] + SW_Volume)

  # Ensure mass-balance, since curve number is empirical
  if ((excess[i])>=totQ[i])
  {
    TM_S[i]<-TM_S[i-1] - baseflow[i] + excess[i] - totQ[i]
    TM_S_heavy_O[i] <- TM_S_heavy_O[i-1] - baseflow_heavy_O[i] + excess_heavy_O[i] - totQ_heavy_O[i]
    TM_S_light_O[i] <- TM_S_light_O[i-1] - baseflow_light_O[i] + excess_light_O[i] - totQ_light_O[i]
    TM_S_heavy_H[i] <- TM_S_heavy_H[i-1] - baseflow_heavy_H[i] + excess_heavy_H[i] - totQ_heavy_H[i]
    TM_S_light_H[i] <- TM_S_light_H[i-1] - baseflow_light_H[i] + excess_light_H[i] - totQ_light_H[i]
  } else 
  {
    TM_S[i]<-TM_S[i-1] - baseflow[i]
    TM_S_heavy_O[i]<-TM_S_heavy_O[i-1] - baseflow_heavy_O[i]
    TM_S_light_O[i]<-TM_S_light_O[i-1] - baseflow_light_O[i]
    TM_S_heavy_H[i]<-TM_S_heavy_H[i-1] - baseflow_heavy_H[i]
    TM_S_light_H[i]<-TM_S_light_H[i-1] - baseflow_light_H[i]
    SoilWater[i]<-SoilWater[i]-(totQ[i]-excess[i])
    SoilWater_heavy_O[i]<-SoilWater_heavy_O[i]-(totQ_heavy_O[i]-excess_heavy_O[i])
    SoilWater_light_O[i]<-SoilWater_light_O[i]-(totQ_light_O[i]-excess_light_O[i])
    SoilWater_heavy_H[i]<-SoilWater_heavy_H[i]-(totQ_heavy_H[i]-excess_heavy_H[i])
    SoilWater_light_H[i]<-SoilWater_light_H[i]-(totQ_light_H[i]-excess_light_H[i])
  }
  
  #Update initial abstraction for next time step
  Se[i] <- Se_min+C1*(AWC-SoilWater[i])
  sigmaS[i,]<- eswsc*Se[i]  ## Amount of storage available in each wetness class [mm]
  Ia[i]<-Ia_coef*Se[i]
  
  #Impervious Runoff
  impervRunoff[i]<-impervPe[i] 

  #Convert R to delta notation for isotopic state variables
  del_Transpiration_O[i] = ((PET_Plant_heavy_O[i]/PET_Plant_light_O[i])/standard-1)*1000
  del_TM_S_O[i] = ((TM_S_heavy_O[i]/TM_S_light_O[i])/standard - 1)*1000
  del_baseflow_O[i] = ((baseflow_heavy_O[i]/baseflow_light_O[i])/standard - 1)*1000
  del_totQ_O[i] = ((totQ_heavy_O[i]/totQ_light_O[i])/standard - 1)*1000
  del_Transpiration_H[i] = ((PET_Plant_heavy_H[i]/PET_Plant_light_H[i])/standard-1)*1000
  del_TM_S_H[i] = ((TM_S_heavy_H[i]/TM_S_light_H[i])/standard - 1)*1000
  del_baseflow_H[i] = ((baseflow_heavy_H[i]/baseflow_light_H[i])/standard - 1)*1000
  del_totQ_H[i] = ((totQ_heavy_H[i]/totQ_light_H[i])/standard - 1)*1000
} 

#Partition CN runoff into direct runoff and shallow interflow (Archibald et al 2014)
# Use the coefficients generated from time of concentration
runoff_breakdown[which(runoff_breakdown < 0.01)] <- 0
OverlandFlow[1:4]<-totQ[1:4]*runoff_breakdown[1]
ShallowInterflow[1]<-0# Assuming no runoff in previous days before this model run
ShallowInterflow[2]<-totQ[1]*runoff_breakdown[2]
ShallowInterflow[3]<-totQ[1]*runoff_breakdown[3] + totQ[2]*runoff_breakdown[2]
ShallowInterflow[4]<-totQ[1]*runoff_breakdown[4] + totQ[2]*runoff_breakdown[3] + totQ[3]*runoff_breakdown[2]

for (i in 5:length(totQ))
  {
    OverlandFlow[i]<- totQ[i]*runoff_breakdown[1]
    ShallowInterflow[i]<- totQ[i-4]*runoff_breakdown[5]+totQ[i-3]*runoff_breakdown[4]+totQ[i-2]*runoff_breakdown[3]+totQ[i-1]*runoff_breakdown[2]
  }

#Post Process Model Outputs
Streamflow <- (OverlandFlow+ShallowInterflow)*(1-0.01*percentImpervious)+impervRunoff*(0.01*percentImpervious)+baseflow#totQ+baseflow
Streamflow_heavy_O = Ri_O*OverlandFlow + (SoilWater_heavy_O[i]/SoilWater[i])*ShallowInterflow + baseflow_heavy_O
Streamflow_light_O = (1-Ri_O)*OverlandFlow + (1 - SoilWater_heavy_O[i]/SoilWater[i])*ShallowInterflow + baseflow_light_O
del_Streamflow_O = ((Streamflow_heavy_O/Streamflow_light_O)/standard-1)*1000
Streamflow_heavy_H = Ri_H*OverlandFlow + (SoilWater_heavy_H[i]/SoilWater[i])*ShallowInterflow + baseflow_heavy_H
Streamflow_light_H = (1-Ri_H)*OverlandFlow + (1 - SoilWater_heavy_H[i]/SoilWater[i])*ShallowInterflow + baseflow_light_H
del_Streamflow_H = ((Streamflow_heavy_H/Streamflow_light_H)/standard-1)*1000

Quickflow <- (OverlandFlow+ShallowInterflow)*(1-0.01*percentImpervious)+impervRunoff*(0.01*percentImpervious)
Quickflow_heavy_O <- Ri_O*(OverlandFlow + ShallowInterflow)
Quickflow_light_O <- (1-Ri_O)*(OverlandFlow + ShallowInterflow)
del_Quickflow_O = ((Quickflow_heavy_O/Quickflow_light_O)/standard-1)*1000
Quickflow_heavy_H <- Ri_H*(OverlandFlow + ShallowInterflow)
Quickflow_light_H <- (1-Ri_H)*(OverlandFlow + ShallowInterflow)
del_Quickflow_H = ((Quickflow_heavy_H/Quickflow_light_H)/standard-1)*1000

#Hydrologic results output table
HydOutput <- data.frame(matrix(nrow=length(P)))
HydOutput$Precip <- P
HydOutput$PET <- PET
HydOutput$SoilWater <- SoilWater
HydOutput$GroundWater <- TM_S
HydOutput$Baseflow <- baseflow
HydOutput$OverlandFlow <- OverlandFlow
HydOutput$ShallowInterflow <- ShallowInterflow
HydOutput$Streamflow <- Streamflow
HydOutput$PET_Evap <- PET_Evap
HydOutput$PET_Plant <- PET_Plant
HydOutput$matrix.nrow...length.P.. <- NULL
HydOutput[is.na(HydOutput)] <- -9999

filename = paste(outpath,'/HydModel_Output_Mix_',mixing_model,"_",mew,'_ET_',E_ET_ratio,'.csv',sep="")
#filename = paste(outpath,'/HydModel_Output_EcoHydro.csv',sep="")
write.table(HydOutput, file=filename, sep=",", col.names=F, row.names = F)

#Oxygen Isotopic results output table
IsoOutput <- data.frame(matrix(nrow=length(P)))
IsoOutput$del_p_O <- del_p_O
IsoOutput$del_a_O <- del_a_O
IsoOutput$del_e_O <- del_e_O
IsoOutput$del_Transpiration_O <- del_Transpiration_O
IsoOutput$del_Unsat_O <- del_SoilWater_O
IsoOutput$del_TM_S_O <- del_TM_S_O
IsoOutput$del_baseflow_O <- del_baseflow_O
IsoOutput$del_Streamflow_O <- del_Streamflow_O
IsoOutput$matrix.nrow...length.P.. <- NULL
IsoOutput[is.na(IsoOutput)] <- -9999

filename = paste(outpath,'/IsotopeModel_Oxygen_Output_Mix_',mixing_model,"_",mew,'_ET_',E_ET_ratio,'_RipMix_',SW_Volume,'.csv',sep="")
#filename = paste(outpath,'/IsotopeModel_Oxygen_Output.csv',sep="")
write.table(IsoOutput, file=filename, sep=",", col.names=F, row.names = F)

#Hydrogen Isotopic results output table
IsoOutput <- data.frame(matrix(nrow=length(P)))
IsoOutput$del_p_H <- del_p_H
IsoOutput$del_a_H <- del_a_H
IsoOutput$del_e_H <- del_e_H
IsoOutput$del_Transpiration_H <- del_Transpiration_H
IsoOutput$del_Unsat_H <- del_SoilWater_H
IsoOutput$del_TM_S_H <- del_TM_S_H
IsoOutput$del_baseflow_H <- del_baseflow_H
IsoOutput$del_Streamflow_H <- del_Streamflow_H
IsoOutput$matrix.nrow...length.P.. <- NULL
IsoOutput[is.na(IsoOutput)] <- -9999

filename = paste(outpath,'/IsotopeModel_Hydrogen_Output_Mix_',mixing_model,"_",mew,'_ET_',E_ET_ratio,'_RipMix_',SW_Volume,'.csv',sep="")
#filename = paste(outpath,'/IsotopeModel_Hydrogen_Output.csv',sep="")
write.table(IsoOutput, file=filename, sep=",", col.names=F, row.names = F)

#Write Snowmelt Model Output
filename = paste(outpath,'/SnowModel.csv',sep="")
write.table(snowmelt, file=filename, sep=",", col.names=T, row.names = F)

#Write Plant Model Output
filename = paste(outpath,'/PlantModel.csv',sep="")
write.table(PlantModel, file=filename, sep=",", col.names=T, row.names = F)