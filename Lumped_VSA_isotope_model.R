#Jo-Flow model modified to include
#  - stable water isotopes (Gibson, 2002; Smith et al., 2016)
#  - plant growth model (citation)
#  - dual crop coefficient partitioning of ET (Allen et al., 2006)
#
#
# Last Modified: JK 05/05/2016

library(EcoHydRology)

#1. Set Model Inputs
#Simulation Settings
outpath = 'C:/Isotope/CaseStudy_OH/output/' #James
outpath= '~/Dropbox/IsotopeWatershedModel/Scratch'#Chelsea
record = 365*5 #Simulation length (days)

#Watershed model inputs
Depth = 1000 # Average soil depth in watershed (mm)
SATper = 0.2 # Porosity (fraction)
AWCper = 0.13 # Available Water Capacity (AWC) (fraction)
no_wet_class = 10 # The number of wetness classes for saturated area designation (not currently used)
percentImpervious = 1 # Percent of the watershed that is impervious
Tp = 4 # Time to peak (hours)
latitudeDegrees = 38.64 #decimal degrees
albedo = 0.23 # Average albedo
StartCond = "avg" # Watershed conditions before first day of run ("wet", "dry", "avg")
SW1 = NULL # Soil water on the first day (depth, mm)
BF1 = 10 #  mm/day can use nearby watershed baseflow, ConvertFlowUnits(cfs,WA=W_Area)
PETcap = 5 # Does not let PET get larger than this cut off (mm)
A1 = 0.002 #Baseflow recession coefficient (Rossman, 2009)
B = 2 #Baseflow recession exponent  (Rossman, 2009)
Se_min = 78 # Minimum effective watershed storage a la Curve Number method (mm)
C1 = 3.1 # Coefficient relating soil water to Curve Number S
Ia_coef = 0.08 # range ~ 0.05 - 0.2
TEW = 20 #total evaporable water, based on Allen et al (2006) Table 19
REW = 10 #readily evaporable water, basde on Allen et al (2006) Table 19

#Isotope model inputs
standard = 2005/1000000 #Vienna standard
h = 0.7 #relative humidity, assumed constant for simplicity
R_Soil_Init = 0.00199 #Initial isotope ratio for unsaturated zone storage
E_ET_ratio = 2 #Ratio of evaporation to total ET, 
#note: use E_ET_ratio = 2 to use dual crop coefficient model
mixing_model = 3 #unsaturated zone mixing model: 
#1 - no mixing (preferrential flow via macro pores)
#2 - complete mixing (uniform infiltration)
#3 - partial mixing (a weighting of the above models)
mew = 0.22 #weighting factor for partial mixing model (0-1)

#Plant growth model inputs - (citation)
Tb = 1 #base temperature for computing growing degree days (degrees C)
GDDmax = 2500 #maximum growing degree days in annual plant cycle
GDD_1 = 0*GDDmax #stage 1 growth
GDD_2 = 0.1*GDDmax #stage 2 growth
GDD_3 = 0.225*GDDmax #stage 3 growth
GDD_4 = 0.9*GDDmax #stage 4 growth
GDD_5 = 1*GDDmax #stage 5 growth
Kcmin = 0.25 #minimum crop coefficient (Kc) for deciduous forest
Kcmax = 1 #max crop coefficient (Kc)

#Constants to calculate
AWC = Depth*AWCper # AWC as a depth (mm)
SAT = Depth*SATper # porosity as a depth (mm)
runoff_breakdown = RunoffBreakdown(Tp, HrPrcDelay = (Tp/2-4)) 
latitude<-latitudeDegrees*pi/180 ## latitude in radians

#2.a Read in Meteorological data
filepath="~/Dropbox/IsotopeWatershedModel/ModelInputData/" #Chelsea
#Reads in distribution of del18O obesrvations for precipitation
del_precip_obs = read.csv(paste(filepath,'Precip_Isotopes.txt', sep=""))
colnames(del_precip_obs)="del_p"
R_precip_obs = (del_precip_obs$del_p/1000 + 1)*standard
#Reads in NCDC data with PRCP, TMIN, TMAX
precip= read.csv(paste(filepath,'NCDC_MetData_Station1.csv', sep="")) #precip input is in inches
precip$PRCP[precip$PRCP == -9999] = 0

#Make sure Tmax > Tmin
#Note: Some NCDC datasets have errors where TMIN > TMAX
for (j in 1:length(precip$TMAX))
{
  if (precip$TMIN[j] == -9999)
  {precip$TMIN[j] = precip$TMIN[j-1]}
  
  if (precip$TMAX[j] == -9999)
  {precip$TMAX[j] = precip$TMAX[j-1]}
  
  if (precip$TMIN[j]>precip$TMAX[j])
  {precip$TMAX[j] = precip$TMAX[j-1]
   precip$TMIN[j] = precip$TMIN[j-1]}
}

#Met Data conversion to SI
precip$PRCP=precip$PRCP*25.4 #precipitation depth is now in mm
precip$TMAX=(precip$TMAX - 32)*5/9 #temperature in degrees C
precip$TMIN=(precip$TMIN - 32)*5/9 #temperature in degrees C

#Create date column for snowmelt model
precip$DATE=as.Date.numeric(precip$Date, origin="2000-01-01")

#3. Snowmelt Model (based on Walter et al (2005))
snowmelt=SnowMelt(Date=precip$DATE, precip_mm=precip$PRCP, Tmax_C=precip$TMAX, Tmin_C=precip$TMIN, lat_deg=latitudeDegrees)

#Divide Precipitation plus smowmelt into two isotopic inputs
P <- (snowmelt$SnowMelt_mm[1:record]+snowmelt$Rain_mm[1:record])
Tmin <- precip$TMIN[1:record]
Tmax <- precip$TMAX[1:record]
Tav <- 0.5*(Tmin + Tmax)
#Since we do not know the exact isotopic signature on each day, we work with probability distributions
#Snowmelt is assigned the same distribution as precip. This is OK since snowmelt is minimal
Ri = sample(R_precip_obs, size=record, replace=TRUE)
Ri = Ri/(1 + Ri)
P_heavy <- snowmelt$SnowMelt_mm[1:record]*Ri + snowmelt$Rain_mm[1:record]*Ri
P_light <- snowmelt$SnowMelt_mm[1:record]*(1-Ri) + snowmelt$Rain_mm[1:record]*(1 - Ri)
del_p = 1000*((P_heavy/P_light)/standard - 1)

#4. PET Model
# Create Julian Date format for use in PET model
x=as.POSIXlt(precip$DATE[1:record], format="Y%b%d")
day=x$yday+1

#Run PET from daily temp model (based on Archibald et al (year)) 
PET<-PET_fromTemp(Jday=day,Tmax_C=Tmax,Tmin_C=Tmin,AvgT=Tav,albedo=albedo,lat_radians=latitude)*1000## mm (Priestley-Taylor)
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
  #Calculate GDD
  PlantModel$Accumulate_GDD[j] = PlantModel$Accumulate_GDD[j-1] + PlantModel$DegreeDays[j]
  if (PlantModel$Year[j] > PlantModel$Year[j-1])
  {PlantModel$Accumulate_GDD[j] = PlantModel$DegreeDays[j]}
  #Calculate GDD as percent of GDDmax
  PlantModel$DegreeDays_percent[j] = 100*min(PlantModel$Accumulate_GDD[j] / GDDmax,1)
}

#Calculate alpha
for (j in 1:nrow(PlantModel))
{
  if (PlantModel$Accumulate_GDD[j] <= GDD_2)
  {PlantModel$alphaK[j] = PlantModel$Accumulate_GDD[j] / GDD_4}
  if (PlantModel$Accumulate_GDD[j] > GDD_2 & PlantModel$Accumulate_GDD[j] <= GDD_3)
  {PlantModel$alphaK[j] = (GDD_2/GDD_4) + ((GDD_4 - GDD_2)/(GDD_3-GDD_2))*(PlantModel$Accumulate_GDD[j] - GDD_2)/GDD_4}
  if (PlantModel$Accumulate_GDD[j] > GDD_3 & PlantModel$Accumulate_GDD[j] <= GDD_4)
  {PlantModel$alphaK[j] = 1}
  if (PlantModel$Accumulate_GDD[j] > GDD_4 & PlantModel$Accumulate_GDD[j] < GDD_5)
  {PlantModel$alphaK[j] = 1 - 0.6*((PlantModel$Accumulate_GDD[j] - GDD_4)/(GDDmax - GDD_4))}
  if (PlantModel$Accumulate_GDD[j] > GDD_5)
  {PlantModel$alphaK[j] = 0}
  
}

#Calculate K and Z
PlantModel$Kc = Kcmin + PlantModel$alphaK*(Kcmax - Kcmin)

#6. Adjust PET initially based on Kc coefficient
ETo = PET*PlantModel$Kc
ET <- ETo    #  Initializing modeled actual ET, equal to ETo when deltaP > 0, but less than ETo on drier days

#7.Watershed Isotope Model
#Effective soil water storage coefficients, eswsc = Se*(2..... see Schneiderman et al 2007)
eswsc <- vector(mode="numeric", length=no_wet_class)
eqn16 <- function(x){(2*(sqrt(1-(1/no_wet_class)*(x-1))-sqrt(1-(1/no_wet_class)*x))/(1/no_wet_class))-1}
x <- seq(1,no_wet_class, by=1)
eswsc <- eqn16(x)

# Initializing vectors for Water Budget Loop, and Runoff estimation
Pe<-vector(length=length(P))# Effective Precipitation (non-impervious areas)

Ratio <- vector(length=length(P)) #mixing model ratio
del_a <- vector(length=length(P)) #del atmospheric
alpha_star <- vector(length=length(P))
epsilon <- vector(length=length(P))
del_e <- vector(length=length(P)) #del evaporation
R_evap<- vector(length=length(P)) #R of evaporation

PET_Evap <- vector(length=length(P))  
PET_Plant <- vector(length=length(P))
PET_Evap_heavy <- vector(length=length(P))
PET_Evap_light <- vector(length=length(P))  
PET_Plant_heavy <- vector(length=length(P))
PET_Plant_light <- vector(length=length(P))  

Evap_light <- vector(length=length(P))
Evap_heavy <- vector(length=length(P))

SoilWater <- vector(length=length(P))##(mm)
SoilWater_heavy<-vector(length=length(P))##(mm)
SoilWater_light <-vector(length=length(P))##(mm)

del_SoilWater <- vector(length=length(P))
del_Transpiration <- vector(length=length(P))
del_TM_S<-vector(length=length(P))
del_baseflow<-vector(length=length(P))
del_totQ<-vector(length=length(P))

baseflow<-vector(length=length(P))
baseflow_heavy <- vector(length=length(P))
baseflow_light <- vector(length=length(P))

excess<-vector(length=length(P))
excess_heavy<-vector(length=length(P))
excess_light<-vector(length=length(P))

TM_S<-vector(length=length(P))##This is the daily time-step T-M storage, used to calc baseflow
TM_S_heavy<-vector(length=length(P))
TM_S_light<-vector(length=length(P))
GW_Ratio<-vector(length=length(P))

totQ<-vector(length=length(P))
totQ_heavy<-vector(length=length(P))
totQ_light<-vector(length=length(P))

Se<-vector(length=length(P))
Q2<-vector(length=length(P))

MaxWetClass<-vector(length=length(P))
sigmaS<-matrix(nrow=length(P), ncol=no_wet_class)  # Storage in each wetness class
runoff_by_class<-matrix(nrow=length(P), ncol=no_wet_class)

impervRunoff<-vector(length=length(P))
OverlandFlow<-vector(length=length(totQ))
OverlandFlow_heavy<-vector(length=length(totQ))
OverlandFlow_light<-vector(length=length(totQ))

ShallowInterflow<-vector(length=length(totQ))
ShallowInterflow_heavy<-vector(length=length(totQ))
ShallowInterflow_light<-vector(length=length(totQ))

#Calculate net precipitation
deltaP <- P - ETo## (mm) neg values indicate net evap
impervIa<-Ia_coef * 5  ##Se = 5 mm for impervious areas (CN=98)
impervPe<-deltaP - impervIa##  Effective precipitation on impervious areas
impervPe[which(impervPe < 0)]<-0

#Model initial conditions
#these values generally do not affect the results, particularly if a spin-up year is used
del_p[1] = -10  #just setting this to something, doesn't affect results

TM_S[1]<- 50
TM_S_heavy[1] = R_Soil_Init*TM_S[1]
TM_S_light[1] = (1-R_Soil_Init)*TM_S[1]

# Initialize soil wetness based on initial conditions ("wet", "dry", "avg") and initial isotopic ratio
SoilWater_heavy[1] <- ifelse(StartCond == "wet", AWC, ifelse(StartCond == "dry", R_Soil_Init*0.2*AWC, R_Soil_Init*0.75*AWC))
SoilWater_light[1] <- ifelse(StartCond == "wet", AWC, ifelse(StartCond == "dry", (1-R_Soil_Init)*0.2*AWC, (1-R_Soil_Init)*0.75*AWC))
SoilWater[1] = SoilWater_heavy[1] + SoilWater_light[1]

excess_heavy[1] <- 0# Assume no runoff from the day before
excess_light[1] <- 0# Assume no runoff from the day before

Se[1]<- Se_min+C1*(AWC-SoilWater[1])
sigmaS[1,] <- eswsc*Se[1]

totQ[1] <- Pe[1]*Pe[1]/(Pe[1]+Se[1])
totQ_heavy[1] = totQ[1]*Ri[1]
totQ_light[1] = totQ[1]*(1-Ri[1])

Ia<-Ia_coef*Se  # Initial abstraction = depth of precip before runoff starts

PlantModel$De[1] = 0
PlantModel$Kr[1] = 0
PlantModel$Ke[1] = 0

#Thornthwaite-Mather Function and Runoff Generation 
for(i in (2):length(deltaP))
{
  
  #Calculate time step runoff with CN method
  if (deltaP[i]-Ia[i-1]>0)
  {Pe[i] <- deltaP[i]-Ia[i-1]} 
  else {Pe[i]<-0}
  
  totQ[i] <- Pe[i]*Pe[i]/(Pe[i]+Se[i-1]) ## Effective storage is from previous day
  totQ_heavy[i] <- totQ[i]*Ri[i]
  totQ_light[i] <- totQ[i]*(1-Ri[i])
  
  #Partition PET into transpiration and evaporation based on Ke and Kc
  PET_Evap[i] = ET[i]*(PlantModel$Ke[i-1] / (PlantModel$Ke[i-1] + PlantModel$Kc[i-1]))
  PET_Plant[i] = ET[i]*(PlantModel$Kc[i-1] / (PlantModel$Ke[i-1] + PlantModel$Kc[i-1]))

  #Determine how much transpiration pulls from each isotope pool
  #Assumes uniform distribution of 18O stored in the soil matrix
  PET_Plant_heavy[i] = max(0,PET_Plant[i]*SoilWater_heavy[i-1]/SoilWater[i-1])
  PET_Plant_light[i] = max(0,PET_Plant[i]*SoilWater_light[i-1]/SoilWater[i-1])
  if (is.nan(PET_Plant_light[i]))  {PET_Plant_light[i] = 0}
  if (is.nan(PET_Plant_heavy[i])) {PET_Plant_heavy[i] = 0}
  
  #Determine evaporation fractionation (from Gibson, 2002)
  Ta = 0.5*(snowmelt$MaxT_C[i] + snowmelt$MinT_C[i])
  del_SoilWater[i-1] = ((SoilWater_heavy[i-1]/SoilWater_light[i-1])/standard - 1)*1000
  #If no rain today, assume precipitation del 18O same as day before
  if (is.nan(del_p[i]))  {del_p[i] = del_p[i-1]}
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
  del_a[i] = (del_p[i] - e_star)/alpha_star[i]
  #Estimate evaporation del 18O
  del_e[i] = (del_SoilWater[i-1] - h*del_a[i] - epsilon[i])/(1-h+e_k/1000)
  R_evap[i] = (del_e[i]/1000 + 1)*standard
  R_evap_total = R_evap[i]/(1 + R_evap[i])
  
  #PET partitioning based on previous time step, dual crop coefficient model (Allen et al., 1998)
  PET_Evap_heavy[i] = PET_Evap[i]*(R_evap_total)
  PET_Evap_light[i] = PET_Evap[i]*(1 - R_evap_total)
  
  # Water Balance
  #a. DRYING CONDITION
  if(deltaP[i]<=0)
  {
      #Determine total water loss
      SoilWater[i] <- SoilWater[i-1]*exp(deltaP[i]/AWC)
      
      #Determine fractionation evaporation
      if (deltaP[i] != 0)
      {
        Losses = SoilWater[i-1] - SoilWater[i]
        Heavy_losses = Losses*(PET_Evap_heavy[i] + PET_Plant_heavy[i])/(PET_Evap_heavy[i] + PET_Plant_heavy[i] + PET_Evap_light[i] + PET_Plant_light[i])
        Light_losses = Losses - Heavy_losses
        SoilWater_light[i] = SoilWater_light[i-1] - Light_losses
        SoilWater_heavy[i] = SoilWater_heavy[i-1] - Heavy_losses
      }
      else
      {
        SoilWater_light[i] = SoilWater_light[i-1]
        SoilWater_heavy[i] = SoilWater_heavy[i-1]
      }

      ET[i] <- SoilWater[i-1]*(1-exp((P[i] - ET[i])/AWC))  # total amount that gets evaporated (mm)
      excess[i] <- 0
      excess_heavy[i] <- 0
      excess_light[i] <- 0
      
      #Cummulative evaporation losses since last precip / snowmelt event (All et al., 2006)
      PlantModel$De[i] = PlantModel$De[i-1] +  ET[i]
  }
  
  #b. WETTING BELOW FIELD CAPACITY
  else if (deltaP[i]+SoilWater[i-1] <= AWC)
  {  
      SoilWater_light[i] <- SoilWater_light[i-1] + P_light[i] - PET_Evap_light[i] - PET_Plant_light[i]
      SoilWater_heavy[i] <- SoilWater_heavy[i-1] + P_heavy[i] - PET_Evap_heavy[i] - PET_Plant_heavy[i]

      SoilWater[i]<- SoilWater_heavy[i] + SoilWater_light[i]
      excess[i] <- 0
      excess_heavy[i] <- 0
      excess_light[i] <- 0
  
      PlantModel$De[i] = 0
  }
  
  #c. WETTING ABOVE FIELD CAPACITY
  else { 
    
    SoilWater[i] <- AWC  ## So overall SW cannot exceed field capacity (AWC)

    #Determine total volume of percolated water
    excess[i] <- deltaP[i]+SoilWater[i-1]-AWC
    
    #Model 1: No Mixing
    if (mixing_model == 1)
    {Ratio[i] = Ri[i]}
    
    #Model 2: Complete Mixing
    if (mixing_model == 2)
    {Ratio[i] = (SoilWater_heavy[i-1] +  P_heavy[i])/(SoilWater[i-1] + P_light[i] + P_heavy[i])}  
    
    #Model 3: Partial Mixing       
    if (mixing_model == 3)
    {
      #Ratio of heavy isotope for no mixing
      Ratio_nomix = Ri[i]
      #Ratio of heavy isotope for complete mixing
      Ratio_complete = (SoilWater_heavy[i-1] +  P_heavy[i])/(SoilWater[i-1] + P_light[i] + P_heavy[i])                
      #Ratio of heavy isotope for partial mixing    
      Ratio[i] = mew*Ratio_complete + (1-mew)*Ratio_nomix
    }
    
    #Compute new unsaturated zone storage and excess from excess loss ratio
    SoilWater_heavy[i] <- SoilWater_heavy[i-1] + (AWC - SoilWater[i-1])*Ratio[i]
    SoilWater_light[i] <- SoilWater_light[i-1] + (AWC - SoilWater[i-1])*(1-Ratio[i])
    excess_heavy[i] <- Ratio[i]*(deltaP[i]+SoilWater[i-1]-AWC)
    excess_light[i] <- (1-Ratio[i])*(deltaP[i]+SoilWater[i-1]-AWC)
    PlantModel$De[i] = 0
  }
  
  #Solve evaporation loss coefficients (Allen et al 1998)
  PlantModel$Kr[i] = 0.5*max(0,(TEW - PlantModel$De[i]))/(TEW - REW)
  PlantModel$Ke[i] = min(PlantModel$Kr[i]*(max(PlantModel$Kc)+0.1 - PlantModel$Kc[i]), max(PlantModel$Kc)+0.1)
  
  #Compute baseflow
  #Baseflow modified to use an exponential, slightly better fit for summer months (JK 05/04/2016)
  baseflow[i] = min(A1*(TM_S[i-1]^B),TM_S[i-1]-1)
  GW_Ratio[i] = TM_S_heavy[i-1] / TM_S[i-1]
  baseflow_heavy[i] <- GW_Ratio[i]*baseflow[i]
  baseflow_light[i] <- (1-GW_Ratio[i])* baseflow[i]
    
  # Ensure mass-balance, since curve number is empirical
  if ((excess[i])>=totQ[i])
  {
    TM_S[i]<-TM_S[i-1] - baseflow[i] + excess[i] - totQ[i]
    TM_S_heavy[i] <- TM_S_heavy[i-1] - baseflow_heavy[i] + excess_heavy[i] - totQ_heavy[i]
    TM_S_light[i] <- TM_S_light[i-1] - baseflow_light[i] + excess_light[i] - totQ_light[i]
  } 
  else 
  {
    TM_S[i]<-TM_S[i-1] - baseflow[i]
    TM_S_heavy[i]<-TM_S_heavy[i-1] - baseflow_heavy[i]
    TM_S_light[i]<-TM_S_light[i-1] - baseflow_light[i]
    SoilWater[i]<-SoilWater[i]-(totQ[i]-excess[i])
    SoilWater_heavy[i]<-SoilWater_heavy[i]-(totQ_heavy[i]-excess_heavy[i])
    SoilWater_light[i]<-SoilWater_light[i]-(totQ_light[i]-excess_light[i])
  }
  
  #Update initial abstraction for next time step
  Se[i] <- Se_min+C1*(AWC-SoilWater[i])
  sigmaS[i,]<- eswsc*Se[i]  ## Amount of storage available in each wetness class [mm]
  Ia[i]<-Ia_coef*Se[i]
  
  #Impervious Runoff
  impervRunoff[i]<-impervPe[i] 

  #8. Post-Processing
  #Convert R to delta notation
  del_Transpiration[i] = ((PET_Plant_heavy[i]/PET_Plant_light[i])/standard-1)*1000
  del_TM_S[i] = ((TM_S_heavy[i]/TM_S_light[i])/standard - 1)*1000
  del_baseflow[i] = ((baseflow_heavy[i]/baseflow_light[i])/standard - 1)*1000
  del_totQ[i] = ((totQ_heavy[i]/totQ_light[i])/standard - 1)*1000
  
} 

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

Streamflow <- (OverlandFlow+ShallowInterflow)*(1-0.01*percentImpervious)+impervRunoff*(0.01*percentImpervious)+baseflow#totQ+baseflow
Streamflow_heavy = Ri*OverlandFlow + (SoilWater_heavy[i]/SoilWater[i])*ShallowInterflow + baseflow_heavy
Streamflow_light = (1-Ri)*OverlandFlow + (1 - SoilWater_heavy[i]/SoilWater[i])*ShallowInterflow + baseflow_light
del_Streamflow = ((Streamflow_heavy/Streamflow_light)/standard-1)*1000

Quickflow <- (OverlandFlow+ShallowInterflow)*(1-0.01*percentImpervious)+impervRunoff*(0.01*percentImpervious)
Quickflow_heavy <- Ri*(OverlandFlow + ShallowInterflow)
Quickflow_light <- (1-Ri)*(OverlandFlow + ShallowInterflow)
del_Quickflow = ((Quickflow_heavy/Quickflow_light)/standard-1)*1000

#CDF Figure
tiff(filename = paste(outpath,'dels',mixing_model,'_EcoHydro_',mew,'.tiff',sep=""), width = 800, height = 600, units = "px", pointsize = 12)
par(mar = c(5, 4, 4, 4) + 0.4) 
plot(ecdf(del_p),xlab="",main="",lwd=3,xlim=c(-30,20),ylab="",col="black",cex=1.0, cex.lab=2.0, cex.axis=2.0)
par(new=TRUE)
plot(ecdf(del_e),xlab="",main="",lwd=3,xlim=c(-30,20),ylab="",col="blue",cex=1.0, cex.lab=2.0, cex.axis=2.0)
par(new=TRUE)
plot(ecdf(del_SoilWater),main="",lwd=3,xlim=c(-30,20),xlab="",ylab="",col="red",cex=1.0, cex.lab=2.0, cex.axis=2.0)
par(new=TRUE)
plot(ecdf(del_TM_S),main="",lwd=3,xlim=c(-30,20),ylab="",xlab="",col="orange",cex=1.0, cex.lab=2.0, cex.axis=2.0)
par(new=TRUE)
plot(ecdf(del_Streamflow),main="",lwd=3,xlim=c(-30,20),ylab="CDF",xlab=expression(paste(delta, "18 O")),col="brown",cex=1.0, cex.lab=2.0, cex.axis=2.0)

legend(-30, 0.95, legend=c("del P", "del E","del UnSat", "del Groundwater", "del_Streamflow"),
       col=c("black","blue","red","orange","brown"), lwd=5, cex=1.5)
dev.off()

#Time Series Figure
tiff(filename = paste(outpath,'dels_TS_EcoHydro_.tiff'), width = 1600, height = 800, units = "px", pointsize = 12)
par(mar = c(5, 4, 4, 4) + 0.4) 
plot(del_p,type="l",ylim=c(-20,0),xlab="",ylab="",col="gray",lwd=2,cex=1.0, cex.lab=2.0, cex.axis=2.0)
par(new=TRUE)
plot(del_SoilWater,type="l",ylim=c(-20,0),xlab="",ylab="",col="red",lwd=2,cex=1.0, cex.lab=2.0, cex.axis=2.0)
par(new=TRUE)
plot(del_TM_S,type="l",ylim=c(-20,0),xlab="",ylab="",col="green",lwd=4,cex=1.0, cex.lab=2.0, cex.axis=2.0)
par(new=TRUE)
plot(del_Streamflow,type="l",ylim=c(-20,0),xlab="Elapsed Time (days)",lwd=2,ylab=expression(paste(delta, "18 O")),col="blue",cex=1.0, cex.lab=2.0, cex.axis=2.0)
legend(0, -15, legend=c("del Precip","del UnSat", "del Groundwater","del Streamflow"),
       col=c("gray","red", "green","blue"), lty=1:2, cex=1.5)
dev.off()

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

filename = paste(outpath,'/HydModel_Output_',mixing_model,'_EcoHydro_',mew,'.csv',sep="")
write.table(HydOutput, file=filename, sep=",", col.names=F, row.names = F)

#Isotopic results output table
IsoOutput <- data.frame(matrix(nrow=length(P)))
IsoOutput$del_p <- del_p
IsoOutput$del_a <- del_a
IsoOutput$del_e <- del_e
IsoOutput$del_Transpiration <- del_Transpiration
IsoOutput$del_Unsat <- del_SoilWater
IsoOutput$del_TM_S <- del_TM_S
IsoOutput$del_baseflow <- del_baseflow
IsoOutput$del_Streamflow <- del_Streamflow
IsoOutput$matrix.nrow...length.P.. <- NULL
IsoOutput[is.na(IsoOutput)] <- -9999

filename = paste(outpath,'/IsotopeModel_StochasticOutput_',mixing_model,'_EcoHydro_',mew,'.csv',sep="")
write.table(IsoOutput, file=filename, sep=",", col.names=F, row.names = F)
