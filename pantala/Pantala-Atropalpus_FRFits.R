#Load required packages
require(emdbook)
require(bbmle)
require(R2WinBUGS)
require(plotrix)
require(plyr)
require(Hmisc)
require(broom)

frdat <- read.csv("TDFR_Pantala-Atro_data.csv")
source("ModelingFunctions.R")


frdat$TempK <- frdat$Temp + 273.15
sv.log <- list(b=log(1.105563e-01), h0 = log(1.99775e-14), Eh = log(0.65))

frfit_pan <-
  mle2(
    Killed ~ dbinom(prob = (
      rogers.sentis.risk.log(N0, b, T0=log(283.15), Tl=log(315.35), h0, M=401.6, bh=0.75, Eh, k=8.62e-5, TempK, t=1)
    ),
    size = N0),
    start = sv.log,
    trace = TRUE,
    data = frdat
  )

exp(coef(frfit_pan))

summary(frfit_pan)

parms_pan <- tidy(frfit_pan)
write.csv(parms_pan,"Parms_Pan_FR.csv",row.names=FALSE)

############# Functional Response Plots

predparms <- read.csv( "Parms_Pan_FR_Wide.csv" )

N0_vec <- seq( 0, 125, by = 1 )
TempK_fac <- c( 293.15, 297.15, 301.15, 305.15, 309.15 )
ddepfr_out.dat <- mat.or.vec( nrow( predparms ) * length( TempK_fac ) * length ( N0_vec ), 4)
colnames( ddepfr_out.dat ) <- c( "Pred", "Temp", "N0", "Killed" )
row <- 0

for(h in 1:nrow( predparms ) ){
  for ( i in 1:length( N0_vec ) ){
    for (j in 1:length( TempK_fac ) ){
      row <- row + 1
      ddepfr_out.dat[ row, "Pred" ] <- predparms[ h, 1 ]
      ddepfr_out.dat[ row,"Temp" ] <- TempK_fac[ j ] - 273.15
      ddepfr_out.dat[ row, "N0" ] <- N0_vec[ i ]
      ddepfr_out.dat[ row, "Killed" ] = rogers.sentis( N0 = N0_vec[ i ], 
                                                       b = exp( predparms[ h, 2 ] ), 
                                                       T0 = predparms[ h, 6 ], 
                                                       Tl = predparms[ h, 7 ], 
                                                       h0 = exp( predparms[ h, 3 ] ), 
                                                       Eh = exp( predparms[ h, 4 ] ), 
                                                       M = predparms[ h, 5 ], 
                                                       bh = 0.75, 
                                                       k = 8.617e-5, 
                                                       TempK = TempK_fac[ j ], 
                                                       t = 1 )
    }
  }
}

ddepfr_out.dat <- as.data.frame( ddepfr_out.dat )
ddepfr_out.dat$Pred <- as.factor(ddepfr_out.dat$Pred)
ddepfr_out.dat$Temp <- as.factor(ddepfr_out.dat$Temp)
ddepfr_out.dat$N0 <- as.numeric(ddepfr_out.dat$N0)
ddepfr_out.dat$Killed <- as.numeric(ddepfr_out.dat$Killed)


frdat$Temp <- as.factor( frdat$Temp )

ddep.plot <- ggplot( data.frame( N0 = c( 0, 120 ) ), aes( N0 ) ) +
  geom_point( data = frdat, 
              aes( x = N0, y= Killed, colour = Temp ) ) +
  geom_line( data = ddepfr_out.dat, aes( x = N0, y = Killed, colour = Temp ), size = 0.75 ) +
  ylab( "Prey Killed" ) +
  xlab( "Prey Density" ) +
  scale_y_continuous( expand = c( 0, 0 ), limits = c( 0, 65 ) ) +
  scale_x_continuous( expand = c( 0, 0 ), limits = c( 0, 125 ), breaks = c(0, 20, 40, 60, 80, 100, 120) ) +
  scale_color_manual( "Temp (°C)", values = c( "blue", "green3", "yellow2", "orange2", "red2" ) ) +
  theme_classic( base_size = 16, base_line_size = 1 ) + 
  theme( panel.border = element_rect( color = "black", fill = NA, size = 0.75 ),
         aspect.ratio = 1 )

ggsave( "FRFig.png",
        plot = ddep.plot,
        units = c( "mm" ),
        dpi = 300 )


############ Model Runs

PreyDev.df <- read.csv( "Parms_PreyDev.csv" )
PreyMort.df <- read.csv( "Parms_PreyMort.csv" )

# Set prediction interval
temp <- seq( 20, 36 )
ti_vec <- rep(20, length(temp))

# Vector of temperature-dependent transition rates
rvec <- ( ( arr.rate( c_dev = PreyDev.df[ 1, 2 ], b_dev = PreyDev.df[ 2, 2 ], temp ) ) * 5 )

# Vector of temperature-dependent mortality rates
mvec <- logisticfun( temp, r_mort = PreyMort.df[ 3, 2 ], x1 = PreyMort.df[ 1, 2 ], y1 = PreyMort.df[ 2, 2 ] )

num_reps <- 3000
Summary <- mat.or.vec( num_reps * length( temp ), 6 ) # Creates a matrix of r_dev and p values for storing the model output, with five columns
colnames( Summary ) <- c( "rep", "b", "h0", "Eh", "temp", "num_survive" ) # Names the columns in the matrix

#Parameter estimates:
Parms <- read.csv( "Parms_Pan_FR.csv" )
Parms$spp <- as.factor(Parms$spp)
Parms$term <- as.factor(Parms$term)
Parms$std.error <- as.numeric(Parms$std.error)


#Some constants:
k <- 8.617e-5
bh <- 0.75
t <- 1
M <- Parms[ 4, 3 ]
T0 <- Parms[ 5, 3 ]
Tl <- Parms[ 6, 3 ]

row <- 0
u <- 1
seednum <- 1

for ( r in 1:num_reps ) {
  set.seed( 1 + r )
  
  b <- rlnorm( 1, meanlog = Parms$estimate[ which( Parms$term == "b" ) ],
               sdlog = Parms$std.error[ which( Parms$term == "b" ) ] )
  h0 <- rlnorm( 1, meanlog = Parms$estimate[ which( Parms$term == "h0" ) ],
                sdlog = Parms$std.error[ which( Parms$term == "h0" ) ] )
  Eh <- rlnorm( 1, meanlog = Parms$estimate[ which( Parms$term == "Eh" ) ],
                sdlog = Parms$std.error[ which( Parms$term == "Eh" ) ] )
  
  for ( l in 1:length( temp ) ) {
    row <- row + 1
    N_vec <- c( 300, 0, 0, 0, 0) #Vector of prey abundance at each size class at the current time step. This starts with 1000 1st instar larvae.
    Summary[ row, "b" ] <- b
    Summary[ row, "Eh" ] <- Eh
    Summary[ row, "h0" ] <- h0
    Summary[ row, "rep" ] = r #rep number
    Summary[ row, "temp" ] = temp[ l ] #Temperature for this loop
    Summary[ row, "num_survive" ] = Stage_pred_devel(
      timesteps = 30,
      N_vec,
      mval = mvec[ l ],
      rval = rvec[ l ],
      bval = b,
      T0val = T0,
      Tlval = Tl,
      h0val = h0,
      bhval = bh,
      Ehval = Eh,
      Mval = M,
      kval = k,
      tval = t,
      tempval = temp[ l ]
    )
  }
}


Pan_Summary.df <- as.data.frame( Summary )

write.csv( Pan_Summary.df, "Output_PanSpp_PredTI.csv", row.names = FALSE)



######## Plots

## Some data wrangling first

pandat <- read.csv( "Output_PanSpp.csv" )
pandat_preyTI <- read.csv( "Output_PanSpp_PreyTI.csv" )
pandat_predTI <- read.csv( "Output_PanSpp_PredTI.csv" )

pandat$Tdep <- c( "Both" )
pandat$Pred <- c( "Pantala" )

pandat_preyTI$Tdep <- c( "Predator Only" )
pandat_preyTI$Pred <- c( "Pantala" )

pandat_predTI$Tdep <- c( "Prey Only" )
pandat_predTI$Pred <- c( "Pantala" )


modelout <- rbind( pandat, pandat_preyTI )
modelout <- rbind( modelout, pandat_predTI )
modelout$Tdep <- factor( modelout$Tdep, levels = c( "Prey Only", "Both", "Predator Only" ) )
modelout <- subset(modelout, Tdep == "Both")


modelfig <- ggplot( modelout, aes( x = temp, y = num_survive, colour = Tdep ) ) +
  geom_smooth( aes( colour = Tdep ), size = 1  ) +
  scale_x_continuous( breaks = c( 20, 24, 28, 32, 36 ) ) +
  scale_y_continuous( expand = c( 0, 0 ) ) +
  coord_cartesian( ylim = c(25, 110)) +
  scale_color_manual( name = "Temp Dependence", 
                      values = c( "green3", "green3", "red2" ), 
                      labels = c( "Prey", "Both", "Predator" ),
                      guide = "legend" ) +
  theme_classic( base_size = 24, base_line_size = 1 ) + 
  theme( legend.position = c( 0.25, 0.86 ),
         legend.text = element_text( size = 14 ),
         legend.title = element_text( size = 16 ),
         legend.key = element_blank(),
         legend.background = element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.line = element_line( colour = "black" ), 
         panel.spacing = unit( .05, "lines" ), 
         panel.border = element_rect( color = "black", fill = NA, size = 1 ), 
         strip.background = element_rect( color = "black", size = 1 ) ) +
  xlab( "Temperature (°C)" ) +
  ylab( "Surviving Prey" ) +
  theme( aspect.ratio = 1 ) +
  theme( plot.margin = unit( c( 5.5, 7.5, 5.5, 5.5 ), unit = "mm" ) )


ggsave( "modelfig2.png", plot = modelfig, units = c( "mm" ), dpi = 300 )
