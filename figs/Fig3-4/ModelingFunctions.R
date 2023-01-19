## Functional Response Functions

rogers.sentis <- function( N0, b, T0, Tl, h0, M, bh, Eh, k, TempK, t ){ # Type II FR, fitted with depletion, with a and h modeled as temp-dependent as described in Sentis et al. 2012
  h = h0 * ( M ^ bh ) * exp( ( Eh ) / ( k * TempK ) )
  a = b *( TempK - T0 ) * ( ( Tl - TempK )^( 1 / 2 ) )
  N0 - lambertW( a * h * N0 * exp( -a * ( t - h * N0 ) ) ) / ( a * h )
}

rogers.sentis.risk <- function( N0, b, T0, Tl, h0, M, bh, Eh, k, TempK, t ){ # The same function, expressed as a risk function
  h = h0 * ( M ^ bh ) * exp( ( Eh ) / ( k * TempK ) )
  a = b * ( TempK - T0 )*( ( Tl - TempK ) ^ ( 1 / 2 ) )
  ( N0 - lambertW( a * h * N0 * exp( -a * ( t - h * N0 ) ) ) / ( a * h ) ) / N0
}

rogers.sentis.risk.log <- function( N0, b, T0, Tl, h0, M, bh, Eh, k, TempK, t ){ # Log-transformed, to help with fitting the FR
  h = exp( h0 ) * ( M ^ bh ) * exp( exp ( Eh ) / ( k * TempK ) )
  a = exp( b )*( TempK - exp( T0 ) )*( ( exp( Tl ) - TempK ) ^ ( 1 / 2 ) )
  ( N0 - lambertW( a * h * N0 * exp( -a * ( t - h * N0 ) ) ) / ( a * h ) ) / N0
}

sentis.attack <- function( b, TempK, T0, Tl ){ # Broken down further, this is the temp-dependence of attack rate described by Sentis et al. 2012
  b * ( TempK - T0 ) * ( ( Tl - TempK ) ^ ( 1 / 2 ) )
}

sentis.handle <- function( h0, M, bh, Eh, k, TempK ){ # And the temp-dependence of handling time
  h0 * (M ^ bh) * exp( Eh / ( k * TempK ) )
}


## Development Rate Functions

arr.rate <- function( c_dev, b_dev, temp ) { # The Arrhenius function, which describes the temp-dependence of prey development
  ( c_dev * exp( b_dev / temp ) )
}

arr.time <- function( c_dev, b_dev, temp ) { # The Arrhenius function, which describes the temp-dependence of prey development
  1 / ( c_dev * exp( b_dev / temp ) )
}

logisticfun <- function( temp, r_mort, x1, y1 ) { # A logistic function used to model prey mortality as a function of temperature for Ae. atropalpus
    1 / ( 1 + exp( -r_mort * (temp - x1 ) ) ) + y1
} 



## Modeling Functions

Stage_pred_devel <-
  function( timesteps,# Number of timesteps for the simulation
           N_vec,# Vector of individuals in each size/stage class at the current time
           mval,# Temperature-dependent prey mortality rate
           rval, # Temperature-dependent prey development rate
           bval,# A fitted constant, part of the equation for attack rate from Sentis et al. 2012
           T0val,# The lower temperature bound for attack rate
           Tlval,# The upper temperature bound for attack rate
           h0val,# A fitted constant, part of the equation for handling time from Sentis et al. 2012
           bhval,# The constant for body mass-scaling, 0.75
           Ehval, # A fitted activation energy, part of the equation for handling time from Sentis et al. 2012
           Mval, # Avg. predator body mass
           kval, # Part of the equation for handling time from Sentis et al. 2012; this is a constant set to 8.617e-5
           tval, # Time to evaluate the functional response for, default is 1
           tempval ) # Current temperature
  {
    emerge <- 0
    new_N_vec <- mat.or.vec( 1, length( N_vec ) ) # Creates a vector to store the prey densities for the next time step
    for ( i in 1:timesteps ){
      # The first loop runs the model for the number of timesteps specified above
      
      for ( j in 2:4 ){ # This loop handles prey development
        #At each timestep a percentage of our larvae in each box moves to the next box according to r, the development rate
        new_N_vec[ j ] <-
          N_vec[ j ] + rval * N_vec[ j - 1 ] - rval * new_N_vec[ j ] # Each cell gains from the cell on the left, and loses to the right.
      }
      
      new_N_vec[ 1 ] <-
        N_vec[ 1 ] - rval * N_vec[ 1 ] # The first cell can only lose individuals w/o reproduction
      emerge <-
        emerge + rval * N_vec[ 5 ] # Larvae that made it to pupation last time step emerge
      new_N_vec[ 5 ] <-
        N_vec[ 5 ] + rval * N_vec[ 4 ] - rval * N_vec[ 5 ] # The last cell gains from cell 4, then loses individuals that emerge
      
      # Calculate p, the proportion of prey eaten, based on current density of prey and temp
      preykilled <- rogers.sentis(
        N0 = sum( N_vec ),
        b = bval,
        T0 = T0val,
        Tl = Tlval,
        h0 = h0val,
        M = Mval,
        bh = bhval,
        Eh = Ehval,
        k = kval,
        TempK = tempval + 273.15,
        t = tval
      )
      
      p <- ifelse( sum( new_N_vec ) == 0, 0, preykilled / sum( new_N_vec ) )
      p <- ifelse( p > 1, 1, p )

      # Predation/mortality
      
      new_N_vec <-
        new_N_vec - ( p * new_N_vec )
      new_N_vec <- 
        new_N_vec - ( mval * new_N_vec )
      N_vec <-
        new_N_vec #At the end of the loop, N_vec (the number of larvae at the most recent time step) is updated with the new values stored in new_N_vec.
    }
    return( emerge )
  }
