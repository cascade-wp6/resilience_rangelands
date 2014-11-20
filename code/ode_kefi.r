# system of ODEs

d_rho_pm <- function(rho, parms, z = 4) { 
  with(parms,                             
       d * (rho[4] - rho[1] - rho[2]) + 
         (rho[5] - rho[3] - rho[1]) * 
         (delta*rho[4] + (z-1)/z * (1 - delta) * 
            (rho[4] - rho[1] - rho[2])/(1 - rho[4] - rho[5]) ) *
         (b - c*rho[4]) - 
         rho[1]*( r + f/z + (z-1)/z * f * rho[1]/rho[5] + m)
  )
}
d_rho_pp <- function(rho, parms, z = 4) { 
  with(parms,                             
       2 * (rho[4] - rho[1] - rho[2]) * 
         (delta*rho[4] + (1-delta)/z +(z-1)/z * (1 - delta) * 
            (rho[4] - rho[1] - rho[2])/(1 - rho[4] - rho[5]) ) *
         (b - c*rho[4]) - 2*rho[2]*m
  )
}
d_rho_mm <- function(rho, parms, z = 4) { 
  with(parms,                             
       2 * (rho[5] - rho[3] - rho[1]) - 
         2*rho[3]*(r+(z-1)/z*f*rho[1]/rho[5])
  )
}
d_rho_p <- function(rho, parms, z = 4) { 
  with(parms,                             
       (delta*rho[4] + (1 - delta) * 
          (rho[4] - rho[1] - rho[2])/(1 - rho[4] - rho[5]) ) *
         (b-c*rho[4])*(1-rho[4]-rho[5])-m*rho[4]
  )
}

d_rho_m <- function(rho, parms, z = 4) { 
  with(parms,                             
       d*(1-rho[4]-rho[5])-(r+f*rho[1]/rho[5])*rho[5]
  )
}

odesys <- function (t, rho, parms = model_parms) {
  list(c(
    d_rho_pm(rho, parms),
    d_rho_pp(rho, parms),
    d_rho_mm(rho, parms),
    d_rho_p(rho, parms),
    d_rho_m(rho, parms)
  ))
}
