

load(file = 'data/z_values.rda')

max(genuine_z_values);

min(genuine_z_values);

mean(genuine_z_values);

sqrt(var(genuine_z_values));

N <- length(genuine_z_values);

###########################################################
# Compute the quantity Ahat from Eq. (55) in Efron (2007) #
###########################################################
z = genuine_z_values

Ahat <- function(x0){
        Y0 <- sum((z < x0) & (z > -x0));
        P0 <- 2 * pnorm(x0) - 1;
        Q0 <- sqrt(2) * x0 * dnorm(x0);
        P0hat <- Y0 / N;
        return((P0-P0hat) / Q0);
}

Ahat(x0 = 1.0)

##########################################
# Plot Ahat(x0) as a function of x0      #
##########################################
Ahatv <- Vectorize(Ahat);

curve(Ahatv, from=0.5, to=3.0);

##########################################
#              Trimming                  #
##########################################
z_ordered <- sort(z)
trimm <- 0.1
N_lower <- ceiling(N*trimm)
N_upper <- floor(N*(1.0-trimm))
z_trimm <- z_ordered[N_lower:N_upper]

sqrt(var(z_trimm))
