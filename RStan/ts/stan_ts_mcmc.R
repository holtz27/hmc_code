library(rstan)
library(stringr)

options(mc.cores = parallel::detectCores())
path = 'Área de Trabalho/Mestrado/Projeto/Projeto II/Simulação/algoritmo/rstan/ts/stan_ts_model.stan'
model  = stan_model(path)

# Running stan code
T = length(y)
M  = 5e3
fit = sampling(model, 
               list(T = T, 
                    y = y,
                    y0 = 0),
               iter = M,
               chains = 3,
               control = list(metric = 'unit_e')
               )

#------------------ Save outputs
#saveRDS(parameters, file = 'output_model_V')
#parameters = readRDS('output_model_V')

#------------------ Plots
parameters = extract(fit)
x = summary(fit)
Y = x$summary
#------------------ Tratando Y
Y = data.frame(Y, row.names = row.names(Y))
#------------------ Avaliando convergência
rows = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v')
cols = c('mean', 'sd', 'X2.5.', 'X97.5.', 'n_eff','Rhat')
data = Y[rows, cols]
data

# h
rhat_h = Y[str_detect(row.names(Y), pattern = '^h'),][, 'Rhat']
plot(rhat_h, main = 'h')
abline(h = 0.95)
abline(h = 1)
abline(h = 1.05)

# l
rhat_l = Y[str_detect(row.names(Y), pattern = '^l'),][, 'Rhat']
plot(rhat_l, main = 'l')
abline(h = 0.95)
abline(h = 1)
abline(h = 1.05)

# mu
mu = parameters$mu
par(mfrow=c(1,3))
plot(mu, type='l', main='', xlab='')
hist(mu, main = round(mean(mu), 3), xlab='', col = 'white', breaks = 40)
abline(v = mean(mu), col = 'blue', lwd=2, lty=2)
abline(v = Y['mu', ]['X2.5.'], col = 'red', lwd=2, lty=2)
abline(v = Y['mu', ]['X97.5.'], col = 'red', lwd=2, lty=2)
plot(acf(mu, lag.max = 100, plot = FALSE)[1:100])
par(mfrow=c(1,1))  

# phi
phi = parameters$phi
par(mfrow=c(1,3))
plot(phi, type='l', main='', xlab='')
hist(phi, main = round(mean(phi), 3), xlab='', col = 'white', breaks = 40)
abline(v = mean(phi), col = 'blue', lwd=2, lty=2)
abline(v = Y['phi', ]['X2.5.'], col = 'red', lwd=2, lty=2)
abline(v = Y['phi', ]['X97.5.'], col = 'red', lwd=2, lty=2)
plot(acf(phi, lag.max = 100, plot = FALSE)[1:100])
par(mfrow=c(1,1))

# sigma
s = parameters$sigma
par(mfrow=c(1,3))
plot(s, type='l', main='', xlab='')
hist(s, main = round(mean(s), 3), xlab='', col = 'white', breaks = 40)
abline(v = mean(s), col = 'blue', lwd=2, lty=2)
abline(v = Y['sigma', ]['X2.5.'], col = 'red', lwd=2, lty=2)
abline(v = Y['sigma', ]['X97.5.'], col = 'red', lwd=2, lty=2)
plot(acf(s, lag.max = 100, plot = FALSE)[1:100])
par(mfrow=c(1,1)) 

################################################################################
############### Análise numérica
theta = matrix( c( mu, phi, s ), ncol = 3 )
mcmcchain_theta = coda::as.mcmc( theta )
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_theta = coda::geweke.diag( mcmcchain_theta )
CD_theta
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_theta = coda::effectiveSize( mcmcchain_theta )
IF_theta = nrow(theta) / N_eff_theta
IF_theta
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_theta = round( apply( theta, MARGIN = 2, FUN = sd) / sqrt( N_eff_theta ), 5 )
mc_error_theta
################################################################################

# b0
b0 = parameters$b0
par(mfrow=c(1,3))
plot(b0, type='l', main='', xlab='')
hist(b0, main=round(mean(b0), 8), xlab='', col = 'white', breaks = 40)
abline(v = mean(b0), col = 'blue', lwd=2, lty=2)
abline(v = Y['b0', ]['X2.5.'], col = 'red', lwd=2, lty=2)
abline(v = Y['b0', ]['X97.5.'], col = 'red', lwd=2, lty=2)
plot(acf(b0, lag.max = 100, plot = FALSE)[1:100])
par(mfrow=c(1,1)) 

# b1
b1 = parameters$b1
par(mfrow=c(1,3))
plot(b1, type='l', main='', xlab='')
hist(b1, main=round(mean(b1), 8), xlab='', col = 'white', breaks = 40)
abline(v = mean(b1), col = 'blue', lwd=2, lty=2)
abline(v = Y['b1', ]['X2.5.'], col = 'red', lwd=2, lty=2)
abline(v = Y['b1', ]['X97.5.'], col = 'red', lwd=2, lty=2)
plot(acf(b1, lag.max = 100, plot = FALSE)[1:100])
par(mfrow=c(1,1)) 

# b2
b2 = parameters$b2
par(mfrow=c(1,3))
plot(b2, type='l', main='', xlab='')
hist(b2, main=round(mean(b2), 8), xlab='', col = 'white', breaks = 40)
abline(v = mean(b2), col = 'blue', lwd=2, lty=2)
abline(v = Y['b2', ]['X2.5.'], col = 'red', lwd=2, lty=2)
abline(v = Y['b2', ]['X97.5.'], col = 'red', lwd=2, lty=2)
plot(acf(b2, lag.max = 100, plot = FALSE)[1:100])
par(mfrow=c(1,1)) 
################################################################################
############### Análise numérica
b = matrix( c( b0, b1, b2 ), ncol = 3 )
mcmcchain_b = coda::as.mcmc( b )
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_b = coda::geweke.diag( mcmcchain_b )
CD_b
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_b = coda::effectiveSize( mcmcchain_b )
IF_b = nrow(b) / N_eff_b
IF_b
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_b = round( apply( b, MARGIN = 2, FUN = sd) / sqrt( N_eff_b ), 5 )
mc_error_b
################################################################################

# h
chain_h = parameters$h
h_hat = apply(chain_h, 2, mean)
lim_inf = apply(chain_h, 2, quantile, probs=c(.025, .975))[1,]
lim_sup = apply(chain_h, 2, quantile, probs=c(.025, .975))[2,]
df.1 = data.frame(x=1:T, ver = h, est = h_hat, 
                lim_inf = lim_inf, lim_sup = lim_sup)
g = ggplot(df[1:250, ]) 
g = g + geom_ribbon(aes(x, ymin = lim_inf, ymax = lim_sup), fill="gray80")
g = g + geom_line(aes(x=x, y = ver), color ='red')
g = g + geom_line(aes(x=x, y=est))
#g = g + geom_line(aes(x=x, y = lim_inf), linetype = 'dashed')
#g = g + geom_line(aes(x=x, y = lim_sup), linetype = 'dashed')
g
############### Análise numérica
mcmcchain_h = coda::as.mcmc( chain_h ) 
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_h = coda::geweke.diag( mcmcchain_h )
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_h = coda::effectiveSize( mcmcchain_h )
IF_h = M / N_eff_h
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_h = round( apply( chain_h, MARGIN = 2, FUN = sd) / sqrt( N_eff_h ), 5 )

# plots
par( mfrow = c(1,3) )
plot( CD_h$z )
abline(h = -1.96)
abline(h = 1.96)
plot( IF_h )
abline(h = 1)
plot( mc_error_h )
par( mfrow = c(1,1) )

################################################################################

# l
chain_l = parameters$l
l_hat = apply(chain_l, 2, mean)
lim_inf = apply(chain_l, 2, quantile, probs=c(.025, .975))[1,]
lim_sup = apply(chain_l, 2, quantile, probs=c(.025, .975))[2,]
df.2 = data.frame(x=1:T, ver = l, est = l_hat, 
                lim_inf = lim_inf, lim_sup = lim_sup)
h = ggplot( df.2[1:100, ] ) 
h = h + geom_ribbon(aes(x, ymin = lim_inf, ymax = lim_sup), fill="gray80")
h = h + geom_line(aes(x=x, y = ver), color ='red')
h = h + geom_line(aes(x=x, y=est))
#g = g + geom_line(aes(x=x, y = lim_inf), linetype = 'dashed')
#g = g + geom_line(aes(x=x, y = lim_sup), linetype = 'dashed')
h

############### Análise numérica
mcmcchain_l = coda::as.mcmc( chain_l ) 
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_l = coda::geweke.diag( mcmcchain_l )
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_l = coda::effectiveSize( mcmcchain_l )
IF_l = M / N_eff_l
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_l = round( apply( chain_l, MARGIN = 2, FUN = sd) / sqrt( N_eff_l ), 5 )

# plots
par( mfrow = c(1,3) )
plot( CD_l$z )
abline(h = -1.96)
abline(h = 1.96)
plot( IF_l )
abline(h = 1)
plot( mc_error_l )
par( mfrow = c(1,1) )
