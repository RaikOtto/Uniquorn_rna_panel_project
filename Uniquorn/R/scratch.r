if(F){

library(MASS)

ints = as.integer( as.character( lib_match$Matches))
white_balls_possible = lib_match$All_variants
ints = sort(ints, decreasing = T)

difs <<- c()
for (i in 1:(length(ints)-1))
    difs = c(difs, ints[i] - ints[i + 1] )
difs_log = log(difs + 1)
difs_log_r = round(difs_log,0)
f = fitdistr(round(difs_log), "Negative Binomial")
qnbinom(.95,size = f$estimate[1], mu =  f$estimate[2])


mean(difs_log)
hist(scale(difs_log))

mu = mean(difs_log)
dif_sc = scale(difs_log)
sdd = sd(difs)

hist(scale(rr))

t.test(x = difs[1],mu = mu)

plot( pnorm(seq(mean(difs)-100,mean(difs)+100,by=.1), mean = mean(difs), sd = sd(difs)) )

x = dnbinom( x= 0:11,size= 100, prob = .5)


rr = rnbinom(size = f$estimate[1], mu =  f$estimate[2], n = length(difs_log))

ks.test(rr,difs_log_r)
plot(cumsum(difs_log_r))
scale(difs_log_r)



}