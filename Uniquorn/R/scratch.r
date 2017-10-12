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

match_t = match_t[match_t$All_variants != 0,]
matches = as.integer(as.character(match_t$Matches))
all_vars = as.integer(as.character(match_t$All_variants))
all_vars = all_vars[matches  > 0]
ccl_names = match_t$CCL[matches > 0]
matches = matches[matches > 0]

quot = as.vector( round( matches / all_vars, 2 ) * 100 )
names(quot ) = ccl_names

quot_s = sort(quot, decreasing = T)
quot_s

hist((log2(quot_s)))

difs <<- c()
for (i in 1:(length(quot_s)-1))
    difs = c(difs, quot_s[i] - quot_s[i + 1] )

hist(log2(difs+1))
ds = log2(difs+1)
mean(ds)
var(ds)
plot(dnorm( seq(0,1,by=.1), mean = mean(ds), sd = sd(ds)) )

mean_match = mean( matching_variants )
max_matche  = max( matching_variants )

len_q = length(g_query$Member_CCLs)

match_t_pos = match_t[ ! is.na(match_t$All_variants) ,]
len_r = as.integer(as.character(match_t_pos$All_variants))

pen <<- c()

for (i in 1:length(match_t_pos$CCL)){
    print(i)  
    pen = c(pen,integrate(
      f = pbeta,
      0,
      1,
      len_q,
      len_r[i],
      stop.on.error = FALSE
    )$value)
}


}
