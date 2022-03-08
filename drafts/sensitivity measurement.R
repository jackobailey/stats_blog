z <- rnorm(1000)
x <- rnorm(1000) + z
y <- rnorm(1000) + 3*z - x
z.observed <- z + rnorm(1000)
cor(z.observed, z)
summary(lm(y~x))
summary(lm(y~x + z.observed))
summary(lm(y~x + z))



socialbackground <- rnorm(1000)
drinking <- rnorm(1000) + socialbackground
partnerbackground <- socialbackground + rnorm(1000)
partnerdrinking <- partnerbackground + rnorm(1000)

childiq <- 2.5*socialbackground - 1 * drinking + rnorm(1000)

summary(lm(childiq~drinking))
summary(lm(childiq~partnerdrinking + drinking))


z <- rnorm(1000)

x <- rnorm(1000) + z

y <- rnorm(1000) + 3*z - x
z.observed <- z + rnorm(1000)
cor(z.observed, z)

summary(lm(y~x))
summary(lm(y~x + z.observed))
summary(lm(y~x + z))





y  <- x * -0.5
y