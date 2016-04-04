
 
setwd('C:\\Users\\felipe.farias\\Projetos\\Forecast')

install.packages('reshape')
install.packages('tseries')
install.packages('TSA')
install.packages('forecast')
#install.packages('fGarch')
install.packages('nortest')
install.packages('R2HTML')
install.packages('gplots')

 
library (lattice)
library (reshape)
library('tseries')
#library('fractal')
library('TSA')
library('fGarch')
library('forecast')
library('nortest') # testes de normalidade
library('R2HTML')
library('tgp'); library('maptree'); # Gaussian 

library('gplots') # Saida em PDF

source('script_boot.R')

#-----------------------------------------------------------



d.df <- read.csv(file="rolling_forecast_r.csv",header=TRUE,sep=";",dec=",");
d.df[5,112] <- 100
str(d.df)
dim(d.df)


#a <- diffinv(runif(1000,-1,2)); plot(a)

d.periodo.forecast <- 15;
d.column.name = 'Produ??o de Gasolina'
d.ylab = 'Gasolina'
d.xlab = 'Tempo'
data(gas)
d <- gas

#-----------------------------------------------------------

d.periodo.forecast <- 15;
d.column.name = 'Teste'
d.ylab = 'unidades'
d <- NULL
d <- ts(scan(dec = ",", sep=";"), end=2011+8/12, frequency=12)





relatorio.pdf <- function () {

tryCatch(sinkplot("cancel"), error = function(e) {} ) #cancela caso tenha ficado a captura de tela em aberto
d.par.old <- par()

d.l = log(d)
d.d = diff(d)
d.dl = diff(log(d))
d.dd = diff(diff(d))
d.ddl = diff(diff(log(d)))
d.dpercent = diff(log(d))*100
d.media = mean(d)
d.sd = sd(d)

# forecast
d.ets <- ets(d)
d.arima <- auto.arima(d, trace=TRUE, stepwise=TRUE, allowdrift=TRUE); 
#d.arfima <- arfima(d, trace=TRUE, stepwise=FALSE, allowdrift=TRUE); 
d.forecast.ets <- forecast(d.ets,h=d.periodo.forecast)
d.forecast.arima <- forecast(d.arima,h=d.periodo.forecast)
#try(d.forecast.arfima <- forecast(d.arfima,h=d.periodo.forecast))


write.table(d.forecast.ets, file= paste(d.column.name,'projecao_ets.csv',sep=''),sep=';',dec=',')
write.table(d.forecast.arima, file= paste(d.column.name,'projecao_arima.csv',sep='') ,sep=';',dec=',')

d.sazonalidade = d.ets$states[, "s1"]
write.table(d.sazonalidade, file= paste(d.column.name,'projecao_ets_sazonalidade.csv',sep=''),sep=';',dec=',')


d.gp.saida <- ts(seq(time(d)[1], time(d)[length(d)] + (d.periodo.forecast/12), by=1/frequency(d)),start=start(d), frequency=frequency(d))

d.gp <- btgpllm(X = time(d), Z = d, XX = d.gp.saida, verb = 1)
#d.gp <- btlm(X = time(d), Z = d, XX = d.gp.saida, verb = 1)
#d.gp <- btgp(X = time(d), Z = d, XX = d.gp.saida, verb = 1)

#d.gp <- btlm(X = time(d), Z = d, XX = XX, verb = 1)
#write.table(d.btgpllm, file='analise_gaussiana.csv',sep=';',dec=',')


#Cria o PDF
pdf(paste('Relatorio_Forecast_',d.column.name,'.pdf',sep=''))

#Runplot
par(mfrow=c(1,1));
plot(d, main=d.column.name, xlab='Tempo', ylab= d.ylab)

#Derivadas
par(mfrow=c(3,2));
#plot.div(3,2)
plot.rootunit(d, main=d.column.name, xlab=d.xlab, ylab= d.ylab)
plot.rootunit(d.l, main='Log', xlab=d.xlab, ylab= paste('log(',d.ylab,')'))
plot.rootunit(d.d, main='Diferen?a', xlab=d.xlab, ylab=paste('diff(',d.ylab,')'))
plot.rootunit(d.dl, main='Diferen?a dos Logs', xlab=d.xlab, ylab=paste('diff(log(',d.ylab,'))'))
plot.rootunit(d.dd, main=expression('Diferen?a'^2~' '), xlab=d.xlab, ylab=paste('diff(diff(',d.ylab,'))'))
plot.rootunit(d.ddl, main=expression('Diferen?a'^2~' dos Logs'), xlab=d.xlab, ylab=paste('diff(diff(log(',d.ylab,')))'))

-----------------------------------------------
--- McLeod-Li test for conditional heteroscedascity (ARCH)
-----------------------------------------------



par(mfrow=c(3,2))
McLeod.Li.test(y=d, main='McLeod-Li test (d)')
McLeod.Li.test(y=d.l, main='McLeod-Li test (d.l)')
McLeod.Li.test(y=d.d, main='McLeod-Li test (d.d)')
McLeod.Li.test(y=d.dl, main='McLeod-Li test (d.dl)')
McLeod.Li.test(y=d.dd, main='McLeod-Li test (d.dd)')
McLeod.Li.test(y=d.ddl, main='McLeod-Li test (d.ddl)')

#Sazonalidade

par(mfrow=c(1,1));
seasonplot(d, main='Gr?fico da Sazonalidade', xlab='M?s', ylab= d_ylab)
x11()
monthplot(d, main='Gr?fico de Crescimento Mensal', xlab='M?s', ylab= d_ylab)

par(mfrow=c(1,1));
monthplot(d, ylab= d.ylab)

# Decomposi??o por Loess (Local Polynomial Regression Fitting)

par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0), mfcol = c(4, 2))
d.stl <- stl(d, s.window = frequency(d), t.window = frequency(d), robust = FALSE)
d.stl.robust <- stl(d, s.window = frequency(d), t.window = frequency(d), robust = TRUE)
write.table(d.stl$time.series, file= paste(d.column.name,'decomposicao_sazonalidade.csv',sep=''),sep=';',dec=',')
write.table(d.stl.robust$time.series, file= paste(d.column.name,'decomposicao_sazonalidade_robusta.csv',sep=''),sep=';',dec=',')
plot(d.stl, set.pars=NULL, labels = NULL, main = "Decomposi??o por Loess (Normal / Robusta )")
plot(d.stl.robust, set.pars=NULL)
# mark the 'outliers' :
(iO <- which(d.stl.robust$weights  < 1e-8))
sts <- d.stl.robust$time.series
points(time(sts)[iO], 0.8* sts[,"remainder"][iO], pch = 4, col = "red")


par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0), mfcol = c(2, 1))
monthplot(d.stl.robust, choice = "seasonal", cex.axis = 0.8)
monthplot(d.stl.robust, choice = "trend", cex.axis = 0.8)

	 
#monthplot(d, ylab = "data", cex.axis = 0.8)	 
#par(mfrow = c(3,1))
#monthplot(d.stl, choice = "seasonal", cex.axis = 0.8)
#monthplot(d.stl, choice = "trend", cex.axis = 0.8)
#monthplot(d.stl, choice = "remainder", type = "h", cex.axis = 0.8)


#GP
par(mfrow = c(1,1))
plot(d.gp, main = "An?lise Gaussiana/Bayesiana - ", layout = "surf")
lines(d, col = 4, lty = 2, lwd = 1)


par(mfrow = c(1,1))
plot(d, col = 4, lty = 3, lwd = 1, main="An?lise Gaussiana/Bayesiana - Forecast", xlab=d.xlab, ylab= d.ylab )
lines(d.gp$XX[,1], d.gp$ZZ.mean)
write.table(cbind("data"=d.gp$XX[,1], "valor"=d.gp$ZZ.mean), file= paste(d.column.name,'_gp.csv',sep=''),sep=';',dec=',')

#plot(d.gp$XX[,1], d.gp$ZZ.mean, type='l')
#lines(d, col = 4, lty = 2, lwd = 1)


par(mfrow = c(1,1))
plot(d.gp, main = "An?lise Gaussiana/Bayesiana - ", layout = "as")

par(mfrow = c(1,1))
tgp.trees(d.gp, nodeinfo=FALSE, heights = "map", cex=0.8)


# ETS
par(d.par)
plot(d.forecast.ets, main=d.column.name, sub=paste("Modelo: ", d.forecast.ets$method, "", sep = ""), col="black");
lines(d.forecast.ets$fitted, col="blue");

# ETS Summary
par(mfrow=c(1,1), mar=c(1.0,1.0,1.0,1.0), omi=c(0, 0, 0, 0));
try(sinkplot());
try(summary(d.forecast.ets));
try(sinkplot("plot",col="black", cex=0.7));
try(title(d.column.name));

# ETS Decomposicao
par(d.par)
try(plot(d.ets));

# ARIMA Grafico
par(d.par)
plot(d.forecast.arima, main="Previs?o Arima", sub=paste("Modelo: ", d.forecast.arima$method, "", sep = ""), col="black")
lines(d.forecast.arima$fitted, col="blue")

# ARIMA Summary
par(mfrow=c(1,1), mar=c(1.0,1.0,1.0,1.0), omi=c(0, 0, 0, 0) );
try(sinkplot());
summary(d.forecast.arima)
try(sinkplot("plot",col="black", cex=0.7));
try(title(d.column.name));

---------------
d.season. = season(d)
d.season.lm  = lm(d~d.season.-1) # sem o intercept
d.season.lm2 = lm(d~d.season.)   # com o intercept
summary(d.season.lm)
summary(d.season.lm2)

x11();par(mfrow=c(2,2))
plot(d.season.lm)
x11();par(mfrow=c(2,2))
plot(d.season.lm2)
---------------

d.lm <- lm(d~time(d))
d.lm.summary <- summary(d.lm)

# LM Summary
par(d.par)
plot(d, main="Regress?o linear", xlab="Tempo", ylab= d.ylab)
abline(d.lm,col="blue")
mtext(paste("Reta ajustada: Y=",d.lm$coefficients[1]," + x*",d.lm$coefficients[2],sep=""), col="blue", adj=0)
mtext(substitute(paste(chi^2, "=", k), list(k = d.lm.summary$r.squared)), col="blue", adj=1)

# LM Summary
par(mfrow=c(1,1), mar=c(1.0,1.0,1.0,1.0), omi=c(0, 0, 0, 0) );
try(sinkplot());
print(d.lm.summary)
try(sinkplot("plot",col="black", cex=0.7));
try(title(d.column.name));

# LM Graficos
par(d.par)
par(mfrow=c(3,2))
plot(d.lm, which = c(1:6))

# ARIMA X ETS
par(mfrow=c(2,2))
plot(d.forecast.ets, main='Previs?o ETS', sub=paste('Modelo: ', d.forecast.ets$method, '', sep = ''), col='black')
lines(d.forecast.ets$fitted, col='blue')

plot(d.forecast.arima, main='Previs?o Arima', sub=paste('Modelo: ', d.forecast.arima$method, '', sep = ''), col='black')
lines(d.forecast.arima$fitted, col='blue')

barplot(rbind(accuracy(d.ets)[1:3],accuracy(d.arima)[1:3]),beside=TRUE, col=c('gray10','gray80'))
barplot(rbind(accuracy(d.ets)[4:6],accuracy(d.arima)[4:6]),beside=TRUE, col=c('gray10','gray80'))
legend('topleft', legend=c('ETS','ARIMA'), text.col=c('gray10','gray50'), bg='white')


# Residuos ETS
par(mfrow=c(2,2))
acf(d.ets$residuals, main="Correla??o dos Residuos - ETS")
pacf(d.ets$residuals, main="Correla??o Parcial dos Residuos - ETS")
d.ets.residuals.std <- residuals(d.ets)/d.ets$sigma2^0.5
plot(d.ets.residuals.std, main = "Standardized Residuals - ETS", ylab="Std.Dev. of Residuals", type = "h")
abline(h=0)
try(qqnorm(d.ets$residuals));
try(qqline(d.ets$residuals));

# Residuos ARIMA
par(mfrow=c(2,2))
acf(d.arima$residuals, main="Correla??o dos Residuos - Arima")
pacf(d.arima$residuals, main="Correla??o Parcial dos Residuos - Arima")
d.arima.residuals.std <- residuals(d.arima)/d.arima$sigma2^0.5
plot(d.arima.residuals.std, main = "Standardized Residuals - Arima", ylab="Std.Dev. of Residuals", type = "h")
abline(h=0)
try(qqnorm(d.arima$residuals));
try(qqline(d.arima$residuals));

dev.off()
d.saida <- cbind("ETS"=as.numeric(d.forecast.ets$mean),"Arima"=as.numeric(d.forecast.arima$mean))
write.table(d.saida, file=paste('projecao_colunas_',d.column.name,'.csv',sep=''),sep=';',dec=',')

}
