#LIBRERIAS Y PAQUETES ----------------------------------------------------------
rm(list=ls(all=TRUE))
library(forecast)
library(lmtest)
library(fANCOVA)
library(nlme)
library(TSA)
library(car)

source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Mipoly.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Mytrigon.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-BP.LB.test-pruebaDW1.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-SuavizamientoEstacional.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Descomp.Loessv02.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencialv02.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexpo.ErrorARMA.R") 
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-SelectModel.R") 
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-interpdeltas.R")

#Leer anex-EMMET-TotalNacional-oct2023-Elaboracion de Bebidas.csv, columna 5: Produccion nominal
Datosx=read.table(file.choose(),header=T,sep=";",skip=14,dec=",",colClasses=c(rep("NULL",7),"numeric",rep("NULL",3)))
Datosx=ts(Datosx,freq=12,start=c(2001,1))
win.graph()
plot(Datosx)

#--------------------------------PUNTO 2a: ANALISIS DESCRIPTIVO DE LA SERIE---------------------------------

#Grafica de la serie
win.graph()
plot(Datosx,ylab="Datosx")

#Grafica ACF estimada con la serie: m = 36 en el caso mensual 
win.graph()
acf(as.numeric(Datosx),lag.max=36,ci.type="ma",col=4,ci.col=2)

#--------------------------------PUNTO 2b: MODELO DE REGRESION GLOBAL-----------------------------------------

m=12
n=length(Datosx)-m
t=1:n
yt=ts(Datosx[t],freq=m,start=c(2001,1))
poli=Mipoly(tiempo=t,grado=4)

#-#-#-#-# #-#-#-#-#   IMPORTANTE    #-#-#-#-#  #-#-#-#-#  

#mes=seasonaldummy(yt)   # Utilizo esta si es con Indicadoras
trigon=Mytrigon(tiempo=t,Frecuencias=c(c(1,2,3,4,5)/12),indicej=c(1,2,3,4,5))  #Utilizo esta si es con trigonometricas

#Matriz de diseño en el ajuste

#-#-#-#-# #-#-#-#-#   IMPORTANTE    #-#-#-#-#  #-#-#-#-#  

#X1=data.frame(poli,mes)   # Utilizo esta si es con Indicadoras
X1=data.frame(poli,trigon)  # Utilizo esta si es con trigonometricas
head(X1) 

#PARA LOS PRONOSTICOS
tnuevo=(n+1):length(Datosx)
ytnuevo=ts(Datosx[tnuevo],freq=m,start=c(2023,8))
polinuevo=Mipoly(t=tnuevo,grado=3) #Cambiar Grado del polinomio

#-#-#-#-# #-#-#-#-#   IMPORTANTE    #-#-#-#-#  #-#-#-#-#  
trigonnuevo=Mytrigon(tiempo=tnuevo,Frecuencias=c(c(1,2,3,4,5)/12),indicej=c(1,2,3,4,5)) #Si es con trigonometricas
#mesnuevo=seasonaldummy(yt,h=m)  #Si es con Indicadoras#Si es con indicadoras

#-#-#-#-# #-#-#-#-#   IMPORTANTE    #-#-#-#-#  #-#-#-#-#  
#Matriz de diseño en los pronosticos
X1nuevo=data.frame(polinuevo,trigonnuevo)  #Si es con trigonometricas
#X1nuevo=data.frame(polinuevo,mesnuevo) #Si es con Indicadoras#Si es con indicadoras
head(X1nuevo) 


#PARA LOS PRONOSTICOS
tnuevo=(n+1):length(Datosx)
ytnuevo=ts(Datosx[tnuevo],freq=m,start=c(2023,8))
polinuevo=Mipoly(t=tnuevo,grado=3) #Cambiar Grado del polinomio
trigonnuevo=Mytrigon(tiempo=tnuevo,Frecuencias=c(c(1,2,3,4,5)/12),indicej=c(1,2,3,4,5))

#Matriz de diseño en los pronosticos
X1nuevo=data.frame(polinuevo,trigonnuevo)
head(X1nuevo) 



# AJUSTE DEL MODELO

modglobal=lm((yt)~.,data=X1)
summary(modglobal)
#Calculo valores ajustados del modelo 1
ythatmod1=ts(fitted(modglobal),freq=m,start=start(yt))

modglobal=lm((yt)~.,data=X1)
summary(modglobal)
tabla.parametros.globales=summary(modglobal)$coefficients
tabla.parametros.globales
write.csv2(tabla.parametros.globales,file="tablamod1ymod2trabajo1.csv",row.names = TRUE)

#Calculo valores ajustados del modelo global
ythatglobal=ts(fitted(modglobal),freq=m,start=start(yt))

#Calculo de los criterios AIC y BIC en modelo global
nparmodglobal=length(coef(modglobal)[coef(modglobal)!=0]);nparmodglobal 
Criteriosglobal= exp.crit.inf.resid(residuales=residuals(modglobal),n.par=nparmodglobal);Criteriosglobal

#Grafico del ajuste
win.graph()
plot(Datosx, ylab="Datosx")
lines(ythatglobal,col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo global"),lty=1,col=c(1,2))

#Pronosticos del modelo global en la escala original
pronosglobal=predict(modglobal,newdata=X1nuevo,interval="prediction",level=0.95)
pronosglobal=ts(pronosglobal,freq=m,start=start(ytnuevo))
pronosglobal
ytpronglobal=pronosglobal[,1]

#precision pronosticos puntuales modelo global
accuracy(ytpronglobal,ytnuevo)
#precision pronosticos por I.P modelo 1
amplcobmodglobal=amplitud.cobertura(real=ytnuevo,LIP=pronosglobal[,2],LSP=pronosglobal[,3]);amplcobmodglobal
#score promedio pronosticos por I.P modelo 1
ScoreIPmodglobal=IntervalScore(real=ytnuevo,LIP=pronosglobal[,2],LSP=pronosglobal[,3],alpha=0.05);ScoreIPmodglobal



#--------------------- PUNTO 3: Evaluación de los supuestos de Ruido Blanco para los errores-------------------------------

#--------------------- PUNTO 3a: Validación de supuestos usando los residuos estructurales

#1) Test Durbin Watson
pruebaDW1(modglobal)

#2) 

# Residuos vs. tiempo
win.graph()
plot.ts(residuals(modglobal),ylab="Residuales")
abline(h=c(-2*summary(modglobal)$sigma,0,2*summary(modglobal)$sigma),lty=2,col=2)
legend("topleft",legend="Modelo global")

#Residuos vs. ajustados
win.graph()
plot(fitted(modglobal),residuals(modglobal),ylab="Residuales")
abline(h=c(-2*summary(modglobal)$sigma,0,2*summary(modglobal)$sigma),lty=2,col=2)
legend("topleft",legend="Modelo global")

#ACF: m = 36 en el caso mensual y m = 24 en el trimestral
win.graph()
acf(residuals(modglobal),ci.type="ma",lag.max=36)

#PACF: m = 36 en el caso mensual y m = 24 en el trimestral
win.graph()
pacf(residuals(modglobal),lag.max=36)

#Test LJUNG-BOX m = 6, 12, 18, 24, 30 y 36 en el caso mensual y m = 6, 12, 18 y 24 en el trimestral
BP.LB.test(residuals(modglobal),maxlag=36,type="Ljung")

#--------------------- PUNTO 3b: IDENTIFICACION DE POSIBLES MODELOS PARA Et--------------------------

#EACF
eacf(residuals(modglobal),ar.max = 36, ma.max = 36)  #Modificar si 24X24

#Identificación de modelos AR(p) con SelectModel --- SOLO SI ACF tienen patron tipo cola --
SelectModel(residuals(modglobal),lag.max=36,Criterion="AIC",ARModel ="AR")
SelectModel(residuals(modglobal),lag.max=36,Criterion="BIC",ARModel ="AR")

#Identificación con auto.arima sobre vector de residuales sin fechas
auto.arima(residuals(modglobal),ic="aic")
auto.arima(residuals(modglobal),ic="bic")

#Identificación con auto.arima sobre vector de residuales con fechas
serieEt=ts(residuals(modglobal),freq=m,start=c(2001,1))
auto.arima(serieEt,ic="aic")
auto.arima(serieEt,ic="bic")

#Identificación con armasubsets
plot(armasubsets(residuals(modglobal),nar=12,nma=12,y.name="AR",ar.method="ml"))

#--------------------- PUNTO 4: Modelos de regresion global con errores estructurales Et-------------

#-----MODELO 1:AR(15)-----------------------------------------------------------

modelo1=Arima(yt,order=c(15,0,0),xreg=as.matrix(X1),method="ML") 
k1=length(coef(modelo1)[coef(modelo1)!=0]);k1 
dfmodelo1=n-k1
coeftest(modelo1,df=dfmodelo1)
summary(modelo1)
ythat1=ts(fitted(modelo1),freq=m,start=start(yt)) 

#-----MODELO 2:ARMA(2,13)-------------------------------------------------------

modelo2=Arima(yt,order=c(2,0,13),xreg=as.matrix(X1),method="ML")
k2=length(coef(modelo2)[coef(modelo2)!=0]);k2 
dfmodelo2=n-k2
coeftest(modelo2,df=dfmodelo2) 
summary(modelo2)
ythat2=ts(fitted(modelo2),freq=m,start=start(yt))

#-----MODELO 3:ARMA(10,2)(1,0)[12]-----------------------------------------------

modelo3=Arima(yt,order=c(10,0,2),seasonal=list(order=c(1,0,0)),xreg=as.matrix(X1),method="ML") 
k3=length(coef(modelo3)[coef(modelo3)!=0]);k3 
dfmodelo3=n-k3;
coeftest(modelo3,df=dfmodelo3)
summary(modelo3)
ythat3=ts(fitted(modelo3),freq=m,start=start(yt))

#-----MODELO 4:ARMA(10,3) reglon4, + teta1------------------------------------------

modelo4=Arima(yt,order=c(10,0,3),fixed=c(0,NA,0,0,NA,0,0,0,0,NA,NA,NA,NA,rep(NA,15)),xreg=as.matrix(X1),method="ML") 
k4=length(coef(modelo4)[coef(modelo4)!=0]);k4  
dfmodelo4=n-k4
coeftest(modelo4,df=dfmodelo4) 
summary(modelo4)
ythat4=ts(fitted(modelo4),freq=m,start=start(yt))


#--GRAFICAS VALORES AJUSTADOS PARA LOS 4 MODELOS -------------------------------

#MODELO 1
win.graph()
plot(Datosx, ylab="Datos13")
lines(ythat1,col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 1"),lty=1,col=c(1,2))

#MODELO 2
win.graph()
plot(Datosx, ylab="Datos13")
lines(ythat2,col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 2"),lty=1,col=c(1,2)) 

#MODELO 3
win.graph()
plot(Datosx, ylab="Datos13")
lines(ythat3,col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 3"),lty=1,col=c(1,2))

#MODELO 4
win.graph()
plot(Datosx, ylab="Datos13")
lines(ythat4,col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 4"),lty=1,col=c(1,2))


#MEDIDAS DE AJUSTE

#Calculo de los criterios AIC y BIC en modelo 1
Res.orig.modelo1=yt-ythat1 
Criteriosmodelo1=exp.crit.inf.resid(residuales=residuals(modelo1),n.par=k1); Criteriosmodelo1
#Calculo de los criterios AIC y BIC en modelo 2
Res.orig.modelo2=yt-ythat2 
Criteriosmodelo2=exp.crit.inf.resid(residuales= Res.orig.modelo2,n.par=k2); Criteriosmodelo2
#Calculo de los criterios AIC y BIC en modelo 3
Res.orig.modelo3=yt-ythat3 
Criteriosmodelo3=exp.crit.inf.resid(residuales=Res.orig.modelo3,n.par=k3); Criteriosmodelo3 
#Calculo de los criterios AIC y BIC en modelo 4
Res.orig.modelo4=yt-ythat4 
Criteriosmodelo4=exp.crit.inf.resid(residuales=Res.orig.modelo4,n.par=k4); Criteriosmodelo4 

#Tabla resumen parametros ajustados 
parametrosmod1=coeftest(modelo1,df=dfmodelo1)
parametrosmod2=coeftest(modelo2,df=dfmodelo2)
parametrosmod3=coeftest(modelo3,df=dfmodelo3)
parametrosmod4=coeftest(modelo4,df=dfmodelo4)
tabla.parametros.ajuste=rbind(parametrosmod1,parametrosmod2,parametrosmod3,parametrosmod4)
tabla.parametros.ajuste
write.csv2(tabla.parametros.ajuste,file="Parametrosajustados.csv",row.names = TRUE)

#------------------------ PUNTO 5: ANALISIS DE RESIDUALES Y VALIDACION DE SUPUESTOS---------------------------------


#Analisis de residuales
win.graph()
plot(residuals(modelo1));abline(h=0)
abline(h=c(-2*sqrt(modelo1$sigma2),2*sqrt(modelo1$sigma2)),lty=2, col=2)
win.graph()
plot(residuals(modelo2));abline(h=0)
abline(h=c(-2*sqrt(modelo2$sigma2),2*sqrt(modelo2$sigma2)),lty=2, col=2)
win.graph()
plot(residuals(modelo3));abline(h=0)
abline(h=c(-2*sqrt(modelo3$sigma2),2*sqrt(modelo3$sigma2)),lty=2, col=2)
win.graph()
plot(residuals(modelo4));abline(h=0)
abline(h=c(-2*sqrt(modelo4$sigma2),2*sqrt(modelo4$sigma2)),lty=2, col=2)


win.graph()
plot(as.numeric(modelo1$fitted),residuals(modelo1));abline(h=0)
abline(h=c(-2*sqrt(modelo1$sigma2),2*sqrt(modelo1$sigma2)),lty=2, col=2)
win.graph()
plot(as.numeric(modelo2$fitted),residuals(modelo2));abline(h=0)
abline(h=c(-2*sqrt(modelo2$sigma2),2*sqrt(modelo2$sigma2)),lty=2, col=2)
win.graph()
plot(as.numeric(modelo3$fitted),residuals(modelo3));abline(h=0)
abline(h=c(-2*sqrt(modelo3$sigma2),2*sqrt(modelo3$sigma2)),lty=2, col=2)
win.graph()
plot(as.numeric(modelo4$fitted),residuals(modelo4));abline(h=0)
abline(h=c(-2*sqrt(modelo4$sigma2),2*sqrt(modelo4$sigma2)),lty=2, col=2)

#validacion de supuestos
win.graph()
acf(as.numeric(residuals(modelo1)),ci.type="ma",lag.max=36,ci.col=2)
win.graph()
acf(as.numeric(residuals(modelo2)),ci.type="ma",lag.max=36,ci.col=2)
win.graph()
acf(as.numeric(residuals(modelo3)),ci.type="ma",lag.max=36,ci.col=2) 
win.graph()
acf(as.numeric(residuals(modelo4)),ci.type="ma",lag.max=36,ci.col=2) 

win.graph()
pacf(as.numeric(residuals(modelo1)),lag.max=36,ci.col=2)
win.graph()
pacf(as.numeric(residuals(modelo2)),lag.max=36,ci.col=2)
win.graph()
pacf(as.numeric(residuals(modelo3)),lag.max=36,ci.col=2)
win.graph()
pacf(as.numeric(residuals(modelo4)),lag.max=36,ci.col=2)

#tabla resumen Test Ljung-box
tabla.LjungBox=cbind(BP.LB.test(residuals(modelo1),maxlag=36,type="Ljung"),BP.LB.test(residuals(modelo2),maxlag=36,type="Ljung"),BP.LB.test(residuals(modelo3),maxlag=36,type="Ljung"),BP.LB.test(residuals(modelo4),maxlag=36,type="Ljung"))[,c(1,3,4,6,7,9,10,12)]
colnames(tabla.LjungBox)=c("QLB-M1","vp-M1","QLB-M2","vp-M2","QLB-M3","vp-M3","QLB-M4","vp-M4")
tabla.LjungBox
write.csv2(tabla.LjungBox,file="LJungBox.csv",row.names = TRUE)

#normalidad
win.graph()
qqnorm(residuals(modelo1))
qqline(residuals(modelo1),col=2) 
win.graph()
qqnorm(residuals(modelo2))
qqline(residuals(modelo2),col=2) 
win.graph()
qqnorm(residuals(modelo3))
qqline(residuals(modelo3),col=2) 
win.graph()
qqnorm(residuals(modelo4))
qqline(residuals(modelo4),col=2) 

shapiro.test(residuals(modelo1))
#tabla resumen test shapiro
tabla.Shapiro=rbind(shapiro.test(residuals(modelo1)),shapiro.test(residuals(modelo2)),shapiro.test(residuals(modelo3)),shapiro.test(residuals(modelo4)))[,c(1,2)]
rownames(tabla.Shapiro)=c("Modelo1","Modelo2","Modelo3","Modelo4")
tabla.Shapiro


#---------------------------------- PUNTO 6: PRONOSTICOS PARA LA VALIDACION CRUZADA ----------------------------------------

#Modelo 1
predmodelo1=as.data.frame(forecast(modelo1,xreg=as.matrix(X1nuevo),level=95))
predmodelo1=ts(predmodelo1,freq=m,start=c(2023,8)) 
predmodelo1
ytpron1=predmodelo1[,1] #Tomando el pronóstico puntual
#Modelo 2
predmodelo2=as.data.frame(forecast(modelo2,xreg=as.matrix(X1nuevo),level=95)) 
predmodelo2=ts(predmodelo2,freq=m,start=c(2023,8)) 
predmodelo2
ytpron2=predmodelo2[,1] #Tomando el pronóstico puntual
#Modelo 3
predmodelo3=as.data.frame(forecast(modelo3,xreg=as.matrix(X1nuevo),level=95))
predmodelo3=ts(predmodelo3,freq=m,start=c(2023,8)) 
predmodelo3
ytpron3=predmodelo3[,1]
#Modelo 4
predmodelo4=as.data.frame(forecast(modelo4,xreg=as.matrix(X1nuevo),level=95))
predmodelo4=ts(predmodelo4,freq=m,start=c(2022,11)) #Fechas y frecuencia según datos en pronósticos
predmodelo4
ytpron4=predmodelo4[,1]

#Metricas de pronostico
#Modelo 1
accuracy(ytpron1,ytnuevo)
amplcobmodelo1=amplitud.cobertura(real=ytnuevo,LIP=predmodelo1[,2],LSP=predmodelo1[,3])
ScoreIP1=IntervalScore(real=ytnuevo,LIP=predmodelo1[,2],LSP=predmodelo1[,3],alpha=0.05);ScoreIP1
#Modelo 2
accuracy(ytpron2,ytnuevo)
amplcobmodelo2=amplitud.cobertura(real=ytnuevo,LIP=predmodelo2[,2],LSP=predmodelo2[,3])
ScoreIP2=IntervalScore(real=ytnuevo,LIP=predmodelo2[,2],LSP=predmodelo2[,3],alpha=0.05);ScoreIP2
#Modelo 3
accuracy(ytpron3,ytnuevo)
amplcobmodelo3=amplitud.cobertura(real=ytnuevo,LIP=predmodelo3[,2],LSP=predmodelo3[,3])
ScoreIP3=IntervalScore(real=ytnuevo,LIP=predmodelo3[,2],LSP=predmodelo3[,3],alpha=0.05);ScoreIP3
#Modelo 4
accuracy(ytpron4,ytnuevo)
amplcobmodelo4=amplitud.cobertura(real=ytnuevo,LIP=predmodelo4[,2],LSP=predmodelo4[,3])
ScoreIP4=IntervalScore(real=ytnuevo,LIP=predmodelo4[,2],LSP=predmodelo4[,3],alpha=0.05);ScoreIP4


#Grafica de los pronosticos
win.graph()
plot(ytnuevo,type="b",ylab="Datos13",col=1,pch=19,ylim=c(min(ytnuevo,ytpron1,ytpron2,ytpron3,ytpron4),max(ytnuevo,ytpron1,ytpron2,ytpron3,ytpron4)),lwd=2,xaxt="n")
axis(1,at=time(ytnuevo),labels=c("ago-23","sep-23","oct-23","nov-23","dic-23","ene-24","feb-24","mar-24","abr-24","may-24","jun-24","jul-24"),cex.axis=0.7)
lines(ytpron1,col=2,pch=1,type="b",lwd=2)
lines(ytpron2,col=3,pch=2,type="b",lwd=2)
lines(ytpron3,col=4,pch=3,type="b",lwd=2)
lines(ytpron4,col=5,pch=4,type="b",lwd=2)
legend("topleft",legend=c("Real","Modelo 1","Modelo 2","Modelo 3","Modelo 4"),pch=c(19,1:4),col=c(1:5),lwd=2)

#tabla resumen de medidas de pronosticos
tabla.pronosticos=cbind(predmodelo1,predmodelo2,predmodelo3,predmodelo4)
colnames(tabla.pronosticos)=c("Pron1","lim.inf","lim.sup","Pron2","lim.inf","lim.sup","Pron3","lim.inf","lim.sup","Pron4","lim.inf","lim.sup")
tabla.pronosticos
write.csv2(tabla.pronosticos,file="Pronosticos.csv",row.names = TRUE)

#tabla resumen de medidas de precision de pronosticos
precision.puntuales=rbind(accuracy(ytpron1,ytnuevo),accuracy(ytpron2,ytnuevo),accuracy(ytpron3,ytnuevo),accuracy(ytpron4,ytnuevo))[,c(2,3,5)]
precision.intervalos=rbind(amplcobmodelo1,amplcobmodelo2,amplcobmodelo3,amplcobmodelo4)
tabla.precision=cbind(precision.puntuales,precision.intervalos)
rownames(tabla.precision)=c("Modelo1","Modelo2","Modelo3","Modelo4")
tabla.precision
write.csv2(tabla.precision,file="Pronosticosoprecision.csv",row.names = TRUE)



#------------ PUNTO 7: ARMA VS MEJOR MODELO LOCAL-------------------------------------------

#El mejor modelo local del trabajo 1 fue el modelo 3
modelocal=SuavizamientoEstacional(yt,seasonal="additive",h=m)
str(modelocal) 

#Grafica de su ajuste
win.graph()
plot(Datosx, ylab="Datosx")
lines(fitted(modelocal),col=2,lwd=2)
legend("topleft",legend=c("Original","SEHW"),lty=1,col=c(1,2))

#Calculo de AIC y BIC 
s=m
parmodelocal=(s-1)+2
Criterioslocal=exp.crit.inf.resid(residuales=residuals(modelocal),n.par=parmodelocal);Criterioslocal
MSE=modelocal$MSE 

#Graficos de sus residuales
win.graph()
plot(residuals(modelocal),ylim=c(min(-2*sqrt(MSE),residuals(modelocal)),max(2*sqrt(MSE),residuals(modelocal))),ylab="Residuales")
abline(h=c(-2*sqrt(MSE),0,2*sqrt(MSE)),lty=1,col=2)

#Residuos vs. ajustados
win.graph()
plot(as.numeric(fitted(modelocal)),residuals(modelocal),ylim=c(min(-2*sqrt(MSE),residuals(modelocal)),max(2*sqrt(MSE),residuals(modelocal))),ylab="Residuales")
abline(h=c(-2*sqrt(MSE),0,2*sqrt(MSE)),lty=1,col=2)


#Pronosticos del modelo local
pronoslocal=modelocal$forecast
pronoslocal
ytpronlocal=pronoslocal[,1] 

#precision pronosticos puntuales modelo local
accuracy(ytpronlocal,ytnuevo) 

#Precision pronosticos por I.P Modelo local
amplcobmodelocal=amplitud.cobertura(real=ytnuevo,LIP=pronoslocal[,2],LSP=pronoslocal[,3]);amplcobmodelocal

#Evaluacion del supuesto de ruido blanco


#ACF 
win.graph()
acf(as.numeric(residuals(modelocal)),ci.type="ma",lag.max=36,main="ACF Modelo local")

#PACF
win.graph()
pacf(as.numeric(residuals(modelocal)),lag.max=36,main="PACF Modelo local")

BP.LB.test(residuals(modelocal),maxlag=36,type="Ljung")


#RESUMEN PROGRAMACION-----------------------------------------------------------

#tabla resumen de medidas de ajuste AIC Y BIC 
tabla.ajuste=rbind(Criteriosmodelo1,Criteriosmodelo2,Criteriosmodelo3,Criteriosmodelo4)
rownames(tabla.ajuste)=c("Modelo1","Modelo2","Modelo3","Modelo4")
tabla.ajuste

#tabla resumen de medidas de precision de pronosticos
precision.puntuales=rbind(accuracy(ytpron1,ytnuevo),accuracy(ytpron2,ytnuevo),accuracy(ytpron3,ytnuevo),accuracy(ytpron4,ytnuevo))[,c(2,3,5)]
precision.intervalos=rbind(amplcobmodelo1,amplcobmodelo2,amplcobmodelo3,amplcobmodelo4)
tabla.precision=cbind(precision.puntuales,precision.intervalos)
rownames(tabla.precision)=c("Modelo1","Modelo2","Modelo3","Modelo4")
tabla.precision
