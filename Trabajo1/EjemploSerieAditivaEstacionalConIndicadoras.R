#LIBRERIAS Y PAQUETES ----------------------------------------------------------
rm(list=ls(all=TRUE))
library(TSA)
library(forecast)
library(fANCOVA)

source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Mipoly.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Mytrigon.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-Criterios.Informacion-Calidad.Intervalos.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencialv02.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-SuavizamientoEstacional.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-Descomp.Loessv02.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-interpdeltas.R")

#Lectura de Datos --------------------------------------------------------------------------------------------------------------------
#Reemplazar este fragmento por el codigo que deja la profesora en el archivo PROGRAMA-R-LECTURA-DATOS ASIGNADOS-022024.R

Datos21=read.table(file.choose(),header=T,sep=";",skip=14,dec=",",colClasses=c(rep("NULL",6),"numeric",rep("NULL",4)))
Datos21=ts(Datos21,freq=12,start=c(2001,1))
Datos21
plot(Datos21)

#--ANALISIS DESCRIPTIVO --------------------------------------------------------

#Graficando la serie 
win.graph()
plot(Datos21, ylab="Datos21",xlab="Time")

#Graficando la tendencia 
win.graph()
Tt=decompose(Datos21)$trend
plot(Tt,ylim=c(min(Datos21),max(Datos21)))

#Graficando boxplot 
win.graph()
boxplot(Datos21~cycle(Datos21),names=month.abb)

#periodograma 
win.graph()
periodogram(diff(Datos21),lwd=3)
abline(h=0)
abline(v=c(1:6)/12,col=2,lty=2)

#------------------------------ DEFINIENDO VARIABLES Y CREACION DE LA MATRIZ NO MOVER -------------------------
m=12
n=length(Datos21)-m

#PARA EL AJUSTE
t=1:n
yt=ts(Datos21[t],freq=m,start=c(2001,1))
poli=Mipoly(tiempo=t,grado=4)
mes=seasonaldummy(yt)

#Matriz de diseño en el ajuste
X1=data.frame(t,mes)  #Grado 1 - Modelo 1 Global
head(X1) 

X2=data.frame(poli,mes) #Grado 4 - Modelo 2 Global
head(X2) 

#PARA LOS PRONOSTICOS

tnuevo=(n+1):length(Datos21) 
polinuevo=Mipoly(tiempo=tnuevo,grado=4)
ytnuevo=ts(Datos21[tnuevo],freq=m,start=c(2023,8)) 
mesnuevo=seasonaldummy(yt,h=m)

#Matriz de diseño en los pronosticos
X1nuevo=data.frame(t=tnuevo,mesnuevo)
head(X1nuevo) 
X2nuevo=data.frame(polinuevo,mesnuevo)
head(X2nuevo) 

#-----MODELO 1:POLINOMIAL-GRADO1 LINEAL------------------------------------------------------

mod1=lm((yt)~.,data=X1)
summary(mod1)
#Calculo valores ajustados del modelo 1
ythatmod1=ts(fitted(mod1),freq=m,start=start(yt))

#-----MODELO 2:POLINOMIAL-GRADO 4------------------------------------------------------

mod2=lm((yt)~.,data=X2)
summary(mod2)
#Calculo valores ajustados del modelo 2
ythatmod2=ts(fitted(mod2),freq=m,start=start(yt))

#-----MODELO 3:SEHW ------------------------------------------------------------

mod3=SuavizamientoEstacional(yt,seasonal="additive",h=m,beta=1e-5)
str(mod3) 

#-----MODELO 4:DLL(GCV) -------------------------------------------------------

mod4=Descomp.Loessv02(serie.ajuste=yt,h=m,tipo.descomp="additive",grado=2,criterio="aicc")
str(mod4) 

#---------------------------- CRITERIOS DE AJUSTE ------------------------------------------#

#Calculo de los criterios AIC y BIC en modelo 1
nparmod1=length(coef(mod1)[coef(mod1)!=0]);nparmod1 
Criterios1= exp.crit.inf.resid(residuales=residuals(mod1),n.par=nparmod1);Criterios1
#Calculo de los criterios AIC y BIC en modelo 2
nparmod2=length(coef(mod2)[coef(mod2)!=0]);nparmod2 
Criterios2= exp.crit.inf.resid(residuales=residuals(mod2),n.par=nparmod2);Criterios2
#Calculo de los criterios AIC y BIC en modelo 3
p3=(m-1)+2
Criterios3=exp.crit.inf.resid(residuales=residuals(mod3),n.par=p3);Criterios3
MSE3=mod3$MSE 
#Calculo AIC y BIC 
Criterios4=exp.crit.inf.resid(residuales=residuals(mod4),n.par=mod4$p);Criterios4
MSE4=mod4$MSE 

#---------------------------- PRONOSTICOS ------------------------------------------#

#Pronosticos del modelo 1 
pronos1=predict(mod1,newdata=X1nuevo,interval="prediction",level=0.95)
pronos1=ts(pronos1,freq=m,start=start(ytnuevo))
pronos1
ytpron1=pronos1[,1]
#Pronosticos del modelo 2 
pronos2=predict(mod2,newdata=X2nuevo,interval="prediction",level=0.95)
pronos2=ts(pronos2,freq=m,start=start(ytnuevo))
pronos2
ytpron2=pronos2[,1]
#Pronosticos del modelo 3 
pronos3=mod3$forecast
pronos3
ytpron3=pronos3[,1] #solo los pronosticos puntuales del suavizamiento
#Pronosticos del modelo 4 
pronos4=mod4$forecast
pronos4
ytpron4=pronos4[,1]#Precision pronosticos puntuales

#---------------------------- METRICAS DE PRONOSTICOS ------------------------------------------#

#precision pronosticos puntuales modelo 1
accuracy(ytpron1,ytnuevo)
amplcobmod1=amplitud.cobertura(real=ytnuevo,LIP=pronos1[,2],LSP=pronos1[,3]);amplcobmod1#precision pronosticos por I.P modelo 1
#precision pronosticos puntuales modelo 2
accuracy(ytpron2,ytnuevo)
amplcobmod2=amplitud.cobertura(real=ytnuevo,LIP=pronos2[,2],LSP=pronos2[,3]);amplcobmod2#precision pronosticos por I.P modelo 1
#precision pronosticos puntuales modelo 3
accuracy(ytpron3,ytnuevo) 
amplcobmod3=amplitud.cobertura(real=ytnuevo,LIP=pronos3[,2],LSP=pronos3[,3]);amplcobmod3#Precision pronosticos por I.P Modelo 3
#precision pronosticos puntuales modelo 4
accuracy(ytpron4,ytnuevo)
amplcobmod4=amplitud.cobertura(real=ytnuevo,LIP=pronos4[,2],LSP=pronos4[,3]);amplcobmod4

#--TABLAS DE LAS MEDIDAS DE AJUSTE Y PRONOSTICOS DE LOS 4 MODELOS --------------

Modelo1=summary(mod1)$coefficients
Modelo2=summary(mod2)$coefficients

tabla.parametros.globales=rbind(Modelo1,Modelo2)
rownames(tabla.parametros.globales)=c("Modelo 1","Modelo 2")
tabla.parametros.globales

#Tabulando medidas de ajuste
tabla1.criterios=rbind(Criterios1,Criterios2,Criterios3,Criterios4)
rownames(tabla1.criterios)=c("Modelo 1","Modelo 2","Modelo 3","Modelo 4")
tabla1.criterios

#tabla resumen de pronosticos 
tabla.pronosticos=cbind(pronos1,pronos2,pronos3,pronos4)
colnames(tabla.pronosticos)=c("Modelo 1","Modelo 2","Modelo 3","Modelo4")
tabla.pronosticos

#Tabulando medidas de pronosticos
precision.puntuales=rbind(accuracy(ytpron1,ytnuevo), accuracy(ytpron2,ytnuevo), accuracy(ytpron3,ytnuevo), accuracy(ytpron4,ytnuevo))[,c(2,3,5)]
precision.intervalos=rbind(amplcobmod1,amplcobmod2,amplcobmod3,c(NA,NA))
tabla2.precision=cbind(precision.puntuales,precision.intervalos)
rownames(tabla2.precision)=c("Modelo 1","Modelo 2","Modelo 3","Modelo 4")
tabla2.precision


#--GRAFICAS VALORES AJUSTADOS PARA LOS 4 MODELOS -------------------------------

#MODELO 1
win.graph()
plot(Datos21, ylab="Datos21")
lines(ythatmod1,col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 1"),lty=1,col=c(1,2))

#MODELO 2
win.graph()
plot(Datos21, ylab="Datos21")
lines(ythatmod2,col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 2"),lty=1,col=c(1,2)) 

#MODELO 3
win.graph()
plot(Datos21, ylab="Datos21")
lines(fitted(mod3),col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 3"),lty=1,col=c(1,2))

#MODELO 4
win.graph()
plot(Datos21, ylab="Datos21")
lines(fitted(mod4),col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 4"),lty=1,col=c(1,2))

#**Graficando la serie desestacionalizada y su ajuste loess para la estacionalidad
win.graph()
plot(mod4$ytd,ylab=NA)
lines(mod4$Tt,col=2,lwd=2)
legend("topleft",legend=c("Serie ajustada estacionalmente","Tendencia LOESS"),col=c(1,2),lty=1)

#Graficando St estimada por el filtro de descomposicion
win.graph()
plot(mod4$St,ylab=expression(hat(S)[t])) 

#------EFECTOS ESTACIONALES ESTIMADOS GLOBALES--------------------------------------

deltas1=interpdeltas(mod1,gradopoly=1,aditivo=TRUE,plotit=FALSE);deltas1$deltasi
deltas2=interpdeltas(mod2,gradopoly=4,aditivo=TRUE,plotit=FALSE);deltas1$deltasi

#Grafico de los efectos estacionales estimados, en un mismo plano
win.graph()
interpdeltas(mod1,gradopoly=1,aditivo=TRUE,plotit=TRUE)
lines(deltas1$periodo,deltas1$deltasi,type="b",lty=1,pch=1,col=1,lwd=3)
lines(deltas2$periodo,deltas2$deltasi,type="b",lty=2,pch=2,col=2,lwd=3)
legend("bottomright",legend=c("Modelo 1","Modelo 2"),col=1:2,pch=1:2,lty=1:2,lwd=3)

#--GRAFICAS RESIDUOS DE AJUSTE PARA LOS 4 MODELOS ------------------------------

#MODELO 1
win.graph()
plot.ts(residuals(mod1),ylab="Residuales")
abline(h=c(-2*summary(mod1)$sigma,0,2*summary(mod1)$sigma),lty=1,col=2)
legend("topleft",legend="Modelo 1")
#Residuos vs. ajustados
win.graph()
plot(fitted(mod1),residuals(mod1),ylab="Residuos vs. ajustados")
abline(h=c(-2*summary(mod1)$sigma,0,2*summary(mod1)$sigma),lty=1,col=2)
legend("topleft",legend="Modelo 1")

#MODELO 2
win.graph()
plot.ts(residuals(mod2),ylab="Residuales")
abline(h=c(-2*summary(mod2)$sigma,0,2*summary(mod2)$sigma),lty=1,col=2)
legend("topleft",legend="Modelo 2")
#Residuos vs. ajustados
win.graph()
plot(fitted(mod2),residuals(mod2),ylab="Residuos vs. ajustados")
abline(h=c(-2*summary(mod2)$sigma,0,2*summary(mod2)$sigma),lty=1,col=2)
legend("topleft",legend="Modelo 2")

#MODELO 3
win.graph()
plot(residuals(mod3),ylim=c(min(-2*sqrt(MSE3),residuals(mod3)),max(2*sqrt(MSE3),residuals(mod3))),ylab="Residuales")
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),lty=1,col=2)
legend("topleft",legend="Modelo 3")
#Residuos vs. ajustados
win.graph()
plot(as.numeric(fitted(mod3)),residuals(mod3),ylim=c(min(-2*sqrt(MSE3),residuals(mod3)),max(2*sqrt(MSE3),residuals(mod3))),ylab="Residuos vs. ajustados")
abline(h=c(-2*sqrt(MSE3),0,2*sqrt(MSE3)),lty=1,col=2)
legend("topleft",legend="Modelo 3")

#MODELO 4 
win.graph()
plot(residuals(mod4),ylim=c(min(-2*sqrt(MSE4),residuals(mod4)),max(2*sqrt(MSE4),residuals(mod4))),ylab="Residuales")
abline(h=c(-2*sqrt(MSE4),0,2*sqrt(MSE4)),lty=1,col=2)
legend("topleft",legend="Modelo 4")
#Residuos vs. ajustados
win.graph()
plot(as.numeric(fitted(mod4)),residuals(mod4),ylim=c(min(-2*sqrt(MSE4),residuals(mod4)),max(2*sqrt(MSE4),residuals(mod4))),ylab="Residuos vs. ajustados")
abline(h=c(-2*sqrt(MSE4),0,2*sqrt(MSE4)),lty=1,col=2)
legend("topleft",legend="Modelo 4")

#--GRAFICA COMPARATIVA DE LOS PRONOSTICOS DE LOS 4 MODELOS ---------------------

win.graph()
plot(ytnuevo,type="b",ylab="Datos21",col=1,pch=19,ylim=c(min(ytnuevo,ytpron1,ytpron2,ytpron3,ytpron4),max(ytnuevo,ytpron1,ytpron2,ytpron3,ytpron4)),lwd=2,xaxt="n")
axis(1,at=time(ytnuevo),labels=c("ago-23","sep-23","oct-23","nov-23","dic-23","ene-24","feb-24","mar-24","abr-24","may-24","jun-24","jul-24"),cex.axis=0.7)
lines(ytpron1,col=2,pch=1,type="b",lwd=2)
lines(ytpron2,col=3,pch=2,type="b",lwd=2)
lines(ytpron3,col=4,pch=3,type="b",lwd=2)
lines(ytpron4,col=5,pch=4,type="b",lwd=2)
legend("topleft",legend=c("Real","Modelo 1","Modelo 2","Modelo 3","Modelo 4"),pch=c(19,1:4),col=c(1:5),lwd=2)

#Exportacion tablas al directorio de trabajo------------------------------------
# SI QUIEREN  QUE GUARDE EL ARCHIVO EN UNA CARPETA EN ESPECIFICO DEBEN CAMBIAR LA RUTA DE file
write.csv2(tabla.parametros.globales,file="tablamod1ymod2trabajo1.csv",row.names = TRUE)
write.csv2(mod3$coefficients,file="tablamod3trabajo1.csv",row.names = TRUE)
write.csv2(mod4$deltasi,file="tablamod4trabajo1.csv",row.names = TRUE)
write.csv2(tabla1.criterios,file="tablacriteriostrabajo1.csv",row.names = TRUE)
write.csv2(tabla.pronosticos,file="tablapronosticostrabajo1.csv",row.names = TRUE)
write.csv2(pronos4,file="tablapronosticosmodelo4trabajo1.csv",row.names = TRUE)
write.csv2(tabla2.precision,file="tablamedidasprecisiondepronosticostrabajo1.csv",row.names = TRUE)
