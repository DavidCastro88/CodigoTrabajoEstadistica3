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
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-predict_expo.R")
#Lectura de Datos --------------------------------------------------------------------------------------------------------------------

#Reemplazar este fragmento por el codigo que deja la profesora en el archivo PROGRAMA-R-LECTURA-DATOS ASIGNADOS-022024.R asociado a sus datos

#Leer Datos
DatosX=read.table(file.choose(),header=T,sep=";",skip=14,dec=",",colClasses=c(rep("NULL",6),"numeric",rep("NULL",4)))
DatosX=ts(DatosX,freq=12,start=c(2001,1))
plot(DatosX)

#--ANALISIS DESCRIPTIVO ------------------------------------------------------------------------------------------------------------------

#Graficando la serie 
win.graph()
plot(DatosX, ylab="DatosX",xlab="Time")

#Graficando la serie en escala logaritmica 
win.graph()
plot(log(DatosX), ylab="DatosX")

#Graficando la tendencia 
win.graph()
Tt.log=decompose(log(DatosX))$trend
plot(Tt.log,ylim=c(min(log(DatosX)),max(log(DatosX))))

#Graficando boxplot 
win.graph()
boxplot(log(DatosX)~cycle(DatosX))

#periodograma 
win.graph()
x=diff(log(DatosX))  
plot(x,ylab=expression(log(Y[t])-log(Y[t-1])));abline(h=mean(x))
periodogram(x);
abline(v=c(1/12,2/12,3/12,4/12,5/12,6/12),col=2,lty=2)


#DEFINIENDO VARIABLES Y CREACION DE LA MATRIZ NO MOVER -----------------------------------------------------------------------------------
m=12 
n=length(DatosX)-m

#PARA EL AJUSTE

t=1:n
yt=ts(DatosX[t],freq=m,start=c(2001,1))
poli=Mipoly(tiempo=t,grado=3) # Grado del polinomio
mes=seasonaldummy(yt)

#Matriz de diseño en el ajuste
X1=data.frame(poli,mes)
head(X1) 

#PARA LOS PRONOSTICOS
tnuevo=(n+1):length(DatosX)
ytnuevo=ts(DatosX[tnuevo],freq=m,start=c(2023,8))
polinuevo=Mipoly(t=tnuevo,grado=3)
mesnuevo=seasonaldummy(yt,h=m)

#Matriz de diseño en los pronosticos
X1nuevo=data.frame(polinuevo,mesnuevo)
head(X1nuevo) 

#-----MODELO 1:COMPLETAMENTE MULTIPLICATIVO----------------------------------------------------------------------------------------------------------------

mod1=lm(log(yt)~.,data=X1)
summary(mod1)
#Calculo valores ajustados del modelo 1
ythatmod1=ts(exp(fitted(mod1))*exp(summary(mod1)$sigma^2/2),freq=m,start=start(yt))

#-----MODELO 2:PARCIALMENTE MULTIPLICATIVO------------------------------------------------------

mod2=regexponencialv02(respuesta=yt,data=X1)
summary(mod2)

#Calculo valores ajustados del modelo 2
ythatmod2=ts(fitted(mod2),freq=m,start=start(yt))

#-----MODELO 3:SEHW ------------------------------------------------------------

mod3=SuavizamientoEstacional(yt,seasonal="multiplicative",h=m)
str(mod3) 

#-----MODELO 4:DLC(aicc) -------------------------------------------------------

mod4=Descomp.Loessv02(serie.ajuste=yt,h=m,tipo.descomp="multiplicative",grado=1,criterio="aicc")
str(mod4) 

#---------------------------- CRITERIOS DE AJUSTE ------------------------------------------#

#Calculo de AIC y BIC en modelo 1
nparmod1=length(coef(mod1)[coef(mod1)!=0]);nparmod1 
res.orig.mod1=yt-ythatmod1
Criterios1= exp.crit.inf.resid(residuales=res.orig.mod1,n.par=nparmod1);Criterios1
#Calculo de AIC y BIC en modelo 2
nparmod2=length(coef(mod2)[coef(mod2)!=0]);nparmod2          
Criterios2=exp.crit.inf.resid(residuales=residuals(mod2),n.par=nparmod2); Criterios2
#Calculo de AIC y BIC en modelo 3
s=m
p3=(s-1)+2
Criterios3=exp.crit.inf.resid(residuales=residuals(mod3),n.par=p3);Criterios3
MSE3=mod3$MSE 
#Calculo AIC y BIC en modelo 4
Criterios4=exp.crit.inf.resid(residuales=residuals(mod4),n.par=mod4$p);Criterios4
MSE4=mod4$MSE 

#---------------------------- PRONOSTICOS ------------------------------------------#

#Pronosticos del modelo 1 en la escala original
pronos1=exp(predict(mod1,newdata=X1nuevo,interval="prediction",level=0.95))*exp(summary(mod1)$sigma^2/2)
pronos1=ts(pronos1,freq=m,start=start(ytnuevo))
pronos1
ytpron1=pronos1[,1]
#Pronosticos del modelo 2 en la escala original, solo son de tipo puntual por ser modelo no lineal
pronos2=predict_expo(mod2,new.data=X1nuevo,interval="prediction",level=0.95)
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
ScoreIP1=IntervalScore(real=ytnuevo,LIP=pronos1[,2],LSP=pronos1[,3],alpha=0.05);ScoreIP1
#precision pronosticos puntuales modelo 2
accuracy(ytpron2,ytnuevo)
amplcobmod2=amplitud.cobertura(real=ytnuevo,LIP=pronos2[,2],LSP=pronos2[,3]);amplcobmod2#precision pronosticos por I.P modelo 1
ScoreIP2=IntervalScore(real=ytnuevo,LIP=pronos2[,2],LSP=pronos2[,3],alpha=0.05);ScoreIP2
#precision pronosticos puntuales modelo 3
accuracy(ytpron3,ytnuevo) 
amplcobmod3=amplitud.cobertura(real=ytnuevo,LIP=pronos3[,2],LSP=pronos3[,3]);amplcobmod3#Precision pronosticos por I.P Modelo 3
ScoreIP3=IntervalScore(real=ytnuevo,LIP=pronos3[,2],LSP=pronos3[,3],alpha=0.05);ScoreIP3
#precision pronosticos puntuales modelo 4
accuracy(ytpron4,ytnuevo)
amplcobmod4=amplitud.cobertura(real=ytnuevo,LIP=pronos4[,2],LSP=pronos4[,3]);amplcobmod4
ScoreIP4=IntervalScore(real=ytnuevo,LIP=pronos4[,2],LSP=pronos4[,3],alpha=0.05);ScoreIP4


#--------------------------------TABLAS DE LAS MEDIDAS DE AJUSTE Y PRONOSTICOS DE LOS 4 MODELOS --------------

Modelo1=summary(mod1)$coefficients
Modelo2=summary(mod2)$coefficients

tabla.parametros.globales=cbind(Modelo1,Modelo2)
tabla.parametros.globales

#Tabulando medidas de ajuste
tabla1.criterios=rbind(Criterios1,Criterios2,Criterios3,Criterios4)
rownames(tabla1.criterios)=c("Modelo 1","Modelo 2","Modelo 3","Modelo 4")
tabla1.criterios

#tabla resumen de pronosticos 
tabla.pronosticos=cbind(pronos1,pronos2,pronos3,pronos4)
tabla.pronosticos

#Tabulando medidas de pronosticos
precision.puntuales=rbind(accuracy(ytpron1,ytnuevo), accuracy(ytpron2,ytnuevo), accuracy(ytpron3,ytnuevo), accuracy(ytpron4,ytnuevo))[,c(2,3,5)]
precision.intervalos=rbind(amplcobmod1,amplcobmod2,amplcobmod3,amplcobmod4)
tabla2.precision=cbind(precision.puntuales,precision.intervalos)
rownames(tabla2.precision)=c("Modelo 1","Modelo 2","Modelo 3","Modelo 4")
tabla2.precision

#--------------------------------GRAFICAS VALORES AJUSTADOS PARA LOS 4 MODELOS -------------------------------

#MODELO 1
win.graph()
plot(DatosX,ylab="DatosX")
lines(ythatmod1,col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 1"),lty=1,col=c(1,2))

#MODELO 2
win.graph()
plot(DatosX,ylab="DatosX")
lines(ythatmod2,col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 2"),lty=1,col=c(1,2)) 

#MODELO 3
win.graph()
plot(DatosX,ylab="DatosX")
lines(fitted(mod3),col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 3"),lty=1,col=c(1,2))

#MODELO 4**
win.graph()
plot(DatosX,ylab="DatosX")
lines(fitted(mod4),col=2,lwd=2)
legend("topleft",legend=c("Original","Modelo 4"),lty=1,col=c(1,2))

#------------------------------ Graficando la serie desestacionalizada y su ajuste loess para la estacionalidad
win.graph()
plot(mod4$ytd,ylab=NA,main="Gráfica del ajuste loess sobre la serie desestacionalizada")
lines(mod4$Tt,col=2,lwd=2)
legend("topleft",legend=c("Serie ajustada estacionalmente","Tendencia LOESS"),col=c(1,2),lty=1)

#Graficando St estimada por el filtro de descomposicion
win.graph()
plot(mod4$St,ylab=expression(hat(S)[t])) 

#------EFECTOS ESTACIONALES ESTIMADOS GLOBALES--------------------------------------
expdeltasi1=interpdeltas(modelo=mod1,gradopoly=3,aditivo=FALSE,plotit=FALSE);expdeltasi1$expdeltasi100
expdeltasi2=interpdeltas(modelo=mod2,gradopoly=3,aditivo=FALSE,plotit=FALSE);expdeltasi2$expdeltasi100

#Grafico de los efectos estacionales estimados, en un mismo plano
win.graph()
interpdeltas(modelo=mod1,gradopoly=1,aditivo=FALSE,plotit=TRUE)
lines(expdeltasi1$periodo,expdeltasi1$expdeltasi100,type="b",lty=1,pch=1,col=1,lwd=3)
lines(expdeltasi2$periodo,expdeltasi2$expdeltasi100,type="b",lty=2,pch=2,col=2,lwd=3)
legend("topleft",legend=c("Modelo 1","Modelo 2"),col=1:2,pch=1:2,lty=1:2,lwd=3)

#------------------------------ GRAFICAS RESIDUOS DE AJUSTE PARA LOS 4 MODELOS ------------------------------

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

#--------------------------------GRAFICA COMPARATIVA DE LOS PRONOSTICOS DE LOS 4 MODELOS ---------------------

win.graph()
plot(ytnuevo,type="b",ylab="DatosX",col=1,pch=19,ylim=c(min(ytnuevo,ytpron1,ytpron2,ytpron3,ytpron4),max(ytnuevo,ytpron1,ytpron2,ytpron3,ytpron4)),lwd=2,xaxt="n")
axis(1,at=time(ytnuevo),labels=c("ago-23","sep-23","oct-23","nov-23","dic-23","ene-24","feb-24","mar-24","abr-24","may-24","jun-24","jul-24"),cex.axis=0.7)
lines(ytpron1,col=2,pch=1,type="b",lwd=2)
lines(ytpron2,col=3,pch=2,type="b",lwd=2)
lines(ytpron3,col=4,pch=3,type="b",lwd=2)
lines(ytpron4,col=5,pch=4,type="b",lwd=2)
legend("topleft",legend=c("Real","Modelo 1","Modelo 2","Modelo 3","Modelo 4"),pch=c(19,1:4),col=c(1:5),lwd=2)

#------------------------------Exportacion tablas al directorio de trabajo------------------------------------
# SI QUIEREN  QUE GUARDE EL ARCHIVO EN UNA CARPETA EN ESPECIFICO DEBEN CAMBIAR LA RUTA DE file
write.csv2(tabla.parametros.globales,file="tablamod1ymod2trabajo1.csv",row.names = TRUE) 
write.csv2(mod3$coefficients,file="tablamod3trabajo1.csv",row.names = TRUE)
write.csv2(mod4$deltasi,file="tablamod4trabajo1.csv",row.names = TRUE)
write.csv2(tabla1.criterios,file="tablacriteriostrabajo1.csv",row.names = TRUE)
write.csv2(tabla.pronosticos,file="tablapronosticostrabajo1.csv",row.names = TRUE)
write.csv2(pronos4,file="tablapronosticosmodelo4trabajo1.csv",row.names = TRUE)
write.csv2(tabla2.precision,file="tablamedidasprecisiondepronosticostrabajo1.csv",row.names = TRUE)
