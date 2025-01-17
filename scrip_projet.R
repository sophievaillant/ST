library(zoo)
library(tseries)
#install.packages("fUnitRoots")
library(fUnitRoots)
library(forecast)

#on charge les donnees
path="C:/Users/Sophie/Desktop/Projet"
setwd(path)
getwd()
list.files()
dta=read.csv("valeurs_mensuelles.csv", sep = ";")
str(dta)

#mise en forme
colnames(dta)[1] <- "date"
colnames(dta)[2] <- "IBPI"
colnames(dta)[3] <- "classe"
dta = dta[-(1:2),]
dta <- dta[order(dta$date),]

## PARTIE 1

#on met la serie en zoo
dta$IBPI=zoo(dta$IBPI,order.by=dta$date)
str(dta$IBPI)

#on la represente
plot(dta$IBPI, main = "Indice brut de production industrielle de la viande de volaille en fonction du temps", xlab = "Dates", ylab = "IBPI")
#elle n'est visiblement pas stationnaire du coup:

#differenciation de la serie
a = as.numeric(dta$IBPI)
dta$DIBPI=zoo(c(NA,diff(a,1)),order.by=dta$date)
str(dta$DIBPI)

#representation graphique de la serie differenciee
plot(dta$DIBPI[2:nrow(dta)],type="l",xlab="Dates",ylab="Diff IBPI",xaxt='n')
#la serie semble maintenant etre stationnaire (centree autour de 0)

T2=length(dta$IBPI)
xm=dta$IBPI[1:nrow(dta)]
str(xm)

#on enleve la tendance si elle existe pour comparer les resultats

#pour retirer la tendance de xm on la regresse sur une tendance lineaire
#et une eventuelle constante

trend = 1:length(xm)
lm = lm(xm ~ trend) 
summary(lm) 
r <- lm$residuals #on garde la serie corrigee de sa tendance
plot(r)
acf(r)
# pas de signe visuel de saisonnalite ?

X_notrend <- xm - lm$coefficients[1] - lm$coefficients[2]*obs 
serie_notrend <- zoo(X_notrend, order.by=dta$date)
plot(serie_notrend,main="Repr�sentation de la s�rie corrig�e de sa tendance")

#On compare la s�rie corrig�e de sa tendance et la s�rie diff�renci�e pour s�lectionner le melleur r�sultat de stationnarit�e
split.screen(1:2)
screen(2) ; plot(serie_notrend)
screen(1) ; plot(dta$DIBPI[2:nrow(dta)],type="l",xlab="Dates",ylab="Diff IBPI",xaxt='n')
#Enlever la tendance ne semble pas stationnariser la s�rie autant que la diff�rentiation le fait 


##Transformation logarithmique de la s�rie et diff�renciation de la s�rie en log ####
#On regarde si le passage au log permet d'avoir des r�sultats plus stationnaires 
a = as.numeric(log(dta$IBPI))
dlogIBPI <- zoo(diff(a,1))
split.screen(1:2)
screen(1) ; plot(a) #pour grandeurs comparables
screen(2) ; plot(dta$dlogIBPI)
#Le passage au log semble utile ici; comparons le log et la diff�renciation

split.screen(1:2)
screen(2) ; plot(dlogIBPI)
screen(1) ; plot(dta$DIBPI[2:nrow(dta)],type="l",xlab="Dates",ylab="Diff IBPI",xaxt='n')
#La serie differenciee en log semble etre plus stationnaire que la differenciee

dta$dlogIBPI <- zoo(c(NA,dlogIBPI), order.by=dta$date)

## Verificiation de la stationnarit� ####

##a. Test de Dickey Fuller ####
#On v�rifie d'abord que la s�rie en niveau n'est pas stationnaire.

adf <- adfTest(dta$IBPI, lag=0, type="ct") #test ADF avec constante et tendance
adf

#il faut d'abord v�rifier la validit� du test, c�d l'absence d'autocorr�lation des r�sidus
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

#on regarde pour la s�rie lm (pas stationnaris�e)
Qtests(adf@test$lm$residuals[0:24], 24, fitdf = length(adf@test$lm$coefficients))
# Toutes les p-value sont superieurs a 10%, on ne peut pas rejeter
# l'hypothese nulle de non correlation des residus
# Le test de Dickey Fuller avec lags = 0 est valide

adf
# On obtient une p-value de 0.01
# On rejette l'hypothese nulle d'existence d'une racine unitaire au seuil de 99%
# donc on rejette l'hypoth�se selon laquelle la s�rie ne serait pas stationnaire

######### revoir ##########
class(dta$date)
#testons donc la racine unitaire pour la s�rie log diff�renci�e

Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
adf
#la racine unitaire est rejet�e au seuil de 95%. La s�rie diff�renci�e est donc I(0).
#locom est donc I(1).


# On realise d'autres tests pour verifier la stationnarite de la serie

# Test Philippe - Peron (marche pas)
pp.test(dta$dlogIBPI)
# On observe une p-value inferieur à 1%
# On rejete à nouveau l'hypothèse nulle d'existence d'une racine unitaire au seuil de 99%

# Test KPSS
kpss.test(dta$dlogIBPI,null = 'Level')
# On obtient une p-value de plus de 10% ("p-value greater than printed p-value"), 
#on ne peut rejeter l'hypoth�se de stationnarit�

# Tout les tests semblent indiquer que la serie diff�renci�e � l'ordre 1 est stationnaire


#3.Repr�sentation de la s�rie choisie avant et apr�s transformation ####
#La comparaison avant/apr�s est la suivante : 

#je copie une comparaison qui est au-dessus
a = as.numeric(log(dta$IBPI))
dlogIBPI <- zoo(diff(a,1))
split.screen(1:2)
screen(1) ; plot(log(dta$IBPI)) #pour grandeurs comparables
screen(2) ; plot(dta$dlogIBPI)

## PARTIE 2

#1.Choix d'un mod�le ARMA
##a.Identification des param�tres ####
#On regarde les ordres max des AR et MA

ggAcf(dlogIBPI) #q* = 25 wesh
ggPacf(dlogIBPI) #p* = aussi 25
