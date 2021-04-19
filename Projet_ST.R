#Packages nécessaire 
require(zoo)
require(tseries)
require(forecast)
require(dplyr)
require(fUnitRoots)
require(normtest)
require(xtable)
require(forecast)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############
## Partie I ##
##############

## Conversion du fichier .csv en dataframe ##

pates_df <- read.table(file = "valeurs_mensuelles.csv",sep = ";",
                       header = FALSE, skip = 3,
                      colClasses = c("character", "character", "NULL"),col.names = c("Date","IPI_pates","Dropped"))

pates_df$IPI_pates <- as.numeric(pates_df$IPI_pates) 
pates_df$Date<-as.Date(paste(pates_df$Date,'01',sep = '-'), format = "%Y-%m-%d")

pates_df <- pates_df[order(pates_df$Date),]
rownames(pates_df) <- 1:nrow(pates_df)

## Création de l'objet TimeSeries ##

ipi_pates_ts <- ts(pates_df$IPI_pates,start =c(1990,1),end = c(2020,3),frequency = 12)


## Visualisation de la série ##
plot(ipi_pates_ts)
# On observe une légère tendance à la baisse avant 2000 et une tendance à la hausse après 2000


## QUESTION 2 ##

## Autocorrélation ##
ggAcf(ipi_pates_ts)
# L'autocorrélation décroit vers 0, ce qui montre que la série n'est pas encore stationnaire

## Décomposition de la série ##
decompose_pates <- decompose(ipi_pates_ts, type = 'additive')
plot(decompose_pates)
#La série admet bien une tendance et n'est pas stationnaire
#La série est désaisonnalisé, l'amplitude et la période des variations saisonnières
#étant constante et très faible
#Le bruit semble homoscedastique

## Différencation de la série ##
# On va d'abord differencié la série à l'ordre 1

diff_pates_1 <- diff(ipi_pates_ts, lag = 1, differences = 1)
plot(diff_pates_1)
ggAcf(diff_pates_1)
# diff_pates_1 semble être stationnaire, l'autocorrélation varie autour de zéro sans être monotone

## Verification de la stationnarité ##

#On va verifier l'absence de constante et de tendance avec une regression
summary(lm(diff_pates_1 ~ pates_df$Date[-1]))
# La regression de la série différencié sur la date donne des p-valeurs très importante pour le coefficient de Date (p = 0.947)
# et pour le coefficient de la constante (p = 0.970)
# On en conclut que les coefficients de la regression ne sont pas significatifs


# On réalise un premier test de Dickey-Fuller augmenté (ADF) avec le package
# fUnitRoots comme dans le TD5 sans constante ni tendance
adf_test_pates <- adfTest(diff_pates_1, lags = 0, type = "nc")

# On utilise de plus la fonction Qtest du TD5 pour vérifier que les résidus du test
# ne sont pas corrélés
Qtests <- function (residuals, max_lag, fitdf = 0){
  lags_list <- seq(1,max_lag)
  pvals <- data.frame('lag' = numeric(), "p_value" = numeric())
  for (lag in lags_list){
    if (lag > fitdf){
      pvalue <- Box.test(residuals, lag = lag, type = "Ljung-Box", fitdf = fitdf)$p.value
      pvals <- rbind(pvals, data.frame("lag" = lag, "p_value" = pvalue))
    }
  }
  return(pvals)
}

Qtests(adf_test_pates@test$lm$residuals[0:24], 24)

# Toutes les p-value sont superieurs à 10%, on ne peut pas rejeter
# l'hypothèse nulle de non corrélation des résidus
# Le test de Dickey Fuller avec lags = 0 est valide

adf_test_pates@test$p.value

# On obtient une p-value de 0.01
# On rejete l'hypothèse nulle d'existence d'une racine unitaire au seuil de 99%

# On réalise d'autres tests pour verifier la stationnarité de la série

# Test PP
pp.test(diff_pates_1)
# On observe une p-value inferieur à 1%
# On rejete à nouveau l'hypothèse nulle d'existence d'une racine unitaire au seuil de 99%

# Test KPSS
kpss.test(diff_pates_1,null = 'Level')
# On obtient une p-value de plus de 10%, on ne peut rejeter l'hypothèse de stationnarité

# Tout les tests semble indiquer que la serie différencié à l'ordre 1 est stationnaire


## Question 3 ##

plot(ipi_pates_ts)
plot(diff_pates_1)


#################
### Partie II ###
#################



# On considère alors la série différencié à l'ordre 1 diff_pates_1

ggAcf(diff_pates_1)

ggPacf(diff_pates_1)

# On lit q* = 1 et p* = 6

# On doit donc tester les modèles ARMA(p,q) avec p<= 6 et q <= 1
# Il y a donc 7x2 = 14 couples de valeurs (p,q) à tester
# Pour qu'un modèle ARMA(p,q) soit retenu, il faut qu'il soit
# Bien ajusté (coefficients statistiquement significatifs)
# Valide (pas d'autocorrélation des résidus -> Test de Ljung-Box)
# On teste alors pour toutes les combinaisons possibles


# On crée d'abord des fonctions coeff_sign qui prend en argument
# un modèle ARIMA et renvoient resepctivement la p-valeur associé à chacund des coefficients
# On utilisera la fonction Qtests du TD5 pour l'autocorrélation des résidus

# On utilise la fonction ARIMA et non ARMA car ARMA ne fonctionne pas.

signif <- function(estim){
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(data.frame(coefficient = coef, variance = se,pvaleur = pval))
  }


modeles_test <- function(timeserie, p_etoile, q_etoile){
  
  resultat <- data.frame('ordre' = factor(),
                         'coefficients' = factor(),
                         'autocorr' = factor())
  
  for (p in 0:p_etoile){
    for (q in 0:q_etoile){
      modele_arima <- arima(timeserie,c(p,0,q), include.mean = FALSE) #Création du modèle
      signif_coeff <- signif(modele_arima) 
      signif_autocorr <- Qtests(modele_arima$residuals[0:25],24,fitdf = 6)
      ordre <- paste('ARMA(',p,',',q,')')
      if (sum(signif_coeff$pvaleur > 0.05) == 0){
        coefficients <- 'significatifs'
      }
      else{
        coefficients <- 'pas significatifs'
        }
      if (sum(signif_autocorr < 0.05) <= 1) #La condition == 0 est trop discriminante
        #et rejette tout les modèles ARMA sauf ARMA(1,1) donc pas interessant de comparer
        #les BIC et AIC après
        {
        autocorr <- 'Pas d\'autocorrélation'
      }
      else{
        autocorr <- 'autocorrélation'
      }
      resultat<- rbind(resultat, data.frame(ordre = ordre, 
                                            coefficients = coefficients,autocorr = autocorr))
    }
  }
  return(resultat)
  }
  
p_etoile <- 6
q_etoile <- 1

df <- modeles_test(diff_pates_1, p_etoile = p_etoile, q_etoile = q_etoile)


# On compe alors 6 modèles bien ajustés et valide
# ARMA(0,1) ,ARMA(1,1) ,ARMA(2,0) ,ARMA(3,0), ARMA(6,0)
# On va utiliser les critères AIC et BIC pour les départager

ar2 <- arima(diff_pates_1, c(2,0,0), include.mean = FALSE)
ar3 <- arima(diff_pates_1, c(3,0,0), include.mean = FALSE)
ar6 <- arima(diff_pates_1, c(6,0,0), include.mean = FALSE)
ma1 <- arima(diff_pates_1, c(0,0,1), include.mean = FALSE)
ar1ma1 <- arima(diff_pates_1, c(1,0,1), include.mean = FALSE)

models <- c("ar2","ar3","ar6","ma1","ar1ma1"); names(models) <- models
criterias <- apply(as.matrix(models),1, function(m) c("AIC"=AIC(get(m)), "BIC"=BIC(get(m))))

# C'est donc le modèle ARMA(1,1) qui minimise les critères AIC et BIC tout en etant le seul
# à n'avoir aucune autocorrélation des residus et dont tous les coefficients sont significatifs

# Finalement on sélectionne le modèle ARMA(1,1) différencié à l'ordre 1 donc
# on selectionne le modèle ARIMA(1,1,1) pour la série initialement choisie

modele <- arima(ipi_pates_ts,c(1,1,1))

modele$coef
plot(modele$residuals)

jarqueberaTest(modele$residuals)
# On obtient une p-value de moins de 2.2e-16, on rejete donc l'hypothèse nulle de normalité des résultats


### Exportation des ressources en TEX ###
xtable(df)
xtable(criterias)

###############
## Partie III #
###############  

#On prédit les deux prochaines valeurs ainsi que leur IC à 95% à l'aide du package forecast
futurVal <- forecast(modele,h=2, level=c(95))
plot(futurVal, xlim = c(2018,2020.5))
