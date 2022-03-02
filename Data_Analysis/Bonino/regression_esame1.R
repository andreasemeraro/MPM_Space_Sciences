
R version 3.6.1 (2019-07-05) -- "Action of the Toes"
Copyright (C) 2019 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64 (64-bit)

R è un software libero ed è rilasciato SENZA ALCUNA GARANZIA.
Siamo ben lieti se potrai redistribuirlo, ma sotto certe condizioni.
Scrivi 'license()' o 'licence()' per dettagli su come distribuirlo.

R è un progetto di collaborazione con molti contributi esterni.
Scrivi 'contributors()' per maggiori informazioni e 'citation()'
per sapere come citare R o i pacchetti di R nelle pubblicazioni.

Scrivi 'demo()' per una dimostrazione, 'help()' per la guida in linea, o
'help.start()' per l'help navigabile con browser HTML.
Scrivi 'q()' per uscire da R.

[Caricato workspace precedentemente salvato]

> hyades <- read.table("Sun_data.dat", header=TRUE, fill=TRUE)
Error in file(file, "rt") : non posso aprire questa connessione
Inoltre: Warning message:
In file(file, "rt") :
  cannot open file 'Sun_data.dat': No such file or directory
> 
Error in file(file, "rt") : non posso aprire questa connessione
Inoltre: Warning message:
In file(file, "rt") :
  cannot open file 'Sun_data.txt': No such file or directory
> hyades <- read.table("Sun_data.txt", header=TRUE, fill=TRUE)
Error in file(file, "rt") : non posso aprire questa connessione
Inoltre: Warning message:
In file(file, "rt") :
  cannot open file 'Sun_data.txt': No such file or directory
> hyades <- read.table("Sun_data.txt", header=TRUE, fill=TRUE)
Error in file(file, "rt") : non posso aprire questa connessione
Inoltre: Warning message:
In file(file, "rt") :
  cannot open file 'Sun_data.txt': No such file or directory
> hyades <- read.table("Sun_data.dat", header=TRUE, fill=TRUE)
> attach(hyades)
> plot(flux_1GHz, flux_2GHz, main="Flux_2GHz vs Flux_1GHz",ylab="flux_2GH", xlab="flux_1GHz")
> cor(hyades)
               time flux_1GHz flux_2GHz
time      1.0000000 0.2934571 0.2132121
flux_1GHz 0.2934571 1.0000000 0.7032788
flux_2GHz 0.2132121 0.7032788 1.0000000
> cor.test(flux_2GHz, flux_1GHz, alternative="two.sided", conf.level=0.95)

        Pearson's product-moment correlation

data:  flux_2GHz and flux_1GHz
t = 59.66, df = 3637, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.6864730 0.7193336
sample estimates:
      cor 
0.7032788 

> model1 = lm(flux_2GHz ~ flux_1GHz)
> summary(model1)

Call:
lm(formula = flux_2GHz ~ flux_1GHz)

Residuals:
    Min      1Q  Median      3Q     Max 
-5.6572 -0.9011 -0.1815  0.5205 14.4760 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept) 10.59862    0.72290   14.66   <2e-16 ***
flux_1GHz    0.95593    0.01602   59.66   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.753 on 3637 degrees of freedom
Multiple R-squared:  0.4946,    Adjusted R-squared:  0.4945 
F-statistic:  3559 on 1 and 3637 DF,  p-value: < 2.2e-16

> abline(model1, col = "red")
> res<- residuals(model1)
> > plot(res, main = "Residuals")
Errore: unexpected '>' in ">"
> > abline(0,0, col="red")
Errore: unexpected '>' in ">"
> plot(res, main = "Residuals")
> abline(0,0, col="red")
> library("lmtest")
Carico il pacchetto richiesto: zoo

Attaching package: ‘zoo’

The following objects are masked from ‘package:base’:

    as.Date, as.Date.numeric

Warning messages:
1: package ‘lmtest’ was built under R version 3.6.3 
2: package ‘zoo’ was built under R version 3.6.3 
> modello <- formula(model1)
> dwtest(modello, alternative="two.sided")

        Durbin-Watson test

data:  modello
DW = 0.2499, p-value < 2.2e-16
alternative hypothesis: true autocorrelation is not 0

> qqnorm(res, main = "Normal Q-Q Plot", xlab = "Theoretical 
+ Quantiles", ylab = "Residuals quantiles")
> qqline(res, distribution=qnorm, col="red")
> bgtest(modello, order=15)

        Breusch-Godfrey test for serial correlation of order up to 15

data:  modello
LM test = 2984.5, df = 15, p-value < 2.2e-16

> f1<-hyades[,2]
> f2<-hyades[,3]
> f1diff<-diff(f1,differences=1)
> f2diff<-diff(f2,differences=1)
> 
> attach(f1diff)
Error in attach(f1diff) : 
  'attach' only works for lists, data frames and environments
> plot(f1diff, f2diff, main="Flux_2GHz vs Flux_1GHz differenciated",ylab="flux_2GHz diff", xlab="flux_1GHz diff")
> plot(f2diff, f1diff, main="Flux_2GHz vs Flux_1GHz differenciated",ylab="flux_2GHz diff", xlab="flux_1GHz diff")
> model1 = lm(flux_2GHz ~ flux_1GHz)
> model2 = lm(f2diff ~ f1diff)
> abline(model2, col = "red")
> res1<- residuals(model2)
> plot(res, main = "Residuals")
> plot(res1, main = "Residuals")
> abline(0,0, col="red")
> cor.test(f2diff, f1diff, alternative="two.sided", conf.level=0.95)

        Pearson's product-moment correlation

data:  f2diff and f1diff
t = 3.3177, df = 3636, p-value = 0.0009165
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.02248064 0.08727864
sample estimates:
       cor 
0.05493748 

> modello2<-formula(model2)
> dwtest(modello, alternative="two.sided")

        Durbin-Watson test

data:  modello
DW = 0.2499, p-value < 2.2e-16
alternative hypothesis: true autocorrelation is not 0

> summary(model2)

Call:
lm(formula = f2diff ~ f1diff)

Residuals:
     Min       1Q   Median       3Q      Max 
-12.3123  -0.3382  -0.0034   0.3442   3.6851 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept) -0.003902   0.011195  -0.349 0.727456    
f1diff       0.059624   0.017971   3.318 0.000917 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.6752 on 3636 degrees of freedom
Multiple R-squared:  0.003018,  Adjusted R-squared:  0.002744 
F-statistic: 11.01 on 1 and 3636 DF,  p-value: 0.0009165

> dwtest(modello2, alternative="two.sided")

        Durbin-Watson test

data:  modello2
DW = 2.8213, p-value < 2.2e-16
alternative hypothesis: true autocorrelation is not 0

> qqnorm(res, main = "Normal Q-Q Plot", xlab = "Theoretical 
+ 
+ 
+ )
+ 
+ 
Errore: unexpected end of input
> qqnorm(res1,main = "Normal Q-Q Plot", xlab = "Theoretical 
+  Quantiles", ylab = "Residuals quantiles")
> qqline(res, distribution=qnorm, col="red")
> 
