###SIMULATIONS

##IMPORTATION DES PACKAGES

library(spatstat)
library(tidyverse)
library(ggplot2)

## DEFINITION DES 4 SOUS ZONES

#On d√©finit 4 sous zones dans la zone d'√©tude, chaque sous zone correspondant
#√† un type de sol : par exemple sable, boue, zoche et gravier
#on attribue une valeur fix√©e (1,2,3 ou 4) a chacune des 4 modalit√©s
#la fonction zonage prend pour argument les coordonn√©es x et y d'un point de la grille
zonage <- function(x,y, xmed = 5, ymed = 5){
  ifelse(x<xmed, 1, 2) + ifelse(y<ymed, 0, 2) 
}


## CALCUL DU CHAMP LATENT

#dans le champ latent, on ne prend pas en compte l'effet spatial al√©atoire,
#et on considere qu'une seule covariable (type de sol)
alpha_s <- -0.5 #on fixe l'intercept √† -0,5
betas_s <- c(0, 1, -1, 2) #on fixe les beta correspondant a chacune des modalit√©s de la
#covariable type de sol
#la fonction du champ latent prend pour argument les coordonn√©es x et y d'un point
#de la grille, les valeurs fix√©es de alpha_s et des beta_s et renvoie la valeur du 
#champ latent en ce point
latent_field_function = function(x,y, alpha_s, betas_s){
  exp(alpha_s +  betas_s[zonage(x,y)])
}

## REPRESENTATION GRAPHIQUE DU CHAMP LATENT DANS LA GRILLE

dta <- data.frame(x = runif(n = 100, 0,10 ), y = runif(100, 0, 10)) %>% 
  mutate(z= latent_field_function(x, y, alpha_s, betas_s))
ggplot(data= dta) + geom_point(aes(x=x, y=y, col = z))


## CALCUL DE L'INTENSITE LAMBDA

alpha_x <- 1
b <- 1  # b peut aussi valoir 0 ou 3

#fonction de calcul du lambda
lambda_function <- function(x, y, alpha_x, b, alpha_s, betas_s){
  lambda <- exp(alpha_x + b * log(latent_field_function(x,y, alpha_s, betas_s)))
}

## GENERATION DES POINTS DE PECHE

# on utilise la fonction rpoispp pour g√©n√©rer les points de peche
# La fonction rpoispp prend comme arguments :
# - lambda : la fonction de densite
# - lmax : pour fixer un maximum pour la densite, eviter qu'elle n'explose
# - win : la fenetre dans laquelle on simule nos points de densite, pour nous c'est 10*10
# - les autres arguments utilises par lambda (sauf x et y qui du coup sont definis par la fenetre win)

pp = rpoispp(lambda_function, lmax=4.5, win=owin(c(0,10),c(0,10)),
             alpha_x=alpha_x, b=b, alpha_s= alpha_s, betas_s = betas_s)
#On repr√©sente graphiquement les points de peche dans les 4 sous zones
plot(pp)



### ZERO INFLATING

Y = rep(0, len = length(pp$x))
# Vecteur qui contient un nombre par point de pÍche
# Pour l'instant ce nombre est fixÈ ‡ 0 mais aprËs il contiendra le nombre de poissons pÍchÈs ‡ ce point

# La fonction peche_attendue permet de calculer uf(xi) = qf * s(xi)
# Equation 8 du papier de Baptiste
# Elle prend en argument la position du point, la capturabilitÈ de la flotte, et retourne la peche attendue
q = 1
# On fixe la capturabilitÈ a 1
peche_attendue = function(x, y, q)
{
  q * latent_field_function(x, y, alpha_s, betas_s)
}

# La fonction p0_binomiale calcule, pour chaque point de peche, la probabilitÈ de ne rien pecher
# Equation 7 du papier de Baptiste
# Elle prend en argument la position du point, ksi (parametre qui controle l'intensitÈ de 0), q
ksi = 0
#  On fixe l'intensitÈ de 0, ksi, ‡ 0.001
p0_binomiale = function(x, y, ksi, q)
{
  exp(-exp(ksi) * peche_attendue(x, y, q))
}

# D'aprËs le papier de Baptiste, dans le cas d'un "succes" de Bernoulli) on n'a rien pÍchÈ
# Donc on fait 1 - proba pour avoir 1 si pÍche successfull, 0 si pÍche ratÈe
# Loi binomiale de cette proba pour dÈterminer le succËs ou l'Èchec de chaque point de pÍche
Y = 1 - rbinom(n=length(Y), size=1, prob=p0_binomiale(pp$x, pp$y, ksi, q))

table = data.frame(x=pp$x, y=pp$y, peche=Y)

# On affiche les pÍches successfull et les pÍches ratÈes
ggplot(table) + geom_point(aes(x=x, y=y, col=as.factor(peche)))

# Maintenant on regarde le nombre de poissons potentiellement pÍchÈs ‡ chaque point
# Pour Áa on suit une loi lognormale de moyenne pecheattendue/(1-psucces) et d'ecart type 1
L = rlnorm(n=length(Y), meanlog=peche_attendue(pp$x, pp$y, q)/(1 - p0_binomiale(pp$x, pp$y, ksi, q)), sdlog=0.1)

# Le nombre de poissons pechÈs = le nombre de poissons potentiellement peches (L) * le succes ou non de la peche (Y)
table$peche = L * Y
ggplot(table) + geom_point(aes(x=x, y=y, col=peche))
hist(L)



### Reallocation uniforme des captures

capturetotale = sum(table$peche)


## Si on considere 1 zone
table$L_unezone = mean(table$peche)
hist(table$L_unezone)

## Si on considere 2 zones

# On coupe selon l'axe des x, ‡ une valeur dÈfinie
xmed = 5
premierezone = table[which(table$x < xmed),]
deuxiemezone = table[which(table$x >= xmed),]
mean_premierezone = mean(premierezone$peche)
mean_deuxiemezone = mean(deuxiemezone$peche)
table$L_deuxzones = rep(0, nrow(table))
table[which(table$x < xmed), "L_deuxzones"] = mean_premierezone
table[which(table$x >= xmed), "L_deuxzones"] = mean_deuxiemezone
hist(table$L_deuxzones)

## Si on considere 4 zones

# On coupe selon l'axe des x et l'axe des y, ‡ des valeurs dÈfinies
xmed = 5
ymed = 5
premierezone = table[which(table$x < xmed & table$y < ymed),]
deuxiemezone = table[which(table$x < xmed & table$y >= ymed),]
troisiemezone = table[which(table$x >= xmed & table$y < ymed),]
quatriemezone = table[which(table$x >= xmed & table$y >= ymed),]
mean_premierezone = mean(premierezone$peche)
mean_deuxiemezone = mean(deuxiemezone$peche)
mean_troisiemezone = mean(troisiemezone$peche)
mean_quatriemezone = mean(quatriemezone$peche)
table$L_quatrezones = rep(0, nrow(table))
table[which(table$x < xmed & table$y < ymed), "L_quatrezones"] = mean_premierezone
table[which(table$x < xmed & table$y >= ymed), "L_quatrezones"] = mean_deuxiemezone
table[which(table$x >= xmed & table$y < ymed), "L_quatrezones"] = mean_troisiemezone
table[which(table$x >= xmed & table$y >= ymed), "L_quatrezones"] = mean_quatriemezone
hist(table$L_quatrezones)
# Si la diff entre les 4 groupes est pas visible, faire :
hist(table$L_quatrezones, breaks = 100)
# Et augmenter le 100 jusqu'‡ ce que ce soit visible
# C pcq certains groupes sont trop proches (genre la ou les peches sont faibles)


# Representation graphique

# Vraies peches
ggplot(table) + geom_point(aes(x=x, y=y, col=peche))

# Separation en une zone
ggplot(table) + geom_point(aes(x=x, y=y, col=as.factor(L_unezone)))

# Separation en deux zones
ggplot(table) + geom_point(aes(x=x, y=y, col=as.factor(L_deuxzones)))

# Separation en quatre zones
ggplot(table) + geom_point(aes(x=x, y=y, col=as.factor(L_quatrezones)))
