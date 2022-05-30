"""
GUEVEL Antoine
MOREAU Quentin 

TP Chimie : Borax

Regression lineaire

"""




import numpy as np
import matplotlib.pyplot as plt
from math import log, sqrt






# régression linéaire



def moyenne(liste):
    somme = 0
    for i in liste:
        somme += i
    return somme/len(liste)

def variance(liste):
    moyenneListe = moyenne(liste)
    somme = 0
    for i in liste:
        somme += (i-moyenneListe)*(i-moyenneListe)
    return somme/len(liste)

def covariance(liste1,liste2):
    moyenne1 = moyenne(liste1)
    moyenne2 = moyenne(liste2)
    numerateur = 0
    for i,j in zip(liste1,liste2):
        numerateur += (i-moyenne1)*(j-moyenne2)
    return numerateur/len(liste1)






constanteDesGazParfaits = 8.314

# mesures des groupes
temperatures = [54.0, 0.5, 45.1, 53.3, 2.9, 35.0, 53.5, 0.5, 40.0, 0.4, 53.6, 29.7, 0.0, 53.5, 36.9, 1.2, 30.0, 52.0, 0.9, 29.6, 50.5, 0.0, 54.0, 37.0, 0.9, 29.6, 0.7, 34.0, 1.0, 43.7, 53.0, 2.1, 35.0, 53.0, 1.0, 45.0, 0.5, 36.5, 0.9, 43.7, 0.9, 34.0]
volumesEquivalences = [31.3, 9.6, 17.9, 28.8, 10.3, 12.1, 27.8, 8.9, 10.0, 9.5, 28.5, 7.0, 10.3, 31.4, 14.4, 10.0, 6.5, 25.5, 12.0, 7.5, 12.0, 10.0, 30.0, 13.0, 10.0, 7.5, 9.7, 9.5, 8.6, 8.9, 27.7, 10.2, 13.1, 27.9, 10.6, 17.8, 9.9, 10.2, 9.4, 9.3, 9.3, 9.8]
volumesPreleves = [5.0, 20.0, 5.0, 5.0, 20.0, 5.0, 5.0, 20.0, 5.0, 20.0, 5.0, 5.0, 20.0, 5.0, 5.0, 20.0, 5.0, 5.0, 20.0, 5.0, 5.0, 20.0, 5.0, 5.0, 20.0, 5.0, 20.0, 5.0, 20.0, 5.0, 5.0, 20.0, 5.0, 5.0, 20.0, 5.0, 20.0, 5.0, 20.0, 5.0, 20.0, 5.0]

pointsAberrant = [20,29,39]
for i in pointsAberrant[::-1]:
    temperatures.pop(i)
    volumesEquivalences.pop(i)
    volumesPreleves.pop(i)


concentration = 0.200

#on enregistre dans ordonnees -Rln(K(T)) et dans abscisses 1/T
ordonnees = [-constanteDesGazParfaits*(log(4)+3*log(0.2*concentration*volumesEquivalences[i]/volumesPreleves[i])) for i in range(len(temperatures))]
abscisses = [1/(temperatures[i]+273.15) for i in range(len(temperatures))]







# affichage



covarianceListe1Liste2 = covariance(ordonnees,abscisses)

coefficientCorrelation = covarianceListe1Liste2/(sqrt(variance(abscisses)*variance(ordonnees)))

coefficientDirecteur = covarianceListe1Liste2/variance(abscisses)
ordonneeOrigine = moyenne(ordonnees)-coefficientDirecteur*moyenne(abscisses)

print(coefficientDirecteur)
print(ordonneeOrigine)

droite = [coefficientDirecteur*abscisses[i]+ordonneeOrigine for i in range(len(abscisses))]

f = plt.figure()
ax = f.add_subplot(111)

plt.text(0.5,0.9,"coefficient de corrélation : " + str(coefficientCorrelation), horizontalalignment = "center", transform = ax.transAxes)
plt.text(0.02,0.4,"coefficient directeur : " + str(coefficientDirecteur) + r" $J.mol^{-1}$", transform = ax.transAxes)
plt.text(0.02,0.2,"ordonnée a l'origine : " + str(ordonneeOrigine) + r" $J.mol^{-1}.K^{-1}$", transform = ax.transAxes)

plt.plot(abscisses, ordonnees, "bx", label="nuage de points")
plt.plot(abscisses, droite, "r", label="régression linéaire")

plt.xlabel(r"$1/T$ (en $K^{-1}$)")
plt.ylabel(r"$-Rln(K(T))$ (en $J.mol^{-1}.K^{-1}$)")
plt.title("Régression linéaire")

plt.legend()
plt.show()


