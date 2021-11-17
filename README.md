https://github.com/marc-l-m-francois/Funambule2

Included files
by Marc Louis Maurice FRANCOIS
marc.francois@univ-nantes.fr

# PROGRAMS

Programs are in capital letters, routines are in lowercase letters. The list of main programs is presented here.

## FUNAMBULE
`FUNAMBULE` is dedicated to curve or silhouette boundaries precise measurements. It uses the VIC method developed by myself (see publications in `Doc/Bibliography` folder). It is basically a DIC algorithm in which the image is compared to a virtual, user chosen, mathematical curve. The precision of the measurement is sub-pixel and better than every image softwares it has been compared to. Presently available curves are:
* B-Splines of various orders, closed or opened
* CAD_IBOX is an example of CAD curve
* Circle
* Ellipse
* Polygon of any number of summits
* Beams adresses any problem of type `y=f(x)`:
	* Cantilever beam with concentrated force at end
	* Chainette curve
	* Cubic polynomial equation
	* PureFlexion (beam loaded by two opposite momentums)
	* Segment (straight line)

## DIGIMCO
DigImCo is a more classical Integrated *Digital Image Correlation software* (IDIC). The available full-fields are (in folder `Champs`):
* Corsol for rigid body motion
* Defhom for an homogenous strain field
* EbHomHPP for an Euler-Bernoulli field of a beam in flexion
* Q4 for a Q4 finite element kinematic field

## YOUNG_POISSON
Using the results of a DigImCo calculus with Defhom, this program extracts Young's modulus and Poisson's ratio from a series of images of a tension or compression test. 


# FOLDERS CONTENT

`Champs`: contains the displacement fields used by DigImCo
`Curves`: contains the curve equations used by Funambule
`Docs`: includes the bibliography on both VIC and DIC and tutorial
`Images`: is the default folder for the images. It actually contains examples for practise.
`Routines`: includes the subroutines of both Digico, Digimco_Q4 and Funambule
`Routines_graphiques`: includes the graphical subroutines of both Digico, Digimco_Q4 and Funambule
`Temporaire`: is the temporary folder in which current images are stored
