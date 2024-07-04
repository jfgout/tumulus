# tumulus
Some R code to check the probability of random alignments with stars in a tumulus

This mini project started after an article in Ciel & Espace suggesting that the tumulus of Concœur contains stones aligned with specific stars.

The people involved in the study described by C&E have published some results from statistical analysis here: https://hal.science/hal-03223357v3

Here, I use random shuffling of the position of stones in a tumulus to compute probabilities of observing certain alignments. With a modest catalog of 13 objects (stars of magnitude < 1.0 + Bellatrix and Pleiades) I find a 55% probability that the Concœur tumulus would have an alignment with both set and rise of at least one object. I also find that random positions of stones produce, on average, alignments to the rise of 4 stars out of 13. This was based on a rough description of the Concœur tumulus and the results will be updated if a better description of the stones location and size becomes available. If the drawing of the tumulus that I used is not correctly to scale (if the stones are drawn bigger than they should), it could impact the results very significantly.

Usage:

Download the two csv files (objects.csv and tumulus_azimuths.csv) and the R code (tumulus_shuffle.R). The only modification you need to do in the R code is the "setwd" command. Change the folder to the location on your computer where you have downloaded the two csv files and run the code.
