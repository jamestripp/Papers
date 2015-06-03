# Distinctiveness and range

A lot of data in the subjective judgement literature shows us that the range position influences judgements. This finding is apparently incompatable with explanations based on ranks.

My analysis shows how a rank based account is consistent with the range effects. If items are weighted by their distinctiveness then range effects can be seen. 

This folder contains the code used to analyse data from 4 previous studies. I am unable to make the data public because it was collected by other researchers. 

This code shows several cool techniques which would could be applied to any other project:

1. A grid search is used. The grid contains equally spaced parameter combinations. At each point a maximum likelihood estimate is carried out using a function minimisation algorithm. In effect, we're running thousands of model fits with different starting parameters. Then we take the best. We do thi because the function mimisation algorithm may get stuck.
2. To speed things up I use the 'forEach' package in R to do parallel processing. This means that on my 8 processor machnine there are 8 models fits running at a given time; one on each processor. The thousands of model fits are completed 8 times faster.
3. The SIMPLE model is implemented using matrices instead of for loops. This code operates on all the distances at one time. Predictions are calculated quicker and each grid point is evaluated faster.

If you have any questions or would like to request the data, then feel free to contact me on james.tripp@warwick.ac.uk.
