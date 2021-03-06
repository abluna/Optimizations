Market Optimization Using an Evolutionary Algorithm
================
Abraham B.
10/5/2019

## R Markdown

Conjoint modelling is one of the most flexible survey-based analysis
that allows the user to isolate the importance of each feature tested.
And even more useful is the ability to create ‘what-if’ scenarios. For
instance, I can simulate changes in the market if Brand A changes or
introduces a new product, and in addition, I can simulate those same
changes in the market if the competitor also chooses to make changes to
their current product portfolio.

Given that flexibility, it is possible then to identify the best (or
optimal) strategy. This objective isn’t trivial–both computational and
strategically (how do we define “optimal”? Is it generating most
interest? maximizing profit? And if profit, long or short term? Or is
market share?). But for the sake of this exercise, we’ll assume demand
is of outmost importance and we will assume that all configurations are
feasible (in reality, a better feature will be more costly to produce)

``` r
## Bring in data
Utilities <- read.csv("C:/Users/abbal/Documents/LRW - DCM Bayesian Optimization/Utilities for Optimization.csv")

## In this simulation, a product must choose 23 features, and here are the number of features within.
## (For instance, a feature (or attribute) can be battery life, and there are three options (levels): 8/16/24hr battery life)
Level_Count <- c(3, 4, 4, 4, 5, 4, 3, 5, 4, 3, 3, 4, 5, 3, 3, 7, 7, 7, 4, 2, 2, 2, 2)

## for the purposes of this analysis, I'm masking all  features tested. But the application is universal among most DCMs.
```

## Objective

The objective of this analysis is to optimize “share of preference”,
which is the shared sum of the exponential values of the levels
selected. In other words, which combination of features will yield the
greatest demand.

The number of combinations to test is the product of `Level_Count`. And
in this case, the total number of combinations equals 8,193,540,096,000.
It would take an unreasonable amount of time to test all combinations.
Fortunately, there are methods that can work towards finding the optimal
combination, leaving behind poorly performing combinations.

We will use a package called `gramEvol` to run an evolutionary
algorithm. And this will require:

  - A scoring function
  - A monitoring function (optional, but is used to monitor the
    performance of the algorithm)
  - The genetic algorithm that will use the two above

<!-- end list -->

``` r
## Scoring function that tests a set of levels activated
Get_Demand <- function(x){

  LevelSumCum <- c(0,cumsum(Level_Count))
  LevelsActivate1 <- LevelSumCum[-length(LevelSumCum)] + x

  Configuration_Matrix <- matrix(0, nrow = 2, ncol = sum(Level_Count)+1)
  Configuration_Matrix[1, LevelsActivate1] <- 1
  Configuration_Matrix[2, ncol(Configuration_Matrix)] <- 1
  
  XB <- t(Configuration_Matrix %*% t(Utilities[,-1]))
  XB_exp <- exp(XB)
  
  XB_Prop <- XB_exp / apply(XB_exp, 1, sum)
  
  SOP <- apply(XB_Prop, 2, mean)
  
  Product_Demand <- SOP[1]
  
  ## must return negative value because the
  ## genetic algorithm seeks to minimize the output
  ## and our goal is to maximize demand
  return(Score=-1*Product_Demand)
  
}
```

``` r
## Monitoring function
monitorFunc1 <- function(result){
  #cat("Best of gen", result$population$currentIteration,":", round(-1*min(result$best$cost),3), "\n")
  points(result$population$currentIteration, -1*min(result$best$cost), cex = 1, col = "dark red", pch = 19)
  
  if(result$population$currentIteration > 1){
    XY <- rbind(c(result$population$currentIteration - 1, -1*result$population$best[result$population$currentIteration - 1]),
                c(result$population$currentIteration, -1*min(result$best$cost)))
    lines(XY)
  }
  
}
```

## Evolution function

Luckily, this function was written to leverage the power of
parallelization. I will use all 8 cores in my machine. The three main
parameters parameters to test are:

Population Size: How many configurations do you want to have at each
step (the more you have the longer it takes to run) Iterations: How many
steps to find the optimal configuration. Mutation: How much of each
winner should change

From my experience, I like to choose a smaller population with more
iterations, that way, I’m only focusing on the best performing
configurations. However, this approach tends to converge around a local
optimum, which is why I also like to increase the mutation probability

``` r
library(doParallel)
library(gramEvol)
cl <- makeCluster(8)
registerDoParallel(cl)
clusterExport(cl,c('monitorFunc1', 'Get_Demand','Utilities','Level_Count'))

GA_Plot <- plot(1, type="n", xlab="Iteration", ylab="Performance", 
                main = "Optimization Performance\nGenetic Algorithm", 
                xlim=c(0, 100), ylim=c(0.7, 0.90))

Results1 <- GeneticAlg.int(genomeLen = length(Level_Count), 
                          genomeMin = rep(1, length(Level_Count)), 
                          genomeMax = Level_Count,
                          popSize = 100,
                          iterations = 100, 
                          monitorFunc = monitorFunc1, 
                          mutationChance = 0.4,
                          evalFunc = Get_Demand, 
                          verbose = FALSE, 
                          plapply = function(...) parLapply(cl, ...))
```

<img src="Optimization-Markdown_files/figure-gfm/Genetic Algorithm-1.png" style="display: block; margin: auto;" />

``` r
## the best level to activate for each attribute
Results1
```

    ## Genetic Algorithm Search Results:
    ##   No. Generations: 100 
    ##   Best Genome:     3 1 3 1 4 2 2 2 3 2 2 3 3 2 1 1 1 1 2 1 2 1 1 
    ##   Best Cost:       -0.843858611111133

## Results

I only let the algorithm run for a few min. But this requires many
iterations to achieve something close to the optimal solution. The
longer it runs, the more difficult it is to see a bump in performance
(think BitCoin mining\!).

Because it is quite easy to identify a configuration with the highest
demand (usually the best performing level from each feature), it doesn’t
take long for the algorithm to reach that performance score. Where it
becomes challenging is when the criteria changes to something more
complex such as: Identify a portfolio of 5 products that will generate
the highest unique demand, while reducing the probability of
cannibalization of products, but that will come later.
