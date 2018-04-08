# HCA Jamboree 2: Where should I more deeply sequence?

The objective of this task is to select sequential experiments such that we find the most accurate estimates of densities across the space of cell states (as defined by transcriptional profiles) with the smallest number of experiments. 

# Definitions for this task

* **Objective function:** Define a function to optimize in the process of iterative experiment design 
  with the fewest number of cells assayed and a limited number of experiments. Possible objective 
  functions include: likelihood of a held out test set, maximizing coverage of the space, such as 
  defining `k` nearest neighborhoods (k-nns) and maximizing the number of k-nn neighborhoods that 
  are covered, maximizing diversity or silhouette methods. It is worthwhile to spend some time 
  thinking of the most appropriate objective here. 
* **Density estimation:** capture the probability of gene expression levels `x_i` across the space of 
  available cell states. In other words, estimate a probability density function for gene 
  expression levels across the space of cell states.
  * Some possible density estimates include mixture models, `K` nearest neighbors or something more 
    sophisticated.
  *	Important methodological challenges include estimates of the means, variance, or covariance 
    terms, high-dimensional density estimates in a low dimensional space; and regularizing these 
    estimates appropriately.
*	**Gate, or Experiment:** to mimic flow cytometry experiments, an Experiment is defined by a set of 
  thresholds (gates) for up to six gene markers or cell surface proteins that may be used to sort 
  the cells, and a value k between 100 and the total remainder that specifies the number of cells 
  to assay.
  *	The simulated experiment will be performed as follows: given the data set processed through 
    MAGIC, we will perform simulated FACS sorting on the chosen genes or cell surface proteins. We 
    will uniformly sample `k` cells (user specifies `k`, within a total budget of cells) from all 
    of the cells that are identified by the experiment. The returned data will be in the cell count 
    matrix format.
* **Sampling Strategy:** Strategy or data driven the function used to select the appropriate next 
  experiment.
* Possible acquisition functions include expected improvement (EI), UCB, uniform sampling and 
    others.

# Task sequence

## Day 1: 

### For the first scRNA-seq dataset (human bone marrow),

1. We will start the task with a random uniform sub-sample of the larger dataset containing 15K of the original cells. Groups will receive this subset of cells as a cells by molecules count matrix. 
2. You will design and a density estimator using these data as well as a function to compare the density between a subsample and full dataset and implement these. 
3. You will design a sampling strategy and acquisition function and use these to sample more cell from the data.   Given a budget of 6 experiments and an additional 25K cells you will:
   1. Propose an experiment using your chosen acquisition function that defines a sampling gate and number of cells to sample from this gate. 
   2. Receive a gene count matrix of cells resulting from that experiment
   3. Re-estimate the transcriptional densities 
   4. Compute the objective function, as well as receive the improvement in the objective function we defined.  

### Validation:

1.  We will evaluate your success on the final set of 15+25K cells, both using our objective function and the one you define on the full dataset and compare it to a baseline of subsampling of 40K cells. 

### Visualization and Reflection:

1. Can you think of ways to visualize the density? 
2. Can you think of metrics and ways to score and visualize regions of phenotypic space that are well sampled, verses regions that contain a degree of uncertainty (require more sampling). 
3. Following the first attempt, can you think of ways to improve on your sample strategy, acquisition function and/or objective function? 

## What you should produce:
* A software tool to perform this task available in GitHub – including 
  * Evaluation of the objective function
  * Density estimate for a collection of cells
  * Density estimate updates for an additional collection of cells
  * A tool that determines what additional experiments to perform, including gates and #cells per gate
  * Evaluation of the acquisition function
  * A tool to estimate regions where density estimation is confident (additional sampling unlikely to change it). 
  * Visualization tools. 
* Documentation
  * Code and pipelines should be clearly documented
  * Alternative approaches should be quantified and evaluated
  * Choices about methodologies should be motivated and described in a write up
  * Ultimately: journal publication with methods developed in your group

## Note:

Additional validation datasets and information will be provided on day 2. 

## Validation

The task will be validated in two ways. We will reveal these metrics towards the end of the Jamboree to avoid biasing the methods the groups use to solve this task. 
