# Include the necessary packages
library(gtools)
library(partitions)

# Initialize digits
(digits <- c(0:9))

# List, allow repetition
# Permutate three out of 10 digits with repetition

# List, allow repetition
# gtools::permutations
(digits_lists <- permutations(10, 3, digits, repeats.allowed = TRUE))
# head(digits_lists, 50)
# the number of possible outcomes
nrow(digits_lists)
# the number of possible outcomes = 10*10*10
10^3

# List, allow repetition
#The first two elements are selected from the set "letters[1:2]", the middle three elements are selected from the set "digits[1:3]", and the last three elements are selected from the set "letters[24:26]".

# List, allow repetition
# base::expand.grid
(sample_lists <- expand.grid(letters[1:2], digits[1:3], letters[24:26]))
# the number of possible outcomes
nrow(sample_lists)
# the number of possible outcomes = 2*3*3
2*3*3

# Permutation, no repetition, order matters
# Permutate three out of 10 digits without repetition

# Permutation, no repetition
# gtools::permutations
(digits_perm <- permutations(10, 3, digits, repeats.allowed = FALSE))
# head(digits_perm, 50)
# the number of possible outcomes
nrow(digits_perm)
# the number of possible outcomes = 10!/(10-3)! = 10*9*8
factorial(10)/factorial(10-3)

# Combination, no repetition, order ignored
# Select three out of 10 digits without repetition

# Combination, no repetition
# gtools::combinations
(digits_comb_wo_repeat <- combinations(10, 3, digits, repeats.allowed = FALSE))
# head(digits_comb_wo_repeat, 50)
# the number of possible outcomes
nrow(digits_comb_wo_repeat)
# the number of possible outcomes = 10!/(10-3)!3! = (10*9*8)/(3*2*1) 
choose(10,3)

# Combination, allow repetition
# Select three out of 10 digits with repetition

# Combinations with repetition
# gtools::combinations
(digits_comb_w_repeat <- combinations(10, 3, digits, repeats.allowed = TRUE))
# head(digits_comb_w_repeat, 50)
# the number of possible outcomes

# [0] [1] [2] [3] [4] [5] [6] [7] [8] [9]
# move from left to right
# two actions: select (S) or shift right (->)
# 000: S S S -> -> -> -> -> -> -> -> ->
# 001: S S -> S -> -> -> -> -> -> -> ->
# ...
# 999: -> -> -> -> -> -> -> -> -> S S S
nrow(digits_comb_w_repeat)
# the number of possible outcomes = (3+9)!/(3!9!)
factorial(3+10-1)/(factorial(3)*factorial(10-1))

# Partition
# Divide 9 children into 3 teams of 3 each. How many different divisions? (3 teams are unordered)

# Partition
# number of group
n_g <- 3
# number of children in each group
n_g1 <- 3
n_g2 <- 3
n_g3 <- 3
# partitions::setparts
(team_unordered_partitions <- setparts(c(n_g1,n_g2,n_g3)))
# the number of possible outcomes
ncol(as.matrix(team_unordered_partitions))
# the number of possible outcomes = (9!/(3!3!3!))/3!
(factorial(n_g1+n_g2+n_g3)/(factorial(n_g1)*factorial(n_g2)*factorial(n_g3)))/factorial(n_g)