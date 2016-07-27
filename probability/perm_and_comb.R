library(gtools)
library(partitions)

digits <- c(0:9)

# List: 10*10*10
# (O, O, O);
(digits_lists <- permutations(10, 3, digits, repeats.allowed = TRUE))
nrow(digits_lists)

# List: 2*3*3
# expand.grid {base}
(sample_lists <- expand.grid(letters[1:2], digits[1:3], letters[24:26]))
nrow(sample_lists)

# Permutations: 10*9*8
# 10!/(10-3)!
(digits_perm <- permutations(10, 3, digits, repeats.allowed = FALSE))
nrow(digits_perm)

# Combinations without repetition: select three out of 10 digits without repetition
# (10*9*8)/(3*2*1)
# 10!/(10-3)!3!
(digits_comb_wo_repeat <- combinations(10, 3, digits, repeats.allowed = FALSE))
nrow(digits_comb_wo_repeat)
choose(10,3)

# Combinations with repetition: select three out of 10 digits with repetition
(digits_comb_w_repeat <- combinations(10, 3, digits, repeats.allowed = TRUE))
nrow(digits_comb_w_repeat)

# [0] [1] [2] [3] [4] [5] [6] [7] [8] [9]
# move from left to right
# two actions: select (S) or shift right (->)
# 000: S S S -> -> -> -> -> -> -> -> ->
# 001: S S -> S -> -> -> -> -> -> -> ->
# ...
# 999: -> -> -> -> -> -> -> -> -> S S S
# (3+9)!/3!9!

# 9 children divided into 3 teams of 3 each. How many different divisions?
# 3 teams are unordered
# (9!/(3!3!3!))/3!
(team_unordered_partitions <- setparts(c(3,3,3)))
ncol(as.matrix(team_unordered_partitions))
