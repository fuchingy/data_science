## LNp.3-2

## define x y z(sample space)
x <- c(1,2,5)
y <- c(5,1,8,9)
z <- c(1,2,3,4,5,6,7,8,9)

## set union
union(x, y)
## intersect(x,y): Intersection of the sets x and y
intersect(x, y)
## setdiff(x,y): Set difference between x and y, consisting of all elements of x that are not in y
setdiff(x, y)
setdiff(y, x)
## x's complement
setdiff(z,x)
## setequal(x,y): Test for equality between x and y
setequal(x, y)
setequal(x,c(1,2,5))

## True for all possible x & y :
setequal( union(x, y),
          c(setdiff(x, y), intersect(x, y), setdiff(y, x)))
## c %in% y: Membership, testing whether c is an element of the set y
## this syntax are more common than "is.elemet(x,y)"
2 %in% x
2 %in% y
## is.element(x, y) is identical to x %in% y.
is.element(2, x) 
is.element(2, y) 

## choose(n,k): Number of possible subsets of size k chosen from a set of size n
choose(5,2)




