#= Things to delete =#

a = [1 2 3]
b = [7 8 ]

c = [ a[i] + b[j] for i in 1:3, j in 1:2 ][:]
