using HarwellBoeing

matrices = ["lock1074.pse" "well1850.rra" "young3c.csa" ]

for matrix in matrices
  M = HarwellBoeingMatrix(matrix)
  print(M)
  @printf("\n")
end
