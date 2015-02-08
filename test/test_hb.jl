using HarwellRutherfordBoeing

hb_matrices = ["lock1074.pse" "well1850.rra" "young3c.csa" ]

for hb_matrix in hb_matrices
  M = HarwellBoeingMatrix(hb_matrix)
  print(M)
  @printf("\n")
end
