using HarwellRutherfordBoeing

hb_matrices = ["lock1074.pse" "well1850.rra" "young3c.csa" ]

for hb_matrix in hb_matrices
  M = HarwellBoeingMatrix(hb_matrix)
  print(M)
  @printf("\n")
end

rb_data = ["commanche_dual.rb"]

for rb_datum in rb_data
  M = RutherfordBoeingData(rb_datum)
  print(M)
  @printf("\n")
end
