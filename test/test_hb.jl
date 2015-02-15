using HarwellRutherfordBoeing

hb_matrices = ["lock1074.pse" "well1850.rra" "young3c.csa" "illc1033.rra" "mahindas.rua"]

for hb_matrix in hb_matrices
  M = HarwellBoeingMatrix(hb_matrix)
  show(STDOUT, M)
  print(M)
  @printf("\n")
end

rb_matrices = ["commanche_dual.rb", "conf5_4-8x8-10.rb"]

for rb_matrix in rb_matrices
  M = RutherfordBoeingData(rb_matrix)
  show(STDOUT, M)
  print(M)
  @printf("\n")
end

# Sample supplementary data files from the technical report.
rb_data = ["ordering.rb", "denserhs.rb", "sparserhs.rb", "estimate.rb", "solution.rb"]

for rb_datum in rb_data
  M = RutherfordBoeingData(rb_datum)
  show(STDOUT, M)
  print(M)
  @printf("\n")
end
