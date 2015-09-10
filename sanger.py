import numpy as np
strg_A = "pGATGATGGCCTGAGCA"
col_A = list(strg_A)
strg_C = "pGATGATGGCCTGAGCATTC"
col_C = list(strg_C)
strg_G = "pGATGATGGCCTGAGCATTCG"
col_G = list(strg_G)
strg_T = "pGATGATGGCCTGAGCATT"
col_T = list(strg_T)
seq_matrix = list([col_A, col_C, col_G, col_T])

outz = []
for col in seq_matrix:
    print col
    for n in col:
        outz.append(n)
#print outz
