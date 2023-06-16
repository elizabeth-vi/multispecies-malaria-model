p_mask = 0.50
mda1 = [270.0]  # lower bound of when MDA occurs
mda2 = [270+30] # upper bound of when MDA occurs (i.e. 30 days of MDA)
c_vec = [[0.3, 0.3]] # coverage scenarios
pP_vec = [0.217]
pG_vec = [[0.0066, 0.000006]]
pN_vec = [[1.0, 0.01]]#0.843]] # treatment scenarios
mda_coverage = [0.5]
pN_mda_vec = [[[1.0, 0.01]]]#0.843]]]
MDA_vec = [False]

FSAT_vec = [False] #focused screening and treatment

    # "General": {
    #     "p_mask": 0.50,
    #     "mda1": [270.0],  # lower bound of when MDA occurs
    #     "mda2": [270+30], # upper bound of when MDA occurs (i.e. 30 days of MDA)
    #     "c_vec": [[0.3, 0.3]], # coverage scenarios
    #     "pP_vec": [0.217],
    #     "pG_vec": [[0.0066, 0.000006]],
    #     "pN_vec": [[1.0, 0.843]], # treatment scenarios
    #     "mda_coverage": [0.5],
    #     "pN_mda_vec": [[[1.0, 0.843]]],
    #     "MDA_vec": [False],
    #     "FSAT_vec": [False] #focused screening and treatment
    # },