############### TESTS ##########################

# Data generation
gen_CRT_binarydata(
    ndatasets = 10,
    n1 = 20,
    n2 = 30,
    var_u0 = 0.6,
    p_intv = 0.7,
    p_ctrl = 0.2,
    batch_size = 10,
    seed = 4
)

# SSD null
# Finding number of clusters
SSD_crt_null_binary(
    p_intv = 0.5,
    p_ctrl = 0.2,
    n1 = 15,
    n2 = 30,
    ndatasets = 15,
    var_u0 = 0.6,
    BF_thresh1 = 2,
    eta1 = 0.8,
    fixed = "n1",
    b_fract = 2,
    max = 100,
    batch_size = 10,
    seed = 15
)

# Finding cluster size
SSD_crt_null_binary(
    p_intv = 0.5,
    p_ctrl = 0.2,
    n1 = 30,
    n2 = 30,
    ndatasets = 15,
    var_u0 = 0.6,
    BF_thresh1 = 2,
    eta1 = 0.7,
    fixed = "n2",
    b_fract = 2,
    max = 100,
    batch_size = 10,
    seed = 5
)

# SSD inf
# Cluster size
SSD_crt_inf_binary(p_intv = 0.8, p_ctrl = 0.5, n2 = 40, ndatasets = 100,
                   var_u0 = 0.1, BF_thresh1 = 4, eta1 = 0.8, fixed = "n2",
                   max = 300, batch_size = 100, seed = 45)

# Number of clusters
SSD_crt_inf_binary(p_intv = 0.8, p_ctrl = 0.5, n1 = 5, ndatasets = 100,
                   var_u0 = 0.1, BF_thresh1 = 4, eta1 = 0.8, fixed = "n1",
                   max = 300, batch_size = 100, seed = 45)
