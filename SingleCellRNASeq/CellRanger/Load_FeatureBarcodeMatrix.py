import scipy.io, csv, gzip, os

mat_path = "/data/CTRL/"
mat = scipy.io.mmread(os.path.join(mat_path, "matrix.mtx.gz"))
