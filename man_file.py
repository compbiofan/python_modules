#!/home/xf2/.conda/envs/my_root/bin/python
# this is a file containing functions to manipulate file IO
def read_csv_f(f, sep, header, col_a):
    fh = open(f, "r")
    a = fh.readline().rstrip("\n")
    ret_a = []
    while(a != ""):
        arr = a.split(sep)
        if header:
            # header is in the file, ignore it
            header = False
            a = fh.readline().rstrip("\n")
            continue
        if len(col_a) == 1:
            ret_a.append(arr[col_a[0]])
        else:
            # return a matrix
            tmp_a = []
            for i in col_a:
                tmp_a.append(arr[i]) 
            ret_a.append(tmp_a)
        a = fh.readline().rstrip("\n")
    fh.close()
    return ret_a

# given a matrix, output it to a file
def print_matrix(m, wHeader, column, row, f):
    fh = open(f, "w+")
    if wHeader:
        row = "\t".join([""] * len(column[0])) + "\t" + "\t".join(row)
        fh.write(row + "\n")
    for i in range(len(m)):
        if wHeader:
            fh.write("\t".join(column[i]) + "\t")
        fh.write("\t".join(m[i]) + "\n")
    fh.close()
