
# str1 and str2 are the columns (0-based) in f1 and f2 that are to be listed, in the same order, f1 first. The columns are separated by ;. 
# ignore the lines starting with "#" for both files
# print to stdout
def combine_files(f1, f2, str1, str2, sep): 
    # separate str1 and str2
    cols1 = str1.split(";")
    cols2 = str2.split(";")
    if sep == "NA":
	# if it s tab
	sep = "\t"
    fh = open(f1, "r")
    a = fh.readline().rstrip("\n")
    mat = []
    while(a != ""):
	if a[0] == "#":
	    continue
	arr = a.split(sep)
	tmp_arr = []
	for i in cols1:
	    tmp_arr.append(a[int(i)])
	mat.append(tmp_arr)
    fh.close()

    fh = open(f2, "r")
    a = fh.readline().rstrip("\n")
    i = 0
    while(a != ""):
	if a[0] == "#":
	    continue
	arr = a.split(sep)
	tmp_arr = mat[i]
	for j in cols2:
	    tmp_arr.append(a[int(j)])
	print sep.join(tmp_arr)
	i = i + 1
	


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
