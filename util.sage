def vector_to_matrix(v,n):
    res = matrix(Fq,n)
    count = 0
    for i in range(n):
        for j in range(n):
            res[i,j] = v[count]
            count += 1
    return res