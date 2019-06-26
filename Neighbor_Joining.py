import numpy as np

def lowest_cell(table):
    min_cell = float("inf")
    x, y = -1, -1
    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j] < min_cell:
                min_cell = table[i][j]
                x, y = i, j
    return x, y



def join_labels(labels, a, b):
    # Swap if the indices are not ordered
    if b < a:
        a, b = b, a

    # Join the labels in the first index
    labels[a] = "(" + labels[a] + "," + labels[b] + ")"

    # Remove the (now redundant) label in the second index
    del labels[b]


# join_table:
#   Joins the entries of a table on the cell (a, b) by averaging their data entries
def join_table(table, a, b):
    # Swap if the indices are not ordered
    if b < a:
        a, b = b, a

    # For the lower index, reconstruct the entire row (A, i), where i < A
    row = []
    for i in range(0, a):
        row.append((table[a][i] + table[b][i] - table[b][a])/2)
    table[a] = row
    
    # Then, reconstruct the entire column (i, A), where i > A
    #   Note: Since the matrix is lower triangular, row b only contains values for indices < b
    for i in range(a+1, b):
        table[i][a] = (table[i][a]+table[b][i] - table[b][a])/2
        
    #   We get the rest of the values from row i
    for i in range(b+1, len(table)):
        table[i][a] = (table[i][a]+table[i][b]-table[b][a])/2
        # Remove the (now redundant) second index column entry
        del table[i][b]

    # Remove the (now redundant) second index row
    del table[b]

def Q_matrix(table):
    dd=np.size(table,1)
    Q=np.zeros(dd*dd).reshape(dd,dd)
    for i in range(dd):
        for ii in range(i+1,dd):
            q=(dd-2)*table[i,ii]-np.sum(table,0)[ii]-np.sum(table,1)[i]
            Q[ii,i]=q
    return Q

def list_to_matrix(l):
    dd=len(l)
    z=np.zeros(dd*dd).reshape(dd,dd)
    for i in range(dd):
        dd1=len(l[i])
        for ii in range(dd1):
            z[i,ii]=l[i][ii]
    z=z+z.T
    return z
        

def Neighbor_joining(table, labels):
    # Until all labels have been joined...   
    final_distance=[]
    table=table.tolist()
    dd=len(labels)
    for i in range(len(table)):
        table[i]=[mm for mm in table[i] if mm!=0]
    while len(labels) > 2:    
        table_matrix=list_to_matrix(table)
#        print(table_matrix)
        Qmatrix=Q_matrix(table_matrix)
#        print(Qmatrix)
        Qmatrix=Qmatrix.tolist()
        for i in range(len(Qmatrix)):
            Qmatrix[i]=[mm for mm in Qmatrix[i] if mm!=0]          
        # Locate lowest cell in the table
        x, y = lowest_cell(Qmatrix)
#        print(table_matrix)
        distance=table_matrix[x,y]/2+0.5*(np.sum(table_matrix,0)[x]-
                             np.sum(table_matrix,1)[y])/(len(table_matrix)-2)
#        print(table_matrix[x,y]-distance)
        text=labels[x]+str(round(distance,2))+'--'+\
        str(round(table_matrix[x,y]-distance,2))+labels[y]
        final_distance.append(text)
        if len(labels)==3:
            vvv=np.array(table_matrix)
            vvv=np.delete(vvv,x,1)
            vvv=np.delete(vvv,y,0)
            ddd=vvv[0,1]-distance
        join_table(table, x, y)
        join_labels(labels, x, y)
#        print(labels,x,y)
    text=labels[x]+str(round(ddd,2))+'--'+labels[y]
    final_distance.append(text)
    table_matrix=list_to_matrix(table)
    Qmatrix=Q_matrix(table_matrix)
    Qmatrix=Qmatrix.tolist()
    for i in range(len(Qmatrix)):
        Qmatrix[i]=[mm for mm in Qmatrix[i] if mm!=0]          
    x, y = lowest_cell(Qmatrix)
    join_table(table, x, y)
    join_labels(labels, x, y)

    
# Return the final label
    return labels[0],final_distance

    

