import numpy as np

A1=np.array([
[1,2,3,4,5,6,7,8,9,0],
[1,2,3,4,5,6,7,8,9,0],
[1,2,3,4,5,6,7,8,9,0],
[1,2,3,4,5,6,7,8,9,0],
[1,2,3,4,5,6,7,8,9,0],
[1,2,3,4,5,6,7,8,9,0],
])

A0=np.array([
[0,2,3,2,0,1,2,3,9,0],
[0,2,3,4,5,6,7,8,9,0],
[0,2,3,4,6,6,7,8,9,0],
[0,2,3,4,5,6,4,6,6,7],
[0,2,3,6,5,6,7,8,9,0],
[0,2,3,4,5,4,0,1,9,0],
])

AR=np.array([
[0,2,2,2,8,0,5,2,2,0],
[1,2,8,8,3,0,4,1,2,6],
[1,2,1,9,2,0,3,6,2,2],
[0,9,2,6,3,0,3,8,2,6],
[0,8,4,2,3,0,2,2,7,9],
[1,5,3,5,4,0,5,2,0,4]
])

whole=np.array(
[[0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0],
[0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0], 
[0, 9, 8, 7, 6, 5, 4, 3, 2, 0, 2, 3, 4, 5, 6, 7, 8, 9, 0],
[0, 9, 8, 7, 6, 5, 4, 3, 2, 0, 2, 3, 4, 5, 6, 7, 8, 9, 0], 
[0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0]] )



whole1=np.array(
[[0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0],
[0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0], 
[0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0],
[0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0], 
[0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0]] )

posOfInt=10
#whole.shape[1]
a=whole[0:,:posOfInt-1]
b=whole[0:, posOfInt:]



rA0 = np.flip(A0, axis=1)
rA1 = np.flip(A1, axis=1)
rAR = np.flip(AR, axis=1)

print ('rA1\n', rA1) 
print ('A1\n', A1 )
print ('--------------') 

print ('rA0\n' , rA0) 
print ('A0\n', A0)
print ('--------------') 

print ('rAR\n', rAR) 
print ('AR\n', AR)
print ('--------------') 




def calc_EHH(haplotypes):
    num_haplotypes = haplotypes.shape[0]
    
    EHH = np.zeros(haplotypes.shape[1])
    
    for i in range(haplotypes.shape[1]):
        homozygous_pairs = 0
        for j in range(num_haplotypes):
            for k in range(j+1, num_haplotypes):
                if np.all(haplotypes[j,:i+1] == haplotypes[k,:i+1]):
                    homozygous_pairs += 1
        
        EHH[i] = round(homozygous_pairs / (num_haplotypes * (num_haplotypes - 1) / 2), 3) 
        
    return EHH

#print('A1', np.flip(calc_EHH(rA1)), calc_EHH(A1)) 
#print('A0', np.flip(calc_EHH(rA0)), calc_EHH(A0)) 
#print('AR', np.flip(calc_EHH(rAR)) ,  calc_EHH(AR)) 


print('A1', np.concatenate((np.flip(calc_EHH(rA1)), calc_EHH(A1)))) 
print('A0', np.concatenate((np.flip(calc_EHH(rA0)), calc_EHH(A0)))) 
print('AR', np.concatenate((np.flip(calc_EHH(rAR)), calc_EHH(AR)))) 
