import argparse
import numpy as np 
import gzip
import sys

def calc_EHH(haplotypes):
    num_haplotypes = haplotypes.shape[0]
    EHH = np.zeros(haplotypes.shape[1])

    for i in range(haplotypes.shape[1]):
        homozygous_pairs = 0
        for j in range(num_haplotypes):
            for k in range(j+1, num_haplotypes):
                if np.all(haplotypes[j,:i+1] == haplotypes[k,:i+1]):
                    homozygous_pairs += 1

        if num_haplotypes < 2:
            return np.full(haplotypes.shape[1], fill_value=500)
        else: 
            EHH[i] = round(homozygous_pairs / (num_haplotypes * (num_haplotypes - 1) / 2), 3) 
    return EHH


def main():
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser()

    # Add command-line arguments
    parser.add_argument("-i", help="Path to the input file, matrix of haplotypes, no header")
    parser.add_argument("-p",  type=int, help="Position of the test SNP in the haplotype window")
    parser.add_argument("-w",  type=int, help="Window size")
    #parser.add_argument("-rwn",  type=str, help="file with row names ")
    parser.add_argument("-refpos",  type=int, help="reference position ")
    parser.add_argument("-o",  type=str, help="outputfile  ")
    # Parse the command-line arguments
    args = parser.parse_args()

    #actions 
    sys.stdout = open(args.o, "w") 

    whole = np.loadtxt(args.i)  # load haplotype matrix 
    testSNP=args.p-1     # 0-based position of test SNP in the window 
    window_size=args.w
    window_name=1 
    colstart=0
    colend=colstart+window_size
    for i in range (whole.shape[1]): 
        #print ('check',  window_name)
        window=whole[0:, colstart:colend]
        window[window != 0] = 1 # replace non zero  non one with one  (see exampple DARC ) 
        #print ('nonzero')
        if window.shape[1]==0: continue   # skip empty windows 
        testAlleles= np.unique(window[:, testSNP])
        refall=window[args.refpos-1 ,testSNP]
        #if len(testAlleles)==1: continue  # skip monomorphic sites 
        for al in testAlleles: 
            #print (al, 'allele')
            sub=window[window[:, testSNP] == al]
            a=sub[0:,:testSNP]  # half haplotypes: start to test SNP --> exclude test SNP itself 
            b=sub[0:, testSNP+1:]   # second half of  halpotypes  test SNP to the end --> exclude test SNP itself 
            rb=np.flip(b, axis=1)   # reverse of second half to calculate EHH 
            ehhvec=np.concatenate((np.flip(calc_EHH(rb)), calc_EHH(b))).tolist()
            # Calculate integral 
            area = np.cumsum(ehhvec)[-1]
            #print( window_name, colstart, colend, al, integral, " ".join(map (str, ehhvec)))  # flip again teh revers to rpvide results in order of position 
            typeal='REF' if  al==refall else  'ALT'
            print (window_name, colstart, colend, al ,typeal,  area , flush=True)
        colstart=colend
        colend=colstart+window_size
        window_name+=1 
if __name__ == "__main__":
    main()
