import numpy as np

from math import ceil
from copy import deepcopy

def kFoggy(k,nMax = 1000000):

    seq = [0]*nMax;

    # The largest element in the sequence is ceil(k/2)+1, but python counts from zero!
    # counter is a list such that counter[i] is the number of digits i in the previous k digits
    counter = [0]*(ceil(k/2)+2); 

    seq[0] = 1;
    counter[1] += 1;
    for n in range(1,nMax):
        seq[n] = counter[seq[n-1]];
        counter[seq[n]] += 1;
        counter[seq[n-k]] += -1; 
    
    return seq;

def periodFind(k):
    # Find the period for the kFoggy sequence for given k, starting with 1.
    # periods may be different if you start with a different sequence, worth considering...
    # This algorithm works by using a reference sequence (seqRef), and going N iterations 
    # forward to look for a repitition of seqRef. If none is found, then a new seqRef 
    # is taken at the new point and N is doubled

    seq = [];
    counter = [0] * (ceil(k/2)+2);
    seq.append(1);
    counter[1] = 1;

    for i in range(2,k+1):
        seq.append(counter[seq[-1]]);
        counter[seq[-1]] += 1;

    seqRef = deepcopy(seq);
    # Iteration counter
    n = 0;

    # this is the cutoff to select a new seqRef, it will be doubled every time it is invoked. 
    # I start it at 64 to avoid the silly small powers of 2
    N = 64;

    loopedFlag = False;

    # This can be made *much* more efficient by not actually popping and appending
    # at each iteration, but using start and end indices and wrapping around in a 
    # fixed array of size k+1. But this should get the point accross.

    while(not loopedFlag):
        n += 1;
        # add the next element to the sequence
        seq.append(counter[seq[-1]])
        # add 1 to the counter for the new elment
        counter[seq[-1]] += 1;
        # Remove the term which is now k+1 terms behind, and decrease the counter
        counter[seq.pop(0)] += -1;

        if(n > N):
            # Reset the reference sequence
            seqRef = deepcopy(seq);
            # Reset the iteration counter
            n = 0;
            # double the cutoff
            N = 2*N;
            if(N > 1e7):
                print("You wanna be doing more than ten million elements?")
                print("There are some complexity improvements which should be implemented first...")
                return -1;

            # obviously the reference and sequence are equal here, cause they 
            # were just copied, so set the flag to false and keep searching
            loopedFlag = False;
        else:
            loopedFlag = (seq == seqRef);

    return n;

def unitTest():
    print("3Foggy seq to 20 places: ",kFoggy(3,20))
    correct = [1,1,2,1,2,2,2,3,1,1,2,1,2,2,2,3,1,1,2,1];
    print("Compared to hardcoded answer (true is correct): ",kFoggy(3,20) == correct);

    print("8Foggy seq to 50 places: ",kFoggy(8,50))
    correct = [1, 1, 2, 1, 3, 1, 4, 1, 5, 1, 4, 2, 1, 4, 3, 1, 3, 2, 2, 3, 3, 4, 1, 2, 3, 3, 4, 2, 2, 3, 3, 4, 2, 3, 3, 4, 2, 2, 3, 3, 4, 2, 3, 3, 4, 2, 2, 3, 3, 4];
    print("Compared to hardcoded answer (true is correct): ",kFoggy(8,50) == correct);
    # This should be [1, 1, 2, 1, 3, 1, 4, 1, 5, 1, 4, 2, 1, 4, 3, 1, 3, 2, 2, 3, 3, 4, 1, 2, 3, 3, 4, 2, 2, 3, 3, 4, 2, 3, 3, 4, 2, 2, 3, 3, 4, 2, 3, 3, 4, 2, 2, 3, 3, 4]
    periods = [3,8,11,25,20,40,9,45,41]
    for i in range(2,11):
        print("Calculated/True period of ",i," Foggy sequnce: ",periodFind(i),"/",periods[i-2]);

if __name__ == '__main__':
    unitTest()