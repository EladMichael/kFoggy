from math import ceil
from copy import deepcopy
import matplotlib.pyplot as plt

def kFoggy(k,nMax = 1000000):

    sequence = [0]*nMax;

    # The largest element in the sequence is ceil(k/2)+1, but python counts from zero!
    # counter is a list such that counter[i] is the number of digits i in the previous k digits
    counter = [0]*(ceil(k/2)+2); 

    sequence[0] = 1;
    counter[1] += 1;
    for n in range(1,nMax):
        sequence[n] = counter[sequence[n-1]];
        counter[sequence[n]] += 1;
        counter[sequence[n-k]] += -1; 
    
    return sequence;

def longestRuns(seq,k):
    elements = list(range(1,k+1));
    longestRuns = [0]*(k);
    runLength = 1;
    runEl = seq[0];
    for a in seq[1:]:
        if(a == runEl):
            runLength += 1;
        else:
            longestRuns[runEl-1] = max(longestRuns[runEl-1],runLength);
            runEl = a;
            runLength = 1;

    plt.bar(elements,longestRuns);
    plt.show();
    return;

def longestRuns(seq,k):
    elements = list(range(1,k+1));
    longestRuns = [0]*(k);
    runLength = 1;
    runEl = seq[0];
    for a in seq[1:]:
        if(a == runEl):
            runLength += 1;
        else:
            longestRuns[runEl-1] = max(longestRuns[runEl-1],runLength);
            runEl = a;
            runLength = 1;

    for i in range(1,k):
        if(longestRuns[-i] > 1):
            return;

    return -1;

# def longestRunSweep(kMin,kMax,N):
#     for 

def periodFind(k):
    # Find the period for the kFoggy sequence for given k, starting with 1.
    # periods may be different if you start with a different sequence, worth considering...
    # This algorithm works by using a reference sequence (seqRef), and going 2^p for some p 
    # iterations forward to look for a repitition of seqRef. If none is found, then a new seqRef 
    # is taken at the new point (2^p steps from the previous) and the process is repeated for 2^(p+1) steps
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
    # fixed array of size k+1. But this should get the point accross. The C++ version
    # does that, but no point trying to hyper-optimize the python version eh? 

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

    periods = [3,8,11,25,20,40,9,45,41]
    for i in range(2,11):
        print("Calculated/True period of ",i," Foggy sequnce: ",periodFind(i),"/",periods[i-2]);

def advanceOne(seq):
    nextSeq = seq[1:];
    nextSeq.append(sum([el==seq[-1] for el in seq]));
    return nextSeq;

def buildGraph(k):
    from pyvis.network import Network
    # In order to find node indices which have already been added to the graph
    nodes = {};

    net = Network(height='100%', width='100%', directed=True);
    
    # this will be the state for each node
    seq = [1]*k;

    # Just an offset for the first iteration
    # I like to have the loop updating conditions at top, for some reason!
    seq[-1] = 0;

    # at each iteration, the state will look like
    # a_{k} a_{k-1} a_{k-2} ... a_{1} a_{0}
    # we increment a_{0} by 1, and if it is above 
    # the upperLim (k if you want all states), then
    # it rolls over to a_{1} and goes to 1. If a_{1} 
    # is over upperLim, it rolls over to a_{2}, etc...
    upperLim = k;

    # this for loop works when upperLim == k, but for other
    # termination conditions you could use while(not all([a == upperLim for a in seq])) 
    # or whatever the right syntax is
    for count in range(upperLim**k):

        # generate the new sequence for this step
        # there is no good reason to increment from the back.
        i = 1;
        seq[-i]+=1;
        while(seq[-i] > upperLim):
            seq[-i] = 1;
            i+=1;
            seq[-i] += 1;

        # Either give this seq a new number, 
        # or get its previously assigned number
        if(str(seq) in nodes):
            thisNum = nodes[str(seq)];
        else:
            thisNum = len(nodes); # new key number
            net.add_node(thisNum,label=str(seq)[1:-1]);
            nodes[str(seq)] = thisNum;

        # where does this sequence lead to in the next step?
        nextStep = advanceOne(seq);

        # get nextStep a node number, or get its existing one
        if(str(nextStep) in nodes):
            nextNum = nodes[str(nextStep)];
        else:
            nextNum = len(nodes); # new key number
            net.add_node(nextNum,label=str(nextStep)[1:-1]);
            nodes[str(nextStep)] = nextNum;

        net.add_edge(thisNum,nextNum);

    net.toggle_physics(True);
    net.show_buttons(filter_=['physics'])
    net.show("%ifoggyGraph.html" % k);
    return

k = 500;
seq = kFoggy(k,100000);
longestRuns(seq,k);
