#include <iostream>
#include <vector>
#include <utility>      // std::pair, std::make_pair
#include <math.h>
#include <stdint.h>
#include <chrono>
#include <map>

// Recommended compilation:
//    g++ -Wall -O4 foggy.cpp -o executableName
// -O4 for speed , -Wall for warnings  

inline std::chrono::high_resolution_clock::time_point tic()
{ return std::chrono::high_resolution_clock::now(); }
inline double toc( const std::chrono::high_resolution_clock::time_point& t2)
{ return std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-t2).count(); }

// the more efficient version of the python code
uint64_t kfoggyPeriod(size_t k){

	// seq will represent the current state of the 
	// sequence, along with the beginning and finish 
	// indices beg and fin
	std::vector<size_t> seq(k+1,0);
	std::vector<size_t> counting(int(k/2) + 2,0);

	// beginning with the canonical 1 and nothing else
	seq[0] = 1;
	// counting contains a running tally of how many of 
	// the digit i-1 are in the state currently
	counting[0] = 1;
	for(size_t i = 1; i < k; i++){
		seq[i] = counting[seq[i-1]-1];
		counting[seq[i]-1] += 1;
	}
	size_t beg = 0;
	size_t fin = k-1;

	// reference state to detect repeated states
	std::vector<size_t> seqRef(seq.begin(),seq.end()-1);

	// p is a counter
	uint64_t p = 0;

	// p_max is the number of steps to take before resetting
	// the reference state seqRef (and doubling p_max threshold)
	uint64_t p_max = 64;

	bool flag = false;
	while(!flag){

		size_t begPrev = beg;
		size_t finPrev = fin;

		// the beginning and finish markers are iterated 1 forward
		// if the go past the end of the alloted k+1 spaces, they 
		// wrap to zero. This means that seq is not always in 
		// order! It must be read from beg to fin.
		if(beg == k){
			beg = 0;
		}else{
			beg++;
		}

		if(fin == k){
			fin = 0;
		}else{
			fin++;
		}

		// the new element is the number of instances of the 
		// previous element in the state
		seq[fin] = counting[seq[finPrev]-1];
		// Add the new element to the counting vector
		counting[seq[fin]-1]++;
		// Decrement the "lost" element from the start of the state
		// this is why the container has to be of size k+1
		counting[seq[begPrev]-1] += -1;

		// Check if the sequence has repeated
		// continue checking elements from beg to fin against seqRef
		// as long as the elements are lining up
		flag = true;
		for(size_t i = 0; (i < k) && flag; i++){
			if(beg+i > k){
				flag = (seqRef[i] == seq[beg+i-k-1]); 
			}else{
				flag = (seqRef[i] == seq[beg+i]);
			} 
		}	

		// counter
		p += 1;

		// check if it's time to reset the counter and threshold
		// unless we just caught a repitition of course
		if((p > p_max) && !flag){

			// reset seqRef to be the current seq
			for(size_t i = 0; i < k; i++){
				if(beg+i > k){
					seqRef[i] = seq[beg+i-k-1]; 
				}else{
					seqRef[i] = seq[beg+i];
				} 
			}

			// reset counter, and double threshold
			p = 0;
			p_max *= 2;
		}
	}

	return p;
}

// a comparison function for a seqRef and seq w/ beg and fin indices, lexicographic comparison
bool isGreater(const std::vector<size_t>& seqRef, const std::vector<size_t>& seq, size_t beg, size_t fin){
	size_t k = seqRef.size();

	for(size_t i = 0; i < k; i++){
		size_t idx = beg+i;
		if(idx > k){
			idx = beg+i-k-1;
		}
		if(seqRef[i] < seq[idx]){
			return false;
		}else if(seqRef[i] > seq[idx]){
			return true;
		}
	}
	return false;
}

// this code finds the cycles as in kfoggyPeriod, but returns the state which is 
// "smallest" within the cycle in order to uniquely identify each cycle (to find them all)
// this also generalizes to beginning with any initial set of k numbers
std::pair<uint64_t,std::vector<size_t>> findCycle(const std::vector<size_t>& seqStart){
	
	size_t k = seqStart.size();


	// seq will represent the current state of the 
	// sequence, along with the beginning and finish 
	// indices beg and fin
	std::vector<size_t> seq(k+1,0);
	std::vector<size_t> counting(k,0);

	// counting contains a running tally of how many of 
	// the digit i-1 are in the state currently
	for(size_t i = 0; i < k; i++){
		seq[i] = seqStart[i];
		counting[seqStart[i]-1] += 1;
	}
	size_t beg = 0;
	size_t fin = k-1;

	// reference state to detect repeated states
	std::vector<size_t> seqRef(seq.begin(),seq.end()-1);

	// p is a counter
	uint64_t p = 0;

	// p_max is the number of steps to take before resetting
	// the reference state seqRef (and doubling p_max threshold)
	uint64_t p_max = 64;

	// this is the same loop as above, to detect a periodic sequence
	bool flag = false;
	while(!flag){

		size_t begPrev = beg;
		size_t finPrev = fin;

		// the beginning and finish markers are iterated 1 forward
		// if the go past the end of the alloted k+1 spaces, they 
		// wrap to zero. This means that seq is not always in 
		// order! It must be read from beg to fin.
		if(beg == k){
			beg = 0;
		}else{
			beg++;
		}

		if(fin == k){
			fin = 0;
		}else{
			fin++;
		}

		// the new element is the number of instances of the 
		// previous element in the state
		seq[fin] = counting[seq[finPrev]-1];
		// Add the new element to the counting vector
		counting[seq[fin]-1]++;
		// Decrement the "lost" element from the start of the state
		// this is why the container has to be of size k+1
		counting[seq[begPrev]-1] += -1;

		// Check if the sequence has repeated
		// continue checking elements from beg to fin against seqRef
		// as long as the elements are lining up
		flag = true;
		for(size_t i = 0; (i < k) && flag; i++){
			if(beg+i > k){
				flag = (seqRef[i] == seq[beg+i-k-1]); 
			}else{
				flag = (seqRef[i] == seq[beg+i]);
			} 
		}	

		// counter
		p += 1;

		// check if it's time to reset the counter and threshold
		// unless we just caught a repitition of course
		if((p > p_max) && !flag){

			// reset seqRef to be the current seq
			for(size_t i = 0; i < k; i++){
				if(beg+i > k){
					seqRef[i] = seq[beg+i-k-1]; 
				}else{
					seqRef[i] = seq[beg+i];
				} 
			}

			// reset counter, and double threshold
			p = 0;
			p_max *= 2;
		}
	}

	//found the cyclic period, and seqRef is a state within it,
	//now I want to find the "smallest" element of the cycle
	//to represent it (just need a unique specifier of some kind)
	for(uint64_t cycle = 0; cycle < p; cycle++){

		//iterate to the next step in the cycle
		size_t begPrev = beg;
		size_t finPrev = fin;

		if(beg == k){
			beg = 0;
		}else{
			beg++;
		}

		if(fin == k){
			fin = 0;
		}else{
			fin++;
		}

		seq[fin] = counting[seq[finPrev]-1];
		counting[seq[fin]-1]++;
		counting[seq[begPrev]-1] += -1;

		// if the state is a lexicographically smaller state
		// than the reference sequence, replace the reference sequence
		if(isGreater(seqRef,seq,beg,fin)){
			for(size_t i = 0; i < k; i++){
				if(beg+i > k){
					seqRef[i] = seq[beg+i-k-1]; 
				}else{
					seqRef[i] = seq[beg+i];
				} 
			}		
		}
	}

	// return the smallest state in the sequence as well as the period
	return std::pair<uint64_t,std::vector<size_t>>(p,seqRef);
}

std::map<std::vector<size_t>, uint64_t> allCycles(size_t k){

	std::map<std::vector<size_t>, uint64_t> foundCycles;
 
	//largest individual element in a cycle (except for the kkkkkk... cycle)
	double upper = ceil(double(k)/2.0)+1;
	double bound = std::pow(upper,k);

	std::vector<size_t> seed(k,1);
	seed[k-1] = 0; //setting up the first iteration

	for(double iter = 0; iter < bound; iter++){
		size_t i = k-1;
		seed[i]++;
		while(seed[i] > upper){
			seed[i] = 1;
			i = i-1;
			seed[i]++;
		}

		std::pair< uint64_t,std::vector<size_t> > cycle = findCycle(seed);

		if(!foundCycles.count(cycle.second)){
			//sequence is not in map
			foundCycles[cycle.second] = cycle.first;
		}
	}

	foundCycles[std::vector<size_t>(k,k)] = 1;
	return foundCycles;
}

std::vector< size_t > unrollOne(std::vector<size_t> seq){
	// unroll a sequence one step back. Returns a vector
	// of possible preceding elements, if empty, it cannot be unrolled
	size_t k = seq.size();
	// this is the frequency it must have had in the last k elements
	size_t freq = seq[k-1];
	// this is the element 
	size_t el = seq[k-2];
	size_t found = 0;
	std::cout<<"the sequence: ";
	for(size_t i = 0; i < k-1; i++){
		std::cout<<seq[i]<<" , ";
		if(seq[i] == el){
			found++;
		}
	}
	std::cout<<seq[k-1]<<" could have been preceeded by: ";

	if(found == freq){
		// the previous element could be anything that's 
		// NOT el
		std::vector<size_t> everythinBut;
		everythinBut.reserve(k-1);
		for(size_t i = 1; i <= k; i++){
			if(i == el){
				continue;
			}
			std::cout<<i<<" , ";
			everythinBut.push_back(i);
		}
		std::cout<<"\n";
		return everythinBut;
	}else if(freq == found + 1){
		// the previous element has to be el
		std::cout<<el<<std::endl;
		return std::vector<size_t>(1,el);
	}else{
		std::cout<<"nothin\n";
		// could not be unrolled
		return std::vector<size_t>();
	}
}

size_t numberToUnRoll(std::vector<size_t> seq){
	// not a good definition of entropy
	
	size_t k = seq.size();

	// this finds all the ways to unroll a state, by 
	// literally just unrolling it one step at a time
	// and enumerating each possible previous state.
	// The possibilities are passed back and forth 
	// between poss1 and poss2 as they are unrolled.
	std::vector< std::vector<size_t> > poss1;
	std::vector< std::vector<size_t> > poss2;

	poss1.push_back(seq);
	bool into2 = true;
	for(size_t step = 0; step < k; step++){
		if(into2){
			//prepare poss2 to be filled
			poss2.clear();
			// for each possible preceding subsequence in poss1
			for(size_t p = 0; p < poss1.size(); p++){
				// the set of possible preceding elements
				std::vector<size_t> oneStep = unrollOne(poss1[p]);
				for(size_t newEl = 0; newEl < oneStep.size(); newEl++){
					// create a new state for the preceeding element
					std::vector<size_t> newPoss(k);
					newPoss[0] = oneStep[newEl];
					for(size_t i = 0; i < k-1; i++){
						newPoss[i+1] = poss1[p][i];
					}
					// store it in poss2
					poss2.push_back(newPoss);
				}
			}
		}else{
			//prepare poss1 to be filled
			poss1.clear();
			// for each possible preceding subsequence in poss2
			for(size_t p = 0; p < poss2.size(); p++){
				// the set of possible preceding elements
				std::vector<size_t> oneStep = unrollOne(poss2[p]);
				for(size_t newEl = 0; newEl < oneStep.size(); newEl++){
					// create a new state for the preceeding element
					std::vector<size_t> newPoss(k);
					newPoss[0] = oneStep[newEl];
					for(size_t i = 0; i < k-1; i++){
						newPoss[i+1] = poss2[p][i];
					}
					// store it in poss1
					poss1.push_back(newPoss);
				}
			}
		}
		into2 = !into2;
	}

	//////////////////////// TO SEE ALL POSSIBLITIES, UNCOMMENT: //////////////////////
	/*
	std::cout<<"For sequence: ";
	for(size_t i = 0; i < k; i++){
		std::cout<<seq[i]<<" , ";
	}
	std::cout<<"\nThere are ";
	if(into2){
		std::cout<<poss1.size()<<" ways to unroll, the following: \n";
		for(size_t p = 0; p < poss1.size(); p++){
			std::cout<<p<<" : ";
			for(size_t i = 0; i < k; i++){
				std::cout<<poss1[p][i]<<" , ";
			}
			std::cout<<std::endl;
		}
	}else{
		std::cout<<poss2.size()<<" ways to unroll, the following: \n";
		for(size_t p = 0; p < poss2.size(); p++){
			std::cout<<p<<" : ";
			for(size_t i = 0; i < k; i++){
				std::cout<<poss2[p][i]<<" , ";
			}
			std::cout<<std::endl;
		}
	}
	*/
	//////////////////////// END VERBOSE BLOCK ////////////////////////////////////////

	if(into2){
		return poss1.size();
	}else{
		return poss2.size();
	}
}


int main(){

	for(int k = 1; k <= 50; k++){
		std::cout<<"------ "<<k<<"-------\n";
		auto t1 = tic();
		uint64_t period = kfoggyPeriod(k);
		auto t2 = toc(t1);
		if(t2 < 2){
			std::cout<<"K = "<<k<<" has period of : "<<period<<" which took "<<1000*t2<<"ms"<<std::endl;
		}else{
			std::cout<<"K = "<<k<<" has period of : "<<period<<" which took "<<t2<<"s"<<std::endl;
		}
		// std::cout<<"For K = "<<k<<" periods are equal: "<<(period1 == period2)<<" and vector is "<<t2D/t2V<<" times the speed!\n";
	}

	// for(size_t k = 1; k < 10; k++){
	// 	auto t1 = tic();
	// 	std::map<std::vector<size_t>, uint64_t> cycles = allCycles(k);
	// 	auto t2 = toc(t1);
	// 	std::cout<<"#####################\n";
	// 	if(t2 < 2){
	// 		std::cout<<"For K = "<<k<<" found "<<cycles.size()<<" cycles in "<<1000*t2<<"ms\n";
	// 	}else{
	// 		std::cout<<"For K = "<<k<<" found "<<cycles.size()<<" cycles in "<<t2<<" seconds\n";
	// 	}
	// 	for(auto it = cycles.begin(); it != cycles.end(); it++){
	// 		std::cout<<"-------------------\n";
	// 		std::cout<<"Element: ";
	// 		for(size_t i = 0; i < k; i++){
	// 			std::cout<<it->first[i]<<" ";
	// 		}
	// 		std::cout<<"\nOf length: "<<it->second<<std::endl;
	// 	}
	// }


	return 1;
}