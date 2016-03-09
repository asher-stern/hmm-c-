/*
 * Hmm.h
 *
 *  Created on: Mar 8, 2016
 *      Author: Asher Stern
 */

#ifndef HMM_H_
#define HMM_H_

#include"data_structures.h"

#include<list>
#include<vector>
#include<stdexcept>
#include<map>

namespace hmm
{



template <class H, class O>
class Hmm
{
public:
	Hmm(H _nullHiddenState, std::list<H> _hiddenStates):nullHiddenState(_nullHiddenState),hiddenStates(_hiddenStates){}
	Hmm(H _nullHiddenState, std::list<H> _hiddenStates, TransitionMatrix<H,H> _hiddenTransitions, TransitionMatrix<H,O> _hiddenToObservedProbabilities):
		nullHiddenState(_nullHiddenState),hiddenStates(_hiddenStates),hiddenTransitions(_hiddenTransitions),hiddenToObservedProbabilities(_hiddenToObservedProbabilities){}
	virtual ~Hmm(){}
	
	double computeProbabilityOfObserved(Vector1<O> observed) const;
	Vector1<H> computeMostLikelyHiddenSequence(Vector1<O> observed) const; // Viterbi algorithm

private:
	TimedProbability<H> computeForward(Vector1<O> observed) const;
	TimedProbability<H> computeBackward(Vector1<O> observed) const;

	H nullHiddenState;
	std::list<H> hiddenStates; // must not include the null state
	TransitionMatrix<H,H> hiddenTransitions; // A_i,j
	TransitionMatrix<H,O> hiddenToObservedProbabilities; // B_i,x
};



template <class H, class O>
TimedProbability<H> Hmm<H,O>::computeForward(Vector1<O> observed) const
{
	TimedProbability<H> ret;
	ret.set(nullHiddenState,0,1.0);
	if (observed.size()==0) return ret;
	
	O firstObserved = observed[1];
	for (const H state : hiddenStates)
	{
		// Probability that at time 1 we are in state "state", given the observations up through observed[1].
		ret.set(state,1,hiddenTransitions(nullHiddenState,state)*hiddenToObservedProbabilities(state,firstObserved));
	}
	
	const Vector1<O>::size_type length = observed.size();
	for (Vector1<O>::size_type index=2; index<=length; ++index) // <= this is Vector1
	{
		for (const H state : hiddenStates)
		{
			double probability = 0.0;
			for (const H previousState : hiddenStates)
			{
				probability += ret(previousState,index-1)*hiddenTransitions(previousState,state)*hiddenToObservedProbabilities(state,observed[index]);
			}
			ret.set(state,index,probability);
		}
	}
	
	return ret;
}

template <class H, class O>
TimedProbability<H> Hmm<H,O>::computeBackward(Vector1<O> observed) const
{
	TimedProbability<H> ret;
	if (observed.size()==0)
	{
		ret.set(nullHiddenState,0,1.0);
		return ret;
	}
	
	const Vector1<O>::size_type length = observed.size();
	for (const H state : hiddenStates)
	{
		ret.set(state,length,1.0);
	}
	
	for (Vector1<O>::size_type index = (length-1); index>0; --index)
	{
		for (const H state : hiddenStates)
		{
			double probability = 0;
			for (const H nextState : hiddenStates)
			{
				probability += hiddenTransitions(state, nextState)*hiddenToObservedProbabilities(nextState,observed[index+1])*ret(nextState,index+1);
			}
			ret.set(state,index,probability);
		}
	}
	
	double probability = 0;
	for (const H nextState : hiddenStates)
	{
		probability += hiddenTransitions(nullHiddenState, nextState)*hiddenToObservedProbabilities(nextState,observed[1])*ret(nextState,1);
	}
	ret.set(nullHiddenState,0,probability);
}

template <class H, class O>
double Hmm<H,O>::computeProbabilityOfObserved(Vector1<O> observed) const
{
	if (observed.size()==0) return 1.0;
	TimedProbability<H> forward = computeForward(observed);
	const Vector1<O>::size_type length = observed.size();
	
	double ret = 0.0;
	bool firstIteration = true;
	for (H state : hiddenStates)
	{
		double probability = forward(state,length);
		if (firstIteration)
		{
			ret = probability;
			firstIteration=false;
		}
		else 
		{
			if (ret<probability)
			{
				ret = probability;
			}
		}
	}
	return ret;
}

template <class H, class O>
Vector1<H> Hmm<H,O>::computeMostLikelyHiddenSequence(Vector1<O> observed) const // Viterbi algorithm
{
	using namespace std;
	Vector1<map<H,H>> backtracking;
	backtracking.reserve(observed.size());
	
	TimedProbability<H> maxProbability;
	maxProbability.set(nullHiddenState,0,1.0);
	
	const Vector1<O>::size_type length = observed.size();
	list<H> nullList;
	nullList.push_back(nullHiddenState);
	
	for (Vector1<O>::size_type index = 1; index<=length; ++index)
	{
		for (const H state : hiddenStates)
		{
			list<H>* previousList;
			if (index==1) {previousList = &nullList;}
			else {previousList = &hiddenStates;}
			
			H* maxPrevious = nullptr;
			double maxProbabilityByPreviousState = 0.0;
			bool firstIteration = true;
			for (const H previousState : (*previousList))
			{
				double probability = TransitionMatrix(previousState,state)*hiddenToObservedProbabilities(state,observed[index])*maxProbability(previousState,index-1);
				
				bool itIsHigherThanOthers = false;
				if (firstIteration)
				{
					itIsHigherThanOthers = true;
					firstIteration = false;
				}
				else
				{
					if (probability>maxProbabilityByPreviousState)
					{
						itIsHigherThanOthers = true;
					}
				}
				
				if (itIsHigherThanOthers)
				{
					maxProbabilityByPreviousState = probability;
					maxPrevious = &previousState;
				}
			} // end for each previous state
			maxProbability.set(state,index,maxProbabilityByPreviousState);
			backtracking[index].insert(pair<H,H>(state, (*maxPrevious) ));
		}
	}
	
	Vector1<H> ret;
	ret.reserve(observed.size());
	H* bestLastState = nullptr;
	double maxLastProbability = 0.0;
	bool firstIteration = true;
	for (const H state : hiddenStates)
	{
		double probability = maxProbability(state, length);
		bool itIsBestSoFar = false;
		if (firstIteration)
		{
			itIsBestSoFar = true;
			firstIteration = false;
		}
		else
		{
			if (probability>maxLastProbability)
			{
				itIsBestSoFar = true;
			}
		}
		if (itIsBestSoFar)
		{
			maxLastProbability = probability;
			bestLastState = &state;
		}
	}
	
	ret[length]=(*bestLastState);
	for (Vector1<O>::size_type index = length-1; index>=1; --index)
	{
		ret[index] = backtracking[index+1].at(ret[index+1]);
	}
	
	return ret;
}

}

#endif /* HMM_H_ */
