/*
 * data_structures.h
 *
 *  Created on: Mar 8, 2016
 *      Author: Asher Stern
 */

#ifndef DATA_STRUCTURES_H_
#define DATA_STRUCTURES_H_

#include<list>
#include<map>
#include<stdexcept>
#include<vector>
#include<utility>

namespace hmm
{

template <class T, class U>
class TransitionMatrix
{
public:
	TransitionMatrix():values(){}
	TransitionMatrix(const TransitionMatrix& other):values(other.values){}
	TransitionMatrix(TransitionMatrix&& other):values(std::move(other.values)){}
	virtual ~TransitionMatrix(){}
	
	double operator()(T t, U u) const;
	void set(T t, U u, double value);
private:
	std::map<T, std::map<U, double>> values;
};



template <class H>
class TimedProbability
{
public:
	TimedProbability():values(){}
	virtual ~TimedProbability(){}

	double operator()(H state, int time) const;
	void set(H state, int time, double probability);
private:
	std::map<int, std::map<H, double>> values;
};



// like std::vector, but index starts at 1, not at 0.
template <class T>
class Vector1
{
public:
	using size_type = std::vector<T>::size_type;
	
	Vector1():internal(){}
	Vector1(Vector1&& other):internal(std::move(other.internal)){}
	Vector1(const Vector1& other):internal(other.internal){}
	~Vector1(){}
	
	size_type size() const {return internal.size();}
	void reserve (size_type n) {internal.reserve(n);}
	T operator[](size_type index);
	const T operator[](size_type index) const;
	void push_back(T t);
	
private:
	std::vector<T> internal;
};





//////////////////// IMPLEMENTATION ////////////////////

template <class T, class U>
double TransitionMatrix<T,U>::operator()(T t, U u) const
{
	if (values.count(t)>0)
	{
		if (values.at(t).count(u)>0)
		{
			return values.at(t).at(u);
		}
	}
	throw std::invalid_argument("value does not exist for the given keys.");
}

template <class T, class U>
void TransitionMatrix<T,U>::set(T t, U u, double value)
{
	if (values.count(t)<=0)
	{
		values.insert(std::pair<T, std::map<U, double>>(t, std::map<U, double>()));
	}
	std::map<U, double>& mapForT = values.at(t);
	if (mapForT.count(u)>0)
	{
		mapForT.at(u)=value;
	}
	else
	{
		mapForT.insert(std::pair<U, double>(u,value));
	}
}



template <class H>
double TimedProbability<H>::operator()(H state, int time) const
{
	if (values.count(time)>0)
	{
		const std::map<H, double> timed = values.at(time);
		if (timed.count(state)>0)
		{
			return timed.at(state);
		}
	}
	throw std::invalid_argument("bad state/time");
}

template <class H>

void TimedProbability<H>::set(H state, int time, double probability)
{
	if (values.count(time)<=0)
	{
		values.insert(std::pair<int, std::map<H, double>>(time,std::map<H, double>()));
	}
	std::map<H, double> timed = values.at(time);
	if (timed.count(state)<=0)
	{
		timed.insert(std::pair<H, double>(state, probability));
	}
	else
	{
		timed.at(state)=probability;
	}
}


template <class T>
T Vector1<T>::operator[](size_type index)
{
	if (index<1) throw std::invalid_argument("index<1");
	return internal[index-1];
	
}

template <class T>
const T Vector1<T>::operator[](size_type index) const
{
	if (index<1) throw std::invalid_argument("index<1");
	return internal[index-1];
}

template <class T>
void Vector1<T>::push_back(T t)
{
	internal.push_back(t);
}


}



#endif /* DATA_STRUCTURES_H_ */
