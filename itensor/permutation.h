//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_PERMUTATION_H
#define __ITENSOR_PERMUTATION_H
#include "global.h"

namespace itensor {

//
// Tell where each index will go, 
// if(p.dest(2) == 1) then 2 -> 1, etc.
//
struct Permutation
    {
    using storage = std::vector<long>;
    using size_type = storage::size_type;

    private: 
    ///////////
    storage store_;
    ///////////
    public:

    Permutation(size_type size = 0);

    long
    size() const { return long(store_.size()); }


    void 
    setFromTo(long from, 
              long to);

    long 
    dest(long j) const { return GET(store_,j); }

    //bool 
    //check(int d);

    const storage& 
    store() const { return store_; }

    };

inline Permutation::
Permutation(size_type size) 
    : 
    store_(size)
    { 
    for(size_type n = 0; n < size; ++n)
        store_[n] = n;
    }

void inline Permutation::
setFromTo(long from, 
          long to) 
    { 
    GET(store_,from) = to; 
    }

Permutation inline
inverse(const Permutation& P)
    {
    Permutation inv(P.size());
    for(long n = 0; n < P.size(); ++n) 
        inv.setFromTo(P.dest(n),n);
    return inv;
    }

bool inline
isTrivial(const Permutation& P)
    {
    for(long n = 0; n < P.size(); ++n) 
        {
        if(P.dest(n) != n) return false;
        }
    return true;
    }

template<typename Iterable>
void
permute(const Permutation& P,
        const Iterable& from,
        Iterable& to)
    {
    using size_type = typename Iterable::size_type;
    for(size_type i = 0; i < from.size(); ++i)
        {
        to[P.dest(i)] = from[i];
        }
    }

//bool inline Permutation::
//check(int d)
//	{
//    for(int i = 1; i <= d; i++)
//        {
//        if(ind_[i] > d || ind_[i] < 1) 
//            {
//            std::cerr << "\nbad Permutation level 1\n\n";
//            return false;
//            }
//        }
//
//    for(int i = 1; i <= d; ++i)
//    for(int j = 1; j <= d; ++j)
//        {
//        if(i == j) continue;
//        if(ind_[i] == ind_[j]) 
//            {
//            std::cerr << "\nbad Permutation level 2\n\n";
//            return false;
//            }
//        }
//
//    return true;
//	}

inline 
std::ostream& 
operator<<(std::ostream& s, const Permutation& P)
    {
    for(long i = 0; i < P.size(); ++i) 
        s << "(" << i << "," << P.dest(i) << ")";
    return s;
    }

}; //namespace itensor

#endif
