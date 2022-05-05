#ifndef Split_H
#define Split_H

#include "Cell.h"

template <class CellType1, class CellType2>
inline void CalcSplit(
    bool& split1, bool& split2, const CellType1& c1, 
    const CellType2& c2, const double d, const double b)
{
    // This function determines whether either input cell needs to be
    // split.  It is written as a template so that the second cell
    // can be either a Cell or a PairCell.  (The same rules apply.)
    // If you already know that c1 needs to be split, then split1 can
    // be input as true, and we only check c2.  (and vice versa)
    // In normal operation, both are input as false, and we check
    // whether they need to be split.  
    //
    // If (s1+s2)/d > b, then we need to split one or both.
    //
    // If s1 > b*d, then it doesn't matter what s2 is -- we should
    // definitely split c1.  
    // Likewise if s2 > b*d
    //
    // If neither of these conditions is true, then we test
    // to see if s1+s2 > b*d
    // If we're close to the threshold, it will generally be quicker
    // to only split the larger one.  But if both s1 and s2 are 
    // reasonably large (compared to b/d) then we will probably end 
    // up needing to split both, so go ahead and do so now.
    // This small vs. large test is quantified by the parameter
    // splitfactor.  I varied split factor with the 2-point
    // correlation code until it ran fastest.  The result is 
    // given above.  I don't know if this value is also best for
    // 3 point uses, but it's probably reasonably close.


    const double splitfactor = 0.585;
    // The split factor helps determine whether to split
    // both cells or just one when the factor (s1+s2)/d 
    // is too large.  
    // If s1+s2 > f*d*b then split both.
    // Otherwise just split the large Cell.
    // The value of f was determined empirically by seeing 
    // when the code ran fastest.  This may be specific
    // to the data I was testing it on, but I would guess
    // that this value is close enough to optimal for most
    // datasets.

    const double s1 = c1.getSize();
    const double s2 = c2.getSize();
    if (s2 > s1) {
        // Make s1 the larger value.
        CalcSplit(split2,split1,c2,c1,d,b);
    } else if (s1 > 2.*s2) {
        // If one cell is more than 2x the size of the other, only split that one.
        if (!split1) {
            split1 = s1 > b*d;
        } // else do nothing.
    } else {
        if (!split1) {
            const double maxs = b*d;
            split1 = s1 > maxs;
            split2 = s2 > maxs;
            if (!split1 && !split2) {
                // If both are small, need to check the sum.
                double sum = c1.getSize()+c2.getSize();
                if (sum > maxs) {
                    double modmax = splitfactor*maxs;
                    split1 = true;
                    split2 = (s2 > modmax);
                }
            }
        } else if (split1) {
            const double s2 = c2.getSize();
            const double maxs = b*d;
            split2 = s2 > maxs;
        } else if (split2) {
            const double s1 = c1.getSize();
            const double maxs = b*d;
            split1 = s1 > maxs;
        } // else do nothing.
    }
}

template <class CellType1, class CellType2>
inline void CalcSplitSq(
    bool& split1, bool& split2, const CellType1& c1, 
    const CellType2& c2, const double dsq, const double bsq)
{
    // The same as above, but when we know the distance squared rather
    // than just the distance.  We get some speed up by saving the 
    // square roots in some parts of the code.
    const double splitfactorsq = 0.3422;
    const double s1sq = c1.getSizeSq();
    const double s2sq = c2.getSizeSq();
    if (s2sq > s1sq) {
        // Make s1 the larger value.
        CalcSplitSq(split2,split1,c2,c1,dsq,bsq);
    } else if (s1sq > 4.*s2sq) {
        // If one cell is more than 2x the size of the other, only split that one.
        if (!split1) {
            split1 = s1sq > bsq*dsq;
        } // else do nothing.
    } else {
        if (!split1) {
            const double maxssq = bsq*dsq;
            split1 = s1sq > maxssq;
            split2 = s2sq > maxssq;
            if (!split1 && !split2) {
                // If both are small, need to check the sum.
                double sumsq = c1.getSize()+c2.getSize();
                sumsq *= sumsq;
                if (sumsq > maxssq) {
                    double modmax = splitfactorsq*maxssq;
                    split1 = true;
                    split2 = (s2sq > modmax);
                }
            }
        } else if (split1) {
            const double s2sq = c2.getSizeSq();
            const double maxssq = bsq*dsq;
            split2 = s2sq > maxssq;
        } else if (split2) {
            const double s1sq = c1.getSizeSq();
            const double maxssq = bsq*dsq;
            split1 = s1sq > maxssq;
        } // else do nothing.
    }
}

template <class CellType1, class CellType2>
inline bool NoSplit(const CellType1& c1, const CellType2& c2, const double d, const double b)
{
    static const double altb = b/(1.-b);
    // A debugging routine.  Usually of the form:
    // XAssert(NoSplit(c1,c2,d,b))
    // This just asserts that the cells obey the non-splitting eqn:
    // (s1 + s2)/d < b
    // Technically we use altb = b/(1-b) which = b for small b.
    if (c1.getSize() + c2.getSize() < d*altb+0.0001) {
        return true;
    } else {
        std::cerr<<c1.getSize()<<" + "<<c2.getSize()<<" > "<<
            d<<" * "<<altb<<std::endl;
        return false;
    }
}

template <class CellType1, class CellType2, class CellType3>
inline bool Check(
    const CellType1& c1, const CellType2& c2, const CellType3& c3,
    const double d1, const double d2, const double d3)
{
    // Checks that d1,d2,d3 are correct for the three Cells given.
    // Used as a debugging check.
    bool ok=true;
    if (Dist(c3.getData()->pos,c2.getData()->pos)-d1 > 0.0001) 
    { std::cerr<<"d1\n"; ok = false; }
    if (Dist(c1.getData()->pos,c3.getData()->pos)-d2 > 0.0001) 
    { std::cerr<<"d2\n"; ok = false; }
    if (Dist(c2.getData()->pos,c1.getData()->pos)-d3 > 0.0001) 
    { std::cerr<<"d3\n"; ok = false; }
    if (d1 > d2+d3+0.0001) { std::cerr<<"sum d1\n"; ok = false; }
    if (d2 > d1+d3+0.0001) { std::cerr<<"sum d2\n"; ok = false; }
    if (d3 > d1+d2+0.0001) { std::cerr<<"sum d3\n"; ok = false; }
    return ok;
}

#endif
