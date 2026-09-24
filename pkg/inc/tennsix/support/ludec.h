//
// Created by Lucius Schoenbaum on 1/23/22.
//

#pragma once

template <typename PP, typename T, size_t size, size_t cosize, size_t prosize>
class Super;

template <typename PP, typename T, size_t size_, size_t cosize_ = size_, size_t prosize_ = 1>
class SuperLUdec {
    friend class TSuper;
public: // fields:

    // static, constant expression-initialized, const:
    static constexpr size_t size = size_;
    static constexpr size_t cosize = cosize_;
    static constexpr size_t shortsize = (size < cosize)? size: cosize;
    static constexpr size_t data_size = size * cosize;
    static constexpr size_t dimension = PP::dimension*data_size;

private: // fields:

    // LU decomposition
    // U in the upper triangle,
    // L in the lower triangle,
    // diagonal has diagonal of U,
    // diagonal of L is understood to be 1's (and not stored).
    Super<PP, T, size_, cosize_, prosize_> data;

    // stores the permutation
    int indx[shortsize*PP::size];
    // todo
    T d;
    // todo
    const Super<PP, T, size, cosize, prosize_> & ref;

public: // ctor:

    SuperLUdec(const Super<PP, T, size, cosize, prosize_> & ain):
            data(ain),
            d(T::unit()),
            ref(ain)
    {
        // set the row-pivoting permutation trivially
        // todo shortsize can be replaced with size (later)
        for (int i = 0; i < shortsize*PP::shortsize; ++i) {
            indx[i] = i;
        }
//        init_doolittle(ain);
        init_crout_with_partial_pivoting(ain);
    }

private: // methods:

    // for debugging: print the array representing the permutation P
    std::string indx_str() {
        std::string out = "(" + std::to_string(indx[0]);
        for (int i = 1; i < shortsize*PP::shortsize; ++i) {
            out += ", " + std::to_string(indx[i]);
        }
        return out + ")";
    }

    // LU decomposition algorithm without pivoting.
    // NOT numerically stable,
    // but at least it removes the hassle of pivoting to reveal a cleaner algorithm.
    // ...called Doolittle's method according to N. Higham's textbook.
    void init_doolittle(const Super<PP, T, size, cosize, prosize_> & ain) {
        int i, j, k;
        for (/* column */ j = 0; j < shortsize; ++j) {
            for (/* row */ i = 0; i <= j; ++i) {
                for (k = 0; k < i; ++k) {
                    // subtractive loss of precision is possible.
                    data[{i,j}] -= data[{i,k}] * data[{k,j}];
                }
            }
            for (/* row */ i = j+1; i < shortsize; ++i) {
                for (k = 0; k < j; ++k) {
                    // subtractive loss of precision is possible.
                    data[{i,j}] -= data[{i,k}] * data[{k,j}];
                }
                // division by zero will occur here if any diagonal element is zero.
                data[{i,j}] /= data[{j,j}];
            }
        }
    }// todo fix after reaching into T from Super

    // LU decomposition algorithm with partial pivoting
    // this version is suboptimal in terms of space usage
    // (I think, I haven't checked really carefully)
    // but it is easier to debug than a memory-sipping revision.
    // ...for example, it could be used as a source of truth/reference for a
    // more optimized implementation.
    void init_crout_with_partial_pivoting(const Super<PP, T, size, cosize, prosize_> & ain) {
        T tmp;
        size_t Bj, jj, Bi, ii, Bj2, jj2;
//        Super<PP, T, shortsize, 1> scaling;
//        populate_scaling(scaling);
        for (/* column */ int j = 0; j < shortsize*PP::shortsize; ++j) {
            pivot_step(j);
            Bj = j/PP::shortsize;
            jj = j%PP::shortsize;
            // a step sometimes useful if the matrix is singular:
            if (data[{Bj, Bj}][{jj, jj}] == T::zero())
                data[{Bj, Bj}][{jj, jj}] = T::tiny();
            // go _down_ by looping over i. For each stop, we:
            // - divide by the multiplier, this builds L.
            for (/* row */ int i = j + 1; i < shortsize*PP::shortsize; ++i) {
                Bi = i/PP::shortsize;
                ii = i%PP::shortsize;
                data[{Bi, Bj}][{ii,jj}] /= data[{Bj, Bj}][{jj, jj}];
                tmp = data[{Bi, Bj}][{ii, jj}];
                for (/* column */ int j2 = j + 1; j2 < shortsize*PP::shortsize; ++j2) {
                    Bj2 = j2/PP::shortsize;
                    jj2 = j2%PP::shortsize;
                    data[{Bi, Bj2}][{ii, jj2}] -= tmp * data[{Bj, Bj2}][{jj, jj2}];
                }
            }
        }
//        cout << "after init:\nLU " << data.str();
//        cout << "P " << indx_str() << endl;
    }

    void pivot_step(int currenti) {
        int pivoti = find_pivot_row(currenti);
        if (/* the current row wasn't the pivot row */ currenti != pivoti)
            interchange_rows(pivoti, currenti);
    }

    // check the value at {currenti, currenti} to see whether it's the best pivot.
    size_t find_pivot_row(size_t currenti) {
        size_t out = currenti;
        T tmp;
        T big = T::zero();
        size_t Bi, ii;
        for (size_t i = currenti; i < shortsize*PP::shortsize; ++i) {
            Bi = i/PP::shortsize;
            ii = i%PP::shortsize;
            tmp = data[{Bi, currenti/PP::shortsize}][{ii, currenti%PP::shortsize}].abs();
            if (tmp > big) {
                big = tmp;
                out = i;
            }
        }
        return out;
    }

    void interchange_rows(int i1, int i2) {
        T tmp;
        size_t itmp;
        size_t Bj, jj;
        size_t Bi1 = i1/PP::shortsize;
        size_t ii1 = i1%PP::shortsize;
        size_t Bi2 = i2/PP::shortsize;
        size_t ii2 = i2%PP::shortsize;
        for (size_t j = 0; j < shortsize*PP::shortsize; ++j) {
            Bj = j/PP::shortsize;
            jj = j%PP::shortsize;
            tmp = data[{Bi1, Bj}][{ii1, jj}];
            data[{Bi1, Bj}][{ii1, jj}] = data[{Bi2, Bj}][{ii2, jj}];
            data[{Bi2,Bj}][{ii2, jj}] = tmp;
        }
        itmp = indx[i1];
        indx[i1] = indx[i2];
        indx[i2] = itmp;
        // bookkeeping for determinants
        d = T::zero() - d;
    }

public: // API:


//    // solve Ax = b; answer is returned in b
//    // tmp: debug only...
//    void solve_tmp(Super<PP, T, size, 1> & b) {
//        int i, ip, j;
//        int i2 = 0;
//        PP sum;
//        for (i = 0; i < shortsize*PP::dimension; ++i) {
//            ip = indx[i];
//            sum = b[ip];
//            b[ip] = b[i];
//            if (i2 != 0) {
//                for (j = i2 - 1; j < i; ++j) {
//                    sum -= data[{i, j}] * b[j];
//                }
//            }
//            else if (sum != PP::zero())
//                i2 = i + 1;
//            b[i] = sum;
//        }
//        for (i = shortsize*PP::dimension-1; i >= 0; --i) {
//            sum = b[i];
//            for (j = i + 1; j < shortsize*PP::dimension; ++j) {
//                sum -= data[{i, j}] * b[j];
//            }
//            b[i] = sum/data[{i,i}];
//        }
//    }

    // solve Ax = b; answer is returned in b
    // ...key routine.
    void solve(Super<PP, T, size, 1, 1> & b) {
        size_t i, j, J;
        size_t Bi, ii, Bj, jj;
        size_t i2 = 0;
        T sum;
        Super<PP, T, size, 1, 1> tmp;
        // LU = PA, so apply P to b and solve LU = P b
        for (i = 0; i < shortsize*PP::shortsize; ++i) {
            Bi = i/PP::shortsize;
            ii = i%PP::shortsize;
            for (jj = 0; jj < PP::cosize; ++jj) {
                tmp[Bi][{ii, jj}] = b[indx[i]/PP::shortsize][{indx[i]%PP::shortsize, jj}];
            }
        }
        b = tmp;
        // PP may have cosize > 1, so loop over J (Jth column in PP)
        for (J = 0; J < PP::cosize; ++J) {
            for (i = 0; i < shortsize*PP::shortsize; ++i) {
                Bi = i/PP::shortsize;
                ii = i%PP::shortsize;
                sum = b[Bi][{ii, J}];
                if (i2 != 0)
                    for (j = i2 - 1; j < i; ++j) {
                        Bj = j / PP::shortsize;
                        jj = j % PP::shortsize;
                        sum -= data[{Bi, Bj}][{ii,jj}] * b[Bj][{jj, J}];
                    }
                else if (sum != T::zero())
                    i2 = i + 1;
                b[Bi][{ii, J}] = sum;
            }
            for (i = shortsize*PP::shortsize-1; i != -1; --i) {
                Bi = i/PP::shortsize;
                ii = i%PP::shortsize;
                sum = b[Bi][{ii, J}];
                for (j = i + 1; j < shortsize*PP::shortsize; ++j) {
                    Bj = j/PP::shortsize;
                    jj = j%PP::shortsize;
                    sum -= data[{Bi, Bj}][{ii, jj}] * b[Bj][{jj, J}];
                }
                b[Bi][{ii, J}] = sum/data[{Bi,Bi}][{ii,ii}];
            }
        }
    }

    // solve Ax = b, preserving b
    void solve(Super<PP, T, size, 1, 1> & x, const Super<PP, T, size, 1, 1> b) {
        x = b;
        solve(x);
    }

    // solve Ax = B; answer is returned in B
    // proceed by straightforward extension of preceding solve() routine.
    template<size_t B_cosize>
    void solve(Super<PP, T, size, B_cosize, 1> & B) {
        size_t i, j;
        Super<PP, T, size, 1, 1> tmp;
        for (j = 0; j < B_cosize; ++j) {
            for (i = 0; i < size; ++i) {
                tmp[i] = B[{i,j}];
            }
            solve(tmp);
            for (i = 0; i < size; ++i) {
                B[{i,j}] = tmp[i];
            }
        }
    }

    // find the inverse of the subject matrix.
    Super<PP, T, size, cosize, 1> inverse() {
        Super<PP, T, size, cosize, 1> out = Super<PP, T, size, cosize, 1>::unit();
        solve<cosize>(out);
        return out;
    }

    // get the "L" piece of the LU,
    // mainly for testing/verification/sanity check purposes.
    // Note PA = LU, P permutation
    Super<PP, T, shortsize, shortsize, 1> L() {
        Super<PP, T, shortsize, shortsize, 1> out = data;
        for (size_t i = 0; i < shortsize; ++i) {
            for (size_t j = i+1; j < shortsize; ++j) {
                out[{i,j}] = T::zero();
            }
            out[{i, i}] = T::unit();
        }
        return out;
    }

    // get the "U" piece of the LU
    // mainly for testing/verification/sanity check purposes.
    // Note PA = LU, P permutation
    Super<PP, T, shortsize, shortsize, 1> U() {
        Super<PP, T, shortsize, shortsize, 1> out = data;
        for (size_t i = 1; i < shortsize; ++i) {
            for (size_t j = 0; j < i; ++j) {
                out[{i,j}] = T::zero();
            }
        }
        return out;
    }

};


