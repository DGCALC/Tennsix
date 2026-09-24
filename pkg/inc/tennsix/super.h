//
// Created by Lucius Schoenbaum on 11/6/21.
//

// Super space
//
// superspace of vectors over a vector space,
// with l2 norm (the usual cartesian norm).
//
// supports operators: +, -,
// scalar multiplication *, /
// dot product %
// as well as norm operation x^1, x^2, x^3, etc.
//
// Also has: zero(), dimension, str()

#pragma once

#include <DGCALC.h>
#include <types/viz.h>
#include <tennsix/support/ludec.h>

template <typename PP, typename T_, size_t size_, size_t cosize_, size_t prosize_ = 1>
class Super {
    friend class TSuper;

public: // types:

    using T = T_;

public: // fields:

    // static, constant expression-initialized, const:
    static constexpr size_t size = size_;
    static constexpr size_t cosize = cosize_;
    static constexpr size_t prosize = prosize_;
    // todo check shortsize usage vs. presence of prosize of 1
    static constexpr size_t shortsize = (size < cosize)? size: cosize;
    static constexpr size_t data_size = size * cosize * prosize;
    static constexpr size_t dimension = PP::dimension*data_size;

protected: // fields:

    PP data[data_size];

public: // ctor:

    explicit Super(const PP & val = PP::zero())
    {
        for (int i = 0; i < data_size; ++i) {
            data[i] = val;
        }
    }

    Super(std::initializer_list<PP> lst)
    {
        if (lst.size() > data_size) assert(0==1);
        size_t i = 0;
        for (auto j: lst) {
            data[i++] = j;
        }
        // fill in any remaining cells
        while (i < data_size) {
            data[i++] = PP::zero();
        }
    }

public: // copy:

    Super(const Super & copyin) {
        if (this == &copyin) return;
        for (size_t i = 0; i < data_size; ++i) {
            data[i] = copyin[i];
        }
    }

    Super & operator=(const Super & rhs) {
        if (this == &rhs) return *this;
        for (size_t i = 0; i < data_size; ++i) {
            data[i] = rhs[i];
        }
        return *this;
    }

public: // dtor:

    ~Super() = default;

public: // artifacts:

    // Equality derived from base type, so you may wish to avoid it
    bool operator==(const Super & rhs) const {
        if (this == &rhs) return true;
        bool out = true;
        for (size_t i = 0; i < data_size; ++i) {
            if (data[i] != rhs.data[i]) {
                out = false;
                break;
            }
        }
        return out;
    }

    bool operator!=(const Super & rhs) const {
        return !(*this == rhs);
    }

    // dictionary order
    bool operator<(const Super & rhs) const {
        int i = 0;
        do {
            if (data[i] < rhs.data[i]) return true;
            if (++i == data_size) break;
        } while (data[i] == rhs.data[i]);
        return false;
    }

    bool operator>(const Super & rhs) const {
        return (rhs < *this);
    }

    bool operator<=(const Super & rhs) const {
        return !(rhs < *this);
    }

    bool operator>=(const Super & rhs) const {
        return !(*this < rhs);
    }

    // ~component-wise addition/subtraction~

    // compound +
    Super & operator+=(const Super & rhs) {
        for (int i = 0; i < data_size; ++i) {
            data[i] += rhs.data[i];
        }
        return *this;
    }

    // compound -
    Super & operator-=(const Super & rhs) {
        for (int i = 0; i < data_size; ++i) {
            data[i] -= rhs.data[i];
        }
        return *this;
    }



    // compound *
    // scalar multiply
    //) aka broadcast.
    //) Super v;
    //+ T a;
    //+ v *= a; // multiply v by scalar a
    // todo examine/test when prosize > 1
    Super & operator*=(const T /* call-by-copy here avoided in case of "big" T */ & rhs) {
        for (int i = 0; i < data_size; ++i) {
            data[i] *= rhs;
        }
        return *this;
    }

    // compound *
    // component-wise multiply
    //) This operation is sometimes denoted circle-times.
    //) A * B is dot product, not derived from this operation.
    Super & operator *=(const Super & rhs) {
        for (auto i = 0; i < data_size; ++i) {
            data[i] *= rhs.data[i];
        }
        return *this;
    }

    // compound /
    // scalar division
    //) If your type PP does not have division,
    //+ you must not use this operator.
    Super & operator/=(const T rhs) {
        for (int i = 0; i < data_size; ++i) {
            data[i] /= rhs;
        }
        return *this;
    }

    // compound /
    // component-wise division
    // todo:
    // ...can't decide what to do.
    // new setup has / and * as matrix multiplication.
    // Not even entirely sure what those methods will do, in general,
    // for a super-vector/super-matrix setup.
    // todo:
    // My sense right now is that this is "ok":
    // It seems crazy to me to divide matrices with /=.
    // So I am comfortable allowing this to be compoenent-wise,
    // even though conventionally, /= and operator/ are the same.
    Super & operator/=(const Super & rhs) {
        for (auto i = 0; i < data_size; ++i) {
            data[i] /= rhs.data[i];
        }
        return *this;
    }

    //-----

    // component-wise addition
    friend Super operator+(const Super & lhs, const Super & rhs) {
        Super out = lhs;
        out += rhs;
        return out;
    }

    friend Super operator-(const Super & lhs, const Super & rhs) {
        Super out = lhs;
        out -= rhs;
        return out;
    }

    // scalar times matrix
    friend Super operator*(const T & lhs, const Super & rhs) {
        Super out = rhs;
        out *= lhs;
        return out;
    }
    // matrix divided by scalar
    friend Super operator/(const Super & lhs, const T rhs) {
        Super out = lhs;
        out /= rhs;
        return out;
    }

    // l2 norm of vector, with some simple options:
    // v^1 is the length of v               << think of |v|^1
    // v^2 is the length of v squared  << think of |v|^2
    // v^3 is |v|^3, etc.
    // if the class is a matrix, it is the l2 matrix norm
    friend T operator^(const Super & ain, int pow) {
        T out = T::zero();
        T tmp;
        if (pow == 2) {
            for (int i = 0; i < data_size; ++i) {
                tmp = ain.data[i]^2;
                out += tmp;
            }
            return out;
        }
        else {
            out = (ain^2).sqrt();
            if (--pow > 0) {
                while (pow-- > 0) {
                    out *= out;
                }
            }
            return out;
        }
    }

    // return a unit instance with same direction as input instance ain.
    //) It will return the input divided by the input's L2 norm.
    //) Does not check for divide-by-zero.
    /*
     * It is tempting to check whether ain is zero
     * before dividing. However, that is left up
     * to the application, as it is expensive
     * and the application may not want to pay
     * that cost in cases where it has a priori knowledge
     * that the input is nonzero.
     * ...So just rely on T class to do "whatever it does"
     * upon divide-by-zero (and all other unusual cases).
     *
     */
    friend Super operator~(const Super & ain) {
        return ain/(ain^1);
    }

    // dot product v%w
    // more generally, the "total dot product",
    //+ the sum of the dot products of all rows (or columns).
    //) DevNote: scalar types must all use % as a renaming of multiplication.
    template <typename PP2>
    friend PP operator%(
            const Super & lhs,
            const Super<PP2, T, size, cosize> & rhs
    ) {
        PP out = PP::zero();
        //\\ the two inputs have the same data_size.
        for (auto i = 0; i < rhs.data_size; ++i) {
            out += lhs.data[i] % rhs[i];
        }
        return out;
    }

    // broadcasting operator%
    // Column Matrix Application:
    // operator% be used for ordinary matrix multiplication of a column matrix
    // (modeled as a vector of vectors, often a convenient format)
    // A, and a vector of scalars b, giving Ab in the usual sense.
    friend Super operator%(
            const Super & lhs,
            const T & rhs
    ) {
        return rhs*lhs;
    }

    // explicit column matrix multiplication (for reference)
    PP colmul(
        const Super<T, T, size, 1> & rhs
    ) {
        PP out = PP::zero();
        for (auto i = 0; i < rhs.data_size; ++i) {
            out += rhs[i]*data[i];
        }
        return out;
    }

    // cross product v & w
    friend Super operator&(const Super & lhs, const Super & rhs) {
        auto out = Super();
        //\\ the two inputs have the same size.
        // if the size is < 3 then it is "UB".
        out.data[0] = lhs.data[1]*rhs.data[2] - lhs.data[2]*rhs.data[1];
        out.data[1] = lhs.data[2]*rhs.data[0] - lhs.data[0]*rhs.data[2];
        out.data[2] = lhs.data[0]*rhs.data[1] - lhs.data[1]*rhs.data[0];
        return out;
    }

    // matrix multiplication A * B
    // method #1
    // v%w is the same as v.t() * w
    // Ax = b is the same as A * x = b
    // ABx = b is the same as A * B * x = b
    template<unsigned long intermediate_size>
    friend Super<PP, T, size, cosize> operator*(
            const Super<PP, T, size, intermediate_size> & lhs,
            const Super<PP, T, intermediate_size, cosize> & rhs
    ) {
        Super out;
        // todo multiplication with LUdec if present???
        for (auto i = 0; i < size; ++i) {
            for (auto k = 0; k < cosize; ++k) {
                // a strided dot product
                for (auto j = 0; j < intermediate_size; ++j) {
                    out.data[cosize*i + k] += lhs[intermediate_size * i + j] * rhs[cosize*j + k];
                }
            }
        }
        return out;
    }

    // matrix multiplication A * B
    // method #2
    // This method is applied whenever the two inputs have the same size and cosize.
    // When A and B are square matrices, it is identical to method #1.
    // When A and B are vectors, it is a component-wise multiplication.
    // Otherwise, the result is most likely not useful. (It's possible the user
    // made a mistake, in that case).
    // The method is implemented anyway, though. The main reason for this
    // is so that the compiler will compile the code, including various
    // templated methods that are present but that are not intended for actual use.
    // (For example, the inverse of a non-square matrix.)
    friend Super operator*(
            const Super & lhs,
            const Super & rhs
    ) {
        Super out;
        for (auto i = 0; i < size; ++i) {
            for (auto k = 0; k < cosize; ++k) {
                // a strided dot product
                for (auto j = 0; j < shortsize; ++j) {
                    out.data[cosize*i + k] += lhs[shortsize * i + j] * rhs[cosize*j + k];
                }
            }
        }
        return out;
    }

    // matrix multiplication A*B
    // method #3
    // When multiplying a block matrix A and a block vector x,
    // the block type for each matrix are slightly unmatched.
    template <typename PPM, size_t intermediate_size>
    friend Super operator*(
            const Super<PPM, T, size, intermediate_size> & lhs,
            const Super<PP, T, intermediate_size, cosize> & rhs
    ) {
        Super out;
        // todo multiplication with LUdec if present???
        for (auto i = 0; i < size; ++i) {
            for (auto k = 0; k < cosize; ++k) {
                // a strided dot product
                for (auto j = 0; j < intermediate_size; ++j) {
                    out.data[cosize*i + k] += lhs[intermediate_size * i + j] * rhs[cosize*j + k];
                }
            }
        }
        return out;
    }

    // division operator (in a non-division ring...)
    friend Super operator/(
            const Super & lhs,
            const Super & rhs
    ) {
        return lhs * rhs.inverse();
    }

    // Bracket notation: mainly intended for vectors,
    // but can be used for matrices and calces when one loops
    // over all cells without needing an index for each dimension.
    PP & operator[](size_t index) {
        return data[index];
    }

    const PP & operator[](size_t index) const {
        return data[index];
    }

    PP & operator[](DGCALC::IPoint<2> index) {
        return data[cosize*index.i[0] + index.i[1]];
    }

    const PP & operator[](DGCALC::IPoint<2> index) const {
        return data[cosize*index.i[0] + index.i[1]];
    }

public: // static methods:

    static Super zero() {
        return Super();
    }

    // returns "the obvious thing", even if size != cosize,
    // in which case there's no mathematical unit.
    static Super unit() {
        Super out = Super();
        for (auto i = 0; i < shortsize; ++i) {
            out.data[i + cosize*i] = PP::unit();
        }
        return out;
    }
	
	static  Super ones(T value = T::unit()) {
		Super out = Super();
		for (auto i = 0; i < data_size; ++i) {
			out.data[i] = PP::ones(value);
		}
	}

    static Super tiny() {
        Super out;
        for (auto i = 0; i < data_size; ++i) {
            out.data[i] = PP::tiny();
        }
        return out;
    }

public: // const methods:

    // sqrt, mainly so that the class itself fulfills the T constraints;
    // perhaps one day the method is used.
    // ...component-wise square root.
    Super sqrt() const {
        Super out;
        for (auto i = 0; i < data_size; ++i) {
            out.data[i] = data[i].sqrt();
        }
        return out;
    }

    Super abs() const {
        Super out;
        for (auto i = 0; i < data_size; ++i) {
            out.data[i] = data[i].abs();
        }
        return out;
    }

    Super neg() const {
        Super out;
        for (auto i = 0; i < data_size; ++i) {
            out.data[i] = data[i].neg();
        }
        return out;
    }

    T max() const {
        T out = data[0].max();
        T tmp;
        for (auto i = 1; i < data_size; ++i) {
            if ((tmp = data[i].max()) > out) out = tmp;
        }
        return out;
    }

    T min() const {
        T out = data[0].max();
        T tmp;
        for (auto i = 1; i < data_size; ++i) {
            if ((tmp = data[i].max()) < out) out = tmp;
        }
        return out;
    }

public: // API:

    // todo normal operations are convenient, but in place methods are
    // todo also necessary.
    // todo The in place methods must be *loudly* declared so...
    // todo such methods should be namespaced as InPlace::t, InPlace::neg, etc.

    // transpose in place
    // overwrites data but technically does not destroy it.
    // only succeeds if size = cosize.
    void t() {
        T tmp;
        for (unsigned long i = 0; i < shortsize; ++i) {
            for (unsigned long j = i+1; j < shortsize; ++j) {
                tmp = data[cosize*i + j];
                data[cosize*i + j] = data[i + cosize*j];
                data[i + cosize*j] = tmp;
            }
        }
    }

    // transpose via producing a copy of the data.
    // For square case, use t().
    Super<PP, T, cosize, size> get_t() {
        Super<PP, T, cosize, size> out;
        for (unsigned long i = 0; i < size; ++i) {
            for (unsigned long j = 0; j < cosize; ++j) {
                // This is just out[{j,i}] = data[{i,j}];
                out.raw()[size*j + i] = data[cosize*i + j];
            }
        }
        return out;
    }

public: // unsafe API:

    /*
     * Raw pointer to data in terms of PP type.
     * (In case of emergency, break glass.)
     *
     */
    PP * raw() {
        return data;
    }

public: // printing:

    std::string str_di() const {
        std::string out;
        if (cosize == 1) {
            out += "(";
            for (auto i = 0; i < size-1; ++i) {
                out += data[i].str() + ", ";
            }
            out += data[size-1].str() + ")";
        }
        else {
            //\\ cosize > 1.
            for (auto /* row */ i = 0; i < size; ++i) {
                out += "[";
                for (auto /* column */ j = 0; j < cosize-1; ++j) {
                    out += data[cosize*i + j].str() + "  ";
                }
                out += data[cosize*(i+1) - 1].str() + "]";
                if (i == size-1) out += ".\n";
                else out += " \n";
            }
        }
        return out;
    }

    // print as data for plotting/analysis
    // method is mainly intended for vectors.
    std::string str_dat(std::string sep = " ") const {
        std::string out;
        if (data_size == 0) return out;
        for (size_t i = 0; i < data_size-1; ++i) {
            out += data[i].str_dat() + sep;
        }
        out += data[data_size-1].str_dat();
        return out;
    }

    // print as data, in so-called "matrix" form:
    // 0th row <--> 0th line
    // 1st row <--> 1st line
    // etc.
    // ...If it helps, what this method prints
    // is exactly what str_di() prints, without the syntax clues.
    [[nodiscard]]
    std::string str_dat_matrix() const {
        std::string out;
        for (auto /* row */ i = 0; i < size; ++i) {
            for (auto /* column */ j = 0; j < cosize-1; ++j) {
                out += data[cosize*i + j].str_dat() + "  ";
            }
            out += data[cosize*(i+1) - 1].str_dat();
            out += " \n";
        }
        return out;
    }

    // print as data, in transposed matrix form.
    // ...This is convenient if the matrix is a discretization of xy space,
    // as in A[{i,j}] = f(x,y) for some f.
    // This method prints the matrix with the size dimension (i; x)
    // going across, and the cosize dimension (j; y) going down,
    // inside the string. NOTE: this is in agreement with Gnuplot convention.
    // ...Similar use case as in the previous one, the user wishes
    // that the output represent an xy-plane in the traditional orientation,
    // that is, origin at bottom left, x increases going right, y increases going up,
    // and A[{i,j}] = f(x,y) is some value of f(x,y).
    // In this case, set y_descends = false.
    [[nodiscard]]
    std::string str_dat_matrix_t(bool y_descends = true) const {
        std::string out;
        if (y_descends) {
            for (auto /* column */ j = 0; j < cosize; ++j) {
                for (auto /* row */ i_stride = 0; i_stride < cosize*(size-1); i_stride += cosize) {
                    out += data[i_stride + j].str_dat() + "  ";
                }
                out += data[cosize*(size-1) + j].str_dat();
                out += " \n";
            }
        } else {
            for (auto /* column */ j = cosize; j > 0; --j) {
                for (auto /* row */ i_stride = 0; i_stride < cosize*(size-1); i_stride += cosize) {
                    out += data[i_stride + (j-1)].str_dat() + "  ";
                }
                out += data[cosize*(size-1) + (j-1)].str_dat();
                out += " \n";
            }
        }
        return out;
    }

    [[nodiscard]]
    std::string str() const {
        return str_di();
    }

public: // viz:

    // visualization method
    //) viz() methods are smoky (generate files).
    // visualize a matrix of scalars, or a vector of scalars.
    void viz(std::string handle, std::string title = "A Default Title") const {
        static_assert(PP::size == 1 && PP::cosize == 1);
        if constexpr (cosize == 1) {
            //\\ It is a vector.
            // create a "distribution visualization"
//            std::cout << "need to implement this viz() method!" << std::endl;
        }
        else {
            //\\ cosize != 1. It is a matrix.
            // create a "matrix visualization"
            T vmax = max();
            T vmin = min();
            // todo write script with subpatterns, and call the subber to sub in the patterns......
        }
    }

    // visualization method
    //) viz() methods are smoky (generate files).
    // visualize the object's sparsity: print as matrix.
    void viz_sparse(std::string handle, std::string title = "A Default Title") const {
        static_assert(PP::size == 1 && PP::cosize == 1);

    }

    // visualization intermediate representation
    //) Utility for higher viz() methods.
    [[nodiscard]]
    Super<PP, viz_t, size, cosize> viz_ir() const {
        Super<PP, viz_t, size, cosize> out;
        for (int i = 0; i < data_size; ++i) {
            out.data[i] = data[i].viz_ir();
        }
        return out;
    }

//-------------------------------------------------------------------
public: // algorithms:

    // Performs Gaussian elimination with partial pivoting
    // The matrix A is not preserved.
    // the answer is returned in b.
    // todo if matrix array is not square?
    void solve_GE(Super<PP, T, size, 1> & b) {
        size_t i, j, k, l;
        T larger, largest, ratio;
        // index ordering
        size_t order[size];
        T c[size];
        for (i = 0; i < size; ++i)
            order[i] = i;
        for (/* column */ j = 0; j < size - 1; ++j) {
            // find the largest element in this column
            largest = data[cosize*order[j] + j].abs();
            k = j;
            for ( /* row */ i = j + 1; i < size; ++i) {
                larger = data[cosize*order[i] + j].abs();
                if (larger > largest) {
                    largest = larger;
                    k = i;
                }
            }
            // update the row order
            l = order[j];
            order[j] = order[k];
            order[k] = l;
            // subtract the pivot row (j) from row i
            for ( /* row */ i = j + 1; i < size; ++i) {
                ratio = T::zero() - data[cosize*order[i] + j] / data[cosize*order[j] + j];
                for(/* column */ l = j + 1; l < size; ++l) {
                    data[cosize*order[i] + l] += ratio * data[cosize*order[j] + l];
                }
                b[order[i]] += ratio * b[order[j]];
            }
        }
        // backward substitution loop
        b[order[size - 1]] /= data[cosize*order[size - 1] + size - 1];
        for (/* row */ i = size - 2; i != -1; --i) {
            for (/* column */ j = i + 1; j < size; ++j){
                b[order[i]] -= data[cosize*order[i] + j] * b[order[j]];
            }
            b[order[i]] /= data[cosize*order[i] + i];
        }
        // reorder the solution vector
        for (i = 0; i < size; ++i)
            c[i] = b[order[i]];
        for (i = 0; i < size; ++i)
            b[i] = c[i];
    };

    // solve Ax = b, write solution into b
    // todo memoize the LUdecomposition!!!!!!!!!!!!!
    void solve_LU(Super<PP, T, size, 1> & b) const {
        SuperLUdec<PP, T, size, cosize> lu(*this);
        lu.solve(b);
    }

    // default solve routine
    void solve(Super<PP, T, size, 1> & b) const {
        solve_LU(b);
    }

    // the inverse (call only for square matrix)
    Super inverse() const {
        SuperLUdec<PP, T, size, cosize> lu(*this);
        return lu.inverse();
    }

    // Performs Thomas algorithm on a tridiagonal
    // matrix array.
    // The matrix array is preserved,
    // the answer is returned in b.
    // ...This method is essentially for testing purposes;
    // its storage policy scales very poorly.
    // todo error log
    void tridiagonal_solve(Super<PP, T, size, 1> & b) const {
        int i;
        PP tmp;
        Super<PP, T, size, 1> vtmp;
        //if (data[0] == T::zero()) assert(0 == 1);
        b[0] = b[0] / data[0];
        tmp = data[0];
        for (i = 1; i < size; ++i) {
            vtmp[i] = data[cosize*(i - 1) + i] / tmp;
            tmp = data[cosize*i + i] - data[cosize*i + (i - 1)] * vtmp[i];
            // if (tmp == T::zero()) assert(0 == 1);
            b[i] = (b[i] - data[cosize*i + (i - 1)]*b[i-1]) / tmp;
        }
        for (i = size-2; i >= 0; --i) {
            b[i] -= vtmp[i + 1] * b[i + 1];
        }
    };

    // Performs modified Thomas algorithm on a band-diagonal
    // matrix array.
    // The matrix array is preserved,
    // the answer is returned in b.
    // ...This method is essentially for testing purposes;
    // its storage policy scales very poorly.
    void banddiagonal_solve(Super<PP, T, size, 1> & b) {

    }

    // todo methods for band-rotated matrix: br_tridiagonal_solve, br_banddiagonal_solve, etc.
    // Reminder: these are rotated 45 degrees for space optimization.

};

