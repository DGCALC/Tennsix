//
// Created by Lucius Schoenbaum on 11/6/21.
//

#pragma once


#include <Tbasis.h>

#include <scalar/types/double754.h>
#include <tennsix/super.h>
#include <tennsix/linear.h>

class TSuper {

public: // API:

    static void Tdryrun() {
        Super<double754, double754, 4, 1> foo{double754(123.0), double754(456.0), double754(789.0), double754(222.0)};
        cout << foo.str() << "\n";
        cout << (foo^2).str() << " is the norm squared.\n";
        cout << (foo^1).str() << " is the norm.\n";
        cout << "v[1] = " << foo[1].str() << "\n";
        Super<double754, double754, 4, 1> bar{double754(1.0), double754(1.0), double754(-1.0), double754(-1.0)};
        cout << bar.str() << "\n";
        cout << (bar^2).str() << " is the norm squared.\n";
        cout << (bar^1).str() << " is the norm.\n";
        cout << "v[1] = " << bar[1].str() << "\n";
        cout << "dimension is " << Vector<double754, 4>::dimension << ".\n";
        Vector<double754, 4> foo2{double754(123.0), double754(456.0), double754(789.0), double754(222.0)};
        cout << foo2.str() << "\n";
        cout << (foo2^2).str() << " is the norm squared.\n";
        cout << (foo2^1).str() << " is the norm.\n";
        cout << "v[1] = " << foo2[1].str() << "\n";
        Vector<double754, 4> bar2{double754(1.0), double754(1.0), double754(-1.0), double754(-1.0)};
        cout << bar2.str() << "\n";
        cout << (bar2^2).str() << " is the norm squared.\n";
        cout << (bar2^1).str() << " is the norm.\n";
        cout << "dimension is " << Vector<double754, 4>::dimension << ".\n";
        cout << "v[1] = " << bar2[1].str() << "\n";
    }
    // test vectors
    static inline Vector<double754, 2> v = {double754(1), double754(2)};
    static inline Vector<double754, 2> v11 = {double754(-2), double754(3)};
    static inline Vector<double754, 2> v2 = {double754(-1), double754(1)};
    static inline Vector<double754, 3> v3 = {double754(1), double754(2), double754(1)};
    static inline Vector<double754, 3> v4 = {double754(-1), double754(1), double754(2)};

    // test matrices
    static inline Matrix<double754, 2, 2> mt = {double754(1), double754(2), double754(3), double754(4)};
    static inline Matrix<double754, 2, 2> mt1 = {double754(4), double754(3), double754(2), double754(1)};
    static inline Matrix<double754, 2, 2> mt2 = {double754(0), double754(1), double754(2), double754(1)};

    // tridiagonal test matrix
    static inline Matrix<double754, 3, 3> tri = {
            double754(1),
            double754(2),
            double754(0),
            //
            double754(1),
            double754(3),
            double754(-1),
            //
            double754(0),
            double754(2),
            double754(3),
    };



    static inline Matrix<double754, 2, 3> mtt = {
            double754(1),
            double754(2),
            double754(3),
            double754(4),
            double754(5),
            double754(6)
    };
    static inline Matrix<double754, 2, 3> mtt2 = {
            double754(1),
            double754(0),
            double754(1),
            double754(0),
            double754(-1),
            double754(0)
    };

public: // API:

    static void Tctor_etc() {
        cout << "0: " << mt.data[0].str();
        cout << "1: " << mt.data[1].str() << endl;
        cout << mt.str() << endl;
        cout << v.str() << endl;
        cout << "empty:\n" << Matrix<double754, 2, 2>().str();
        cout << "unit:\n" << Matrix<double754, 2, 2>::unit().str();
        cout << "off-unit:\n" << Matrix<double754, 3, 2>::unit().str();
        cout << "uniform:\n" << Matrix<double754, 2, 2>(double754(3)).str();
    }

    static void Tcopy() {
        Matrix<double754, 2, 2> cp = mt;
        cout << mt.str();
        cout << cp.str();
        if (cp == mt) {
            cout << "they are equal" << endl;
        }
    }

    static void Tbracket() {
        cout << "test of bracket operation\n";
        for (int i = 0; i < mtt.data_size; ++i) {
            cout << mtt.data[i].str() << endl;
            cout << mtt[i].str() << endl;
        }
        cout << "...\n";
        for (int i = 0; i < mtt.data_size; ++i) {
            mtt[i] = double754(10 + i);
            cout << mtt.str() << endl;
        }
        cout << "...\n";
        for (size_t i = 0; i < 2UL; ++i) {
            for (size_t j = 0; j < 3UL; ++j) {
                cout << mtt[{i, j}].str() << endl;
            }
        }
        cout << "...\n";
        for (size_t i = 0; i < 2UL; ++i) {
            for (size_t j = 0; j < 3UL; ++j) {
                mtt[{i,j}] = double754(i + j);
                cout << mtt.str() << endl;
            }
        }
        cout << "bye\n";
    }

    static void Ttranspose() {
        cout << mt.str() << endl;
        mt.t();
        cout << mt.str() << endl;
        cout << mtt.str() << endl;
        for (int i = 0; i < mtt.data_size; ++i) {
            cout << mtt.data[i].str() << endl;
        }
        mtt.t();
        cout << mtt.str() << endl;
        for (int i = 0; i < mtt.data_size; ++i) {
            cout << mtt.data[i].str() << endl;
        }
    }

    static void Tstr_dat_matrix() {
        cout << mt.str() << endl;
        cout << mt.str_dat_matrix() << endl;
        cout << mt.str_dat_matrix_t() << endl;
        cout << mt.str_dat_matrix_t(false) << endl;
        cout << mtt.str() << endl;
        cout << mtt.str_dat_matrix() << endl;
        cout << mtt.str_dat_matrix_t() << endl;
        cout << mtt.str_dat_matrix_t(false) << endl;
    }

    static void Tdotproduct() {
        cout << mt.str() << endl;
        cout << mt2.str() << endl;
        cout << (mt % mt2).str() << endl;
    }

    static void Tmultiply() {
        cout << "A = \n" << mt.str();
        cout << "v = " << v.str() << endl;
        cout << "A*v = \n" << (mt * v).str() << endl;
        cout << "A should not be overwritten:\n";
        cout << "A = \n" << mt.str();
        cout << "B = \n" << mtt.str();
        cout << "A*B = \n" << (mt * mtt).str() << endl;
    }

    static void Tunit_length() {
        cout << mt.str();
        cout << "length = " << (mt^1).str() << endl;
        cout << "length squared = " << (mt^2).str() << endl;
        mt = ~mt;
        cout << "unit = " << mt.str();
        cout << "length = " << (mt^1).str() << endl;
        cout << "length squared = " << (mt^2).str() << endl;
    }

    static void Tcross_product() {
        cout << v3.str() << endl;
        cout << v4.str() << endl;
        cout << (v3 & v4).str() << endl;
    }

    static void Tsolve() {
        cout << "A =\n";
        cout << mt.str();
        cout << "b =\n";
        cout << v.str() << endl;
        mt.solve(v);
        cout << "...solve...\n";
        cout << "A =\n";
        cout << mt.str();
        cout << "x =\n";
        cout << v.str() << endl;
        cout << "Ax = \n";
        cout << (mt * v).str() << "\n";
    }

    static void Ttridiagonalsolve() {
        cout << "A = \n";
        cout << tri.str();
        cout << "b = ";
        cout << v3.str() << "\n";
        tri.tridiagonal_solve(v3);
        cout << "solve Ax = b. \nx = \n";
//        cout << tri.str();
        cout << v3.str() << "\n";;
        cout << "Check: A*x = \n";
        cout << (tri * v3).str() << "\n";
        cout << "Let b = x. Solve Ax = b.\nx = \n";
        tri.tridiagonal_solve(v3);
        cout << v3.str() << "\n";
        cout << "Check: A*x = \n";
        cout << (tri * v3).str() << "\n";
        cout << "A is preserved: A = \n" << tri.str();
        cout << "bye";
    }


    // note: for mt, LU is [[1,2],[3,-2]]
    static void TSuperLUdec() {
        Super<double754, double754, 2, 2> L, U, LxU;
        cout << "\n\nTEST 1:\n";
        cout << mt1.str();
        SuperLUdec<double754, double754, 2>  lu1(mt1);
        cout << "this is lu?\n";
        cout << lu1.data.str();
        cout << "this is U?\n";
        cout << lu1.U().str();
        cout << "this is L?\n";
        cout << lu1.L().str();
        cout << "this is original (possibly w/ rows permuted)?\n";
        L = lu1.L();
        U = lu1.U();
        LxU = L * U;
        cout << LxU.str();
        cout << "\n\nTEST 2:\n";
        cout << mt.str();
        SuperLUdec<double754, double754, 2>  lu(mt);
        cout << "this is lu? should be [[1, 2] [3, -2]].\n";
        cout << lu.data.str();
        cout << "this is U?\n";
        cout << lu.U().str();
        cout << "this is L?\n";
        cout << lu.L().str();
        cout << "this is original (possibly w/ rows permuted)?\n";
        L = lu.L();
        U = lu.U();
        LxU = L * U;
        cout << LxU.str();
        cout << "\n\nTEST 3:\n";
        cout << mt2.str();
        SuperLUdec<double754, double754, 2>  lu2(mt2);
        cout << "this is lu?\n";
        cout << lu2.data.str();
        cout << "this is U?\n";
        cout << lu2.U().str();
        cout << "this is L?\n";
        cout << lu2.L().str();
        cout << "this is original (possibly w/ rows permuted)?\n";
        cout << (lu2.L() * lu2.U()).str();
        cout << "\n\nTEST 4:\n";
        cout << tri.str();
        SuperLUdec<double754, double754, 3>  lu3(tri);
        cout << "this is lu?\n";
        cout << lu3.data.str();
        cout << "this is U?\n";
        cout << lu3.U().str();
        cout << "this is L?\n";
        cout << lu3.L().str();
        cout << "this is original (possibly w/ rows permuted)?\n";
        cout << (lu3.L() * lu3.U()).str();
    };


    static void inverse_v(const Matrix<double754, 2, 2> & m) {
        Matrix<double754, 2, 2> minv = m.inverse();
        cout << "inverse:" << endl;
        cout << minv.str();
        cout << "check:" << endl;
        Matrix<double754, 2, 2> prod = m*minv;
        cout << prod.str();
    }

    static void Tinverse_LU() {
        cout << "\n\n\ntest inverse via LU\n";
        cout << mt1.str();
        inverse_v(mt1);
        cout << "\n\n\ntest inverse via LU\n";
        cout << mt.str();
        inverse_v(mt);
        cout << "correct inverse: [ -2 1][3/2 -1/2].\n";
        cout << "\n\n\ntest inverse via LU\n";
        cout << mt2.str();
        inverse_v(mt2);
    }

//    [[1.000000  0.000000]
//    [0.000000  1.000000].
    static inline Super<double754, double754, 2, 2> J11{
        double754(1),
        double754(0),
        double754(0),
        double754(1)
    };
//    [0.000000  0.000000]
//    [0.000000  0.000000].
    static inline Super<double754, double754, 2, 2> J12{
            double754(0),
            double754(0),
            double754(0),
            double754(0)
    };
//    [0.000000  0.000000]
//    [0.000000  0.000000].
    static inline Super<double754, double754, 2, 2> J13{
            double754(0),
            double754(0),
            double754(0),
            double754(0)
    };
//    [0.000000  0.000000]
//    [0.000000  0.000000].
//    ]
    static inline Super<double754, double754, 2, 2> J14{
            double754(0),
            double754(0),
            double754(0),
            double754(0)
    };
//    [[1.000000  0.333333]
//    [0.314986  1.000000].
    static inline Super<double754, double754, 2, 2> J21{
            double754(1.0),
            double754(0.333),
            double754(0.314986),
            double754(1.0)
    };
//    [-1.000000  0.333333]
//    [0.314986  -1.000000].
    static inline Super<double754, double754, 2, 2> J22{
            double754(-1.0),
            double754(0.3333),
            double754(0.314986),
            double754(-1.0)
    };
//    [0.000000  0.000000]
//    [0.000000  0.000000].
    static inline Super<double754, double754, 2, 2> J23{
            double754(0.0),
            double754(0.0),
            double754(0.0),
            double754(0.0)
    };
//    [0.000000  0.000000]
//    [0.000000  0.000000].
//    ]
    static inline Super<double754, double754, 2, 2> J24{
            double754(0.0),
            double754(0.0),
            double754(0.0),
            double754(0.0)
    };
//    [[0.000000  0.000000]
//    [0.000000  0.000000].
    static inline Super<double754, double754, 2, 2> J31{
            double754(0.0),
            double754(0.0),
            double754(0.0),
            double754(0.0)
    };
//    [1.000000  0.333333]
//    [0.180101  1.000000].
    static inline Super<double754, double754, 2, 2> J32{
            double754(1.0),
            double754(0.3333),
            double754(0.180101),
            double754(1.0)
    };
//    [-1.000000  0.333333]
//    [0.180101  -1.000000].
    static inline Super<double754, double754, 2, 2> J33{
            double754(-1.0),
            double754(0.3333),
            double754(0.180101),
            double754(-1.0)
    };
//    [0.000000  0.000000]
//    [0.000000  0.000000].
//    ]
    static inline Super<double754, double754, 2, 2> J34{
            double754(0.0),
            double754(0.0),
            double754(0.0),
            double754(0.0)
    };
//    [[0.000000  0.000000]
//    [0.000000  0.000000].
    static inline Super<double754, double754, 2, 2> J41{
            double754(0.0),
            double754(0.0),
            double754(0.0),
            double754(0.0)
    };
//    [0.000000  0.000000]
//    [0.000000  0.000000].
    static inline Super<double754, double754, 2, 2> J42{
            double754(0.0),
            double754(0.0),
            double754(0.0),
            double754(0.0)
    };
//    [1.000000  0.333333]
//    [-0.031908  1.000000].
    static inline Super<double754, double754, 2, 2> J43{
            double754(1.0),
            double754(0.3333),
            double754(-0.031908),
            double754(1.0)
    };
//    [-1.000000  0.333333]
//    [-0.031908  -1.000000].
//    ].
    static inline Super<double754, double754, 2, 2> J44{
            double754(-1.0),
            double754(0.3333),
            double754(-0.031908),
            double754(-1.0)
    };
    static inline Super<Super<double754, double754, 2, 2>, double754, 4, 4> Jt1{
        J11,
        J12,
        J13,
        J14,
        J21,
        J22,
        J23,
        J24,
        J31,
        J32,
        J33,
        J34,
        J41,
        J42,
        J43,
        J44
    };

    static void Tblockmultiply() {
        cout << "Test of block multiply:\n";
        cout << "A IS " << Jt1.str();
        Super<Super<double754, double754, 2, 2>, double754, 4, 4> B;
        Super<Super<double754, double754, 2, 2>, double754, 4, 4> AxB;
        B = Super<Super<double754, double754, 2, 2>, double754, 4, 4>::unit();
        cout << "B IS " << B.str();
        AxB = Jt1*B;
        cout << "A*B IS " << AxB.str();
    }

    static void Tblockinverse() {
//        cout << "Test of block inverse:\n";
//        Super<Super<double754, double754, 2, 2>, double754, 2, 2> A {
//                Super<double754, double754, 2, 2>::unit(), Super<double754, double754, 2, 2>::zero(), Super<double754, double754, 2, 2>::zero(), Super<double754, double754, 2, 2>::unit()
//        };
//        cout << "A IS " << A.str();
//        Super<Super<double754, double754, 2, 2>, double754, 2, 2> Ainv;
//        Ainv = A.inverse();
//        cout << "INVERSE IS " << Ainv.str();
//        cout << "PRODUCT IS " << (A * Ainv).str();

//        cout << "\n\n\n" << endl;
//        Super<Super<double754, double754, 2, 2>, double754, 2, 2> A2 {
//                Super<double754, double754, 2, 2>::zero(),
//                Super<double754, double754, 2, 2>::unit(),
//                Super<double754, double754, 2, 2>::unit(),
//                Super<double754, double754, 2, 2>::zero()
//        };
//        cout << "A IS " << A2.str();
//        Super<Super<double754, double754, 2, 2>, double754, 2, 2> A2inv;
//        A2inv = A2.inverse();
//        cout << "INVERSE IS " << A2inv.str();
//        cout << "PRODUCT IS " << (A2 * A2inv).str();


        cout << "\n\n\n" << endl;
        cout << "A IS " << Jt1.str();
        Super<Super<double754, double754, 2, 2>, double754, 4, 4> Jt1inv;
        Jt1inv = Jt1.inverse();
        cout << "INVERSE IS " << Jt1inv.str();
        cout << "PRODUCT IS " << (Jt1 * Jt1inv).str();
    }

    static void Tblocksolve() {
        cout << "\n\n\n" << endl;
//        Super<Super<double754, double754, 2, 1>, double754, 2, 1> b {
//                {1.0, 2.0},
//                {3.0, 4.0}
//        };
//        cout << "A IS " << Jt1.str();
//        cout << "b IS " << b.str() << endl;
//        Jt1.solve(b);
//        cout << "Ax = b" << endl;
//        cout << "x IS " << b.str() << endl;
    }


    static void Tcolumnmatrix() {
        ColumnMatrix<double754, 2, 2> colmat;
        colmat[0] = v11;
        colmat[1] = v2;
        cout << "colmat =\n" << colmat.str() << endl;
        cout << "v = " << v.str() << endl;
        cout << "v2 = " << v2.str() << endl;
        cout << "colmat*v = " << (colmat.colmul(v)).str() << endl;
        cout << "colmat*v2 = " << colmat.colmul(v2).str() << endl;
        cout << "colmat % v = " << (colmat % v).str() << endl;
        cout << "colmat % v2 = " << (colmat % v2).str() << endl;
    }

    static void Tgeneral_subscript() {
       using Arr = Matrix<Matrix<double754, 5, 4>, 3, 2>;
       Arr arr;
       arr[{1,2,1,2}] = double754(42);
       cout << arr.str() << endl;
    }

};
