/*
    Title:      Lab 1 – main.cpp
    Purpose:    Test Polynomial Class Methods
    Author:     Rami Isaac, Prof. Carlos Arias Arevalo
    Date:       April 19, 2020
*/

/*
#include "polynomial.h"

#include <iostream>
#include <sstream>

using std::cout;
using std::cin;
using std::endl;
using std::stringstream;

int main(int argc, char* argv[]){

    stringstream data[5];
    data[0].str("3 -1 2 0 -2.5"); // a
    data[1].str("4 -1 -2 0 0 3"); // b
    data[2].str("2 -1 0 4");      // c
    data[3].str("3 -6 -5 2 1");   // d
    data[4].str("1 -2 1");        // e

    stringstream answers[6];
    answers[0].str("4 -2 0 0 -2.5 3"); // a + b
    answers[1].str("7 1 0 -4 2.5 2 6 0 -7.5"); // a * b
    answers[2].str("2 3 4 1"); // d / e
    answers[3].str("2 1 0 -4"); // -c
    answers[4].str("4 0 -4 0 2.5 3"); // b - a
    answers[5].str("2 2 0 -7.5"); // d/dx a

    Polynomial a(0), b(0), c(0), d(3), e(1), r(0);
    a.Read(data[0]); b.Read(data[1]);
    c.Read(data[2]); d.Read(data[3]); e.Read(data[4]);

    Polynomial s(0), t(0), u(0), v(0), w(0), x(0);
    s.Read(answers[0]); t.Read(answers[1]); u.Read(answers[2]);
    v.Read(answers[3]); w.Read(answers[4]); x.Read(answers[5]);

    //--------------------------------------------------------------------------
    //                       REQUIRED METHODS - UNIT TESTS
    //--------------------------------------------------------------------------

    // SUM
    cout << endl << endl;
    if (a.Sum(b).Equals(s)){ cout << "SUM PASSED " << endl; }
    else { cout << "SUM FAILED" << endl; }
    cout << " IN: " << a.ToString() << " | "; a.Write(cout); cout << endl;
    cout << " IN: " << b.ToString() << " | "; b.Write(cout); cout << endl;
    cout << "OUT: " << s.ToString() << " | "; s.Write(cout); cout << endl;
    cout << "ACT: " << a.Sum(b).ToString() << " | "; (a.Sum(b)).Write(cout);

    // SUBTRACT
    cout << endl << endl;
    if (b.Subtract(a).Equals(w)){ cout << "SUBTRACT PASSED " << endl;
    }else{ cout << "SUBTRACT FAILED" << endl; }
    cout << " IN: " << b.ToString() << " | "; b.Write(cout); cout << endl;
    cout << " IN: " << a.ToString() << " | "; a.Write(cout); cout << endl;
    cout << "OUT: " << w.ToString() << " | "; w.Write(cout); cout << endl;
    cout << "ACT: " << b.Subtract(a).ToString() << " | "; (b.Subtract(a)).Write(cout);

    // MULTIPLY
    cout << endl << endl;
    if (a.Multiply(b).Equals(t)){ cout << "MULTIPLY PASSED " << endl;
    }else{ cout << "MULTIPLY FAILED" << endl; }
    cout << " IN: " << a.ToString() << " | "; a.Write(cout); cout << endl;
    cout << " IN: " << b.ToString() << " | "; b.Write(cout); cout << endl;
    cout << "OUT: " << t.ToString() << " | "; t.Write(cout); cout << endl;
    cout << "ACT: " << a.Multiply(b).ToString() << " | "; (a.Multiply(b)).Write(cout);

    // DERIVE
    cout << endl << endl;
	if (a.Derive().Equals(x)){ cout << "DERIVE PASSED " << endl;
	}else{ cout << "DERIVE FAILED" << endl;}
    cout << " IN: " << a.ToString() << " | "; a.Write(cout); cout << endl;
    cout << "OUT: " << x.ToString() << " | "; x.Write(cout); cout << endl;
    cout << "ACT: " << a.Derive().ToString() << " | "; (a.Derive()).Write(cout);

    // EVALUATE
    cout << endl << endl;
	if (b.Evaluate(3) == 236.0){ cout << "EVALUATE PASSED " << endl;
	}else{ cout << "EVALUATE FAILED" << endl;}
    cout << " IN: " << b.ToString() << " | "; b.Write(cout); cout << endl;
    cout << " IN: " << "3" << endl;
    cout << "OUT: " << "236" << endl;
    cout << "ACT: " << b.Evaluate(3) << endl;

    //--------------------------------------------------------------------------
    //                       EXTRA CREDIT - UNIT TESTS
    //--------------------------------------------------------------------------

    cout << endl << endl << "EXTRA CREDIT" << endl;

    // INTEGRATE
    cout << endl << endl;
	if (a.Integrate(1, 2) == -7.375){ cout << "INTEGRATE PASSED " << endl;
	}else{ cout << "INTEGRATE FAILED" << endl;}
    cout << " IN: " << a.ToString() << " | "; a.Write(cout); cout << endl;
    cout << " IN: " << "1, 2" << endl;
    cout << "OUT: " << "-7.375" << endl;
    cout << "ACT: " << a.Integrate(1,2) << endl;

    // DIVIDE
    cout << endl << endl;
    if (d.Divide(e).Equals(u)){cout << "DIVIDE PASSED " << endl;
    }else{cout << "DIVIDE FAILED" << endl;}
    cout << " IN: " << d.ToString() << " | "; d.Write(cout); cout << endl;
    cout << " IN: " << e.ToString() << " | "; e.Write(cout); cout << endl;
    cout << "OUT: " << u.ToString() << " | "; u.Write(cout); cout << endl;
    cout << "ACT: " << "n/a" << endl;

    return 0;
}
*/