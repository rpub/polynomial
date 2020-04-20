/*
    Title:      Lab 1 – polynomial.cpp
    Purpose:    Implement Polynomial Class Methods
    Author:     Rami Isaac, Prof. Carlos Arias Arevalo
    Date:       April 19, 2020
*/

#include "polynomial.h"

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cfloat>
#include <algorithm>
#include <math.h>

using namespace std;

using std::istream;
using std::ostream;
using std::string;
using std::stringstream;
using std::fixed;
using std::setprecision;
using std::showpos;


/**
* Polynomial
* Creates a polynomial of degree 'degree'.
 * Dynamically allocates array of size 'degree' + 1.
 * Sets coefficients to 0.0
* @param element the degree of the polynomial
* @returns n/a
*/
Polynomial::Polynomial(size_t degree) : _degree(degree){
	_coefficients = new float[_degree + 1];
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = 0.0;
	}
}

/**
* Polynomial
* Creates a polynomial of degree degree.
 * Assigns values of coefficients to 'coefficients' array
* @param element the degree of the polynomial and the array of coefficients
* @returns n/a
*/
Polynomial::Polynomial(size_t degree, const float* coefficients): _degree(degree){
	_coefficients = new float[_degree + 1];
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = coefficients[i];
	}
}

/**
* Polynomial
* Copy constructor
 * Creates deep copy of other polynomials
* @param element the degree of the polynomial
* @returns n/a
*/
Polynomial::Polynomial(const Polynomial& polynomial): _degree(polynomial._degree){
	_coefficients = new float[_degree + 1];
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = polynomial._coefficients[i];
	}
}

/**
* ~Polynomial
* Destructor
 * Deallocates allocated memory
* @param element n/a
* @returns n/a
*/
Polynomial::~Polynomial(){
     delete[] _coefficients;
}

/**
* Sum
* Returns polynomial with sum of this and rhs
* @param element the two polynomials you want to add
* @returns 'polySum', the sum of '*this' and 'rhs'
*/
const Polynomial Polynomial::Sum(const Polynomial& rhs)const {

    /*
     * Concept:
     * The addition of two polynomials can be achieved algorithmically through...
      * 1. Use of an 'if' statement to determine which polynomial has a greater degree
      * 2. Use of a 'for' loop limited in iterations to the smaller degree + 1
      *    -> Make a copy of the longer polynomial (this is also the return value)
      *    -> Add values of shorter polynomial to values of longer polynomial
      * 3. Return polynomial 'polySum' with summed value
     */

    // Bool to determine if degree of this is greater than that of rhs
    bool thisIsGreater = (_degree > rhs._degree) ? true : false;

    size_t i; // Loop counter

    // If thisIsGreater...
    if(thisIsGreater) {
        Polynomial polySum(*this); // Create polynomial of degree 'this'
        for (i = 0; i < rhs._degree + 1; i++) { // Loop until less than rhs degree + 1
            polySum._coefficients[i] += rhs._coefficients[i]; // Add, then assign values at corresponding index [i]
        }
        return polySum; // Returns new, summed polynomial
    }

    // If degree of rhs is greater...
    if(!thisIsGreater) {
        Polynomial polySum(rhs); // Create polynomial of degree 'rhs'
        for (i = 0; i < _degree + 1; i++){ // Loop until less than degree + 1
            polySum._coefficients[i] += _coefficients[i]; // Add, then assign values at corresponding index [i]
        }
        return polySum; // Returns new, summed polynomial
    }
}

/**
* Subtract
* Returns polynomial with difference of this and rhs
* @param elements the two polynomials you want to substract
* @returns 'polySub', the difference of '*this' and 'rhs'
 */
const Polynomial Polynomial::Subtract(const Polynomial& rhs)const{

    /*
     * Concept:
     * 5 - 2 = 3
     * ... is the same as 5 + (-2) = 3
     *
     * Equation: P1 + (-P2)
     *
     * Therefore, the Sum method can be used to add two polynomials if...
      * 1. [ P1 ] The 1st Polynomial is copied to a new Polynomial 'polySub'
      * 2. [ P2 ] The Minus method is used to find the additive inverse of the 2nd Polynomial
      * 3. [ P1 + (-P2) ] P1 is added (Sum method) to the additive inverse of P2
     */

    // P1, Copy to 'polySub'
    Polynomial polySub(*this); // copy of P1

    // Return P1 + (-P2)
    return polySub.Sum(rhs.Minus());
}

/**
* Minus
* Returns polynomial with addidtive inverse of this
* @param element the polynomial you want to find the addititve inverse of
* @returns 'retval', additive inverse of '*this'
*/
const Polynomial Polynomial::Minus()const{
	Polynomial retVal(*this);
	for (size_t i = 0; i < _degree + 1; i++) {
		retVal._coefficients[i] *= -1;
	}
	return retVal;
}

/**
* Multiply
* Returns polynomial with product of 'this' and 'rhs'
* @param element the two polynomials you want to multiply
* @returns 'polyMult', the product of '*this' and 'rhs'
*/
const Polynomial Polynomial::Multiply(const Polynomial& rhs)const{

    /*
     * Concept:
     * The multiplication of two polynomials can be achieved algorithmically through...
      * 1. The declaration of a new polynomial with a 'degree' = '_degree' + 'rhs._degree'
      * 2. The implementation of a double-nested 'for' loop where...
      *    -> Loop 'i' iterates through the values of the first polynomial array
      *    -> Loop 'x' iterates though the values of the second polynomial array
      * 3. The coefficients of corresponding index [i] and [x] and multiplied
      * 4. The resulting coefficient is then assigned to the new polynomial 'polyMult', then returned
     */

    // New polynomial (and return value) declared with the combined degree of both polynomials
    Polynomial polyMult(_degree + rhs._degree);

    size_t i, x; // Loop counter

    // Double-nested 'for' loop to iterate through both polynomials
    for (i = 0; i < _degree + 1; i++) { // Polynomial 1
        for (x = 0; x < rhs._degree + 1; x++) { // Polynomial 2
            // Assign value of P1[i] * P2[x] to 'polyMult' [i+x]
            polyMult._coefficients[i + x] += _coefficients[i] * rhs._coefficients[x];
        }
    }
    return polyMult; // return value
}

/**
* Divide
* Returns polynomial with division of this by rhs
* @param element n/a
* @returns n/a
*/
const Polynomial Polynomial::Divide(const Polynomial& rhs)const{
    return Polynomial(0);
}

/**
* Derive
* Returns polynomial with derivative of this
* @param element the polynomial you want to derive
* @returns 'polyDer', the derivative of '*this' and 'rhs'
*/
const Polynomial Polynomial::Derive()const {

    /*
     * Concept:
     * The derivative of 5x^3 = 15x^2
     * ...which can be broken down into 5(3)x^(3-1)
     * ... which is the same as coefficient*exponent*x^(exponent-1)
     * ... where the coefficient is multiplied by the degree, then the degree is reduced by 1
     *
     * Equation: c*e*x^(e-1)
     *
     * The derivative of a polynomial can be found algorithmically through the following:
      * 1. [ e - 1 ] Define new polynomial of 1 degree less than the original polynomial
      * 2. [ c ] Assign the values at index [i + 1] of the original polynomial to [i] of the new polynomial
      *    -> This accounts for the lesser degree and skips the old polynomial's first value
      * 3. [ c * e ] Multiply the values at index [i] by (counter i + 1)
      *    -> This multiplies the coefficient by the proper exponent
     */

    // The polynomial returned, 1 degree less original
    Polynomial polyDer(_degree - 1);

    size_t i; // Loop counter

    for(i = 0; i < _degree; i++){
        // Assign value of index [i] in new coefficient array = index [i+1] in '_coefficients' array
        polyDer._coefficients[i] = _coefficients[i + 1];
        // Multiply value of [i] by (counter 'i' + 1)
        polyDer._coefficients[i] *= (i + 1);
    }
    return polyDer;
}

/**
* Evaluate
* Returns value of evaluating polynomial with value x
* @param element the value you want to substitute for 'x'
* @returns 'solution', the value evaluated from substituting x
*/
float Polynomial::Evaluate(float x)const{

    /*
     * Concept:
     * f(3) = 5x^3 can be solved via 5 * (3)^3
     * ... which can be written as the coefficient (c) multiplied by the input (x) to the power of exponent (e)
     *
     * Equation: c * (x^e)
     *
     * Therefore, f(x) can be evaluated algorithmically through...
      * 1. [ x ] Substitute input for 'x'
      *    [ x^e ] Solve individually for a temporary array of x^exponent
      * 2. [ c * (x^e ] Multiply individual values of '_coefficients' to corresponding values in 'tempExpo' array
      * 3. Returning float 'solution' with final answer
    */

    double inputVal = x; // Input value
    double solution = 0; // Return value

    // Values of exponents to the power of 'z'
    double tempExpo [_degree];

    size_t z, y; // Loop counters

    // Solve individually for 'x' to the power of 'z'
    for(z = 0; z < _degree + 1; z++)
        tempExpo[z] = pow(inputVal, z);

    // Solve for coefficients * 'tempExpo' (or x^z)
    for(y = 0; y < _degree + 1; y++)
        solution += (_coefficients[y] * tempExpo[y]);

    return solution; // return value
}

/**
* Integrate
* Returns value of definite integral evaluated from a to b
* @param element the start and end values for the definite integral
* @returns 'solution', the value of the evaluated definite integral
*/
float Polynomial::Integrate(float start, float end)const{

    /*
     * Concept:
     * You can solve for a definite integral by substituting the start and end
     * values in the antiderivative of the original polynomial, then subtracting.
     *
     * Fundamental Theorem of Calculus...
     * f(x) dx = F(b) - F(a)
     *
      * 1. Find antiderivative of polynomial
      * 2. F(a) -> Find antiderivative
      *           -> Substitute 'startVal' into antiderivative
      *             -> Store result in 'tempStart'
      *           -> Multiply 'tempStart[y]' with corresponding array values in '_coefficients[y]'
      *         -> Assign product to 'startValResult' array (this holds F(a))
      *    F(b) -> Find antiderivative
      *           -> Same as F(a) above, processed seperately
      * 3. Solution = F(b) - F(a)
      * 4. Return 'solution'
      *
      * Note* The algorithms in this method use similar concepts from the Evaluate method...
     */

    // Input start/end values
    double startVal = start;
    double endVal = end;

    double solution; // Return value, F(b) - F(a)

    // Result holders for F(a) and F(b)
    double startValResult = 0;
    double endValResult = 0;

    // Temp. arrays used to store values of antiderivative with user input evaluations
    double tempStart [_degree];
    double tempEnd [_degree];

    size_t z, y; // Loop counters

    // --- F(a):
    // Find antiderivative, substitute 'startVal'
    for(z = 0; z < _degree + 1; z++)
        tempStart[z] = pow(startVal, z + 1)/(z + 1);

    // -- Multiply coeffecients by evaluation directly above
    for(y = 0; y < _degree + 1; y++)
        startValResult += (_coefficients[y] * tempStart[y]);

    // --- F(b):
    // Find antiderivative, substitute 'endVal'
    for(z = 0; z < _degree + 1; z++)
        tempEnd[z] = pow(endVal, z + 1)/(z + 1);

    // -- Multiply coeffecients by evaluation directly above
    for(y = 0; y < _degree + 1; y++)
        endValResult += (_coefficients[y] * tempEnd[y]);

    // --- F(b) - F(a):
    solution = endValResult - startValResult;

    return solution; // return value
}

/**
* Operator =
* Assigns polynomial rhs to this
* @param -
* @returns -
*/
const Polynomial& Polynomial::operator=(const Polynomial& rhs){
	if (&rhs == this){
		return *this;
	}
	if (_degree != rhs._degree){
		if (_coefficients){
			delete[] _coefficients;
		}
		_degree = rhs._degree;
		_coefficients = new float[_degree + 1];
	}
	for (size_t i = 0; i < _degree + 1; i++) {
		_coefficients[i] = rhs._coefficients[i];
	}
	return *this;
}

/**
* Equals
* Returns true if this is equal to rhs
* @param -
* @returns -
*/
bool Polynomial::Equals(const Polynomial& rhs)const{
	if (_degree != rhs._degree){
		return false;
	}
	for (size_t i=0; i < _degree; i++){
		if (abs(_coefficients[i] - rhs._coefficients[i]) > 0.0001){
			return false;
		}
	}
	return true;
}

/**
* ToString
* Returns a string representation of the polynomial
* @param -
* @returns -
*/
string Polynomial::ToString()const{
	stringstream ss;
	for (size_t i = _degree; i > 0; i--) {
		ss << showpos << fixed << setprecision(2) << _coefficients[i] << "x^" << i << " ";
	}
	ss << showpos << fixed << setprecision(2) << _coefficients[0];
	return ss.str();
}

/**
* Write
* Writes a polynomial from the output stream
* @param -
* @returns -
*/
ostream& Polynomial::Write(ostream& output)const{
	output << _degree << " ";
	for (size_t i = 0; i < _degree + 1; i++) {
		output << _coefficients[i] << " ";
	}
	return output;
}

/**
* Read
* Reads a polynomial from the input stream
* @param -
* @returns -
*/
istream& Polynomial::Read(istream& input){
	size_t degree;
	input >> degree;
	if (input.fail()){
		return input;
	}
	float* coefficients = new float[degree + 1];
	for (size_t i = 0; i < degree + 1; i++) {
		input >> coefficients[i];
		if (input.fail()){
			delete[] coefficients;
			return input;
		}
	}

	if (degree != _degree){
		if (_coefficients){
			delete[] _coefficients;
		}
		_degree = degree;
		_coefficients = coefficients;
	}else{
		for (size_t i = 0; i < _degree + 1; i++) {
			_coefficients[i] = coefficients[i];
		}
		delete[] coefficients;
	}
	return input;
}