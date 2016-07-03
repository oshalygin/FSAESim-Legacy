


//Start of Rounding Function

function roundNumber(num, dec) {
    var result = Math.round(num * Math.pow(10, dec)) / Math.pow(10, dec);
    return result;
}

//End of Rounding Function

//Start of Matrix Multiplication Function

function matrixMultiplication(a, b) {

	var productMatrix = [];
	var rowsa;
	var colsa;
	var rowsb;
	var colsb;


	if (a instanceof Array) {

		if (a[0].length) {
			var tempArrayA = a;
			var tempArrayB = b;
			rowsa = a[0].length;
			colsa = a.length;
			rowsb = b[0].length;
			colsb = b.length;
			var sum = 0;

			//Building new Array
			var tempArray = new Array(colsb);
			for (var i = 0; i < colsb; i++) {
				tempArray[i] = new Array(rowsa);
			}

			for (i = 0; i < rowsa; i++) {
				for (var k = 0; k < colsb; k++) {
					for (var l = 0; l < rowsb; l++) {
						a = tempArrayA[l][i];
						b = tempArrayB[k][l];
						sum += a * b;

					}
					tempArray[k][i] = sum;
					sum = 0;

				}
			}



			productMatrix = tempArray;
			return productMatrix;

		}
	}



	else {
		tempArray = b;
		//Scalar Multiplication for Large Array
		if (b[0].length) {
			rowsb = b[0].length;
			colsb = b.length;
			for (i = 0; i < colsb; i++) {
				for (var j = 0; j < rowsb; j++) {
					tempArray[i][j] = a * tempArray[i][j];
				}
			}
			productMatrix = tempArray;
			return productMatrix;

		}
		//Scalar Multiplication for 1xn array
		else {
			for (i = 0; i < b.length; i++) {
				tempArray[i] = a * b[i];

			}
			productMatrix = tempArray;
			return productMatrix;
		}
	}
}

//End of Matrix Multiplication Function

//Start of Mean Value Function
function meanValue(arr) {
	var mean;
	var sum = 0;

	if (arr[0].length) {
		for (var j = 0; j < arr[0].length; j++) {
			sum += arr[0][j];
		}
		mean = sum / arr[0].length;

	}
	else {
		for (var i = 0; i < arr.length; i++) {
			sum += arr[i];
		}
		mean = sum / arr.length;


	}

	return mean;


}

//Start of Standard Deviation Function
function standardDeviation(arr) {
	var mean = meanValue(arr);
	var difference;
	var sum = 0;

	var stdev;

	if (arr[0].length) {
		for (var j = 0; j < arr[0].length; j++) {
			difference = (arr[0][j] - mean);
			sum += Math.pow(difference, 2);


		}
		stdev = (sum) / (arr[0].length - 1);
		stdev = Math.sqrt(stdev);
	}
	else {
		for (var i = 0; i < arr.length; i++) {
			difference = (arr[i] - mean);
			sum += Math.pow(difference, 2);
		}
		stdev = (sum) / (arr.length - 1);
		stdev = Math.sqrt(stdev);
	}



	return stdev;

}



//End of Standard Deviation Function


//Start of Matrix Transpose Function
function matrixTranspose(arr) {
	var transposedArray = [];
	var rows;
	var cols;
	if (arr.length <= 0) {
		return transposedArray;
	}


	//Builds the Array
	if (arr[0].length) {
		cols = arr.length;
		rows = arr[0].length;
	}
	else {
		rows = 0;
		cols = arr.length;


		var tempArray = new Array(1);
		tempArray[0] = new Array(cols);

		for (i = 0; i < arr.length; i++) {
			tempArray[0][i] = arr[i];
		}
		transposedArray = tempArray;
		return transposedArray;


	}


	tempArray = new Array(rows);
	for (i = 0; i < tempArray.length; i++) {
		tempArray[i] = new Array(cols);
	}
	transposedArray = tempArray;

	for (i = 0; i < tempArray.length; i++) {
		for (j = 0; j < tempArray[i].length; j++) {

			transposedArray[i][j] = arr[j][i];
		}
	}


	return transposedArray;
}

//End of Matrix Transpose Function

//Start of Matrix Dot Product
function dotProduct(vectorX, vectorY, n) {

	var sum = 0;
	for (var i = 0; i < n; i++) {
		sum += vectorX[i] * vectorY[i];
	}
	return sum;

}

//End of Matrix Dot Product

//Start of Matrix Norm
function norm(vectorX, n) {
	var normValue;
	normValue = Math.sqrt(dotProduct(vectorX, vectorX, n));
	return normValue;

}

//End of Matrix Norm

//Start of Matrix Scalar Multiplication
function scalarMatrixMultiplication(scalarValue, vectorX, n) {

	var result = new Array(n);
	for (var i = 0; i < n; i++) {
		result[i] = scalarValue * vectorX[i];
	}
	return result;

}

//End of Matrix Scalar Multiplication


//Start of Vector Subtraction
function vectorSubtraction(vectorX, vectorY, n) {

	var result = new Array(n);
	for (var i = 0; i < n; i++) {
		result[i] = vectorX[i] - vectorY[i];
	}
	return result;


}
//End of Vector Subtraction


//Start of Zero's Function(Fills all elements in the MxN matrix with zero's)
function zeros(m, n) {

	var zeroArray = new Array(n);
	for (var i = 0; i < n; i++) {
		zeroArray[i] = new Array(m);
	}

	for (var j = 0; j < n; j++) {
		for (var k = 0; k < m; k++) {
			zeroArray[j][k] = 0;
		}

	}


	return zeroArray;

}

//End of Zero's Function(Fills all elements in the MxN matrix with zero's)

//Start of Ones Matrix
function ones(m, n) {
	//Constructs a matrix that is MxN

	//Contruct matrixOnes
	var matrixOnes = new Array(n);
	for (var k = 0; k < n; k++) {
		matrixOnes[k] = new Array(m)
	}


	//Assigns the Value of 1 to each element in the array.
	for (var i = 0; i < n; i++) {
		for (var j = 0; j < m; j++) {
			matrixOnes[i][j] = 1;

		}

	}
	return matrixOnes;

}


//End of Ones Matrix


//Start of Vector Matrix
function vector(matrixInput, number, rowCol) {
	//rowCol=1 means that the function will determine the row vector.  If rowCol=2 then the function will determine the column vector.
	var m = matrixInput[0].length;
	var n = matrixInput.length;
	number = number;

	if (rowCol === 1) {
		var vectorOutput = new Array(n);

		for (var i = 0; i < matrixInput.length; i++) {
			vectorOutput[i] = matrixInput[i][number];
		}
	}
	else {
		var vectorOutput = new Array(m);

		for (i = 0; i < matrixInput[0].length; i++) {
			vectorOutput[i] = matrixInput[number][i];

		}


	}

	return vectorOutput;

}
//End of Vector Matrix

//Start of Matrix Library by Sylvester
// === Sylvester ===
// Vector and Matrix mathematics modules for JavaScript
// Copyright (c) 2007 James Coglan
// 
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

function Vector() { }
Vector.prototype = {

	// Returns element i of the vector
	e: function (i) {
		return (i < 1 || i > this.elements.length) ? null : this.elements[i - 1];
	},

	// Returns the number of elements the vector has
	dimensions: function () {
		return this.elements.length;
	},

	// Returns the modulus ('length') of the vector
	modulus: function () {
		return Math.sqrt(this.dot(this));
	},

	// Returns true iff the vector is equal to the argument
	eql: function (vector) {
		var n = this.elements.length;
		var V = vector.elements || vector;
		if (n != V.length) { return false; }
		do {
			if (Math.abs(this.elements[n - 1] - V[n - 1]) > Sylvester.precision) { return false; }
		} while (--n);
		return true;
	},

	// Returns a copy of the vector
	dup: function () {
		return Vector.create(this.elements);
	},

	// Maps the vector to another vector according to the given function
	map: function (fn) {
		var elements = [];
		this.each(function (x, i) {
			elements.push(fn(x, i));
		});
		return Vector.create(elements);
	},

	// Calls the iterator for each element of the vector in turn
	each: function (fn) {
		var n = this.elements.length, k = n, i;
		do {
			i = k - n;
			fn(this.elements[i], i + 1);
		} while (--n);
	},

	// Returns a new vector created by normalizing the receiver
	toUnitVector: function () {
		var r = this.modulus();
		if (r === 0) { return this.dup(); }
		return this.map(function (x) { return x / r; });
	},

	// Returns the angle between the vector and the argument (also a vector)
	angleFrom: function (vector) {
		var V = vector.elements || vector;
		var n = this.elements.length, k = n, i;
		if (n != V.length) { return null; }
		var dot = 0, mod1 = 0, mod2 = 0;
		// Work things out in parallel to save time
		this.each(function (x, i) {
			dot += x * V[i - 1];
			mod1 += x * x;
			mod2 += V[i - 1] * V[i - 1];
		});
		mod1 = Math.sqrt(mod1); mod2 = Math.sqrt(mod2);
		if (mod1 * mod2 === 0) { return null; }
		var theta = dot / (mod1 * mod2);
		if (theta < -1) { theta = -1; }
		if (theta > 1) { theta = 1; }
		return Math.acos(theta);
	},

	// Returns true iff the vector is parallel to the argument
	isParallelTo: function (vector) {
		var angle = this.angleFrom(vector);
		return (angle === null) ? null : (angle <= Sylvester.precision);
	},

	// Returns true iff the vector is antiparallel to the argument
	isAntiparallelTo: function (vector) {
		var angle = this.angleFrom(vector);
		return (angle === null) ? null : (Math.abs(angle - Math.PI) <= Sylvester.precision);
	},

	// Returns true iff the vector is perpendicular to the argument
	isPerpendicularTo: function (vector) {
		var dot = this.dot(vector);
		return (dot === null) ? null : (Math.abs(dot) <= Sylvester.precision);
	},

	// Returns the result of adding the argument to the vector
	add: function (vector) {
		var V = vector.elements || vector;
		if (this.elements.length != V.length) { return null; }
		return this.map(function (x, i) { return x + V[i - 1]; });
	},

	// Returns the result of subtracting the argument from the vector
	subtract: function (vector) {
		var V = vector.elements || vector;
		if (this.elements.length != V.length) { return null; }
		return this.map(function (x, i) { return x - V[i - 1]; });
	},

	// Returns the result of multiplying the elements of the vector by the argument
	multiply: function (k) {
		return this.map(function (x) { return x * k; });
	},

	x: function (k) { return this.multiply(k); },

	// Returns the scalar product of the vector with the argument
	// Both vectors must have equal dimensionality
	dot: function (vector) {
		var V = vector.elements || vector;
		var i, product = 0, n = this.elements.length;
		if (n != V.length) { return null; }
		do { product += this.elements[n - 1] * V[n - 1]; } while (--n);
		return product;
	},

	// Returns the vector product of the vector with the argument
	// Both vectors must have dimensionality 3
	cross: function (vector) {
		var B = vector.elements || vector;
		if (this.elements.length != 3 || B.length != 3) { return null; }
		var A = this.elements;
		return Vector.create([
      (A[1] * B[2]) - (A[2] * B[1]),
      (A[2] * B[0]) - (A[0] * B[2]),
      (A[0] * B[1]) - (A[1] * B[0])
    ]);
	},

	// Returns the (absolute) largest element of the vector
	max: function () {
		var m = 0, n = this.elements.length, k = n, i;
		do {
			i = k - n;
			if (Math.abs(this.elements[i]) > Math.abs(m)) { m = this.elements[i]; }
		} while (--n);
		return m;
	},

	// Returns the index of the first match found
	indexOf: function (x) {
		var index = null, n = this.elements.length, k = n, i;
		do {
			i = k - n;
			if (index === null && this.elements[i] == x) {
				index = i + 1;
			}
		} while (--n);
		return index;
	},

	// Returns a diagonal matrix with the vector's elements as its diagonal elements
	toDiagonalMatrix: function () {
		return Matrix.Diagonal(this.elements);
	},

	// Returns the result of rounding the elements of the vector
	round: function () {
		return this.map(function (x) { return Math.round(x); });
	},

	// Returns a copy of the vector with elements set to the given value if they
	// differ from it by less than Sylvester.precision
	snapTo: function (x) {
		return this.map(function (y) {
			return (Math.abs(y - x) <= Sylvester.precision) ? x : y;
		});
	},

	// Returns the vector's distance from the argument, when considered as a point in space
	distanceFrom: function (obj) {
		if (obj.anchor) { return obj.distanceFrom(this); }
		var V = obj.elements || obj;
		if (V.length != this.elements.length) { return null; }
		var sum = 0, part;
		this.each(function (x, i) {
			part = x - V[i - 1];
			sum += part * part;
		});
		return Math.sqrt(sum);
	},

	// Returns true if the vector is point on the given line
	liesOn: function (line) {
		return line.contains(this);
	},

	// Return true iff the vector is a point in the given plane
	liesIn: function (plane) {
		return plane.contains(this);
	},

	// Rotates the vector about the given object. The object should be a 
	// point if the vector is 2D, and a line if it is 3D. Be careful with line directions!
	rotate: function (t, obj) {
		var V, R, x, y, z;
		switch (this.elements.length) {
			case 2:
				V = obj.elements || obj;
				if (V.length != 2) { return null; }
				R = Matrix.Rotation(t).elements;
				x = this.elements[0] - V[0];
				y = this.elements[1] - V[1];
				return Vector.create([
          V[0] + R[0][0] * x + R[0][1] * y,
          V[1] + R[1][0] * x + R[1][1] * y
        ]);
				break;
			case 3:
				if (!obj.direction) { return null; }
				var C = obj.pointClosestTo(this).elements;
				R = Matrix.Rotation(t, obj.direction).elements;
				x = this.elements[0] - C[0];
				y = this.elements[1] - C[1];
				z = this.elements[2] - C[2];
				return Vector.create([
          C[0] + R[0][0] * x + R[0][1] * y + R[0][2] * z,
          C[1] + R[1][0] * x + R[1][1] * y + R[1][2] * z,
          C[2] + R[2][0] * x + R[2][1] * y + R[2][2] * z
        ]);
				break;
			default:
				return null;
		}
	},

	// Returns the result of reflecting the point in the given point, line or plane
	reflectionIn: function (obj) {
		if (obj.anchor) {
			// obj is a plane or line
			var P = this.elements.slice();
			var C = obj.pointClosestTo(P).elements;
			return Vector.create([C[0] + (C[0] - P[0]), C[1] + (C[1] - P[1]), C[2] + (C[2] - (P[2] || 0))]);
		} else {
			// obj is a point
			var Q = obj.elements || obj;
			if (this.elements.length != Q.length) { return null; }
			return this.map(function (x, i) { return Q[i - 1] + (Q[i - 1] - x); });
		}
	},

	// Utility to make sure vectors are 3D. If they are 2D, a zero z-component is added
	to3D: function () {
		var V = this.dup();
		switch (V.elements.length) {
			case 3: break;
			case 2: V.elements.push(0); break;
			default: return null;
		}
		return V;
	},

	// Returns a string representation of the vector
	inspect: function () {
		return '[' + this.elements.join(', ') + ']';
	},

	// Set vector's elements from an array
	setElements: function (els) {
		this.elements = (els.elements || els).slice();
		return this;
	}
};

// Constructor function
Vector.create = function (elements) {
	var V = new Vector();
	return V.setElements(elements);
};

// i, j, k unit vectors
Vector.i = Vector.create([1, 0, 0]);
Vector.j = Vector.create([0, 1, 0]);
Vector.k = Vector.create([0, 0, 1]);

// Random vector of size n
Vector.Random = function (n) {
	var elements = [];
	do {
		elements.push(Math.random());
	} while (--n);
	return Vector.create(elements);
};

// Vector filled with zeros
Vector.Zero = function (n) {
	var elements = [];
	do {
		elements.push(0);
	} while (--n);
	return Vector.create(elements);
};



function Matrix() { }
Matrix.prototype = {

	// Returns element (i,j) of the matrix
	e: function (i, j) {
		if (i < 1 || i > this.elements.length || j < 1 || j > this.elements[0].length) { return null; }
		return this.elements[i - 1][j - 1];
	},

	// Returns row k of the matrix as a vector
	row: function (i) {
		if (i > this.elements.length) { return null; }
		return Vector.create(this.elements[i - 1]);
	},

	// Returns column k of the matrix as a vector
	col: function (j) {
		if (j > this.elements[0].length) { return null; }
		var col = [], n = this.elements.length, k = n, i;
		do {
			i = k - n;
			col.push(this.elements[i][j - 1]);
		} while (--n);
		return Vector.create(col);
	},

	// Returns the number of rows/columns the matrix has
	dimensions: function () {
		return { rows: this.elements.length, cols: this.elements[0].length };
	},

	// Returns the number of rows in the matrix
	rows: function () {
		return this.elements.length;
	},

	// Returns the number of columns in the matrix
	cols: function () {
		return this.elements[0].length;
	},

	// Returns true iff the matrix is equal to the argument. You can supply
	// a vector as the argument, in which case the receiver must be a
	// one-column matrix equal to the vector.
	eql: function (matrix) {
		var M = matrix.elements || matrix;
		if (typeof (M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
		if (this.elements.length != M.length ||
        this.elements[0].length != M[0].length) { return false; }
		var ni = this.elements.length, ki = ni, i, nj, kj = this.elements[0].length, j;
		do {
			i = ki - ni;
			nj = kj;
			do {
				j = kj - nj;
				if (Math.abs(this.elements[i][j] - M[i][j]) > Sylvester.precision) { return false; }
			} while (--nj);
		} while (--ni);
		return true;
	},

	// Returns a copy of the matrix
	dup: function () {
		return Matrix.create(this.elements);
	},

	// Maps the matrix to another matrix (of the same dimensions) according to the given function
	map: function (fn) {
		var els = [], ni = this.elements.length, ki = ni, i, nj, kj = this.elements[0].length, j;
		do {
			i = ki - ni;
			nj = kj;
			els[i] = [];
			do {
				j = kj - nj;
				els[i][j] = fn(this.elements[i][j], i + 1, j + 1);
			} while (--nj);
		} while (--ni);
		return Matrix.create(els);
	},

	// Returns true iff the argument has the same dimensions as the matrix
	isSameSizeAs: function (matrix) {
		var M = matrix.elements || matrix;
		if (typeof (M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
		return (this.elements.length == M.length &&
        this.elements[0].length == M[0].length);
	},

	// Returns the result of adding the argument to the matrix
	add: function (matrix) {
		var M = matrix.elements || matrix;
		if (typeof (M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
		if (!this.isSameSizeAs(M)) { return null; }
		return this.map(function (x, i, j) { return x + M[i - 1][j - 1]; });
	},

	// Returns the result of subtracting the argument from the matrix
	subtract: function (matrix) {
		var M = matrix.elements || matrix;
		if (typeof (M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
		if (!this.isSameSizeAs(M)) { return null; }
		return this.map(function (x, i, j) { return x - M[i - 1][j - 1]; });
	},

	// Returns true iff the matrix can multiply the argument from the left
	canMultiplyFromLeft: function (matrix) {
		var M = matrix.elements || matrix;
		if (typeof (M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
		// this.columns should equal matrix.rows
		return (this.elements[0].length == M.length);
	},

	// Returns the result of multiplying the matrix from the right by the argument.
	// If the argument is a scalar then just multiply all the elements. If the argument is
	// a vector, a vector is returned, which saves you having to remember calling
	// col(1) on the result.
	multiply: function (matrix) {
		if (!matrix.elements) {
			return this.map(function (x) { return x * matrix; });
		}
		var returnVector = matrix.modulus ? true : false;
		var M = matrix.elements || matrix;
		if (typeof (M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
		if (!this.canMultiplyFromLeft(M)) { return null; }
		var ni = this.elements.length, ki = ni, i, nj, kj = M[0].length, j;
		var cols = this.elements[0].length, elements = [], sum, nc, c;
		do {
			i = ki - ni;
			elements[i] = [];
			nj = kj;
			do {
				j = kj - nj;
				sum = 0;
				nc = cols;
				do {
					c = cols - nc;
					sum += this.elements[i][c] * M[c][j];
				} while (--nc);
				elements[i][j] = sum;
			} while (--nj);
		} while (--ni);
		var M = Matrix.create(elements);
		return returnVector ? M.col(1) : M;
	},

	x: function (matrix) { return this.multiply(matrix); },

	// Returns a submatrix taken from the matrix
	// Argument order is: start row, start col, nrows, ncols
	// Element selection wraps if the required index is outside the matrix's bounds, so you could
	// use this to perform row/column cycling or copy-augmenting.
	minor: function (a, b, c, d) {
		var elements = [], ni = c, i, nj, j;
		var rows = this.elements.length, cols = this.elements[0].length;
		do {
			i = c - ni;
			elements[i] = [];
			nj = d;
			do {
				j = d - nj;
				elements[i][j] = this.elements[(a + i - 1) % rows][(b + j - 1) % cols];
			} while (--nj);
		} while (--ni);
		return Matrix.create(elements);
	},

	// Returns the transpose of the matrix
	transpose: function () {
		var rows = this.elements.length, cols = this.elements[0].length;
		var elements = [], ni = cols, i, nj, j;
		do {
			i = cols - ni;
			elements[i] = [];
			nj = rows;
			do {
				j = rows - nj;
				elements[i][j] = this.elements[j][i];
			} while (--nj);
		} while (--ni);
		return Matrix.create(elements);
	},

	// Returns true iff the matrix is square
	isSquare: function () {
		return (this.elements.length == this.elements[0].length);
	},

	// Returns the (absolute) largest element of the matrix
	max: function () {
		var m = 0, ni = this.elements.length, ki = ni, i, nj, kj = this.elements[0].length, j;
		do {
			i = ki - ni;
			nj = kj;
			do {
				j = kj - nj;
				if (Math.abs(this.elements[i][j]) > Math.abs(m)) { m = this.elements[i][j]; }
			} while (--nj);
		} while (--ni);
		return m;
	},

	// Returns the indeces of the first match found by reading row-by-row from left to right
	indexOf: function (x) {
		var index = null, ni = this.elements.length, ki = ni, i, nj, kj = this.elements[0].length, j;
		do {
			i = ki - ni;
			nj = kj;
			do {
				j = kj - nj;
				if (this.elements[i][j] == x) { return { i: i + 1, j: j + 1 }; }
			} while (--nj);
		} while (--ni);
		return null;
	},

	// If the matrix is square, returns the diagonal elements as a vector.
	// Otherwise, returns null.
	diagonal: function () {
		if (!this.isSquare) { return null; }
		var els = [], n = this.elements.length, k = n, i;
		do {
			i = k - n;
			els.push(this.elements[i][i]);
		} while (--n);
		return Vector.create(els);
	},

	// Make the matrix upper (right) triangular by Gaussian elimination.
	// This method only adds multiples of rows to other rows. No rows are
	// scaled up or switched, and the determinant is preserved.
	toRightTriangular: function () {
		var M = this.dup(), els;
		var n = this.elements.length, k = n, i, np, kp = this.elements[0].length, p;
		do {
			i = k - n;
			if (M.elements[i][i] == 0) {
				for (j = i + 1; j < k; j++) {
					if (M.elements[j][i] != 0) {
						els = []; np = kp;
						do {
							p = kp - np;
							els.push(M.elements[i][p] + M.elements[j][p]);
						} while (--np);
						M.elements[i] = els;
						break;
					}
				}
			}
			if (M.elements[i][i] != 0) {
				for (j = i + 1; j < k; j++) {
					var multiplier = M.elements[j][i] / M.elements[i][i];
					els = []; np = kp;
					do {
						p = kp - np;
						// Elements with column numbers up to an including the number
						// of the row that we're subtracting can safely be set straight to
						// zero, since that's the point of this routine and it avoids having
						// to loop over and correct rounding errors later
						els.push(p <= i ? 0 : M.elements[j][p] - M.elements[i][p] * multiplier);
					} while (--np);
					M.elements[j] = els;
				}
			}
		} while (--n);
		return M;
	},

	toUpperTriangular: function () { return this.toRightTriangular(); },

	// Returns the determinant for square matrices
	determinant: function () {
		if (!this.isSquare()) { return null; }
		var M = this.toRightTriangular();
		var det = M.elements[0][0], n = M.elements.length - 1, k = n, i;
		do {
			i = k - n + 1;
			det = det * M.elements[i][i];
		} while (--n);
		return det;
	},

	det: function () { return this.determinant(); },

	// Returns true iff the matrix is singular
	isSingular: function () {
		return (this.isSquare() && this.determinant() === 0);
	},

	// Returns the trace for square matrices
	trace: function () {
		if (!this.isSquare()) { return null; }
		var tr = this.elements[0][0], n = this.elements.length - 1, k = n, i;
		do {
			i = k - n + 1;
			tr += this.elements[i][i];
		} while (--n);
		return tr;
	},

	tr: function () { return this.trace(); },

	// Returns the rank of the matrix
	rank: function () {
		var M = this.toRightTriangular(), rank = 0;
		var ni = this.elements.length, ki = ni, i, nj, kj = this.elements[0].length, j;
		do {
			i = ki - ni;
			nj = kj;
			do {
				j = kj - nj;
				if (Math.abs(M.elements[i][j]) > Sylvester.precision) { rank++; break; }
			} while (--nj);
		} while (--ni);
		return rank;
	},

	rk: function () { return this.rank(); },

	// Returns the result of attaching the given argument to the right-hand side of the matrix
	augment: function (matrix) {
		var M = matrix.elements || matrix;
		if (typeof (M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
		var T = this.dup(), cols = T.elements[0].length;
		var ni = T.elements.length, ki = ni, i, nj, kj = M[0].length, j;
		if (ni != M.length) { return null; }
		do {
			i = ki - ni;
			nj = kj;
			do {
				j = kj - nj;
				T.elements[i][cols + j] = M[i][j];
			} while (--nj);
		} while (--ni);
		return T;
	},

	// Returns the inverse (if one exists) using Gauss-Jordan
	inverse: function () {
		if (!this.isSquare() || this.isSingular()) { return null; }
		var ni = this.elements.length, ki = ni, i, j;
		var M = this.augment(Matrix.I(ni)).toRightTriangular();
		var np, kp = M.elements[0].length, p, els, divisor;
		var inverse_elements = [], new_element;
		// Matrix is non-singular so there will be no zeros on the diagonal
		// Cycle through rows from last to first
		do {
			i = ni - 1;
			// First, normalise diagonal elements to 1
			els = []; np = kp;
			inverse_elements[i] = [];
			divisor = M.elements[i][i];
			do {
				p = kp - np;
				new_element = M.elements[i][p] / divisor;
				els.push(new_element);
				// Shuffle of the current row of the right hand side into the results
				// array as it will not be modified by later runs through this loop
				if (p >= ki) { inverse_elements[i].push(new_element); }
			} while (--np);
			M.elements[i] = els;
			// Then, subtract this row from those above it to
			// give the identity matrix on the left hand side
			for (j = 0; j < i; j++) {
				els = []; np = kp;
				do {
					p = kp - np;
					els.push(M.elements[j][p] - M.elements[i][p] * M.elements[j][i]);
				} while (--np);
				M.elements[j] = els;
			}
		} while (--ni);
		return Matrix.create(inverse_elements);
	},

	inv: function () { return this.inverse(); },

	// Returns the result of rounding all the elements
	round: function () {
		return this.map(function (x) { return Math.round(x); });
	},

	// Returns a copy of the matrix with elements set to the given value if they
	// differ from it by less than Sylvester.precision
	snapTo: function (x) {
		return this.map(function (p) {
			return (Math.abs(p - x) <= Sylvester.precision) ? x : p;
		});
	},

	// Returns a string representation of the matrix
	inspect: function () {
		var matrix_rows = [];
		var n = this.elements.length, k = n, i;
		do {
			i = k - n;
			matrix_rows.push(Vector.create(this.elements[i]).inspect());
		} while (--n);
		return matrix_rows.join('\n');
	},

	// Set the matrix's elements from an array. If the argument passed
	// is a vector, the resulting matrix will be a single column.
	setElements: function (els) {
		var i, elements = els.elements || els;
		if (typeof (elements[0][0]) != 'undefined') {
			var ni = elements.length, ki = ni, nj, kj, j;
			this.elements = [];
			do {
				i = ki - ni;
				nj = elements[i].length; kj = nj;
				this.elements[i] = [];
				do {
					j = kj - nj;
					this.elements[i][j] = elements[i][j];
				} while (--nj);
			} while (--ni);
			return this;
		}
		var n = elements.length, k = n;
		this.elements = [];
		do {
			i = k - n;
			this.elements.push([elements[i]]);
		} while (--n);
		return this;
	}
};

// Constructor function
Matrix.create = function (elements) {
	var M = new Matrix();
	return M.setElements(elements);
};

// Identity matrix of size n
Matrix.I = function (n) {
	var els = [], k = n, i, nj, j;
	do {
		i = k - n;
		els[i] = []; nj = k;
		do {
			j = k - nj;
			els[i][j] = (i == j) ? 1 : 0;
		} while (--nj);
	} while (--n);
	return Matrix.create(els);
};

// Diagonal matrix - all off-diagonal elements are zero
Matrix.Diagonal = function (elements) {
	var n = elements.length, k = n, i;
	var M = Matrix.I(n);
	do {
		i = k - n;
		M.elements[i][i] = elements[i];
	} while (--n);
	return M;
};

// Rotation matrix about some axis. If no axis is
// supplied, assume we're after a 2D transform
Matrix.Rotation = function (theta, a) {
	if (!a) {
		return Matrix.create([
      [Math.cos(theta), -Math.sin(theta)],
      [Math.sin(theta), Math.cos(theta)]
    ]);
	}
	var axis = a.dup();
	if (axis.elements.length != 3) { return null; }
	var mod = axis.modulus();
	var x = axis.elements[0] / mod, y = axis.elements[1] / mod, z = axis.elements[2] / mod;
	var s = Math.sin(theta), c = Math.cos(theta), t = 1 - c;
	// Formula derived here: http://www.gamedev.net/reference/articles/article1199.asp
	// That proof rotates the co-ordinate system so theta
	// becomes -theta and sin becomes -sin here.
	return Matrix.create([
    [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
    [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
    [t * x * z - s * y, t * y * z + s * x, t * z * z + c]
  ]);
};

// Special case rotations
Matrix.RotationX = function (t) {
	var c = Math.cos(t), s = Math.sin(t);
	return Matrix.create([
    [1, 0, 0],
    [0, c, -s],
    [0, s, c]
  ]);
};
Matrix.RotationY = function (t) {
	var c = Math.cos(t), s = Math.sin(t);
	return Matrix.create([
    [c, 0, s],
    [0, 1, 0],
    [-s, 0, c]
  ]);
};
Matrix.RotationZ = function (t) {
	var c = Math.cos(t), s = Math.sin(t);
	return Matrix.create([
    [c, -s, 0],
    [s, c, 0],
    [0, 0, 1]
  ]);
};

// Random matrix of n rows, m columns
Matrix.Random = function (n, m) {
	return Matrix.Zero(n, m).map(
    function () { return Math.random(); }
  );
};

// Matrix filled with zeros
Matrix.Zero = function (n, m) {
	var els = [], ni = n, i, nj, j;
	do {
		i = n - ni;
		els[i] = [];
		nj = m;
		do {
			j = m - nj;
			els[i][j] = 0;
		} while (--nj);
	} while (--ni);
	return Matrix.create(els);
};



function Line() { }
Line.prototype = {

	// Returns true if the argument occupies the same space as the line
	eql: function (line) {
		return (this.isParallelTo(line) && this.contains(line.anchor));
	},

	// Returns a copy of the line
	dup: function () {
		return Line.create(this.anchor, this.direction);
	},

	// Returns the result of translating the line by the given vector/array
	translate: function (vector) {
		var V = vector.elements || vector;
		return Line.create([
      this.anchor.elements[0] + V[0],
      this.anchor.elements[1] + V[1],
      this.anchor.elements[2] + (V[2] || 0)
    ], this.direction);
	},

	// Returns true if the line is parallel to the argument. Here, 'parallel to'
	// means that the argument's direction is either parallel or antiparallel to
	// the line's own direction. A line is parallel to a plane if the two do not
	// have a unique intersection.
	isParallelTo: function (obj) {
		if (obj.normal) { return obj.isParallelTo(this); }
		var theta = this.direction.angleFrom(obj.direction);
		return (Math.abs(theta) <= Sylvester.precision || Math.abs(theta - Math.PI) <= Sylvester.precision);
	},

	// Returns the line's perpendicular distance from the argument,
	// which can be a point, a line or a plane
	distanceFrom: function (obj) {
		if (obj.normal) { return obj.distanceFrom(this); }
		if (obj.direction) {
			// obj is a line
			if (this.isParallelTo(obj)) { return this.distanceFrom(obj.anchor); }
			var N = this.direction.cross(obj.direction).toUnitVector().elements;
			var A = this.anchor.elements, B = obj.anchor.elements;
			return Math.abs((A[0] - B[0]) * N[0] + (A[1] - B[1]) * N[1] + (A[2] - B[2]) * N[2]);
		} else {
			// obj is a point
			var P = obj.elements || obj;
			var A = this.anchor.elements, D = this.direction.elements;
			var PA1 = P[0] - A[0], PA2 = P[1] - A[1], PA3 = (P[2] || 0) - A[2];
			var modPA = Math.sqrt(PA1 * PA1 + PA2 * PA2 + PA3 * PA3);
			if (modPA === 0) return 0;
			// Assumes direction vector is normalized
			var cosTheta = (PA1 * D[0] + PA2 * D[1] + PA3 * D[2]) / modPA;
			var sin2 = 1 - cosTheta * cosTheta;
			return Math.abs(modPA * Math.sqrt(sin2 < 0 ? 0 : sin2));
		}
	},

	// Returns true iff the argument is a point on the line
	contains: function (point) {
		var dist = this.distanceFrom(point);
		return (dist !== null && dist <= Sylvester.precision);
	},

	// Returns true iff the line lies in the given plane
	liesIn: function (plane) {
		return plane.contains(this);
	},

	// Returns true iff the line has a unique point of intersection with the argument
	intersects: function (obj) {
		if (obj.normal) { return obj.intersects(this); }
		return (!this.isParallelTo(obj) && this.distanceFrom(obj) <= Sylvester.precision);
	},

	// Returns the unique intersection point with the argument, if one exists
	intersectionWith: function (obj) {
		if (obj.normal) { return obj.intersectionWith(this); }
		if (!this.intersects(obj)) { return null; }
		var P = this.anchor.elements, X = this.direction.elements,
        Q = obj.anchor.elements, Y = obj.direction.elements;
		var X1 = X[0], X2 = X[1], X3 = X[2], Y1 = Y[0], Y2 = Y[1], Y3 = Y[2];
		var PsubQ1 = P[0] - Q[0], PsubQ2 = P[1] - Q[1], PsubQ3 = P[2] - Q[2];
		var XdotQsubP = -X1 * PsubQ1 - X2 * PsubQ2 - X3 * PsubQ3;
		var YdotPsubQ = Y1 * PsubQ1 + Y2 * PsubQ2 + Y3 * PsubQ3;
		var XdotX = X1 * X1 + X2 * X2 + X3 * X3;
		var YdotY = Y1 * Y1 + Y2 * Y2 + Y3 * Y3;
		var XdotY = X1 * Y1 + X2 * Y2 + X3 * Y3;
		var k = (XdotQsubP * YdotY / XdotX + XdotY * YdotPsubQ) / (YdotY - XdotY * XdotY);
		return Vector.create([P[0] + k * X1, P[1] + k * X2, P[2] + k * X3]);
	},

	// Returns the point on the line that is closest to the given point or line
	pointClosestTo: function (obj) {
		if (obj.direction) {
			// obj is a line
			if (this.intersects(obj)) { return this.intersectionWith(obj); }
			if (this.isParallelTo(obj)) { return null; }
			var D = this.direction.elements, E = obj.direction.elements;
			var D1 = D[0], D2 = D[1], D3 = D[2], E1 = E[0], E2 = E[1], E3 = E[2];
			// Create plane containing obj and the shared normal and intersect this with it
			// Thank you: http://www.cgafaq.info/wiki/Line-line_distance
			var x = (D3 * E1 - D1 * E3), y = (D1 * E2 - D2 * E1), z = (D2 * E3 - D3 * E2);
			var N = Vector.create([x * E3 - y * E2, y * E1 - z * E3, z * E2 - x * E1]);
			var P = Plane.create(obj.anchor, N);
			return P.intersectionWith(this);
		} else {
			// obj is a point
			var P = obj.elements || obj;
			if (this.contains(P)) { return Vector.create(P); }
			var A = this.anchor.elements, D = this.direction.elements;
			var D1 = D[0], D2 = D[1], D3 = D[2], A1 = A[0], A2 = A[1], A3 = A[2];
			var x = D1 * (P[1] - A2) - D2 * (P[0] - A1), y = D2 * ((P[2] || 0) - A3) - D3 * (P[1] - A2),
          z = D3 * (P[0] - A1) - D1 * ((P[2] || 0) - A3);
			var V = Vector.create([D2 * x - D3 * z, D3 * y - D1 * x, D1 * z - D2 * y]);
			var k = this.distanceFrom(P) / V.modulus();
			return Vector.create([
        P[0] + V.elements[0] * k,
        P[1] + V.elements[1] * k,
        (P[2] || 0) + V.elements[2] * k
      ]);
		}
	},

	// Returns a copy of the line rotated by t radians about the given line. Works by
	// finding the argument's closest point to this line's anchor point (call this C) and
	// rotating the anchor about C. Also rotates the line's direction about the argument's.
	// Be careful with this - the rotation axis' direction affects the outcome!
	rotate: function (t, line) {
		// If we're working in 2D
		if (typeof (line.direction) == 'undefined') { line = Line.create(line.to3D(), Vector.k); }
		var R = Matrix.Rotation(t, line.direction).elements;
		var C = line.pointClosestTo(this.anchor).elements;
		var A = this.anchor.elements, D = this.direction.elements;
		var C1 = C[0], C2 = C[1], C3 = C[2], A1 = A[0], A2 = A[1], A3 = A[2];
		var x = A1 - C1, y = A2 - C2, z = A3 - C3;
		return Line.create([
      C1 + R[0][0] * x + R[0][1] * y + R[0][2] * z,
      C2 + R[1][0] * x + R[1][1] * y + R[1][2] * z,
      C3 + R[2][0] * x + R[2][1] * y + R[2][2] * z
    ], [
      R[0][0] * D[0] + R[0][1] * D[1] + R[0][2] * D[2],
      R[1][0] * D[0] + R[1][1] * D[1] + R[1][2] * D[2],
      R[2][0] * D[0] + R[2][1] * D[1] + R[2][2] * D[2]
    ]);
	},

	// Returns the line's reflection in the given point or line
	reflectionIn: function (obj) {
		if (obj.normal) {
			// obj is a plane
			var A = this.anchor.elements, D = this.direction.elements;
			var A1 = A[0], A2 = A[1], A3 = A[2], D1 = D[0], D2 = D[1], D3 = D[2];
			var newA = this.anchor.reflectionIn(obj).elements;
			// Add the line's direction vector to its anchor, then mirror that in the plane
			var AD1 = A1 + D1, AD2 = A2 + D2, AD3 = A3 + D3;
			var Q = obj.pointClosestTo([AD1, AD2, AD3]).elements;
			var newD = [Q[0] + (Q[0] - AD1) - newA[0], Q[1] + (Q[1] - AD2) - newA[1], Q[2] + (Q[2] - AD3) - newA[2]];
			return Line.create(newA, newD);
		} else if (obj.direction) {
			// obj is a line - reflection obtained by rotating PI radians about obj
			return this.rotate(Math.PI, obj);
		} else {
			// obj is a point - just reflect the line's anchor in it
			var P = obj.elements || obj;
			return Line.create(this.anchor.reflectionIn([P[0], P[1], (P[2] || 0)]), this.direction);
		}
	},

	// Set the line's anchor point and direction.
	setVectors: function (anchor, direction) {
		// Need to do this so that line's properties are not
		// references to the arguments passed in
		anchor = Vector.create(anchor);
		direction = Vector.create(direction);
		if (anchor.elements.length == 2) { anchor.elements.push(0); }
		if (direction.elements.length == 2) { direction.elements.push(0); }
		if (anchor.elements.length > 3 || direction.elements.length > 3) { return null; }
		var mod = direction.modulus();
		if (mod === 0) { return null; }
		this.anchor = anchor;
		this.direction = Vector.create([
      direction.elements[0] / mod,
      direction.elements[1] / mod,
      direction.elements[2] / mod
    ]);
		return this;
	}
};


// Constructor function
Line.create = function (anchor, direction) {
	var L = new Line();
	return L.setVectors(anchor, direction);
};

// Axes
Line.X = Line.create(Vector.Zero(3), Vector.i);
Line.Y = Line.create(Vector.Zero(3), Vector.j);
Line.Z = Line.create(Vector.Zero(3), Vector.k);



function Plane() { }
Plane.prototype = {

	// Returns true iff the plane occupies the same space as the argument
	eql: function (plane) {
		return (this.contains(plane.anchor) && this.isParallelTo(plane));
	},

	// Returns a copy of the plane
	dup: function () {
		return Plane.create(this.anchor, this.normal);
	},

	// Returns the result of translating the plane by the given vector
	translate: function (vector) {
		var V = vector.elements || vector;
		return Plane.create([
      this.anchor.elements[0] + V[0],
      this.anchor.elements[1] + V[1],
      this.anchor.elements[2] + (V[2] || 0)
    ], this.normal);
	},

	// Returns true iff the plane is parallel to the argument. Will return true
	// if the planes are equal, or if you give a line and it lies in the plane.
	isParallelTo: function (obj) {
		var theta;
		if (obj.normal) {
			// obj is a plane
			theta = this.normal.angleFrom(obj.normal);
			return (Math.abs(theta) <= Sylvester.precision || Math.abs(Math.PI - theta) <= Sylvester.precision);
		} else if (obj.direction) {
			// obj is a line
			return this.normal.isPerpendicularTo(obj.direction);
		}
		return null;
	},

	// Returns true iff the receiver is perpendicular to the argument
	isPerpendicularTo: function (plane) {
		var theta = this.normal.angleFrom(plane.normal);
		return (Math.abs(Math.PI / 2 - theta) <= Sylvester.precision);
	},

	// Returns the plane's distance from the given object (point, line or plane)
	distanceFrom: function (obj) {
		if (this.intersects(obj) || this.contains(obj)) { return 0; }
		if (obj.anchor) {
			// obj is a plane or line
			var A = this.anchor.elements, B = obj.anchor.elements, N = this.normal.elements;
			return Math.abs((A[0] - B[0]) * N[0] + (A[1] - B[1]) * N[1] + (A[2] - B[2]) * N[2]);
		} else {
			// obj is a point
			var P = obj.elements || obj;
			var A = this.anchor.elements, N = this.normal.elements;
			return Math.abs((A[0] - P[0]) * N[0] + (A[1] - P[1]) * N[1] + (A[2] - (P[2] || 0)) * N[2]);
		}
	},

	// Returns true iff the plane contains the given point or line
	contains: function (obj) {
		if (obj.normal) { return null; }
		if (obj.direction) {
			return (this.contains(obj.anchor) && this.contains(obj.anchor.add(obj.direction)));
		} else {
			var P = obj.elements || obj;
			var A = this.anchor.elements, N = this.normal.elements;
			var diff = Math.abs(N[0] * (A[0] - P[0]) + N[1] * (A[1] - P[1]) + N[2] * (A[2] - (P[2] || 0)));
			return (diff <= Sylvester.precision);
		}
	},

	// Returns true iff the plane has a unique point/line of intersection with the argument
	intersects: function (obj) {
		if (typeof (obj.direction) == 'undefined' && typeof (obj.normal) == 'undefined') { return null; }
		return !this.isParallelTo(obj);
	},

	// Returns the unique intersection with the argument, if one exists. The result
	// will be a vector if a line is supplied, and a line if a plane is supplied.
	intersectionWith: function (obj) {
		if (!this.intersects(obj)) { return null; }
		if (obj.direction) {
			// obj is a line
			var A = obj.anchor.elements, D = obj.direction.elements,
          P = this.anchor.elements, N = this.normal.elements;
			var multiplier = (N[0] * (P[0] - A[0]) + N[1] * (P[1] - A[1]) + N[2] * (P[2] - A[2])) / (N[0] * D[0] + N[1] * D[1] + N[2] * D[2]);
			return Vector.create([A[0] + D[0] * multiplier, A[1] + D[1] * multiplier, A[2] + D[2] * multiplier]);
		} else if (obj.normal) {
			// obj is a plane
			var direction = this.normal.cross(obj.normal).toUnitVector();
			// To find an anchor point, we find one co-ordinate that has a value
			// of zero somewhere on the intersection, and remember which one we picked
			var N = this.normal.elements, A = this.anchor.elements,
          O = obj.normal.elements, B = obj.anchor.elements;
			var solver = Matrix.Zero(2, 2), i = 0;
			while (solver.isSingular()) {
				i++;
				solver = Matrix.create([
          [N[i % 3], N[(i + 1) % 3]],
          [O[i % 3], O[(i + 1) % 3]]
        ]);
			}
			// Then we solve the simultaneous equations in the remaining dimensions
			var inverse = solver.inverse().elements;
			var x = N[0] * A[0] + N[1] * A[1] + N[2] * A[2];
			var y = O[0] * B[0] + O[1] * B[1] + O[2] * B[2];
			var intersection = [
        inverse[0][0] * x + inverse[0][1] * y,
        inverse[1][0] * x + inverse[1][1] * y
      ];
			var anchor = [];
			for (var j = 1; j <= 3; j++) {
				// This formula picks the right element from intersection by
				// cycling depending on which element we set to zero above
				anchor.push((i == j) ? 0 : intersection[(j + (5 - i) % 3) % 3]);
			}
			return Line.create(anchor, direction);
		}
	},

	// Returns the point in the plane closest to the given point
	pointClosestTo: function (point) {
		var P = point.elements || point;
		var A = this.anchor.elements, N = this.normal.elements;
		var dot = (A[0] - P[0]) * N[0] + (A[1] - P[1]) * N[1] + (A[2] - (P[2] || 0)) * N[2];
		return Vector.create([P[0] + N[0] * dot, P[1] + N[1] * dot, (P[2] || 0) + N[2] * dot]);
	},

	// Returns a copy of the plane, rotated by t radians about the given line
	// See notes on Line#rotate.
	rotate: function (t, line) {
		var R = Matrix.Rotation(t, line.direction).elements;
		var C = line.pointClosestTo(this.anchor).elements;
		var A = this.anchor.elements, N = this.normal.elements;
		var C1 = C[0], C2 = C[1], C3 = C[2], A1 = A[0], A2 = A[1], A3 = A[2];
		var x = A1 - C1, y = A2 - C2, z = A3 - C3;
		return Plane.create([
      C1 + R[0][0] * x + R[0][1] * y + R[0][2] * z,
      C2 + R[1][0] * x + R[1][1] * y + R[1][2] * z,
      C3 + R[2][0] * x + R[2][1] * y + R[2][2] * z
    ], [
      R[0][0] * N[0] + R[0][1] * N[1] + R[0][2] * N[2],
      R[1][0] * N[0] + R[1][1] * N[1] + R[1][2] * N[2],
      R[2][0] * N[0] + R[2][1] * N[1] + R[2][2] * N[2]
    ]);
	},

	// Returns the reflection of the plane in the given point, line or plane.
	reflectionIn: function (obj) {
		if (obj.normal) {
			// obj is a plane
			var A = this.anchor.elements, N = this.normal.elements;
			var A1 = A[0], A2 = A[1], A3 = A[2], N1 = N[0], N2 = N[1], N3 = N[2];
			var newA = this.anchor.reflectionIn(obj).elements;
			// Add the plane's normal to its anchor, then mirror that in the other plane
			var AN1 = A1 + N1, AN2 = A2 + N2, AN3 = A3 + N3;
			var Q = obj.pointClosestTo([AN1, AN2, AN3]).elements;
			var newN = [Q[0] + (Q[0] - AN1) - newA[0], Q[1] + (Q[1] - AN2) - newA[1], Q[2] + (Q[2] - AN3) - newA[2]];
			return Plane.create(newA, newN);
		} else if (obj.direction) {
			// obj is a line
			return this.rotate(Math.PI, obj);
		} else {
			// obj is a point
			var P = obj.elements || obj;
			return Plane.create(this.anchor.reflectionIn([P[0], P[1], (P[2] || 0)]), this.normal);
		}
	},

	// Sets the anchor point and normal to the plane. If three arguments are specified,
	// the normal is calculated by assuming the three points should lie in the same plane.
	// If only two are sepcified, the second is taken to be the normal. Normal vector is
	// normalised before storage.
	setVectors: function (anchor, v1, v2) {
		anchor = Vector.create(anchor);
		anchor = anchor.to3D(); if (anchor === null) { return null; }
		v1 = Vector.create(v1);
		v1 = v1.to3D(); if (v1 === null) { return null; }
		if (typeof (v2) == 'undefined') {
			v2 = null;
		} else {
			v2 = Vector.create(v2);
			v2 = v2.to3D(); if (v2 === null) { return null; }
		}
		var A1 = anchor.elements[0], A2 = anchor.elements[1], A3 = anchor.elements[2];
		var v11 = v1.elements[0], v12 = v1.elements[1], v13 = v1.elements[2];
		var normal, mod;
		if (v2 !== null) {
			var v21 = v2.elements[0], v22 = v2.elements[1], v23 = v2.elements[2];
			normal = Vector.create([
        (v12 - A2) * (v23 - A3) - (v13 - A3) * (v22 - A2),
        (v13 - A3) * (v21 - A1) - (v11 - A1) * (v23 - A3),
        (v11 - A1) * (v22 - A2) - (v12 - A2) * (v21 - A1)
      ]);
			mod = normal.modulus();
			if (mod === 0) { return null; }
			normal = Vector.create([normal.elements[0] / mod, normal.elements[1] / mod, normal.elements[2] / mod]);
		} else {
			mod = Math.sqrt(v11 * v11 + v12 * v12 + v13 * v13);
			if (mod === 0) { return null; }
			normal = Vector.create([v1.elements[0] / mod, v1.elements[1] / mod, v1.elements[2] / mod]);
		}
		this.anchor = anchor;
		this.normal = normal;
		return this;
	}
};

// Constructor function
Plane.create = function (anchor, v1, v2) {
	var P = new Plane();
	return P.setVectors(anchor, v1, v2);
};

//End of Matrix Library by Sylvester

//Start of Modified Gram Schmidt Orthogonalization

function gram(matrixA) {

	var m = matrixA[0].length;
	var n = matrixA.length;

	var matrixQ = zeros(m, n);
	var matrixR = zeros(n, n);

	var tempArray = matrixA;
	matrixQ = tempArray;

	var tempQVariableOne = new Array(m);
	var tempQVariableTwo = new Array(m);
	var tempQVariableThree = new Array(m);
	var tempSubtractVect = new Array(m);
	var tempNormVariable = new Array(m);

	for (var j = 0; j < n; j++) {
		tempNormVariable = vector(matrixQ, j, 2);
		matrixR[j][j] = norm(tempNormVariable, tempNormVariable.length);

		for (var i = 0; i < tempNormVariable.length; i++) {
			tempQVariableOne[i] = tempNormVariable[i] / matrixR[j][j];
		}

		for (var l = 0; l < tempNormVariable.length; l++) {
			matrixQ[j][l] = tempQVariableOne[l];
		}

		for (k = j + 1; k < n; k++) {
			tempQVariableTwo = vector(matrixQ, k, 2);
			matrixR[k][j] = dotProduct(tempQVariableOne, tempQVariableTwo, tempQVariableOne.length);
			tempSubtractVect = scalarMatrixMultiplication(matrixR[k][j], tempQVariableOne, tempQVariableOne.length);
			tempQVariableThree = vectorSubtraction(tempQVariableTwo, tempSubtractVect, tempSubtractVect.length);

			for (s = 0; s < m; s++) {
				matrixQ[k][s] = tempQVariableThree[s];
			}

		}


	}
	return [matrixQ, matrixR];

}



//End of Modified Gram Schmidt Orthogonalization




//Start of Vandermonde Matrix Construction
function vandermonde(x, orderOfPolynomial) {
	var n = orderOfPolynomial + 1;

	//Builds initial Vandermonde matrix
	var matrixV = new Array(n)
	for (var k = 0; k < n; k++) {
		matrixV[k] = new Array(x[0].length);
	}

	for (var i = 0; i < matrixV.length; i++) {
		for (var j = 0; j < matrixV[0].length; j++) {
			matrixV[i][j] = 0;
			if (i === (n - 1)) {
				matrixV[i][j] = 1;
			}
		}
	}

	for (var p = orderOfPolynomial - 1; p >= 0; p--) {
		for (var w = 0; w < matrixV[0].length; w++) {
			matrixV[p][w] = (x[0][w]) * (matrixV[p + 1][w]);
		}
	}


	return matrixV;

}

//End of Vandermonde Matrix Construction

//Start of POLYFIT 
function polyfit(xArray, yArray, order) {

	var x = matrixTranspose(xArray);


	//To be used later if deemed necessary
	/*
	var meanX = meanValue(x);
	var stdevX = standardDeviation(x);
        
	for (var i = 0; i < x[0].length; i++) {
	x[0][i] = (x[0][i] - meanX) / (stdevX);
	}
	*/
	var y = Matrix.create(yArray);
	y = y.transpose();

	var vMatrix = vandermonde(x, order);
	var qrMatrix = gram(vMatrix);


	var qMatrix = qrMatrix[0];
	var rMatrix = qrMatrix[1];

	var rMatrixInverse = Matrix.create(rMatrix);
	rMatrixInverse = rMatrixInverse.inverse();

	var qTransposed = Matrix.create(qMatrix);
	var qTransposed = qTransposed.transpose();

	var matrixToBeMultiplied = y.multiply(qTransposed);

	var p = matrixToBeMultiplied.multiply(rMatrixInverse)

	var matrixP = Matrix.create(p);
	matrixP = matrixP.transpose();


	return matrixP;


}
//End of POLYFIT

//Start of Polyval
function polyval(funcArray, x) {
	var polynomialOrder = funcArray.rows() - 1;
	var size = polynomialOrder + 1;
	var result = 0;

	for (var i = 0; i < size; i++) {
		result += funcArray.elements[i] * (Math.pow(x, polynomialOrder));
		polynomialOrder -= 1;


	}
	return result;

}

//End of Polyval

//Start of SUM Calculation
function sum(arr) {

	var result = 0;

	for (var i = 0; i < arr.length; i++) {

		result += arr[i];

	}
	return result;

}

//End of SUM Calculation


//Start of Vehicle Velocity Factor Calculations
function velocityFactorCalculation(nt, nf, np, wheelRadius) {

	//Builds the two Arrays
	var rpmIteration = new Array();
	var velocityIteration = new Array();
	var pi = 3.1416;
	var velocityFactorArray = new Array();
	var velocityFactor = 0;

	rpmIteration[0] = 0;
	velocityIteration[0] = 0;

	for (var i = 0; i <= 90; i++) {
		rpmIteration[i + 1] = rpmIteration[i] + 10;
		velocityIteration[i + 1] = (((rpmIteration[i] * 2 * pi) / 60) / (nf * np * nt)) * wheelRadius;
	}
	velocityFactorArray = polyfit(rpmIteration, velocityIteration, 1);
	var velocityFactorWorking = velocityFactorArray.elements[0];
	velocityFactor = 1 / (velocityFactorWorking[0]);

	return velocityFactor;

}

//End of Vehicle Velocity Factor Calculations
//Start of Instantaneous Torque Calculation

function torqueLocator(currentRpm, rpm, tor) {
	var outputTorque;
    

	for (var i = 1; i < rpm.length; i++) {

		if (rpm[i] > currentRpm) {
			outputTorque = (tor[i - 1] + (currentRpm - rpm[i - 1]) * (tor[i] - tor[i - 1]) / (rpm[i] - rpm[i - 1]));
			break;
		}


	}
	return outputTorque;

}

//End of Instantaneous Torque Calculation


//Start of Instantaneous Drag Force Calculation

function dragForceCalculation(rho, velocity, coefficientOfDrag, frontalArea, downforce, liftCoefficient) {

    var outputDrag;
    outputDrag = 0.5 * rho * (Math.pow(velocity, 2)) * coefficientOfDrag * frontalArea + downforce * liftCoefficient;
	return outputDrag;
}

//End of Instantaneous Drag Force Calculation

//Start of Instantaneous Rolling Resistance Force Calculation
function rollingResistanceForceCalculation(f_0, f_s, velocity, vehicleWeight, downforce) {

	var outputRollingResistance;
	vehicleWeight = vehicleWeight + downforce;

	outputRollingResistance = (f_0 + 3.25 * f_s * (Math.pow((velocity / 100), 2.5))) * vehicleWeight;
	return outputRollingResistance;

}

//End of Instantaneous Rolling Resistance Force Calculation
//Start of Instantaneous Longitudinal Acceleration Calculation

function longitudinalAccelerationCalculation(instantaneousTorque, drivetrainEfficiency, wheelRadius, instantaneousDragForce, instantaneousRollingResistanceForce, vehicleWeight, ntf) {

	var outputAcceleration;
	var g = 32.17405;
	outputAcceleration = ((instantaneousTorque * ntf * drivetrainEfficiency) / (wheelRadius) - instantaneousRollingResistanceForce - instantaneousDragForce) / (vehicleWeight / g);

	if (outputAcceleration < 0) {
	    outputAcceleration = 0.5;
    
    }
    
    return outputAcceleration;

}



//End of Instantaneous Longitudinal Acceleration Calculation



//Start of Instantaneous Friction Caclulation

function frictionCalculation(weightOnTire, tireChoice) {

	/*
	Tire Selection.  The equations that follow are based on tabulated Excel data.  The maximum coefficient of friction was assumed in each load distribution.
	A graph was created and the equation extrapolated to be used in the following code.  

	The correction factor is 0.67 for the data, which is multiplied by the instantaneous coefficients of friction to obtain a real world value.
	Used values other than 0.67 to achieve more realistic acceleration results as 0.67 yields a longitudinal acceleration that is sometimes not feasible.

    The absolute minimum coefficient of friction is 0.8 which is set as the bar.

	Tire 1 = Goodyear 13"
	Tire 2 = Hoosier 13"
	Tire 3 = Smaller Version of Hoosier 13"
	Tire 4 = Michelin 13"

	*/

	var coefficientOfFrictionCorrection = 0.55;
	var frictionInstant;
	var lowestFrictionValue = 0.8;

	if (tireChoice === 1) {

	    frictionInstant = coefficientOfFrictionCorrection * ((-2.8103) * Math.pow(10, -13) * Math.pow(weightOnTire, 4) + 1.0733 * Math.pow(10, -9) * Math.pow(weightOnTire, 3) - 9.0420 * Math.pow(10, -7) * Math.pow(weightOnTire, 2) - 6.1908 * Math.pow(10, -4) * weightOnTire + 2.5359);

		if (weightOnTire > 2000) {
		    frictionInstant = 0.8;
        }
	}
	else if (tireChoice === 2) {

		frictionInstant = coefficientOfFrictionCorrection * ((-3.2015) * Math.pow(10, -13) * Math.pow(weightOnTire, 4) + 1.2594 * Math.pow(10, -9) * Math.pow(weightOnTire, 3) - 1.1920 * Math.pow(10, -6) * Math.pow(weightOnTire, 2) - 4.7533 * Math.pow(10, -4) * weightOnTire + 2.6820);

		if (weightOnTire > 2000) {
		    frictionInstant = 0.8;
		}
	}
	else if (tireChoice === 3) {

	frictionInstant = coefficientOfFrictionCorrection * ((-4.3383) * Math.pow(10, -13) * Math.pow(weightOnTire, 4) + 1.8600 * Math.pow(10, -9) * Math.pow(weightOnTire, 3) - 2.3000 * Math.pow(10, -6) * Math.pow(weightOnTire, 2) + 3.4083 * Math.pow(10, -4) * weightOnTire + 2.5170);

		if (weightOnTire > 2000) {
		    frictionInstant = 0.8;
		}
    }
    else if (tireChoice === 5) {

    frictionInstant = 1.3;

         if (weightOnTire > 2000) {
             frictionInstant = 0.8;
        }
    }
	else {

		frictionInstant = coefficientOfFrictionCorrection * (2.1891 * Math.pow(10, -13) * Math.pow(weightOnTire, 4) - 1.3599 * Math.pow(10, -9) * Math.pow(weightOnTire, 3) + 3.0851 * Math.pow(10, -6) * Math.pow(weightOnTire, 2) - 3.0166 * Math.pow(10, -3) * weightOnTire + 2.2424);

		if (weightOnTire > 550) {
		    frictionInstant = 0.8;
		}
	}

	return frictionInstant;

}


//End of Instantaneous Friction Caclulation
//Start of Downforce Calculation
function downforceCalculation(constant, currentVelocity) {

    var downforce = constant * Math.pow(currentVelocity, 2);
    return downforce;

}




//End of Downforce Calculation



//END OF BUILT-IN FUNCTIONS
//END OF BUILT-IN FUNCTIONS
//END OF BUILT-IN FUNCTIONS
//END OF BUILT-IN FUNCTIONS





//End of Mean Value Function
function Main(inVehicleWeight, inTireChoice, inWheelBase, inWheelRadius, inTrackWidth, inHcg, inWeightDistribution, inShiftRpm, inNf, inDownforceAtVelocityInput, inLiftDragCoefficient, inEngineTorque, inForcedInduction) {


	




	//Variable Declaration
	var bsfc;
	var cornerTotal = [];
	var prevCorners;
	var numberOfLaps;
	var straightSections;
	var dist = [];
	var numberOfCorners = [];
	var cornerRadius = [];
	var corneringDist = [];
	var vehicleWeight;
	var coefficientOfDrag;
	var frontalArea;
	var rho;
	var f_0;
	var f_s;
	var drivetrainEfficiency;
	var g;
	var pi;
	var wheelRadius;
	var iterationTime;
	var wheelbase;
	var trackWidth;
	var trackA;
	var hcg;
	var weightDistribution;
	var cLength;
	var bLength;
	var nt_1;
	var nt_2;
	var nt_3;
	var nt_4;
	var nt_5;
	var nt_6;
	var nf;
	var np;
	var ntf;
	var shiftTime;
	var shiftRpm;
	var shifts;
	var totalShifts;
	var tireChoice;
	var competitionAcceleration;
	var competitionSkidpad;
	var competitionAutocross;
	var competitionFuelEconomy;
	var competitionEndurance;
	var competitionStaticTotal;
	var p = [];
	var s = [];
	var n;
	var velocityFactorOne;
	var velocityFactorTwo;
	var velocityFactorThree;
	var velocityFactorFour;
	var velocityFactorFive;
	var velocityFactorSix;
	var currentRpm;
	var instantaneousTorque;
	var instantaneousDragForce;
	var instantaneousRollingResistanceForce;
	var instantaneousLongitudinalAcceleration;
	var weightOnTire;
	var timeLaunching;
	var t;
	var launchingTime;
	var launchingDistance;
	var time;
	var instantaneousCoefficientOfFriction;
	var maximumTractiveForce;
	var maximumLongitudinalAcceleration;
	var corneringAcceleration;
	var corneringCounter;
	var difference;
	var tireOne;
	var tireTwo;
	var tireThree;
	var tireFour;
	var corneringForce;
	var oldCorneringAcceleration;
	var corneringWeightDistribution;
	var corneringComparison;
	var corneringTime;
	var corneringVelocity = [];
	var brakingAcceleration;
	var corneringBrakeTo;
	var velocityFactorCurrent;
	var x;
	var count;
	var velocityStore = [];
	var xStore = [];
	var instantaneousLongitudinalAccelerationStore = [];
	var timeStraightSection;
	var accelerationCurve = [];
	var accelerating;
	var decelerating;
	var decelAccelPoint;
	var totalAccelDecel;
	var brakingToBeSqrt;
	var initialVelocity;
	var horsepowerStore = [];
	var averageLongitudinalAcceleration = [];
	var averageHorsepower = [];
	var decelerationVelocity;
	var brakingDistance = [];
	var brakingTime;
	var totalStraightSectionTime = [];
	var totalNumberOfShifts = [];
	var totalDistance;
	var totalTime;
	var fuelConsumptionRate;
	var totalFuel;
	var fuelUsed;
	var brakingTimeStore = [];
	var brakingTracker;
	var totalStraightSectionTimeOutput;
	var totalNumberOfShiftsOutput;
	var rpm = [];
	var tor = [];
	var timeOpenThrottle;
	var totalTimeWideOpenThrottle = [];
	var outputTotalWideOpenThrottle;
	var percentTimeWideOpenThrottle;
	var timeTractionLimited;
	var totalTimeTractionLimited = [];
	var outputTimeTractionLimited;
	var timeStraightSectionAcceleration;
	var accelerationShifts;
	var accelerationLongitudinalAccelerationStore = [];
	var averageAccelerationLongitudinalAccel;
	var accelerationTrapSpeed;
	var accelerationTimeTractionLimited;
	var accelerationPoints;
	var accelerationWideOpenThrottleTime;
	var accelerationWideOpenThrottlePercentage;
	var skidpadConstantRadius;
	var skidpadLength;
	var skidpadTime;
	var skidpadLateralAccel;
	var skidpadVelocity;
	var skidpadRadius;
	var skidpadPoints;
	var autocrossTime;
	var autocrossPoints;
	var autocrossTotalShifts;
	var autocrossTractionLimitedTime;
	var autocrossAverageVelocity;
	var autocrossVelocityStore = [];
	var autocrossAverageVelocityTotal = [];
	var meanAutocrossVelocity;
	var maxHorsepower;
	var maxHorsepowerCalc;
	var downforceConstant;
	var downforceVelocityInput;
	var instantaneousDownforce;
	var engineChoice;
	var torqueCurveChoice;
	var engineTorqueChoice;
	var accelerationDistance;
	var launchingShifts;
	var forcedInductionIncrease;
	



	//IN PROGRESS

	/* Engine Choices

    1: WR450 Yamaha
    2: R6 Yamaha
    3: CBR600 Honda


    */


	/* Torque Curve Choices

             1:Rennteam Uni Stuttgart --- 2005 Honda CBR600 RR
             2:TU Graz Racing Team --- 2006 Yamaha R6
             3:Joanneum Racing Graz --- Rotax 540
             4:Metropolia Motorsport(Helsinki) --- Yamaha R6
             5:TU Darmstadt Racing Team e.V. --- Suzuki GSX-R 600
             5:Oxford Brookes University --- KTM EXC-R 530
             6:KA-RaceIng(Karlsruhe) --- 2003 Honda CBR600
             7:FaSTTUBe(Berlin) --- Suzuki GSX-R 600
             8:Race UP Team(Padova) --- 2006 Kawasaki NinjaZX6-RR Supersport
             9:Unicorn Race Engineering(Aalborg) --- Honda CBR600-RR
            10:Global Formula Racing (Ravensburg) --- Yamaha YZF R6
            11:Technikum Mittweida Motorsport --- Honda CBR600 RR PC37
            12:Team wob-racing(Wolfebüttel) --- BLANK/NA
            13:University Racing Eindhoven --- Suzuki GSX-R 600
            14:Raceyard Kiel --- Honda CBR600 RR PC37
            15:ISAT Formula Student(Nevers) --- BLANK/NA
            16:Kaiserslautern Racing Team --- Suzuki GSX-R 600 K4
            17:Scuderia Mensa HS RheinMain Racing--- Suzuki GSX-R 600 K3
            18:Ignition Racing Team(Osnabruck) --- BLANK/NA
            19:Eleven O Six Racing Team(Hamburg) --- 2002 Honda CBR 600 F (PC35)
            20:Racetech Racing Team TU Bergakademie Freiberg e.V. --- Honda CBR600
            21:Tokyo Denki University Formula SAE Project --- HONDA CRF450X
            22:Montreal ETS --- BLANK/NA
            23:BA Motors(Berlin) --- KTM LC4 640 turbo
            24:Formula Racing Cologne(Koln) --- 2005 YZF R6
            25:WHZ Racing Team(Zwickau) --- Honda CBR 600


    */



    //Engine Choices
	//engineChoice=1  Honda CBR600
	//engineChoice=2  Honda CRF450
	//engineChoice=3  Kawasaki NinjaZX6-RR
	//engineChoice=4  KTM EXC-R 530
	//engineChoice=5  KTM LC4 640 turbo
	//engineChoice=6  Rotax 540
	//engineChoice=7  Suzuki GSXR600
	//engineChoice=8  Yamaha R6
	//engineChoice=9  Yamaha WR450




	//Torque Curve Choices

	//engineChoice=1  Honda CBR600
	//Honda CBR600
	//torqueCurveChoice=1  Rennteam Uni Stuttgart
	//torqueCurveChoice=2  KA-RaceIng(Karlsruhe)
	//torqueCurveChoice=3  Unicorn Race Engineering(Aalborg)
	//torqueCurveChoice=4  Technikum Mittweida Motorsport
	//torqueCurveChoice=5  Raceyard Kiel
	//torqueCurveChoice=6  Eleven O Six Racing Team(Hamburg)
	//torqueCurveChoice=7  Racetech Racing Team TU Bergakademie Freiberg
	//torqueCurveChoice=8  WHZ Racing Team(Zwickau)

	//engineChoice=2  Honda CRF450
	//Honda CRF450
	//torqueCurveChoice=1  Tokyo Denki University Formula SAE Project

	//engineChoice=3  Kawasaki NinjaZX6-RR
	//Kawasaki NinjaZX6-RR
	//torqueCurveChoice=1  Race UP Team(Padova)

	//engineChoice=4  KTM EXC-R 530
	//KTM EXC-R 530
	//torqueCurveChoice=1  Oxford Brookes University

	//engineChoice=5  KTM LC4 640 turbo
	//KTM LC4 640 turbo
	//torqueCurveChoice=1  BA Motors(Berlin)

	//engineChoice=6  Rotax 540
	//Rotax 540
	//torqueCurveChoice=1  Joanneum Racing Graz

	//engineChoice=7  Suzuki GSXR600
	//Suzuki GSXR600
	//torqueCurveChoice=1  TU Darmstadt Racing Team e.V
	//torqueCurveChoice=2  FaSTTUBe(Berlin)
	//torqueCurveChoice=3  University Racing Eindhoven
	//torqueCurveChoice=4  Scuderia Mensa HS RheinMain Racing
	//torqueCurveChoice=5  Kaiserslautern Racing Team

	//engineChoice=8  Yamaha R6
	//Yamaha R6
	//torqueCurveChoice=1  TU Graz Racing Team
	//torqueCurveChoice=2  Metropolia Motorsport(Helsinki)
	//torqueCurveChoice=3  Global Formula Racing (Ravensburg)
	//torqueCurveChoice=4  Formula Racing Cologne(Koln)

	//engineChoice=9  Yamaha WR450
	//torqueCurveChoice=1  Delft University

	engineTorqueChoice = parseInt(inEngineTorque);

	if (engineTorqueChoice === 1) {
	    engineChoice = 1;
	    torqueCurveChoice = 1;
    }
	else if (engineTorqueChoice === 2) {
	    engineChoice = 1;
	    torqueCurveChoice = 2;
	}
	else if (engineTorqueChoice === 3) {
	    engineChoice = 1;
	    torqueCurveChoice = 3;
	}
	else if (engineTorqueChoice === 4) {
	    engineChoice = 1;
	    torqueCurveChoice = 4;
	}
	else if (engineTorqueChoice === 5) {
	    engineChoice = 1;
	    torqueCurveChoice = 5;
	}
	else if (engineTorqueChoice === 6) {
	    engineChoice = 1;
	    torqueCurveChoice = 6;
	}
	else if (engineTorqueChoice === 7) {
	    engineChoice = 1;
	    torqueCurveChoice = 7;
	}
	else if (engineTorqueChoice === 8) {
	    engineChoice = 1;
	    torqueCurveChoice = 8;
	}
	else if (engineTorqueChoice === 9) {
	    engineChoice = 2;
	    torqueCurveChoice = 1;
	}
	else if (engineTorqueChoice === 10) {
	    engineChoice = 3;
	    torqueCurveChoice = 1;
	}
	else if (engineTorqueChoice === 11) {
	    engineChoice = 4;
	    torqueCurveChoice = 1;
	}
	else if (engineTorqueChoice === 12) {
	    engineChoice = 5;
	    torqueCurveChoice = 1;
	}
	else if (engineTorqueChoice === 13) {
	    engineChoice = 6;
	    torqueCurveChoice = 1;
	}
	else if (engineTorqueChoice === 14) {
	    engineChoice = 7;
	    torqueCurveChoice = 1;
	}
	else if (engineTorqueChoice === 15) {
	    engineChoice = 7;
	    torqueCurveChoice = 2;
	}
	else if (engineTorqueChoice === 16) {
	    engineChoice = 7;
	    torqueCurveChoice = 3;
	}
	else if (engineTorqueChoice === 17) {
	    engineChoice = 7;
	    torqueCurveChoice = 4;
	}
	else if (engineTorqueChoice === 18) {
	    engineChoice = 7;
	    torqueCurveChoice = 5;
	}
	else if (engineTorqueChoice === 19) {
	    engineChoice = 8;
	    torqueCurveChoice = 1;
	}
	else if (engineTorqueChoice === 20) {
	    engineChoice = 8;
	    torqueCurveChoice = 2;
	}
	else if (engineTorqueChoice === 21) {
	    engineChoice = 8;
	    torqueCurveChoice = 3;
	}
	else if (engineTorqueChoice === 22) {
	    engineChoice = 8;
	    torqueCurveChoice = 4;
	}
	else if (engineTorqueChoice === 23) {
	    engineChoice = 9;
	    torqueCurveChoice = 1;
	}

	if (engineChoice === 1) {

	    //Transmission and shifting parameters

	    nt_1 = 33 / 12;                 //1st Gear Ratio
	    nt_2 = 32 / 16;                 //2nd Gear Ratio
	    nt_3 = 30 / 18;                //3rd Gear Ratio
	    nt_4 = 26 / 18;                 //4th Gear Ratio
	    nt_5 = 30 / 23;                //5th Gear Ratio
	    nt_6 = 29 / 24;                  //6th Gear Ratio
	    nf = parseFloat(inNf);          //Final Gear Ratio
	    np = 76 / 36;                   //Primary Reduction Ratio

	    if (torqueCurveChoice === 1) {
	        rpm = [1, 500, 7245, 7457, 7664, 7875, 8078, 8287, 8474, 8665, 8868, 9054, 9248, 9454, 9669, 9889, 10114, 10333, 10540, 10759, 10966, 11162, 11371, 11562, 11759, 11950, 12133, 12320, 12490, 12665, 12830, 12986, 13142, 13269, 13354, 13410, 13442, 13455, 13460];
	        tor = [0.1, 2, 39.5333, 39.60709, 40.71343, 42.1148, 42.1148, 41.30348, 40.49216, 39.60709, 39.45957, 39.38582, 39.38582, 39.90211, 40.86094, 42.48358, 43.14739, 43.73744, 43.9587, 43.66368, 42.99987, 42.48358, 41.74602, 41.15597, 40.56592, 39.97587, 39.23831, 38.27948, 37.32064, 36.14055, 35.25547, 34.3704, 33.41157, 31.27264, 27.43731, 21.24179, 16.07885, 10.62089, 6.490547];
	    }
	    if (torqueCurveChoice === 2) {
	        rpm = [1, 500, 6537, 6677, 6814, 6965, 7107, 7265, 7427, 7579, 7731, 7895, 8060, 8246, 8437, 8634, 8825, 9011, 9209, 9392, 9577, 9766, 9954, 10126, 10308, 10486, 10666, 10841, 11011, 11184, 11341, 11495, 11620, 11713, 11785, 11830, 11862, 11892];
	        tor = [0.1, 2, 34.223, 33.7804, 34.0016, 34.29664, 34.444, 35.2554, 37.3944, 38.35323, 38.13196, 37.54191, 38.72201, 41.45099, 43.58992, 45.43383, 46.54017, 46.90895, 46.90895, 46.76144, 46.54017, 46.17139, 45.50758, 44.99129, 44.54875, 44.32749, 44.10622, 43.81119, 43.2949, 42.70485, 41.59851, 40.27089, 38.05821, 33.26405, 28.02736, 21.53681, 16.52139, 12.02226];
	    }
	    if (torqueCurveChoice === 3) {
	        rpm = [1, 500, 5108, 5194, 5291, 5372, 5466, 5557, 5644, 5741, 5834, 5930, 6018, 6113, 6211, 6305, 6418, 6520, 6623, 6733, 6841, 6948, 7062, 7176, 7289, 7408, 7531, 7652, 7772, 7890, 8005, 8128, 8241, 8356, 8469, 8573, 8681, 8790, 8882, 8996, 9095, 9195, 9294, 9391, 9486, 9581, 9673, 9765, 9856, 9941, 10023, 10118, 10205, 10289, 10374, 10455, 10540, 10611, 10696, 10776, 10850, 10919, 11001, 11072, 11145, 11212, 11283, 11348, 11395, 11419, 11443, 11456, 11474, 11482];
	        tor = [0.1, 2, 41.37724, 41.59851, 41.96729, 42.40982, 42.70485, 42.99987, 43.2949, 44.10622, 44.62251, 44.91753, 45.1388, 45.43383, 46.02388, 47.13022, 48.08905, 48.82661, 49.1954, 50.00671, 50.67052, 51.55559, 52.2194, 52.95696, 53.62077, 54.65336, 55.68594, 56.64477, 57.23482, 57.6036, 57.67736, 57.6036, 57.45609, 57.01355, 56.12848, 54.87462, 53.32574, 52.29316, 51.33433, 50.81803, 50.15423, 49.41666, 48.6791, 47.86778, 47.27773, 46.54017, 46.02388, 45.43383, 44.69627, 44.10622, 43.66368, 43.14739, 42.55734, 42.1148, 41.52475, 41.15597, 40.56592, 39.97587, 39.31206, 38.79577, 38.27948, 37.54191, 36.80435, 36.2143, 35.5505, 35.0342, 34.3704, 33.78035, 31.05137, 22.34813, 15.85759, 12.02226, 10.1046, 8.703233];
	    }
	    if (torqueCurveChoice === 4) {
            rpm = [1, 500, 5027, 5089, 5158, 5223, 5292, 5372, 5439, 5517, 5591, 5665, 5741, 5814, 5887, 5967, 6047, 6121, 6194, 6275, 6354, 6427, 6505, 6585, 6664, 6740, 6815, 6892, 6970, 7047, 7124, 7202, 7278, 7357, 7440, 7518, 7594, 7672, 7747, 7818, 7890, 7959, 8034, 8106, 8179, 8253, 8321, 8391, 8461, 8531, 8599, 8667, 8741, 8808, 8874, 8945, 9016, 9083, 9153, 9223, 9291, 9360, 9427, 9494, 9559, 9627, 9697, 9764, 9836, 9904, 9974, 10045, 10119, 10195, 10269, 10342, 10414, 10490, 10568, 10644, 10716, 10792, 10867, 10938, 11010, 11085, 11157, 11223, 11297, 11367, 11434, 11508, 11574, 11640, 11707, 11767, 11829, 11890, 11950, 12011, 12067, 12123, 12174, 12222, 12270, 12302];
	        tor = [0.1, 2, 22.27, 23.01, 23.97, 24.92, 26.10, 26.92, 27.65, 28.17, 28.54, 28.83, 28.98, 29.28, 29.50, 29.79, 29.94, 29.87, 29.79, 29.87, 30.24, 30.46, 30.46, 30.24, 30.09, 30.09, 30.09, 30.01, 30.09, 30.16, 30.38, 30.53, 30.60, 30.68, 30.75, 30.68, 30.60, 30.46, 30.166, 29.87, 29.50, 29.13, 28.83, 28.61, 28.46, 28.39, 28.32, 28.17, 27.95, 27.73, 27.58, 27.36, 27.28, 27.21, 27.21, 27.14, 27.06, 27.14, 27.06, 27.06, 26.99, 26.99, 26.99, 26.84, 26.77, 26.77, 26.84, 26.84, 27.06, 27.21, 27.43, 27.73, 28.02, 28.46, 28.83, 29.13, 29.42, 29.57, 29.79, 29.79, 29.72, 29.57, 29.42, 29.28, 29.13, 28.98, 28.69, 28.39, 28.17, 28.10, 28.02, 27.87, 27.36, 26.92, 26.55, 26.109, 24.92, 24.48, 24.63, 24.19, 23.74, 23.159, 21.83, 21.02, 20.13, 17.99];
	    }
	    if (torqueCurveChoice === 5) {
	        rpm = [1, 500, 7478, 7660, 7842, 8014, 8192, 8357, 8521, 8669, 8823, 8981, 9138, 9297, 9452, 9620, 9786, 9942, 10099, 10260, 10418, 10578, 10735, 10883, 11039, 11195, 11344, 11489, 11623, 11758, 11888, 12019, 12141, 12255, 12371, 12490, 12591, 12669, 12729, 12780, 12836, 12885, 12930, 12968, 12995];
	        tor = [0.1, 2, 50.22, 49.85, 49.19, 48.45, 48.08, 46.90, 45.50, 43.95, 43.22, 42.70, 42.48, 42.40, 42.40, 42.77, 43.51, 44.03, 44.03, 43.81, 43.44, 43.14, 42.63, 42.18, 41.89, 41.74, 41.59, 40.93, 39.90, 38.72, 37.32, 36.06, 35.10, 34.07, 33.11, 32.37, 31.19, 28.10, 23.60, 19.47, 17.18, 15.63, 14.38, 13.27, 11.80];
	    }
	    if (torqueCurveChoice === 6) {
	        rpm = [1, 500, 6703, 6784, 6871, 6968, 7057, 7156, 7254, 7353, 7440, 7537, 7619, 7710, 7799, 7888, 7974, 8051, 8121, 8211, 8285, 8364, 8461, 8553, 8654, 8750, 8851, 8952, 9057, 9161, 9253, 9353, 9436, 9533, 9617, 9705, 9787, 9868, 9960, 10039, 10132, 10222, 10309, 10397, 10478, 10565, 10636, 10727, 10811, 10889, 10953, 11027, 11090, 11165, 11244, 11319, 11392, 11451, 11520, 11589, 11652, 11728, 11789, 11852, 11914, 11962, 12021, 12080, 12137, 12196, 12249, 12298, 12347, 12394, 12442, 12476, 12519, 12567, 12604, 12649, 12691, 12728, 12764, 12803, 12829, 12863, 12894, 12916, 12940, 12955, 12968];
	        tor = [0.1, 2, 13.64489976, 14.75124299, 16.44763593, 18.2915413, 19.39788453, 19.69290939, 20.28295911, 21.09427747, 21.53681476, 21.02052125, 20.65174018, 19.47164074, 19.17661588, 19.84042182, 19.7666656, 18.51280995, 17.84900401, 17.92276023, 17.99651644, 17.7752478, 18.07027266, 19.47164074, 20.79925261, 21.53681476, 21.97935205, 21.83183962, 22.42188934, 22.79067041, 22.20062069, 21.97935205, 20.57798396, 20.13544668, 20.06169046, 19.84042182, 19.39788453, 19.02910345, 19.02910345, 19.02910345, 19.02910345, 19.2503721, 19.54539696, 19.54539696, 19.39788453, 19.02910345, 18.43905373, 18.07027266, 18.58656616, 18.73407859, 17.7752478, 16.96392943, 16.22636728, 15.63631756, 15.85758621, 16.59514836, 16.66890457, 16.22636728, 15.48880513, 15.48880513, 15.48880513, 15.12002406, 15.04626785, 14.97251163, 14.52997434, 13.49738733, 12.8335814, 13.05485004, 13.05485004, 13.12860626, 13.20236247, 12.68606897, 11.65348196, 11.21094467, 11.06343224, 10.54713873, 9.367039296, 9.662064156, 10.10460145, 9.73582037, 9.73582037, 9.588307941, 9.145770651, 8.555720932, 8.260696072, 7.670646352, 7.523133923, 7.228109063, 6.638059343, 5.457959905, 4.572885326];
	    }
	    if (torqueCurveChoice === 7) {
	        rpm = [1, 500, 5170, 5305, 5436, 5572, 5696, 5832, 5966, 6101, 6236, 6370, 6498, 6623, 6746, 6870, 7006, 7133, 7253, 7385, 7504, 7631, 7758, 7877, 8002, 8120, 8246, 8366, 8488, 8610, 8735, 8867, 8996, 9129, 9264, 9394, 9527, 9655, 9784, 9908, 10038, 10160, 10286, 10411, 10522, 10638, 10753, 10876, 10986, 11092, 11201, 11305, 11415, 11519, 11621, 11717, 11816, 11910, 12004, 12095, 12150, 12183, 12205, 12212];
	        tor = [0.1, 2, 24.708332, 28.10111789, 29.72375462, 31.27263513, 31.93644106, 32.15770971, 32.45273457, 32.52649078, 32.52649078, 32.52649078, 32.52649078, 32.30522214, 31.86268485, 31.56765999, 31.42014756, 31.49390377, 31.34639134, 31.05136648, 30.75634162, 30.60882919, 30.46131677, 30.46131677, 30.38756055, 30.24004812, 30.09253569, 29.94502326, 29.87126705, 29.94502326, 30.24004812, 30.75634162, 31.34639134, 31.78892863, 32.15770971, 32.45273457, 32.52649078, 32.37897835, 32.08395349, 31.78892863, 31.56765999, 31.34639134, 30.97761027, 30.60882919, 30.16629191, 29.72375462, 29.35497354, 28.91243625, 28.46989896, 27.95360546, 27.51106817, 27.06853088, 26.77350602, 26.47848116, 25.88843144, 25.37213794, 24.78208822, 24.33955093, 23.82325742, 23.30696392, 21.90559583, 17.11144186, 11.8747506, 8.113183642];
	    }
	    if (torqueCurveChoice === 8) {
	        rpm = [1, 500, 6131, 6269, 6410, 6537, 6671, 6823, 6947, 7086, 7200, 7335, 7479, 7597, 7745, 7876, 8043, 8202, 8366, 8519, 8685, 8863, 9026, 9175, 9306, 9456, 9620, 9749, 9901, 10058, 10214, 10387, 10530, 10679, 10826, 10983, 11106, 11245, 11374, 11497, 11609, 11641, 11656, 11665, 11668];
	        tor = [0.1, 2, 36.87810746, 37.5419134, 38.50074419, 39.09079391, 38.79576905, 38.35323176, 37.68942583, 36.80435125, 36.5830826, 36.5830826, 36.65683882, 36.73059503, 36.87810746, 37.68942583, 39.38581877, 41.524749, 43.58992302, 45.06504732, 45.58134083, 45.72885326, 45.4338284, 44.84377868, 43.66367924, 42.70484844, 41.89353008, 41.524749, 41.59850522, 41.96728629, 42.4835798, 42.9998733, 43.07362952, 42.85236087, 42.18855494, 41.30348036, 40.41840578, 39.38581877, 38.13196312, 37.02561989, 35.25547074, 28.69116761, 19.39788453, 10.91591981, 5.605472335];
	    }


	}

	if (engineChoice === 2) {

	    //Transmission and shifting parameters

	    nt_1 = 27 / 15;                 //1st Gear Ratio
	    nt_2 = 25 / 17;                 //2nd Gear Ratio
	    nt_3 = 21 / 17;                 //3rd Gear Ratio
	    nt_4 = 21 / 20;                 //4th Gear Ratio
	    nt_5 = 20 / 22;                 //5th Gear Ratio
	    nt_6 = 20 / 22;                 //6th Gear Ratio
	    nf = parseFloat(inNf);          //Final Gear Ratio
	    np = 63 / 23;                   //Primary Reduction Ratio

	    if (torqueCurveChoice === 1) {
	        rpm = [1, 500, 5910, 6074, 6254, 6438, 6615, 6804, 6991, 7158, 7332, 7498, 7646, 7816, 7973, 8122, 8274, 8415, 8566, 8705, 8847, 8981, 9101, 9219, 9314, 9391, 9452, 9502];
	        tor = [0.1, 2, 32.74775943, 32.67400321, 32.74775943, 33.1165405, 33.78034644, 34.22288373, 34.73917723, 34.5916648, 33.92785887, 33.19029672, 32.15770971, 31.1251227, 30.46131677, 29.87126705, 29.20746111, 28.69116761, 28.39614275, 28.1748741, 27.5848243, 26.6997498, 25.66716279, 24.11828228, 21.97935205, 19.10285967, 16.00509864, 12.90733761];
	    }

	}

	if (engineChoice === 3) {

	    //Transmission and shifting parameters

	    nt_1 = 38 / 13;                 //1st Gear Ratio
	    nt_2 = 33 / 16;                 //2nd Gear Ratio
	    nt_3 = 31 / 19;                 //3rd Gear Ratio
	    nt_4 = 29 / 21;                 //4th Gear Ratio
	    nt_5 = 28 / 23;                 //5th Gear Ratio
	    nt_6 = 26 / 24;                 //6th Gear Ratio
	    nf = parseFloat(inNf);          //Final Gear Ratio
	    np = 89 / 44;                    //Primary Reduction Ratio

	    if (torqueCurveChoice === 1) {
	        rpm = [1, 500, 7205, 7390, 7578, 7775, 7998, 8225, 8454, 8661, 8860, 9068, 9265, 9483, 9699, 9917, 10143, 10347, 10531, 10700, 10855, 10983, 11094, 11196, 11276, 11361, 11440, 11520, 11605, 11700, 11795, 11885, 11989, 12098, 12216, 12351, 12491, 12625, 12722, 12792, 12838, 12870, 12893, 12909, 12921, 12932, 12938, 12939, 12943, 12944];
	        tor = [0.1, 2, 35.5504956, 35.99303288, 36.65683882, 38.50074419, 40.04962471, 42.33606737, 44.4749976, 44.69626625, 43.51616681, 41.59850522, 40.492162, 40.71343064, 42.85236087, 44.17997274, 44.4749976, 44.10621653, 42.63109223, 40.63967443, 37.46815718, 34.5916648, 30.83009784, 26.99477466, 23.45447635, 21.09427747, 19.02910345, 17.7752478, 17.18519808, 17.33271051, 17.99651644, 18.43905373, 18.88159102, 19.91417803, 21.38930233, 23.01193906, 25.00335686, 26.33096873, 25.88843144, 22.05310826, 17.92276023, 13.20236247, 9.957089015, 7.080596633, 5.38420369, 4.056591821, 3.245273457, 2.360198878, 1.475124299, 1.180099439];

	    }

	}

	if (engineChoice === 4) {

	    //Transmission and shifting parameters

	    nt_1 = 36 / 14;                 //1st Gear Ratio
	    nt_2 = 32 / 17;              //2nd Gear Ratio
	    nt_3 = 28 / 19;                 //3rd Gear Ratio
	    nt_4 = 26 / 22;                //4th Gear Ratio
	    nt_5 = 23 / 24;                 //5th Gear Ratio
	    nt_6 = 21 / 26;                 //6th Gear Ratio
	    nf = parseFloat(inNf);          //Final Gear Ratio
	    np = 76 / 32;                    //Primary Reduction Ratio

	    if (torqueCurveChoice === 1) {
	        rpm = [1, 500, 5632, 5752, 5869, 5992, 6107, 6239, 6356, 6478, 6587, 6705, 6816, 6931, 7040, 7149, 7255, 7365, 7468, 7569, 7674, 7782, 7881, 7979, 8070, 8168, 8260, 8345, 8446, 8529, 8626, 8714, 8799, 8887, 8970, 9055, 9132, 9219, 9297, 9358, 9414, 9465, 9494, 9517, 9530, 9559, 9588, 9593, 9601, 9611, 9615];
	        tor = [0.1, 2, 38.35323176, 38.50074419, 38.79576905, 39.0170377, 39.23830634, 39.23830634, 39.09079391, 38.86952527, 38.35323176, 37.76318204, 37.32064475, 36.73059503, 36.28805774, 35.84552045, 35.47673938, 34.96044588, 34.5916648, 34.14912751, 33.78034644, 33.41156536, 33.1165405, 32.74775943, 32.37897835, 31.86268485, 31.34639134, 30.75634162, 30.38756055, 30.01877948, 29.72375462, 29.35497354, 29.05994868, 28.69116761, 28.32238653, 27.95360546, 27.6585806, 27.43731195, 27.06853088, 26.10970008, 23.60198878, 20.72549639, 18.21778509, 13.86616841, 10.17835766, 8.260696072, 8.776989576, 8.481964717, 6.195522054, 4.425372896, 3.835323176];

	    }

	}

	if (engineChoice === 5) {

	    //Transmission and shifting parameters

	    nt_1 = 35 / 14;                 //1st Gear Ratio
	    nt_2 = 24 / 15;               //2nd Gear Ratio
	    nt_3 = 21 / 18;                 //3rd Gear Ratio
	    nt_4 = 19 / 20;                 //4th Gear Ratio
	    nt_5 = 18 / 22;                  //5th Gear Ratio
	    nt_6 = 18 / 22;          //6th Gear Ratio
	    nf = parseFloat(inNf);          //Final Gear Ratio
	    np = 81 / 30;                    //Primary Reduction Ratio

	    if (torqueCurveChoice === 1) {
	        rpm = [1, 500, 6713, 6894, 7063, 7230, 7388, 7553, 7722, 7887, 8053, 8230, 8410, 8581, 8766, 8954, 9137, 9325, 9507, 9686, 9856, 10029, 10190, 10355, 10516, 10671, 10826, 10981, 11132, 11279, 11420, 11569, 11699, 11831, 11955, 12072, 12193, 12299, 12415, 12520, 12618, 12716, 12803, 12889, 12978, 13065, 13159, 13251, 13340, 13427, 13514, 13595, 13672, 13743, 13819, 13843, 13846, 13867];
	        tor = [0.1, 2, 35.62425181, 37.32064475, 38.50074419, 38.86952527, 38.42698798, 37.91069447, 37.5419134, 37.76318204, 38.35323176, 38.64825662, 39.16455013, 39.75459985, 40.56591821, 41.524749, 42.26231115, 42.40982358, 42.40982358, 42.11479872, 41.67226143, 40.78718686, 39.97586849, 38.94328148, 38.0582069, 37.32064475, 36.73059503, 36.21430153, 35.47673938, 34.88668966, 34.37039616, 33.70659022, 32.89527186, 31.78892863, 30.75634162, 29.72375462, 28.54365518, 27.58482438, 26.40472494, 25.51965036, 24.78208822, 23.67574499, 22.64315798, 21.68432719, 21.16803368, 20.57798396, 20.35671532, 21.02052125, 21.2417899, 20.87300882, 20.57798396, 19.91417803, 19.10285967, 18.51280995, 17.70149158, 14.60373056, 9.219526866, 4.277860466];
	    }

	}

	if (engineChoice === 6) {

	    //Transmission and shifting parameters

	    nt_1 = 29 / 12;                 //1st Gear Ratio
	    nt_2 = 26 / 15;                 //2nd Gear Ratio
	    nt_3 = 21 / 16;                 //3rd Gear Ratio
	    nt_4 = 21 / 20;                 //4th Gear Ratio
	    nt_5 = 21 / 25;                 //5th Gear Ratio
	    nt_6 = 21 / 35;                 //6th Gear Ratio
	    nf = parseFloat(inNf);          //Final Gear Ratio
	    np = 62 / 22;                   //Primary Reduction Ratio

	    if (torqueCurveChoice === 1) {
	        rpm = [1, 500, 6031, 6227, 6415, 6591, 6772, 6956, 7142, 7316, 7492, 7660, 7832, 8019, 8183, 8358, 8518, 8668, 8815, 8965, 9110, 9244, 9393, 9513, 9584, 9623, 9653, 9671];
	        tor = [0.1, 2, 53.54701204, 56.27599199, 56.49726063, 55.8334547, 55.02213634, 54.72711148, 54.57959905, 54.06330554, 53.03071853, 52.44066881, 51.77686288, 52.44066881, 52.14564395, 51.92437531, 51.18681316, 49.85920129, 48.08905213, 46.4664154, 45.21255975, 44.03246031, 42.92611709, 42.40982358, 38.50074419, 30.53507298, 22.05310826, 15.12002406];
	    }

	}

	if (engineChoice === 7) {

	    //Transmission and shifting parameters

	    nt_1 = 39 / 14;                 //1st Gear Ratio
	    nt_2 = 39 / 19;                 //2nd Gear Ratio
	    nt_3 = 37 / 22;                 //3rd Gear Ratio
	    nt_4 = 29 / 20;                 //4th Gear Ratio
	    nt_5 = 30 / 23;                 //5th Gear Ratio
	    nt_6 = 26 / 22;                 //6th Gear Ratio
	    nf = parseFloat(inNf);          //Final Gear Ratio
	    np = 79 / 41;                   //Primary Reduction Ratio

	    if (torqueCurveChoice === 1) {
	        rpm = [1, 500, 6119, 6256, 6396, 6533, 6676, 6816, 6965, 7094, 7228, 7360, 7496, 7629, 7751, 7903, 8049, 8194, 8328, 8491, 8643, 8789, 8941, 9091, 9245, 9397, 9561, 9720, 9886, 10053, 10220, 10384, 10549, 10711, 10879, 11036, 11194, 11351, 11501, 11632, 11786, 11923, 12059, 12191, 12315, 12437, 12541, 12670, 12777, 12878, 12962, 13066, 13151, 13243, , 13416, 13489, 13555, 13575, 13598, 13615, 13635, 13651, 13687, 13739, 13769, 13797, 13808, 13829, 13848, 13856, 13870];
	        tor = [0.1, 2, 35.91927667, 36.36181396, 37.02561989, 37.76318204, 38.27947555, 38.50074419, 38.20571933, 37.91069447, 37.24688854, 36.73059503, 36.5830826, 36.80435125, 37.17313232, 37.46815718, 38.0582069, 39.09079391, 39.97586849, 39.97586849, 40.04962471, 40.12338092, 40.34464957, 40.56591821, 41.0084555, 41.524749, 42.18855494, 42.85236087, 43.51616681, 44.10621653, 44.62251003, 45.06504732, 45.28631597, 45.13880354, 44.84377868, 44.4749976, 44.03246031, 43.51616681, 42.70484844, 41.81977386, 41.0084555, 39.90211228, 38.57450041, 37.24688854, 36.21430153, 35.25547074, 34.0753713, 32.600247, 31.6414162, 30.68258541, 28.10111789, 26.1834563, 26.77350602, 25.22462551, 25.00335686, 24.26579471, 23.45447635, 21.68432719, 16.816417, 10.54713873, 7.744402567, 6.416790699, 5.900497194, 6.269278269, 10.3996263, 11.80099439, 9.8833328, 7.375621493, 5.826740979, 5.974253409, 5.310447475, 4.277860466];
	    }
	    if (torqueCurveChoice === 2) {
	        rpm = [1, 500, 5876, 6012, 6132, 6257, 6377, 6508, 6642, 6786, 6931, 7076, 7207, 7339, 7472, 7599, 7733, 7877, 8030, 8187, 8353, 8510, 8671, 8826, 8988, 9155, 9301, 9454, 9605, 9749, 9888, 10029, 10167, 10302, 10438, 10581, 10718, 10856, 10999, 11134, 11277, 11415, 11552, 11684, 11817, 11951, 12085, 12211, 12336, 12463, 12585, 12701, 12816, 12932, 13038, 13149, 13252, 13348, 13424, 13475, 13516, 13549, 13563, 13573];
	        tor = [0.1, 2, 41.89353008, 41.22972414, 39.82835606, 38.57450041, 37.61566961, 37.39440097, 38.13196312, 40.41840578, 42.18855494, 43.22114195, 42.70484844, 41.89353008, 40.56591821, 39.45957499, 38.64825662, 40.27089335, 43.14738573, 45.28631597, 46.09763433, 46.90895269, 47.4252462, 47.57275863, 47.79402727, 48.08905213, 48.01529592, 47.35148998, 46.61392783, 45.4338284, 44.40124139, 43.44241059, 42.4835798, 41.81977386, 41.59850522, 41.74601765, 41.96728629, 41.96728629, 42.04104251, 42.11479872, 41.96728629, 41.89353008, 41.59850522, 41.15596793, 40.93469928, 40.63967443, 40.27089335, 39.75459985, 39.09079391, 38.57450041, 37.76318204, 37.02561989, 36.21430153, 35.47673938, 34.51790859, 33.19029672, 32.74775943, 31.6414162, 28.24863032, 22.79067041, 18.14402887, 14.3087057, 10.76840738, 6.711815558];
	    }
	    if (torqueCurveChoice === 3) {
	        rpm = [1, 500, 5267, 5419, 5592, 5752, 5908, 6069, 6241, 6404, 6572, 6739, 6905, 7083, 7250, 7398, 7562, 7739, 7921, 8094, 8269, 8443, 8613, 8777, 8939, 9102, 9265, 9432, 9589, 9738, 9897, 10050, 10198, 10350, 10505, 10646, 10794, 10939, 11093, 11233, 11373, 11510, 11646, 11783, 11911, 12042, 12164, 12281, 12398, 12505, 12615, 12723, 12826, 12911, 12978, 13034, 13084, 13117, 13141, 13159, 13175];
	        tor = [0.1, 2, 18.2915413, 22.79067041, 29.57624219, 33.78034644, 35.62425181, 37.32064475, 38.27947555, 38.86952527, 39.45957499, 39.97586849, 40.27089335, 40.34464957, 40.34464957, 39.82835606, 39.16455013, 39.97586849, 41.15596793, 42.11479872, 41.89353008, 41.37723657, 40.93469928, 40.492162, 40.12338092, 39.75459985, 39.68084363, 39.38581877, 39.0170377, 38.50074419, 38.13196312, 37.61566961, 37.17313232, 36.87810746, 36.5830826, 36.21430153, 35.99303288, 35.77176424, 35.47673938, 35.10795831, 34.66542102, 34.14912751, 33.48532158, 32.96902807, 32.45273457, 32.45273457, 31.78892863, 30.90385405, 30.16629191, 29.28121733, 28.39614275, 27.51106817, 26.84726223, 26.03594387, 24.63457579, 22.05310826, 18.73407859, 15.71007378, 12.90733761, 9.73582037, 7.523133923, 5.974253409];
	    }
	    if (torqueCurveChoice === 4) {
	        rpm = [1, 500, 6828, 7032, 7238, 7447, 7677, 7916, 8139, 8358, 8588, 8800, 9008, 9197, 9380, 9545, 9715, 9868, 10028, 10194, 10362, 10530, 10707, 10891, 11087, 11263, 11453, 11648, 11832, 12004, 12177, 12345, 12498, 12648, 12783, 12918, 13046, 13163, 13273, 13374, 13472, 13559, 13621, 13677, 13695];
	        tor = [0.1, 2, 24.708332, 28.91243625, 36.65683882, 39.75459985, 43.22114195, 46.68768405, 47.9415397, 48.38407699, 48.23656456, 47.72027106, 46.31890297, 44.10621653, 41.08221171, 38.27947555, 37.17313232, 36.28805774, 35.18171452, 34.73917723, 34.96044588, 35.40298317, 35.99303288, 36.73059503, 37.76318204, 39.0170377, 39.82835606, 40.04962471, 39.82835606, 39.38581877, 38.72201284, 37.5419134, 36.21430153, 34.22288373, 32.52649078, 30.68258541, 29.28121733, 28.10111789, 26.55223737, 24.708332, 22.93818284, 21.61057097, 19.02910345, 15.78382999, 12.75982518];
	    }
	    if (torqueCurveChoice === 5) {
	        rpm = [1, 500, 6206, 6303, 6406, 6515, 6631, 6748, 6876, 6999, 7109, 7246, 7356, 7480, 7599, 7708, 7826, 7946, 8062, 8178, 8297, 8417, 8540, 8665, 8786, 8909, 9040, 9172, 9309, 9442, 9577, 9709, 9851, 9984, 10126, 10250, 10390, 10543, 10676, 10806, 10940, 11067, 11196, 11320, 11443, 11563, , 11795, 11914, 12028, 12140, 12245, 12348, 12453, 12557, 12659, 12747, 12851, 12944, 13031, 13118, 13200, 13278, 13355, 13433, 13508, 13585, 13656, 13728, 13796, 13872, 13941, 14010, 14094, 14166, 14201, 14215, 14220];
	        tor = [0.1, 2, 23.89701364, 24.708332, 25.51965036, 26.47848116, 27.80609303, 29.1337049, 30.68258541, 31.49390377, 31.93644106, 31.78892863, 31.6414162, 31.27263513, 30.90385405, 30.68258541, 30.46131677, 30.38756055, 30.46131677, 30.46131677, 30.53507298, 30.97761027, 31.42014756, 31.78892863, 32.08395349, 32.30522214, 32.600247, 33.04278429, 33.70659022, 34.22288373, 34.66542102, 34.88668966, 35.10795831, 35.40298317, 35.62425181, 35.84552045, 35.99303288, 36.0667891, 35.99303288, 35.69800803, 35.25547074, 34.81293345, 34.14912751, 33.70659022, 33.26405293, 32.82151564, 32.23146592, 31.71517242, 31.05136648, 30.53507298, 29.94502326, 29.42872976, 28.83868004, 28.32238653, 27.87984924, 27.36355574, 26.77350602, 26.10970008, 25.51965036, 24.63457579, 23.82325742, 23.2332077, 22.20062069, 21.31554611, 20.72549639, 20.79925261, 20.43047153, 19.91417803, 19.61915317, 19.47164074, 19.2503721, 19.02910345, 18.95534724, 19.32412831, 20.50422775, 17.40646672, 10.98967602, 6.933084203];
	    }

	}

	if (engineChoice === 8) {

	    //Transmission and shifting parameters

	    nt_1 = 31 / 12;                 //1st Gear Ratio
	    nt_2 = 32 / 16;                 //2nd Gear Ratio
	    nt_3 = 30 / 18;                 //3rd Gear Ratio
	    nt_4 = 26 / 18;                 //4th Gear Ratio
	    nt_5 = 27 / 21;                 //5th Gear Ratio
	    nt_6 = 23 / 20;                 //6th Gear Ratio
	    nf = parseFloat(inNf);          //Final Gear Ratio
	    np = 85 / 41;                   //Primary Reduction Ratio

	    if (torqueCurveChoice === 1) {
	        rpm = [1, 500, 5956, 6114, 6268, 6436, 6603, 6770, 6955, 7132, 7313, 7497, 7678, 7862, 8042, 8221, 8406, 8575, 8748, 8929, 9106, 9264, 9432, 9584, 9741, 9895, 10044, 10199, 10342, 10490, 10629, 10764, 10849, 10922, 10955, 10969, 10979, 10986, 10987];
	        tor = [0.1, 2, 44.4749976, 44.77002246, 45.36007218, 46.68768405, 47.57275863, 49.12163914, 50.52300723, 51.85061909, 52.66193746, 53.39949961, 53.91579311, 54.06330554, 54.06330554, 53.69452447, 53.47325582, 53.17823096, 52.80944989, 52.58818124, 52.21940017, 51.55559423, 50.3754948, 49.19539536, 48.08905213, 47.13022134, 46.31890297, 45.65509704, 44.91753489, 44.17997274, 43.44241059, 42.40982358, 39.38581877, 31.78892863, 25.07711308, 16.37387971, 9.957089015, 6.859327988, 4.277860466];
	    }
	    if (torqueCurveChoice === 2) {
	        rpm = [1, 500, 4630, 4724, 4829, 4933, 5037, 5150, 5263, 5372, 5476, 5574, 5677, 5784, 5881, 6008, 6129, 6251, 6383, 6511, 6646, 6804, 6950, 7092, 7240, 7390, 7534, 7681, 7824, 7971, 8124, 8274, 8417, 8567, 8692, 8844, 8982, 9122, 9245, 9385, 9513, 9641, 9783, 9906, 10035, 10156, 10286, 10411, 10536, 10654, 10775, 10895, 11013, 11119, 11247, 11360, 11469, 11576, 11686, 11789, 11890, 11988, 12088, 12182, 12257, 12288, 12291];
	        tor = [0.1, 2, 24.56081957, 26.99477466, 28.76492382, 30.01877948, 30.75634162, 31.93644106, 33.19029672, 33.92785887, 34.00161508, 33.63283401, 33.1165405, 32.600247, 32.74775943, 33.41156536, 34.44415237, 35.69800803, 37.46815718, 38.72201284, 40.41840578, 43.07362952, 44.99129111, 45.58134083, 45.87636569, 45.87636569, 45.58134083, 45.58134083, 45.58134083, 45.65509704, 45.80260947, 46.02387811, 46.17139054, 46.02387811, 45.58134083, 45.13880354, 44.54875382, 43.88494788, 43.29489816, 42.77860466, 42.18855494, 41.59850522, 41.15596793, 40.71343064, 40.34464957, 40.04962471, 39.68084363, 39.45957499, 39.16455013, 38.79576905, 38.27947555, 37.91069447, 37.61566961, 37.39440097, 36.95186368, 36.36181396, 35.69800803, 35.18171452, 34.44415237, 33.48532158, 32.82151564, 32.37897835, 31.71517242, 31.05136648, 29.50248597, 23.67574499, 14.23494948];
	    }
	    if (torqueCurveChoice === 3) {
	        rpm = [1, 500, 5062, 5151, 5242, 5326, 5404, 5471, 5539, 5595, 5657, 5700, 5760, 5821, 5872, 5939, 6015, 6090, 6168, 6254, 6334, 6422, 6508, 6588, 6658, 6740, 6829, 6917, 6997, 7086, 7169, 7255, 7355, 7452, 7546, 7645, 7753, 7870, 7987, 8119, 8256, 8402, 8549, 8708, 8863, 9027, 9185, 9340, 9492, 9654, 9807, 9972, 10118, 10274, 10426, 10573, 10712, 10845, 10975, 11108, 11235, 11363, 11481, 11590, 11689, 11761, 11811, 11846, 11871, 11878, 11888, 11890];
	        tor = [0.1, 2, 26.1834563, 25.59340658, 24.63457579, 23.38072013, 22.64315798, 21.38930233, 20.28295911, 18.36529752, 16.74266079, 15.85758621, 14.8249992, 14.16119327, 14.08743705, 14.52997434, 15.26753649, 16.59514836, 17.70149158, 18.80783481, 19.7666656, 20.50422775, 20.87300882, 21.02052125, 20.50422775, 19.7666656, 20.43047153, 20.94676504, 21.09427747, 21.16803368, 21.16803368, 21.38930233, 22.27437691, 23.15945149, 23.45447635, 23.89701364, 24.92960065, 26.1834563, 27.80609303, 29.87126705, 32.30522214, 33.85410265, 35.77176424, 37.17313232, 38.50074419, 39.0170377, 39.45957499, 39.45957499, 39.23830634, 39.16455013, 39.45957499, 39.68084363, 39.38581877, 38.79576905, 38.35323176, 37.83693826, 37.02561989, 35.99303288, 34.96044588, 34.29663994, 33.63283401, 32.89527186, 31.86268485, 30.38756055, 28.83868004, 25.59340658, 21.68432719, 16.66890457, 12.24353168, 9.072014436, 4.499129111, 1.253855654];

	    }
	    if (torqueCurveChoice === 4) {
	        rpm = [1, 500, 3742, 3790, 3838, 3877, 3923, 3973, 4023, 4074, 4129, 4180, 4235, 4289, 4346, 4400, 4452, 4507, 4561, 4614, 4673, 4738, 4810, 4878, 4953, 5029, 5111, 5195, 5276, 5362, 5450, 5540, 5619, 5706, 5794, 5872, 5962, 6054, 6144, 6236, 6325, 6413, 6508, 6600, 6691, 6782, 6873, 6964, 7052, 7141, 7226, 7313, 7402, 7486, 7577, 7667, 7745, 7829, 7912, 7999, 8083, 8166, 8260, 8351, 8449, 8552, 8653, 8758, 8859, 8964, 9065, 9166, 9263, 9363, 9462, 9561, 9661, 9755, 9854, 9954, 10048, 10140, 10237, 10332, 10421, 10516, 10612, 10702, 10795, 10886, 10976, 11069, 11155, 11245, 11332, 11416, 11501, 11588, 11668, 11750, 11830, 11909, 11988, 12067, 12144, 12221, 12297, 12378, 12428, 12455, 12468, 12486, 12503];
	        tor = [0.1, 2, 18.66032238, 19.17661588, 19.54539696, 19.91417803, 19.98793425, 20.28295911, 20.65174018, 21.53681476, 22.42188934, 23.01193906, 23.52823256, 23.74950121, 24.33955093, 24.708332, 24.85584443, 24.48706336, 24.1920385, 24.26579471, 24.708332, 25.96218765, 27.6585806, 29.28121733, 30.83009784, 32.37897835, 34.00161508, 34.88668966, 35.91927667, 37.17313232, 38.20571933, 39.38581877, 39.68084363, 39.0170377, 38.27947555, 38.0582069, 38.42698798, 39.23830634, 40.19713714, 41.0084555, 41.22972414, 41.08221171, 41.0084555, 41.22972414, 41.45099279, 41.59850522, 41.59850522, 41.45099279, 41.15596793, 40.78718686, 40.34464957, 39.75459985, 39.38581877, 39.38581877, 39.5333312, 39.75459985, 39.60708742, 39.0170377, 38.57450041, 38.27947555, 38.0582069, 38.20571933, 38.79576905, 40.04962471, 41.89353008, 43.44241059, 45.06504732, 46.09763433, 46.68768405, 46.68768405, 46.31890297, 45.87636569, 45.58134083, 45.28631597, 45.21255975, 45.13880354, 45.13880354, 44.99129111, 44.69626625, 44.62251003, 44.25372896, 43.81119167, 43.51616681, 43.22114195, 43.07362952, 42.85236087, 42.55733601, 42.33606737, 42.26231115, 42.18855494, 41.96728629, 41.74601765, 41.22972414, 40.86094307, 40.41840578, 39.90211228, 39.38581877, 39.0170377, 38.57450041, 38.0582069, 37.46815718, 37.02561989, 36.80435125, 36.65683882, 36.50932639, 36.28805774, 36.14054531, 35.99303288, 35.32922695, 31.42014756, 24.78208822, 18.66032238, 14.60373056];
	    }

	}

	if (engineChoice === 9) {

	    //Transmission and shifting parameters

	    nt_1 = 29 / 12;                 //1st Gear Ratio
	    nt_2 = 26 / 15;                 //2nd Gear Ratio
	    nt_3 = 21 / 16;                 //3rd Gear Ratio
	    nt_4 = 21 / 20;                 //4th Gear Ratio
	    nt_5 = 21 / 25;                 //5th Gear Ratio
	    nt_6 = 21 / 35;                 //6th Gear Ratio
	    nf = parseFloat(inNf);          //Final Gear Ratio
	    np = 62 / 22;                   //Primary Reduction Ratio

	    if (torqueCurveChoice === 1) {
        //Done
	        rpm = [1, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4090, 4155, 4226, 4284, 4349, 4405, 4472, 4531, 4594, 4657, 4723, 4786, 4844, 4916, 4980, 5046, 5112, 5179, 5246, 5313, 5386, 5459, 5527, 5603, 5681, 5760, 5839, 5925, 6002, 6088, 6167, 6246, 6329, 6410, 6486, 6565, 6648, 6732, 6814, 6897, 6977, 7064, 7142, 7234, 7322, 7402, 7489, 7573, 7657, 7739, 7824, 7902, 7988, 8065, 8144, 8220, 8299, 8377, 8446, 8521, 8595, 8667, 8743, 8813, 8880, 8952, 9017, 9084, 9142, 9210, 9270, 9332, 9391, 9452, 9509, 9567, 9620, 9683, 9739, 9794, 9840, 9894, 9941, 9987, 10034, 10057, 10096, 10106, 10111, 10112];
	        tor = [.25, 2.0054, 2.3754, 2.7454, 3.1154, 3.4854, 3.8554, 4.2254, 4.5954, 4.9654, 5.3354, 5.7054, 6.0754, 6.4454, 6.8154, 7.1854, 7.5554, 7.9254, 8.2954, 8.6654, 9.0354, 9.4054, 9.7754, 10.1454, 10.5154, 10.8854, 11.2554, 11.6254, 11.9954, 12.3654, 12.7354, 13.1054, 13.4754, 13.8454, 14.2154, 14.5854, 14.9554, 15.244, 18.722, 22.274, 23.754, 24.642, 24.79, 24.79, 24.864, 24.938, 25.012, 25.086, 25.308, 25.53, 25.678, 25.9, 26.122, 26.344, 26.418, 26.566, 26.714, 27.01, 27.454, 27.972, 28.638, 29.674, 30.266, 31.006, 31.894, 32.264, 32.56, 32.634, 32.56, 32.412, 32.19, 32.116, 32.116, 32.264, 32.56, 32.782, 32.856, 33.152, 33.340, 33.892, 34.04, 34.114, 34.114, 34.114, 34.114, 34.04, 33.744, 33.596, 33.374, 33.004, 32.708, 32.486, 32.042, 31.672, 31.228, 30.932, 30.562, 30.266, 29.97, 29.674, 29.23, 28.934, 28.49, 28.046, 27.676, 27.084, 26.566, 25.9, 25.382, 25.086, 24.642, 24.272, 23.902, 23.606, 23.458, 23.384, 23.014, 22.422, 21.386, 20.572, 19.832, 19.462, 18.056, 15.244, 13.394, 9.028, 4.144];
    
        }

	}

    //Linearly increasing the horsepower via Forced induction
	forcedInductionIncrease = parseFloat(inForcedInduction);

	if (forcedInductionIncrease > 0) {
	    if (forcedInductionIncrease === 2) {
	        forcedInductionIncrease = 5;
	    }
	    if (forcedInductionIncrease === 3) {
	        forcedInductionIncrease = 10;
	    }
	    if (forcedInductionIncrease === 4) {
	        forcedInductionIncrease = 20;
	    }
	    if (forcedInductionIncrease === 5) {
	        forcedInductionIncrease = 30;
	    }
	    if (forcedInductionIncrease === 6) {
	        forcedInductionIncrease = 10;
	    }
	    if (forcedInductionIncrease === 7) {
	        forcedInductionIncrease = 20;
	    }


	    for (var uu = 0; uu < tor.length; uu++) {

	        tor[uu] = (((tor[uu] * rpm[uu] / 5252) + forcedInductionIncrease) * 5252) / rpm[uu];

	    }
	    tor[0] = 0.25;
	    tor[1] = 1;
	}

	//Max Horsepower Calculation
	maxHorsepower = 0;
	for (var u = 0; u < tor.length; u++) {
	maxHorsepowerCalc = (tor[u] * rpm[u]) / 5252;
	if (maxHorsepower < maxHorsepowerCalc) {

	maxHorsepower = maxHorsepowerCalc;

	}
    
	}
	maxHorsepower = roundNumber(maxHorsepower, 2);



	//Downforce Calculations based on User Input
	//downforceVelocityInput is the velocity at which the user is inputting the downforce
	//40mph = (40 * 1.46666667) ft/s

	var LiftDragCoefficient = parseFloat(inLiftDragCoefficient);
	var downforceAtVelocity = parseFloat(inDownforceAtVelocityInput);
	downforceVelocityInput = 40 * 1.46666667;
	downforceConstant = downforceAtVelocity / (Math.pow(downforceVelocityInput, 2));



	/*

	Calculating BSFC:
	BSFC=r/P
	P=power produced, P=ax*vehicle velocity
	r is the fuel consumption rate(g-s^-1)
	BSFC is a constant for most engines
	r*total time = fuel used.
	// http://en.wikipedia.org/wiki/Brake_specific_fuel_consumption 
	Link states that the typical BSFC for gasoline engines is 0.45 to 0.37
	Taking the average of 0.45 and 0.37 yields a BSFC of .41
	BSFC is in units of lb/hp*h
	ax=ft/s2*ft/s

	*/

	bsfc = .41 / 3600;


	/*

	FSAE West 2007 Endurance Track Map
	Endurance Event is 13.66 miles, total distance/track distance will equal
	number of laps.

	*/

	numberOfLaps = 20;

	//Straight Section Distances

	straightSections = 18;
	dist[0] = 48.33;
	dist[1] = 15.41;
	dist[2] = 14.75;
	dist[3] = 193.33;
	dist[4] = 169.25;
	dist[5] = 25;
	dist[6] = 51.6667;
	dist[7] = 25.833;
	dist[8] = 92.5;
	dist[9] = 13.1667;
	dist[10] = 12.33;
	dist[11] = 92.5;
	dist[12] = 41.08;
	dist[13] = 31.08;
	dist[14] = 82.667;
	dist[15] = 30.416;
	dist[16] = 136.5;
	dist[17] = 74.1667;


	//Associated Corners

	numberOfCorners[0] = 4;
	cornerRadius[0] = 36.25;
	corneringDist[0] = 81.083;
	cornerRadius[1] = 30.833;
	corneringDist[1] = 42.91667;
	cornerRadius[2] = 57.5;
	corneringDist[2] = 47.1667;
	cornerRadius[3] = 105.833;
	corneringDist[3] = 82.75;


	numberOfCorners[1] = 1;
	cornerRadius[4] = 47.5;
	corneringDist[4] = 58.41667;

	numberOfCorners[2] = 3;
	cornerRadius[5] = 37.5;
	corneringDist[5] = 89.1667;
	cornerRadius[6] = 38.33;
	corneringDist[6] = 78.41667;
	cornerRadius[7] = 53.333;
	corneringDist[7] = 46.75;

	numberOfCorners[3] = 7;
	cornerRadius[8] = 56.25;
	corneringDist[8] = 50.41667;
	cornerRadius[9] = 56.25;
	corneringDist[9] = 50.41667;
	cornerRadius[10] = 19.1667;
	corneringDist[10] = 33.91667;
	cornerRadius[11] = 19.16667;
	corneringDist[11] = 33.91667;
	cornerRadius[12] = 19.16667;
	corneringDist[12] = 33.91667;
	cornerRadius[13] = 19.1667;
	corneringDist[13] = 33.91667;
	cornerRadius[14] = 19.1667;
	corneringDist[14] = 15.41667;

	numberOfCorners[4] = 3;
	cornerRadius[15] = 42.5;
	corneringDist[15] = 52.41667;
	cornerRadius[16] = 42.5;
	corneringDist[16] = 96.5833;
	cornerRadius[17] = 33.333;
	corneringDist[17] = 72.5833;

	numberOfCorners[5] = 3;
	cornerRadius[18] = 105;
	corneringDist[18] = 134;
	cornerRadius[19] = 66.6667;
	corneringDist[19] = 94.0833;
	cornerRadius[20] = 29.16667;
	corneringDist[20] = 64.6667;

	numberOfCorners[6] = 1;
	cornerRadius[21] = 47.9167;
	corneringDist[21] = 121.33;

	numberOfCorners[7] = 1;
	cornerRadius[22] = 51.667;
	corneringDist[22] = 130.8333;

	numberOfCorners[8] = 2;
	cornerRadius[23] = 38.33;
	corneringDist[23] = 36.833;
	cornerRadius[24] = 38.33;
	corneringDist[24] = 36.75;

	numberOfCorners[9] = 2;
	cornerRadius[25] = 37.5;
	corneringDist[25] = 37.08;
	cornerRadius[26] = 37.5;
	corneringDist[26] = 31.833;

	numberOfCorners[10] = 2;
	cornerRadius[21] = 54.1667;
	corneringDist[21] = 49.0833;
	cornerRadius[21] = 54.1667;
	corneringDist[21] = 47;

	numberOfCorners[11] = 6;
	cornerRadius[22] = 80;
	corneringDist[22] = 35.75;
	cornerRadius[23] = 29.1667;
	corneringDist[23] = 38.41667;
	cornerRadius[24] = 29.1667;
	corneringDist[24] = 45.75;
	cornerRadius[25] = 29.1667;
	corneringDist[25] = 46.5833;
	cornerRadius[26] = 29.1667;
	corneringDist[26] = 50.833;
	cornerRadius[27] = 29.1667;
	corneringDist[27] = 32.75;

	numberOfCorners[12] = 1;
	cornerRadius[28] = 56.667;
	corneringDist[28] = 52.0833;

	numberOfCorners[13] = 3;
	cornerRadius[29] = 57.5;
	corneringDist[29] = 67.0833;
	cornerRadius[30] = 77.5;
	corneringDist[30] = 33.083;
	cornerRadius[31] = 77.5;
	corneringDist[31] = 38.25;

	numberOfCorners[14] = 2;
	cornerRadius[32] = 72.5;
	corneringDist[32] = 42;
	cornerRadius[33] = 72.5;
	corneringDist[33] = 38.75;

	numberOfCorners[15] = 2;
	cornerRadius[34] = 24.1667;
	corneringDist[34] = 56;
	cornerRadius[35] = 24.1667;
	corneringDist[35] = 57;

	numberOfCorners[16] = 1;
	cornerRadius[36] = 41.667;
	corneringDist[36] = 80.833;

	numberOfCorners[17] = 2;
	cornerRadius[37] = 62.0833;
	corneringDist[37] = 58.16667;
	cornerRadius[38] = 81.25;
	corneringDist[38] = 75.41667;
	cornerRadius[39] = 243.75;
	corneringDist[39] = 91.75;

	//Saves all of the corner radii into one array.

	for (var i = 0; i < 40; i++) {
		cornerTotal[i] = cornerRadius[i];
	}


	//Car Parameters

    vehicleWeight = parseInt(inVehicleWeight);
	coefficientOfDrag = 1;
	frontalArea = 9.24021;
	rho = 0.00236;
	f_0 = 0.01417;
	f_s = 0.01;
	drivetrainEfficiency = 0.85;
	g = 32.17405;
	pi = 3.1416;
	wheelRadius = parseFloat(inWheelRadius);
	wheelRadius = wheelRadius / 12;
	iterationTime = 0.01;
	wheelbase = parseFloat(inWheelBase);
	trackWidth = parseFloat(inTrackWidth);
	trackA = trackWidth / 2;

	hcg = parseFloat(inHcg);
	weightDistribution = parseFloat(inWeightDistribution); //Percent to the front of the car
	cLength = wheelbase * weightDistribution;
	bLength = wheelbase * (1 - weightDistribution);


	shiftTime = 0.01;   //Time required to shift
	shiftRpm = parseInt(inShiftRpm);     //Shift RPM
	

	if (shiftRpm > rpm[rpm.length - 1]) {
	    shiftRpm = rpm[rpm.length - 1] - 50;
    }



	shifts = 0;
	totalShifts = 0;    //Keeps track of the shifts

	/*
	Tire 1=Goodyear 13"
	Tire 2=Hoosier 13"
	Tire 3=Hoosier Small 13"
	Tire 4=Michelin 13"
	*/

	tireChoice = parseInt(inTireChoice);

	/*
	Competition Points Breakdown - FSAE 2011
	Not Currently Used.
	IN PROGRESS
	*/

	competitionAcceleration = 75;
	competitionSkidpad = 50;
	competitionAutocross = 150;
	competitionFuelEconomy = 100;
	competitionEndurance = 300;
	competitionStaticTotal = 325;


	// Vehicle Velocity Calculations

	velocityFactorOne = velocityFactorCalculation(nt_1, nf, np, wheelRadius);
	velocityFactorTwo = velocityFactorCalculation(nt_2, nf, np, wheelRadius);
	velocityFactorThree = velocityFactorCalculation(nt_3, nf, np, wheelRadius);
	velocityFactorFour = velocityFactorCalculation(nt_4, nf, np, wheelRadius);
	velocityFactorFive = velocityFactorCalculation(nt_5, nf, np, wheelRadius);
	velocityFactorSix = velocityFactorCalculation(nt_6, nf, np, wheelRadius);

	/*
	For 0.5 seconds assume that the maximum acceleration that the tires are able to produce with 50% static weight distribution.  
	Then weight transfer will occur.
	*/

	velocity = 0;  //Starting Velocity
	timeLaunching = 0; //Launching Time
	t = iterationTime;
	launchingTime = 0; //Time for maximum tractive acceleration
	launchingDistance = 0;
	ntf = nt_1 * np * nf; //Setup in First Gear
	time = 0;
	timeStraightSection = 0;
	timeOpenThrottle = 0;
	timeTractionLimited = 0;
	brakingTime = 0;
	brakingTracker = 0;
	accelerationTimeTractionLimited = 0;
	accelerationWideOpenThrottleTime = 0;


	/*
	Maximum tire traction for 0.5 seconds without weight transfer.  After 0.5 seconds, there is weight transfer until 15mph(22ft/s).  The car launches from a stop
	and accelerates until it reaches 15mph.  At that point, the software will jump onto the regular sequence and calculate accordingly.  The following
	while loop will record the distances traveled and the time taken to get to 15mph.

	*/


	while (launchingTime <= 0.5) {


	    weightOnTire = vehicleWeight / 2;

		instantaneousCoefficientOfFriction = frictionCalculation(weightOnTire, tireChoice);
		maximumTractiveForce = instantaneousCoefficientOfFriction * weightOnTire;
		maximumLongitudinalAcceleration = maximumTractiveForce / (vehicleWeight / g);

		launchingDistance += velocity * t + 0.5 * maximumLongitudinalAcceleration * Math.pow(t, 2);
		velocity += maximumLongitudinalAcceleration * t;
		launchingTime += t;


	}

        
        launchingShifts = 0;

        velocityFactorCurrent = velocityFactorOne;
        currentRpm = velocityFactorCurrent * velocity;
        ntf = nt_1 * np * nf;

        if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorTwo)) {

            velocityFactorCurrent = velocityFactorTwo;
            currentRpm = velocityFactorCurrent * velocity;
            launchingTime += shiftTime;
            launchingShifts += 1;
            ntf = nt_2 * np * nf;

        }

        if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorThree)) {

            velocityFactorCurrent = velocityFactorThree;
            currentRpm = velocityFactorCurrent * velocity;
            launchingTime += shiftTime;
            launchingShifts += 1;
            ntf = nt_3 * np * nf;

        }

        if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorFour)) {

            velocityFactorCurrent = velocityFactorFour;
            currentRpm = velocityFactorCurrent * velocity;
            launchingTime += shiftTime;
            launchingShifts += 1;
            ntf = nt_4 * np * nf;

        }

        if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorFive)) {

            velocityFactorCurrent = velocityFactorFive;
            currentRpm = velocityFactorCurrent * velocity;
            launchingTime += shiftTime;
            launchingShifts += 1;
            ntf = nt_5 * np * nf;

        }
        if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorSix)) {

            velocityFactorCurrent = velocityFactorSix;
            currentRpm = velocityFactorCurrent * velocity;
            launchingTime += shiftTime;
            launchingShifts += 1;
            ntf = nt_6 * np * nf;

        }


        //This part of the code is where the engine is bouncing off the rev limiter and will not
        //accelerate any further
        if ((currentRpm > shiftRpm) && (velocityFactorCurrent === velocityFactorSix)) {

            currentRpm = shiftRpm;
        }
        






	while (velocity < 22) {
	    instantaneousTorque = torqueLocator(currentRpm, rpm, tor);
	    instantaneousDownforce = downforceCalculation(downforceConstant, velocity);
	    instantaneousRollingResistanceForce = rollingResistanceForceCalculation(f_0, f_s, velocity, vehicleWeight, instantaneousDownforce);
	    instantaneousDragForce = dragForceCalculation(rho, velocity, coefficientOfDrag, frontalArea, instantaneousDownforce, LiftDragCoefficient);

		instantaneousLongitudinalAcceleration = longitudinalAccelerationCalculation(instantaneousTorque, drivetrainEfficiency, wheelRadius, instantaneousDragForce, instantaneousRollingResistanceForce, vehicleWeight, ntf);


		//Weight Distribution and Instantaneous Friction Between Tires and Ground
		weightOnTire = (vehicleWeight * bLength) / (bLength + cLength - hcg * instantaneousLongitudinalAcceleration / g);
		weightOnTire = weightOnTire + instantaneousDownforce / 2;
		instantaneousCoefficientOfFriction = frictionCalculation(weightOnTire, tireChoice);

		//Maximum Tractive Force Calculation

		maximumTractiveForce = instantaneousCoefficientOfFriction * weightOnTire;
		maximumLongitudinalAcceleration = maximumTractiveForce / (vehicleWeight / g);

		if (instantaneousLongitudinalAcceleration > maximumLongitudinalAcceleration) {

		    instantaneousLongitudinalAcceleration = maximumLongitudinalAcceleration;
		    accelerationTimeTractionLimited += t;
		    accelerationWideOpenThrottleTime -= t;

		}

		launchingDistance = launchingDistance + velocity * t + 0.5 * instantaneousLongitudinalAcceleration * Math.pow(t, 2);
		velocity += instantaneousLongitudinalAcceleration * t;
		launchingTime += t;
		accelerationWideOpenThrottleTime += t;

		currentRpm = velocityFactorCurrent * velocity;

		//Max Vehicle Velocity Reached as the current RPM is at Redline and the vehicle is in the last gear
		if ((currentRpm > shiftRpm) && (velocityFactorCurrent === velocityFactorSix)) {

		    velocity -= instantaneousLongitudinalAcceleration * t;
		}

		//This part of the code is where the engine is bouncing off the rev limiter and will not
		//accelerate any further
		if ((currentRpm > shiftRpm) && (velocityFactorCurrent === velocityFactorSix)) {

		    currentRpm = shiftRpm;
		}

		if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorTwo)) {

		    velocityFactorCurrent = velocityFactorTwo;
		    currentRpm = velocityFactorCurrent * velocity;
		    launchingTime += shiftTime;
		    launchingShifts += 1;
		    ntf = nt_2 * np * nf;

		}

		if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorThree)) {

		    velocityFactorCurrent = velocityFactorThree;
		    currentRpm = velocityFactorCurrent * velocity;
		    launchingTime += shiftTime;
		    launchingShifts += 1;
		    ntf = nt_3 * np * nf;

		}

		if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorFour)) {

		    velocityFactorCurrent = velocityFactorFour;
		    currentRpm = velocityFactorCurrent * velocity;
		    launchingTime += shiftTime;
		    launchingShifts += 1;
		    ntf = nt_4 * np * nf;

		}

		if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorFive)) {

		    velocityFactorCurrent = velocityFactorFive;
		    currentRpm = velocityFactorCurrent * velocity;
		    launchingTime += shiftTime;
		    launchingShifts += 1;
		    ntf = nt_5 * np * nf;

		}
		if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorSix)) {

		    velocityFactorCurrent = velocityFactorSix;
		    currentRpm = velocityFactorCurrent * velocity;
		    launchingTime += shiftTime;
		    launchingShifts += 1;
		    ntf = nt_6 * np * nf;

		}

	}

	initialVelocity = velocity;
	currentRpm = velocityFactorCurrent * velocity;

    /*  Acceleration Event Simulation


    */
	//The Acceleration event is 75 meters, which equals about 246.062992 feet.

	instantaneousTorque = torqueLocator(currentRpm, rpm, tor);
	x = 0;
	t = iterationTime;
	count = -1;
	timeStraightSectionAcceleration = 0;
	accelerationShifts = 0;
	accelerationDistance = 246.06-launchingDistance;

	while (x < accelerationDistance) {
			count += 1;


			instantaneousDownforce = downforceCalculation(downforceConstant, velocity);
			instantaneousRollingResistanceForce = rollingResistanceForceCalculation(f_0, f_s, velocity, vehicleWeight, instantaneousDownforce);
			instantaneousDragForce = dragForceCalculation(rho, velocity, coefficientOfDrag, frontalArea, instantaneousDownforce, LiftDragCoefficient);

			//Acceleration Calculation
			instantaneousLongitudinalAcceleration = longitudinalAccelerationCalculation(instantaneousTorque, drivetrainEfficiency, wheelRadius, instantaneousDragForce, instantaneousRollingResistanceForce, vehicleWeight, ntf);


			//Weight Distribution and Instantaneous Friction Between Tires and Ground
			weightOnTire = (vehicleWeight * bLength) / (bLength + cLength - hcg * instantaneousLongitudinalAcceleration / g);
			weightOnTire = weightOnTire + instantaneousDownforce / 2;
			instantaneousCoefficientOfFriction = frictionCalculation(weightOnTire, tireChoice);

			//Maximum Tractive Force Calculation

			maximumTractiveForce = instantaneousCoefficientOfFriction * weightOnTire;
			maximumLongitudinalAcceleration = maximumTractiveForce / (vehicleWeight / g);

			if (instantaneousLongitudinalAcceleration > maximumLongitudinalAcceleration) {

				instantaneousLongitudinalAcceleration = maximumLongitudinalAcceleration;
				accelerationTimeTractionLimited += t;
				accelerationWideOpenThrottleTime -= t;
			}

			x = x + velocity * t + 0.5 * instantaneousLongitudinalAcceleration * Math.pow(t, 2);
			velocity = velocity + instantaneousLongitudinalAcceleration * t;
			accelerationLongitudinalAccelerationStore[count] = instantaneousLongitudinalAcceleration;

			currentRpm = velocityFactorCurrent * velocity;

			//Max Vehicle Velocity Reached as the current RPM is at Redline and the vehicle is in the last gear
			if ((currentRpm > shiftRpm) && (velocityFactorCurrent === velocityFactorSix)) {

			    velocity -= instantaneousLongitudinalAcceleration * t;
			}

			//This part of the code is where the engine is bouncing off the rev limiter and will not
			//accelerate any further
			if ((currentRpm > shiftRpm) && (velocityFactorCurrent === velocityFactorSix)) {

			    currentRpm = shiftRpm;
			}

			//Shifting Tracker
			if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorTwo)) {

				velocityFactorCurrent = velocityFactorTwo;
				currentRpm = velocityFactorCurrent * velocity;
				timeStraightSectionAcceleration += shiftTime;
				accelerationShifts += 1;
				ntf = nt_2 * np * nf;

			}

			if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorThree)) {

				velocityFactorCurrent = velocityFactorThree;
				currentRpm = velocityFactorCurrent * velocity;
				timeStraightSectionAcceleration += shiftTime;
				accelerationShifts += 1;
				ntf = nt_3 * np * nf;

			}

			if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorFour)) {

				velocityFactorCurrent = velocityFactorFour;
				currentRpm = velocityFactorCurrent * velocity;
				timeStraightSectionAcceleration += shiftTime;
				accelerationShifts += 1;
				ntf = nt_4 * np * nf;

			}

			if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorFive)) {

				velocityFactorCurrent = velocityFactorFive;
				currentRpm = velocityFactorCurrent * velocity;
				timeStraightSectionAcceleration += shiftTime;
				accelerationShifts += 1;
				ntf = nt_5 * np * nf;

            }
            if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorSix)) {

                velocityFactorCurrent = velocityFactorSix;
                currentRpm = velocityFactorCurrent * velocity;
                timeStraightSectionAcceleration += shiftTime;
                accelerationShifts += 1;
                ntf = nt_6 * np * nf;

            }

  
			timeStraightSectionAcceleration += t;
			accelerationWideOpenThrottleTime += t;
			instantaneousTorque = torqueLocator(currentRpm, rpm, tor);

		}
        accelerationTrapSpeed = velocity / 1.4666667;
        accelerationTrapSpeed = roundNumber(accelerationTrapSpeed, 2);

		velocity = initialVelocity;
		
		averageAccelerationLongitudinalAccel = meanValue(accelerationLongitudinalAccelerationStore);
		averageAccelerationLongitudinalAccel = roundNumber(averageAccelerationLongitudinalAccel, 3);


		timeStraightSectionAcceleration = timeStraightSectionAcceleration + launchingTime;
		timeStraightSectionAcceleration = roundNumber(timeStraightSectionAcceleration, 3);

		accelerationTimeTractionLimited = roundNumber(accelerationTimeTractionLimited, 3);

		accelerationPoints = 75;

		accelerationWideOpenThrottlePercentage = (accelerationWideOpenThrottleTime / timeStraightSectionAcceleration) * 100;
		accelerationWideOpenThrottlePercentage = roundNumber(accelerationWideOpenThrottlePercentage, 2);

		accelerationShifts = accelerationShifts + launchingShifts;
		





	/*
	Cornering Simulation
	IN PROGRESS

	The Cornering Acceleration is first assumed and later iterated.  The higher acceleration, the higher the weight transfer
	to the oustide wheels.

	Tire 1 is the left front tire
	Tire 2 is the right front tire
	Tire 3 is the left rear tire
	Tire 4 is the right rear tire

	*/
	corneringAcceleration = 1 * g;
	corneringCounter = 0;
	difference = 100;
	tireOne = 100;
	tireTwo = 100;
	var corneringDownforce = downforceAtVelocity / 4;
	while (difference > 0.1 || corneringCounter < 10) {

		corneringForce = 0;
		oldCorneringAcceleration = corneringAcceleration;

		corneringWeightDistribution = vehicleWeight * ((trackA / trackWidth) + (corneringAcceleration / g) * (hcg / trackWidth));

		if (corneringWeightDistribution >= vehicleWeight) {
			corneringWeightDistribution = vehicleWeight;
		}

        tireOne = (vehicleWeight - corneringWeightDistribution) / 2 + corneringDownforce;
        tireThree = tireOne + corneringDownforce;
        tireTwo = (corneringWeightDistribution) / 2 + corneringDownforce;
        tireFour = tireTwo + corneringDownforce;

		weightOnTire = tireOne;
		instantaneousCoefficientOfFriction = frictionCalculation(weightOnTire, tireChoice);
		corneringForce = instantaneousCoefficientOfFriction * tireOne + corneringForce;

		weightOnTire = tireTwo;
		instantaneousCoefficientOfFriction = frictionCalculation(weightOnTire, tireChoice);
		corneringForce = instantaneousCoefficientOfFriction * tireTwo + corneringForce;

		weightOnTire = tireThree;
		instantaneousCoefficientOfFriction = frictionCalculation(weightOnTire, tireChoice);
		corneringForce = instantaneousCoefficientOfFriction * tireThree + corneringForce;

		weightOnTire = tireFour;
		instantaneousCoefficientOfFriction = frictionCalculation(weightOnTire, tireChoice);
		corneringForce = instantaneousCoefficientOfFriction * tireFour + corneringForce;

		corneringAcceleration = (corneringForce / vehicleWeight) * g;
		corneringComparison = oldCorneringAcceleration - corneringAcceleration;

		difference = (Math.abs(corneringComparison) / ((oldCorneringAcceleration + corneringAcceleration) * 0.5)) * 100;
		corneringCounter += 1;



	}

	corneringTime = 0;
	for (i = 0; i < cornerTotal.length; i++) {

		corneringVelocity[i] = Math.sqrt(corneringAcceleration * cornerTotal[i]);
		corneringTime += corneringDist[i] / corneringVelocity[i];



	}


//Skidpad Simulation

//The inner circles are 15.25 (50.03 feet) in diameter,
//and the outer circles are (69.72 feet) in diameter.
//The diameter of the driving line is then 59.875 ft in diameter, with the radius of the corner at 29.9375.
//The full circumference of the circle is the arc length, 2*pi*radius = 188.103 feet.

    skidpadTime = 0;
    skidpadLength = 188.103;
    skidpadVelocity = 0;
    skidpadLateralAccel = 0;
    skidpadRadius = 29.9375;
    skidpadPoints = 50;


    skidpadVelocity = Math.sqrt(corneringAcceleration * skidpadRadius);
    skidpadTime = skidpadLength / skidpadVelocity;


    skidpadVelocity = roundNumber(skidpadVelocity, 2);
    skidpadTime = roundNumber(skidpadTime, 3);

    skidpadLateralAccel = corneringAcceleration / g;
    skidpadLateralAccel = roundNumber(skidpadLateralAccel, 2);


	brakingAcceleration = (-1.5) * g;


	for (var o = 0; o < straightSections; o++) {

		if (o === 0) {
			corneringBrakeTo = corneringVelocity[0];
		}
		else {
			corneringBrakeTo = corneringVelocity[1 + numberOfCorners[o - 1]];
		}

		velocityFactorCurrent = velocityFactorOne;
		currentRpm = velocityFactorOne * velocity;
		ntf = nt_1 * np * nf;

		//Shift to 2nd Gear.  Doesnt count as a shift because it is set-up after a corner
		if (currentRpm > shiftRpm) {
			velocityFactorCurrent = velocityFactorTwo;
			currentRpm = velocityFactorTwo * velocity;
			ntf = nt_2 * np * nf;
		}
		//Shift to 3rd Gear.  Doesnt count as a shift because it is set-up after a corner
		if (currentRpm > shiftRpm) {
			velocityFactorCurrent = velocityFactorThree;
			currentRpm = velocityFactorThree * velocity;
			ntf = nt_3 * np * nf;
		}
		//Shift to 4th Gear.  Doesnt count as a shift because it is set-up after a corner
		if (currentRpm > shiftRpm) {
			velocityFactorCurrent = velocityFactorFour;
			currentRpm = velocityFactorFour * velocity;
			ntf = nt_4 * np * nf;
		}
		//Shift to 5th Gear.  Doesnt count as a shift because it is set-up after a corner
		if (currentRpm > shiftRpm) {
			velocityFactorCurrent = velocityFactorFive;
			currentRpm = velocityFactorFive * velocity;
			ntf = nt_5 * np * nf;
        }        //Shift to 6th Gear.  Doesnt count as a shift because it is set-up after a corner        if (currentRpm > shiftRpm) {            velocityFactorCurrent = velocityFactorSix;            currentRpm = velocityFactorSix * velocity;            ntf = nt_6 * np * nf;        }
        if (currentRpm > shiftRpm) {
            currentRpm = shiftRpm;
        }


		instantaneousTorque = torqueLocator(currentRpm, rpm, tor);
		x = 0;
		t = iterationTime;
		count = -1;

		while (x < dist[o]) {

			count += 1;

			instantaneousDownforce = downforceCalculation(downforceConstant, velocity);
			instantaneousRollingResistanceForce = rollingResistanceForceCalculation(f_0, f_s, velocity, vehicleWeight, instantaneousDownforce);
			instantaneousDragForce = dragForceCalculation(rho, velocity, coefficientOfDrag, frontalArea, instantaneousDownforce, LiftDragCoefficient);

			//Acceleration Calculation
			instantaneousLongitudinalAcceleration = longitudinalAccelerationCalculation(instantaneousTorque, drivetrainEfficiency, wheelRadius, instantaneousDragForce, instantaneousRollingResistanceForce, vehicleWeight, ntf);


			//Weight Distribution and Instantaneous Friction Between Tires and Ground
			weightOnTire = (vehicleWeight * bLength) / (bLength + cLength - hcg * instantaneousLongitudinalAcceleration / g);
			weightOnTire = weightOnTire + instantaneousDownforce / 2;
			instantaneousCoefficientOfFriction = frictionCalculation(weightOnTire, tireChoice);

			//Maximum Tractive Force Calculation

			maximumTractiveForce = instantaneousCoefficientOfFriction * weightOnTire;
			maximumLongitudinalAcceleration = maximumTractiveForce / (vehicleWeight / g);

			if (instantaneousLongitudinalAcceleration > maximumLongitudinalAcceleration) {

				instantaneousLongitudinalAcceleration = maximumLongitudinalAcceleration;
				
			}

			x = x + velocity * t + 0.5 * instantaneousLongitudinalAcceleration * Math.pow(t, 2);
			velocity = velocity + instantaneousLongitudinalAcceleration * t;

			currentRpm = velocityFactorCurrent * velocity;

			//Max Vehicle Velocity Reached as the current RPM is at Redline and the vehicle is in the last gear
			if ((currentRpm > shiftRpm) && (velocityFactorCurrent === velocityFactorSix)) {

			    velocity -= instantaneousLongitudinalAcceleration * t;
			}

			//Stores the values for velocity, distance, and acceleration to be used in future
			//calculations and analysis.
			velocityStore[count] = velocity;
			xStore[count] = x;
			instantaneousLongitudinalAccelerationStore[count] = instantaneousLongitudinalAcceleration;


			//This part of the code is where the engine is bouncing off the rev limiter and will not
			//accelerate any further
			if ((currentRpm > shiftRpm) && (velocityFactorCurrent === velocityFactorSix)) {

			    currentRpm = shiftRpm;
			}


			//Shifting Tracker
			if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorTwo)) {

				velocityFactorCurrent = velocityFactorTwo;
				currentRpm = velocityFactorCurrent * velocity;
				timeStraightSection += shiftTime;
				shifts += 1;
				ntf = nt_2 * np * nf;

			}

			if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorThree)) {

				velocityFactorCurrent = velocityFactorThree;
				currentRpm = velocityFactorCurrent * velocity;
				timeStraightSection += shiftTime;
				shifts += 1;
				ntf = nt_3 * np * nf;

			}

			if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorFour)) {

				velocityFactorCurrent = velocityFactorFour;
				currentRpm = velocityFactorCurrent * velocity;
				timeStraightSection += shiftTime;
				shifts += 1;
				ntf = nt_4 * np * nf;

			}

			if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorFive)) {

				velocityFactorCurrent = velocityFactorFive;
				currentRpm = velocityFactorCurrent * velocity;
				timeStraightSection += shiftTime;
				shifts += 1;
				ntf = nt_5 * np * nf;

            }

            if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorSix)) {

                velocityFactorCurrent = velocityFactorSix;
                currentRpm = velocityFactorCurrent * velocity;
                timeStraightSection += shiftTime;
                shifts += 1;
                ntf = nt_6 * np * nf;

            }

			instantaneousTorque = torqueLocator(currentRpm, rpm, tor);

		}

		accelerationCurve = polyfit(xStore, velocityStore, 8);

		decelAccelPoint = 0;
		accelerating = 0;
		decelerating = 0;
		totalAccelDecel = 1;

		//If our vehicle did not have enough time to achieve the maximum speed necessary to take the corner
		//then it will corner at the current velocity and not the maximum that is possible for the corner.


		if (corneringBrakeTo < velocity) {

			while (totalAccelDecel > 0.1) {

				decelAccelPoint += 0.01;
				accelerating = polyval(accelerationCurve, decelAccelPoint);
				brakingToBeSqrt = Math.pow(corneringBrakeTo, 2) - 2 * brakingAcceleration * (dist[o] - decelAccelPoint);
				decelerating = Math.sqrt(brakingToBeSqrt);
				totalAccelDecel = Math.abs(accelerating - decelerating);

			}
		}


		if (corneringBrakeTo > velocity) {
			corneringBrakeTo = velocity;
			decelAccelPoint = dist[o]
			//The variable brakingTracker will keep track of how many times the vehicle did not achieve maximum cornering speed
			brakingTracker += 1;
		}

		//End of the calculation for braking point of the First Straight

		//Recalculation of the Acceleration and Deceleration

		velocity = initialVelocity;

		x = 0;
		t = iterationTime;
		timeStraightSection = 0;
		timeOpenThrottle = 0;
		timeTractionLimited = 0;
		count = -1;
		shifts = 0;

		var velocityStore = [];
		var xStore = [];
		var instantaneousLongitudinalAccelerationStore = [];

		velocityFactorCurrent = velocityFactorOne;
		currentRpm = velocityFactorOne * velocity;
		ntf = nt_1 * np * nf;

		//Shift to 2nd Gear.  Doesnt count as a shift because it is set-up after a corner
		if (currentRpm > shiftRpm) {
			velocityFactorCurrent = velocityFactorTwo;
			currentRpm = velocityFactorTwo * velocity;
			ntf = nt_2 * np * nf;
		}
		//Shift to 3rd Gear.  Doesnt count as a shift because it is set-up after a corner
		if (currentRpm > shiftRpm) {
			velocityFactorCurrent = velocityFactorThree;
			currentRpm = velocityFactorThree * velocity;
			ntf = nt_3 * np * nf;
		}
		//Shift to 4th Gear.  Doesnt count as a shift because it is set-up after a corner
		if (currentRpm > shiftRpm) {
			velocityFactorCurrent = velocityFactorFour;
			currentRpm = velocityFactorFour * velocity;
			ntf = nt_4 * np * nf;
		}
		//Shift to 5th Gear.  Doesnt count as a shift because it is set-up after a corner
		if (currentRpm > shiftRpm) {
			velocityFactorCurrent = velocityFactorFive;
			currentRpm = velocityFactorFive * velocity;
			ntf = nt_5 * np * nf;
        }        //Shift to 6th Gear.  Doesnt count as a shift because it is set-up after a corner        if (currentRpm > shiftRpm) {            velocityFactorCurrent = velocityFactorSix;            currentRpm = velocityFactorSix * velocity;            ntf = nt_6 * np * nf;        }

		instantaneousTorque = torqueLocator(currentRpm, rpm, tor);

		while (x < decelAccelPoint) {
			count += 1;


			instantaneousDownforce = downforceCalculation(downforceConstant, velocity);
			instantaneousRollingResistanceForce = rollingResistanceForceCalculation(f_0, f_s, velocity, vehicleWeight, instantaneousDownforce);
			instantaneousDragForce = dragForceCalculation(rho, velocity, coefficientOfDrag, frontalArea, instantaneousDownforce, LiftDragCoefficient);

			//Acceleration Calculation
			instantaneousLongitudinalAcceleration = longitudinalAccelerationCalculation(instantaneousTorque, drivetrainEfficiency, wheelRadius, instantaneousDragForce, instantaneousRollingResistanceForce, vehicleWeight, ntf);


			//Weight Distribution and Instantaneous Friction Between Tires and Ground
			weightOnTire = (vehicleWeight * bLength) / (bLength + cLength - hcg * instantaneousLongitudinalAcceleration / g);
			weightOnTire = weightOnTire + instantaneousDownforce / 2;
			instantaneousCoefficientOfFriction = frictionCalculation(weightOnTire, tireChoice);

			//Maximum Tractive Force Calculation

			maximumTractiveForce = instantaneousCoefficientOfFriction * weightOnTire;
			maximumLongitudinalAcceleration = maximumTractiveForce / (vehicleWeight / g);

			if (instantaneousLongitudinalAcceleration > maximumLongitudinalAcceleration) {

				instantaneousLongitudinalAcceleration = maximumLongitudinalAcceleration;
				timeTractionLimited += t;
				timeOpenThrottle -= t;
			}

			horsepowerStore[count] = instantaneousTorque * currentRpm / 5252;
			x = x + velocity * t + 0.5 * instantaneousLongitudinalAcceleration * Math.pow(t, 2);
			velocity = velocity + instantaneousLongitudinalAcceleration * t;

			currentRpm = velocityFactorCurrent * velocity;

			//Max Vehicle Velocity Reached as the current RPM is at Redline and the vehicle is in the last gear
			if ((currentRpm > shiftRpm) && (velocityFactorCurrent === velocityFactorSix)) {

			    velocity -= instantaneousLongitudinalAcceleration * t;
			}

			instantaneousLongitudinalAccelerationStore[count] = instantaneousLongitudinalAcceleration;
			autocrossVelocityStore[count] = velocity;
			

			//This part of the code is where the engine is bouncing off the rev limiter and will not
			//accelerate any further
			if ((currentRpm > shiftRpm) && (velocityFactorCurrent === velocityFactorSix)) {

			    currentRpm = shiftRpm;
			}


			//Shifting Tracker
			if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorTwo)) {

				velocityFactorCurrent = velocityFactorTwo;
				currentRpm = velocityFactorCurrent * velocity;
				timeStraightSection += shiftTime;
				shifts += 1;
				ntf = nt_2 * np * nf;

			}

			if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorThree)) {

				velocityFactorCurrent = velocityFactorThree;
				currentRpm = velocityFactorCurrent * velocity;
				timeStraightSection += shiftTime;
				shifts += 1;
				ntf = nt_3 * np * nf;

			}

			if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorFour)) {

				velocityFactorCurrent = velocityFactorFour;
				currentRpm = velocityFactorCurrent * velocity;
				timeStraightSection += shiftTime;
				shifts += 1;
				ntf = nt_4 * np * nf;

			}

			if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorFive)) {

				velocityFactorCurrent = velocityFactorFive;
				currentRpm = velocityFactorCurrent * velocity;
				timeStraightSection += shiftTime;
				shifts += 1;
				ntf = nt_5 * np * nf;

            }
            if ((currentRpm > shiftRpm) && (velocityFactorCurrent > velocityFactorSix)) {

                velocityFactorCurrent = velocityFactorSix;
                currentRpm = velocityFactorCurrent * velocity;
                timeStraightSection += shiftTime;
                shifts += 1;
                ntf = nt_6 * np * nf;

            }

			timeStraightSection += t;
			timeOpenThrottle += t;
			instantaneousTorque = torqueLocator(currentRpm, rpm, tor);

		}


		averageLongitudinalAcceleration[o] = meanValue(instantaneousLongitudinalAccelerationStore);
		averageHorsepower[o] = meanValue(horsepowerStore);
		autocrossAverageVelocityTotal[o] = meanValue(autocrossVelocityStore);

		//Clears all the arrays
		var horsepowerStore = [];
		var instantaneousLongitudinalAccelerationStore = [];
		var autocrossVelocityStore = [];
		var velocityStore = [];
		var xStore = [];

		//Time to brake for the first corner
		decelerationVelocity = velocity;
		brakingDistance[o] = dist[o] - decelAccelPoint;
		timeStraightSection += 2 * brakingDistance[o] / (decelerationVelocity + corneringBrakeTo) + shiftTime;


		//Total time storage
		totalStraightSectionTime[o] = timeStraightSection;

		//Total Shifts storage
		totalNumberOfShifts[o] = shifts;

		//Total time at Wide Open Throttle
		totalTimeWideOpenThrottle[o] = timeOpenThrottle;

		//Total time traction limited Acceleration
		totalTimeTractionLimited[o] = timeTractionLimited;

		//Set-up for the next Straight Section Calculation

		velocity = initialVelocity;

		x = 0;
		t = iterationTime;
		timeStraightSection = 0;
		timeOpenThrottle = 0;
		timeTractionLimited = 0;
		brakingTime = 0;
		count = -1;
		shifts = 0;

		velocityFactorCurrent = velocityFactorOne;
		currentRpm = velocityFactorOne * velocity;
		ntf = nt_1 * np * nf;

		//Shift to 2nd Gear.  Doesnt count as a shift because it is set-up after a corner
		if (currentRpm > shiftRpm) {
			velocityFactorCurrent = velocityFactorTwo;
			currentRpm = velocityFactorTwo * velocity;
			ntf = nt_2 * np * nf;
		}
		//Shift to 3rd Gear.  Doesnt count as a shift because it is set-up after a corner
		if (currentRpm > shiftRpm) {
			velocityFactorCurrent = velocityFactorThree;
			currentRpm = velocityFactorThree * velocity;
			ntf = nt_3 * np * nf;
		}
		//Shift to 4th Gear.  Doesnt count as a shift because it is set-up after a corner
		if (currentRpm > shiftRpm) {
			velocityFactorCurrent = velocityFactorFour;
			currentRpm = velocityFactorFour * velocity;
			ntf = nt_4 * np * nf;
		}
		//Shift to 5th Gear.  Doesnt count as a shift because it is set-up after a corner
		if (currentRpm > shiftRpm) {
			velocityFactorCurrent = velocityFactorFive;
			currentRpm = velocityFactorFive * velocity;
			ntf = nt_5 * np * nf;
        }
        //Shift to 5th Gear.  Doesnt count as a shift because it is set-up after a corner
        if (currentRpm > shiftRpm) {            velocityFactorCurrent = velocityFactorSix;            currentRpm = velocityFactorSix * velocity;            ntf = nt_6 * np * nf;        }

		instantaneousTorque = torqueLocator(currentRpm, rpm, tor);



	}


	//Total Distance Traveled Calculation
	totalDistance = (sum(dist) + sum(corneringDist)) * numberOfLaps;

    
	//Total Time
	totalTime = (sum(totalStraightSectionTime) + corneringTime) * numberOfLaps + launchingTime;
	totalTime = roundNumber(totalTime, 3);

	autocrossTime = sum(totalStraightSectionTime) + corneringTime + launchingTime;
	autocrossTime = roundNumber(autocrossTime, 3);


	totalStraightSectionTimeOutput = (sum(totalStraightSectionTime)) * numberOfLaps;
	totalStraightSectionTimeOutput = Math.round(totalStraightSectionTimeOutput);

	outputTotalWideOpenThrottle = (sum(totalTimeWideOpenThrottle)) * numberOfLaps;
	outputTotalWideOpenThrottle = Math.round(outputTotalWideOpenThrottle);

	outputTimeTractionLimited = (sum(totalTimeTractionLimited)) * numberOfLaps;
	outputTimeTractionLimited = roundNumber(outputTimeTractionLimited, 3);

	autocrossTractionLimitedTime = sum(totalTimeTractionLimited);
	autocrossTractionLimitedTime = roundNumber(autocrossTractionLimitedTime, 3);

	//Fuel Consumption Calculation
	fuelConsumptionRate = bsfc * meanValue(averageHorsepower);
	totalFuel = fuelConsumptionRate * totalTime;

	//The weight density of gasoline is 6.073 lb/US gal
	fuelUsed = totalFuel / 6.073;
	fuelUsed = roundNumber(fuelUsed, 2);


	//Shifting Calculation
	totalNumberOfShiftsOutput = (sum(totalNumberOfShifts)) * numberOfLaps;
	totalNumberOfShiftsOutput = roundNumber(totalNumberOfShiftsOutput, 1);

	autocrossTotalShifts = sum(totalNumberOfShifts);

	var startingVelocity = initialVelocity;
	var startingVelocityMPH = startingVelocity / 1.4666667;
	startingVelocityMPH = roundNumber(startingVelocityMPH, 3);

	var totalNumberOfShifting = sum(totalNumberOfShifts) * numberOfLaps;
	totalNumberOfShifting = Math.round(totalNumberOfShifting);

	var meanLongAccel = meanValue(averageLongitudinalAcceleration);
	meanLongAccel = roundNumber(meanLongAccel, 2);

	meanAutocrossVelocity = (meanValue(autocrossAverageVelocityTotal) + meanValue(corneringVelocity)) / 2;
	meanAutocrossVelocity = meanAutocrossVelocity / 1.4666667;
	meanAutocrossVelocity = roundNumber(meanAutocrossVelocity, 2);


	var meanHorsepower = meanValue(averageHorsepower);
	meanHorsepower = roundNumber(meanHorsepower, 2);

	var corneringTimeTotal = corneringTime * numberOfLaps;
	corneringTimeTotal = roundNumber(corneringTimeTotal, 3);

	percentTimeWideOpenThrottle = (outputTotalWideOpenThrottle / totalTime) * 100;
	percentTimeWideOpenThrottle = roundNumber(percentTimeWideOpenThrottle, 2);
    
    var totalPoints = 300;
    autocrossPoints = 150;



    return { output1: totalTime, output2: totalPoints, output3: fuelUsed, output4: percentTimeWideOpenThrottle, output5: corneringTimeTotal, output6: totalStraightSectionTimeOutput, output7: totalNumberOfShiftsOutput, output8: outputTimeTractionLimited, output9: meanLongAccel, output10: meanHorsepower, output11: timeStraightSectionAcceleration, output12: averageAccelerationLongitudinalAccel, output13: accelerationTimeTractionLimited, output14: accelerationPoints, output15: accelerationWideOpenThrottlePercentage, output16: accelerationShifts, output17: accelerationTrapSpeed, output18: skidpadTime, output19: skidpadPoints, output20: skidpadVelocity, output21: skidpadLateralAccel, output22: autocrossTime, output23: autocrossPoints, output24: autocrossTotalShifts, output25: autocrossTractionLimitedTime, output26: meanAutocrossVelocity, output27: maxHorsepower };

}