/*
Matrix3d
A port of Webkit's source code of 4x4 Matrices operations.
Implements methods from both SkMatrix44 (types, concat, invert, ...) and transform_util (decomposition, composition, ...)

Differences between implementations:
 - All operations are non desctructive, applied on a new Matrix3d instance.
 - Convenience methods toString, clone, round, ... mostly specific to JavaScript.
 - Some speed optimizations have been removed for simplicity and brevity.
*/ 'use strict';

var prime = require('prime');

var Vector3 = require('./Vector3');
var Vector4 = require('./Vector4');

var stringify = function(array, places) {
  if (places == null || places > 20) places = 20;

  var strings = [];
  for (var i = 0; i < array.length; i++) strings[i] = array[i].toFixed(10).replace(/\.?0+$/, '');
  return strings;
};

var TypeMask = {
    Identity: 0,
    Translate: 0x01,  //!< set if the matrix has translation
    Scale: 0x02,  //!< set if the matrix has any scale != 1
    Affine: 0x04,  //!< set if the matrix skews or rotates
    Perspective: 0x08,   //!< set if the matrix is in perspective
    All: 0xF,
    Unknown: 0x80
};

var Matrix3d = prime({

  constructor: function Matrix3d() {

    // m.m11, m.m12, m.m13, m.m14,
    // m.m21, m.m22, m.m23, m.m24,
    // m.m31, m.m32, m.m33, m.m34,
    // m.m41, m.m42, m.m43, m.m44

    var values = arguments;

    if (values.length === 1) values = values[0]; // matrix as list

    if (!values.length) values = [
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1
    ];

    var i = 0, j, k = 0;

    if (values.length === 6) { // 2d matrix

      var a = values[0];
      var b = values[1];
      var c = values[2];
      var d = values[3];
      var e = values[4];
      var f = values[5];

      values = [
        a, b, 0, 0,
        c, d, 0, 0,
        0, 0, 1, 0,
        e, f, 0, 1
      ];

    }

    if (values.length !== 16) throw new Error('invalid matrix');

    // always 16

    for (i = 0; i < 4; i++) {
      var col = this[i] = [];
      for (j = 0; j < 4; j++) {
        col[j] = values[k++];
      }
    }

  },

  // get 2x3

  get a() { return this.m11; },
  get b() { return this.m12; },
  get c() { return this.m21; },
  get d() { return this.m22; },
  get e() { return this.m41; },
  get f() { return this.m42; },

  // set 2x3

  set a(value) { this.m11 = value; },
  set b(value) { this.m12 = value; },
  set c(value) { this.m21 = value; },
  set d(value) { this.m22 = value; },
  set e(value) { this.m41 = value; },
  set f(value) { this.m42 = value; },

  // get 4x4

  get m11() { return this[0][0]; },
  get m12() { return this[0][1]; },
  get m13() { return this[0][2]; },
  get m14() { return this[0][3]; },
  get m21() { return this[1][0]; },
  get m22() { return this[1][1]; },
  get m23() { return this[1][2]; },
  get m24() { return this[1][3]; },
  get m31() { return this[2][0]; },
  get m32() { return this[2][1]; },
  get m33() { return this[2][2]; },
  get m34() { return this[2][3]; },
  get m41() { return this[3][0]; },
  get m42() { return this[3][1]; },
  get m43() { return this[3][2]; },
  get m44() { return this[3][3]; },

  // set 4x4

  set m11(value) { this[0][0] = value; },
  set m12(value) { this[0][1] = value; },
  set m13(value) { this[0][2] = value; },
  set m14(value) { this[0][3] = value; },
  set m21(value) { this[1][0] = value; },
  set m22(value) { this[1][1] = value; },
  set m23(value) { this[1][2] = value; },
  set m24(value) { this[1][3] = value; },
  set m31(value) { this[2][0] = value; },
  set m32(value) { this[2][1] = value; },
  set m33(value) { this[2][2] = value; },
  set m34(value) { this[2][3] = value; },
  set m41(value) { this[3][0] = value; },
  set m42(value) { this[3][1] = value; },
  set m43(value) { this[3][2] = value; },
  set m44(value) { this[3][3] = value; },

  // get shortcuts

  get transX() { return this[3][0]; },
  get transY() { return this[3][1]; },
  get transZ() { return this[3][2]; },
  get scaleX() { return this[0][0]; },
  get scaleY() { return this[1][1]; },
  get scaleZ() { return this[2][2]; },
  get perspX() { return this[0][3]; },
  get perspY() { return this[1][3]; },
  get perspZ() { return this[2][3]; },

  // set shortcuts

  set transX(value) { this[3][0] = value; },
  set transY(value) { this[3][1] = value; },
  set transZ(value) { this[3][2] = value; },
  set scaleX(value) { this[0][0] = value; },
  set scaleY(value) { this[1][1] = value; },
  set scaleZ(value) { this[2][2] = value; },
  set perspX(value) { this[0][3] = value; },
  set perspY(value) { this[1][3] = value; },
  set perspZ(value) { this[2][3] = value; },

  // type getter

  get type() {
    var m = this;
    var mask = 0;

    if (0 !== m.perspX || 0 !== m.perspY || 0 !== m.perspZ || 1 !== m[3][3]) {
      return TypeMask.Translate | TypeMask.Scale | TypeMask.Affine | TypeMask.Perspective;
    }

    if (0 !== m.transX || 0 !== m.transY || 0 !== m.transZ) {
      mask |= TypeMask.Translate;
    }

    if (1 !== m.scaleX || 1 !== m.scaleY || 1 !== m.scaleZ) {
      mask |= TypeMask.Scale;
    }

    if (0 !== m[1][0] || 0 !== m[0][1] || 0 !== m[0][2] ||
        0 !== m[2][0] || 0 !== m[1][2] || 0 !== m[2][1]) {
          mask |= TypeMask.Affine;
    }

    return mask;
  },

  // W3C
  is2d: function() {
    var m = this;

    return m.m31 === 0 && m.m32 === 0 &&
           m.m13 === 0 && m.m14 === 0 &&
           m.m23 === 0 && m.m24 === 0 &&
           m.m33 === 1 && m.m34 === 0 &&
           m.m43 === 0 && m.m44 === 1;
  },

  equals: function(m2) {
    var m1 = this;

    return
      m1.m11 === m2.m11 && m1.m12 === m2.m12 && m1.m13 === m2.m13 && m1.m14 === m2.m14 &&
      m1.m21 === m2.m21 && m1.m22 === m2.m22 && m1.m23 === m2.m23 && m1.m24 === m2.m24 &&
      m1.m31 === m2.m31 && m1.m32 === m2.m32 && m1.m33 === m2.m33 && m1.m34 === m2.m34 &&
      m1.m41 === m2.m41 && m1.m42 === m2.m42 && m1.m43 === m2.m43 && m1.m44 === m2.m44;
  },

  clone: function() {
    var m = this;

    return new Matrix3d(
      m.m11, m.m12, m.m13, m.m14,
      m.m21, m.m22, m.m23, m.m24,
      m.m31, m.m32, m.m33, m.m34,
      m.m41, m.m42, m.m43, m.m44
    );
  },

  /**
   *  Return true if the matrix is identity.
   */
  isIdentity: function() {
    return this.type === TypeMask.Identity;
  },

  /**
   *  Return true if the matrix contains translate or is identity.
   */
  isTranslate: function() {
    return !(this.type & ~TypeMask.Translate);
  },

  /**
   *  Return true if the matrix only contains scale or translate or is identity.
   */
  isScaleTranslate: function() {
    return !(this.type & ~(TypeMask.Scale | TypeMask.Translate));
  },

  concat: function(m2) {
    if (this.isIdentity()) return m2.clone();
    if (m2.isIdentity()) return this.clone();

    var m = new Matrix3d;

    for (var j = 0; j < 4; j++) {
      for (var i = 0; i < 4; i++) {
        var value = 0;
        for (var k = 0; k < 4; k++) {
          value += this[k][i] * m2[j][k];
        }
        m[j][i] = value;
      }
    }

    return m;
  },

  translate: function(v3) {
    var translationMatrix = new Matrix3d;
    translationMatrix.m41 = v3[0];
    translationMatrix.m42 = v3[1];
    translationMatrix.m43 = v3[2];
    return this.concat(translationMatrix);
  },

  scale: function(v3) {
    var m = new Matrix3d;
    m.m11 = v3[0];
    m.m22 = v3[1];
    m.m33 = v3[2];
    return this.concat(m);
  },

  rotate: function(v4q) {
    var rotationMatrix = new Matrix3d;

    var x = v4q[0];
    var y = v4q[1];
    var z = v4q[2];
    var w = v4q[3];

    rotationMatrix.m11 = 1 - 2 * (y * y + z * z);
    rotationMatrix.m21 = 2 * (x * y - z * w);
    rotationMatrix.m31 = 2 * (x * z + y * w);
    rotationMatrix.m12 = 2 * (x * y + z * w);
    rotationMatrix.m22 = 1 - 2 * (x * x + z * z);
    rotationMatrix.m32 = 2 * (y * z - x * w);
    rotationMatrix.m13 = 2 * (x * z - y * w);
    rotationMatrix.m23 = 2 * (y * z + x * w);
    rotationMatrix.m33 = 1 - 2 * (x * x + y * y);

    return this.concat(rotationMatrix);
  },

  skew: function(v3) {
    var skewMatrix = new Matrix3d;

    skewMatrix[1][0] = v3[0];
    skewMatrix[2][0] = v3[1];
    skewMatrix[2][1] = v3[2];

    return this.concat(skewMatrix);
  },

  perspective: function(v4) {
    var perspectiveMatrix = new Matrix3d;

    perspectiveMatrix.m14 = v4[0];
    perspectiveMatrix.m24 = v4[1];
    perspectiveMatrix.m34 = v4[2];
    perspectiveMatrix.m44 = v4[3];

    return this.concat(perspectiveMatrix);
  },

  map: function(v4) {
    var result = new Vector4;

    for (var i = 0; i < 4; i++) {
      var value = 0;
      for (var j = 0; j < 4; j++) {
        value += this[j][i] * v4[j];
      }
      result[i] = value;
    }

    return result;
  },

  determinant: function() {
    if (this.isIdentity()) return 1;
    if (this.isScaleTranslate()) return this[0][0] * this[1][1] * this[2][2] * this[3][3];

    var a00 = this[0][0];
    var a01 = this[0][1];
    var a02 = this[0][2];
    var a03 = this[0][3];
    var a10 = this[1][0];
    var a11 = this[1][1];
    var a12 = this[1][2];
    var a13 = this[1][3];
    var a20 = this[2][0];
    var a21 = this[2][1];
    var a22 = this[2][2];
    var a23 = this[2][3];
    var a30 = this[3][0];
    var a31 = this[3][1];
    var a32 = this[3][2];
    var a33 = this[3][3];

    var b00 = a00 * a11 - a01 * a10;
    var b01 = a00 * a12 - a02 * a10;
    var b02 = a00 * a13 - a03 * a10;
    var b03 = a01 * a12 - a02 * a11;
    var b04 = a01 * a13 - a03 * a11;
    var b05 = a02 * a13 - a03 * a12;
    var b06 = a20 * a31 - a21 * a30;
    var b07 = a20 * a32 - a22 * a30;
    var b08 = a20 * a33 - a23 * a30;
    var b09 = a21 * a32 - a22 * a31;
    var b10 = a21 * a33 - a23 * a31;
    var b11 = a22 * a33 - a23 * a32;

    // Calculate the determinant
    return b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
  },

  normalize: function() {
    var m44 = this.m44;
    // Cannot normalize.
    if (m44 === 0) return false;

    var normalizedMatrix = new Matrix3d;

    var scale = 1 / m44;

    for (var i = 0; i < 4; i++)
      for (var j = 0; j < 4; j++)
        normalizedMatrix[j][i] = this[j][i] * scale;

    return normalizedMatrix;
  },

  decompose: function() {
    // We'll operate on a copy of the matrix.
    var matrix = this.normalize();

    // If we cannot normalize the matrix, then bail early as we cannot decompose.
    if (!matrix) return false;

    var perspectiveMatrix = matrix.clone();

    var i, j;

    for (i = 0; i < 3; i++) perspectiveMatrix[i][3] = 0;
    perspectiveMatrix[3][3] = 1;

    // If the perspective matrix is not invertible, we are also unable to
    // decompose, so we'll bail early. Constant taken from SkMatrix44::invert.
    if (Math.abs(perspectiveMatrix.determinant()) < 1e-8) return false;

    var perspective;

    if (matrix[0][3] !== 0 || matrix[1][3] !== 0 || matrix[2][3] !== 0) {
      // rhs is the right hand side of the equation.
      var rightHandSide = new Vector4(
        matrix[0][3],
        matrix[1][3],
        matrix[2][3],
        matrix[3][3]
      );

      // Solve the equation by inverting perspectiveMatrix and multiplying
      // rightHandSide by the inverse.
      var inversePerspectiveMatrix = perspectiveMatrix.invert();
      if (!inversePerspectiveMatrix) return false;

      var transposedInversePerspectiveMatrix = inversePerspectiveMatrix.transpose();

      perspective = transposedInversePerspectiveMatrix.map(rightHandSide);

    } else {
      // No perspective.
      perspective = new Vector4(0, 0, 0, 1);
    }

    var translate = new Vector3;
    for (i = 0; i < 3; i++) translate[i] = matrix[3][i];

    var row = [];

    for (i = 0; i < 3; i++) {
      var v3 = row[i] = new Vector3;
      for (j = 0; j < 3; ++j)
        v3[j] = matrix[i][j];
    }

    // Compute X scale factor and normalize first row.
    var scale = new Vector3;
    scale[0] = row[0].length();
    row[0] = row[0].normalize();

    // Compute XY shear factor and make 2nd row orthogonal to 1st.
    var skew = new Vector3;
    skew[0] = row[0].dot(row[1]);
    row[1] = row[1].combine(row[0], 1.0, -skew[0]);

    // Now, compute Y scale and normalize 2nd row.
    scale[1] = row[1].length();
    row[1] = row[1].normalize();

    skew[0] /= scale[1];

    // Compute XZ and YZ shears, orthogonalize 3rd row
    skew[1] = row[0].dot(row[2]);
    row[2] = row[2].combine(row[0], 1.0, -skew[1]);
    skew[2] = row[1].dot(row[2]);
    row[2] = row[2].combine(row[1], 1.0, -skew[2]);

    // Next, get Z scale and normalize 3rd row.
    scale[2] = row[2].length();
    row[2] = row[2].normalize();
    skew[1] /= scale[2];
    skew[2] /= scale[2];

    // At this point, the matrix (in rows) is orthonormal.
    // Check for a coordinate system flip.  If the determinant
    // is -1, then negate the matrix and the scaling factors.
    var pdum3 = row[1].cross(row[2]);
    if (row[0].dot(pdum3) < 0) {
      for (i = 0; i < 3; i++) {
        scale[i] *= -1;
        for (j = 0; j < 3; ++j)
          row[i][j] *= -1;
      }
    }

    var quaternion = new Vector4(
      0.5 * Math.sqrt(Math.max(1 + row[0][0] - row[1][1] - row[2][2], 0)),
      0.5 * Math.sqrt(Math.max(1 - row[0][0] + row[1][1] - row[2][2], 0)),
      0.5 * Math.sqrt(Math.max(1 - row[0][0] - row[1][1] + row[2][2], 0)),
      0.5 * Math.sqrt(Math.max(1 + row[0][0] + row[1][1] + row[2][2], 0))
    );

    if (row[2][1] > row[1][2]) quaternion[0] = -quaternion[0];
    if (row[0][2] > row[2][0]) quaternion[1] = -quaternion[1];
    if (row[1][0] > row[0][1]) quaternion[2] = -quaternion[2];

    return new DecomposedMatrix(perspective, translate, quaternion, skew, scale);
  },

  invert: function() {
    var a00 = this[0][0];
    var a01 = this[0][1];
    var a02 = this[0][2];
    var a03 = this[0][3];
    var a10 = this[1][0];
    var a11 = this[1][1];
    var a12 = this[1][2];
    var a13 = this[1][3];
    var a20 = this[2][0];
    var a21 = this[2][1];
    var a22 = this[2][2];
    var a23 = this[2][3];
    var a30 = this[3][0];
    var a31 = this[3][1];
    var a32 = this[3][2];
    var a33 = this[3][3];

    var b00 = a00 * a11 - a01 * a10;
    var b01 = a00 * a12 - a02 * a10;
    var b02 = a00 * a13 - a03 * a10;
    var b03 = a01 * a12 - a02 * a11;
    var b04 = a01 * a13 - a03 * a11;
    var b05 = a02 * a13 - a03 * a12;
    var b06 = a20 * a31 - a21 * a30;
    var b07 = a20 * a32 - a22 * a30;
    var b08 = a20 * a33 - a23 * a30;
    var b09 = a21 * a32 - a22 * a31;
    var b10 = a21 * a33 - a23 * a31;
    var b11 = a22 * a33 - a23 * a32;

    // Calculate the determinant
    var det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

    // If det is zero, we want to return false. However, we also want to return false
    // if 1/det overflows to infinity (i.e. det is denormalized). Both of these are
    // handled by checking that 1/det is finite.
    if (det === 0 || !isFinite(det)) return false;

    var invdet = 1.0 / det;

    b00 *= invdet;
    b01 *= invdet;
    b02 *= invdet;
    b03 *= invdet;
    b04 *= invdet;
    b05 *= invdet;
    b06 *= invdet;
    b07 *= invdet;
    b08 *= invdet;
    b09 *= invdet;
    b10 *= invdet;
    b11 *= invdet;

    return new Matrix3d(
      a11 * b11 - a12 * b10 + a13 * b09,
      a02 * b10 - a01 * b11 - a03 * b09,
      a31 * b05 - a32 * b04 + a33 * b03,
      a22 * b04 - a21 * b05 - a23 * b03,
      a12 * b08 - a10 * b11 - a13 * b07,
      a00 * b11 - a02 * b08 + a03 * b07,
      a32 * b02 - a30 * b05 - a33 * b01,
      a20 * b05 - a22 * b02 + a23 * b01,
      a10 * b10 - a11 * b08 + a13 * b06,
      a01 * b08 - a00 * b10 - a03 * b06,
      a30 * b04 - a31 * b02 + a33 * b00,
      a21 * b02 - a20 * b04 - a23 * b00,
      a11 * b07 - a10 * b09 - a12 * b06,
      a00 * b09 - a01 * b07 + a02 * b06,
      a31 * b01 - a30 * b03 - a32 * b00,
      a20 * b03 - a21 * b01 + a22 * b00
    );
  },

  // W3C
  transpose: function() {
    var m = this;

    return new Matrix3d(
      m.m11, m.m21, m.m31, m.m41,
      m.m12, m.m22, m.m32, m.m42,
      m.m13, m.m23, m.m33, m.m43,
      m.m14, m.m24, m.m34, m.m44
    );
  },

  interpolation: function(matrix) {
    return new MatrixInterpolation(this, matrix);
  },

  toArray: function() {
    return this.is2d() ? this.toArray2d() : this.toArray3d();
  },

  toArray3d: function() {
    var m = this;

    return [
      m.m11, m.m12, m.m13, m.m14,
      m.m21, m.m22, m.m23, m.m24,
      m.m31, m.m32, m.m33, m.m34,
      m.m41, m.m42, m.m43, m.m44
    ];
  },

  toArray2d: function() {
    var m = this;

    return  [
      m.a, m.b,
      m.c, m.d,
      m.e, m.f
    ];
  },

  toString: function(places) {
    return this.is2d() ? this.toString2d(places) : this.toString3d(places);
  },

  toString3d: function(places) {
    return 'matrix3d(' + stringify(this.toArray3d()).join(', ') + ')';
  },

  toString2d: function(places) {
    return  'matrix(' + stringify(this.toArray2d()).join(', ') + ')';
  }

});

var DecomposedMatrix = prime({

  constructor: function DecomposedMatrix(perspective, translate, quaternion, skew, scale) {
    this.perspective = perspective;
    this.translate = translate;
    this.quaternion = quaternion;
    this.skew = skew;
    this.scale = scale;
  },

  interpolate: function(to, delta) {
    var from = this;

    var perspective = from.perspective.lerp(to.perspective, delta);
    var translate = from.translate.lerp(to.translate, delta);
    var quaternion = from.quaternion.slerp(to.quaternion, delta);
    var skew = from.skew.lerp(to.skew, delta);
    var scale = from.scale.lerp(to.scale, delta);
    return new DecomposedMatrix(perspective, translate, quaternion, skew, scale);
  },

  compose: function() {
    return new Matrix3d()
      .perspective(this.perspective)
      .translate(this.translate)
      .rotate(this.quaternion)
      .skew(this.skew)
      .scale(this.scale);
  }

});

var MatrixInterpolation = prime({

  constructor: function MatrixInterpolation(from, to) {
    this.matrix1 = from;
    this.matrix2 = to;
    this.from = from.decompose();
    this.to = to.decompose();
  },

  step: function(delta) {
    if (delta === 0) return this.matrix1;
    if (delta === 1) return this.matrix2;
    return this.from.interpolate(this.to, delta).compose();
  }

});

Matrix3d.Decomposed = DecomposedMatrix;
Matrix3d.Interpolation = MatrixInterpolation;

module.exports = Matrix3d;
