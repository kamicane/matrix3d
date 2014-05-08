/*
Vector4
This is mostly from three.js / euclideanspace.com / transform_util source code.
*/ 'use strict';

var prime = require('prime');

var degToRad = function(degrees) {
  return degrees * Math.PI / 180;
};

var radToDeg = function(radians) {
  return radians * 180 / Math.PI;
};

var Vector4 = prime({

  constructor: function Vector4(x, y, z, w) {
    this[0] = x || 0;
    this[1] = y || 0;
    this[2] = z || 0;
    this[3] = w || 0;
  },

  clone: function() {
    return new Vector4(this[0], this[1], this[2], this[3]);
  },

  get x() {
    return this[0];
  },

  get y() {
    return this[1];
  },

  get z() {
    return this[2];
  },

  get w() {
    return this[2];
  },

  equals: function(v3) {
    return this[0] === v3[0] && this[1] === v3[1] && this[2] === v3[2] && this[3] === v3[3];
  },

  // three.js
  length: function() {
    return Math.sqrt(this[0] * this[0] + this[1] * this[1] + this[2] * this[2] + this[3] * this[3]);
  },

  // three.js
  dot: function(v4) {
    return this[0] * v4[0] + this[1] * v4[1] + this[2] * v4[2] + this[3] * v4[3];
  },

  // three.js
  normalize: function() {
    var length = this.length();
    if (length === 0) return new Vector4(0, 0, 0, 1);

    var inv = 1 / length;
    return new Vector4(this[0] * inv, this[1] * inv, this[2] * inv, this[3] * inv);
  },

  // three.js
  // http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/index.htm
  // deg indicates you want the angle in degrees in the resulting Vector4
  quaternionToAngle: function(deg) {
    var v4 = this;
    // if w>1 acos and sqrt will produce errors, this cant happen if quaternion is normalised
    if (v4[3] > 1) v4 = v4.normalize();
    var w = 2 * Math.acos(v4[3]);
    var s = Math.sqrt(1 - v4[3] * v4[3]); // assuming quaternion normalised then w is less than 1, so term always positive.

    if (s < 1e-4) { // test to avoid divide by zero, s is always positive due to sqrt
      // if s close to zero then direction of axis not important
      return new Vector4(v4[0], v4[1], v4[2], w);
    } else {
      // normalise axis
      return new Vector4(v4[0] / s, v4[1] / s, v4[2] / s, deg ? radToDeg(w) : w);
    }

  },

  // three.js
  // http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm
  // deg indicates the Vector4 contains the angle in degrees
  angleToQuaternion: function(deg) {
    var angle = deg ? degToRad(this[3]) : this[3];
    var half = angle / 2, s = Math.sin(half);
    return new Vector4(this[0] * s, this[1] * s, this[2] * s, Math.cos(half));
  },

  // transform_util
  combine: function(v4, scale1, scale2) {
    var result = new Vector4;
    for (var i = 0; i < 4; i++) result[i] = this[i] * scale1 + v4[i] * scale2;
    return result;
  },

  lerp: function (v4, delta) {
    var scale1 = delta;
    var scale2 = 1 - delta;
    return v4.combine(this, scale1, scale2);
  },

  // transform_util
  slerp: function(v4q, delta) {
    var interpolated = new Vector4;

    var product = this.dot(v4q);

    // Clamp product to -1.0 <= product <= 1.0.
    product = Math.min(Math.max(product, -1), 1);

    // Interpolate angles along the shortest path. For example, to interpolate
    // between a 175 degree angle and a 185 degree angle, interpolate along the
    // 10 degree path from 175 to 185, rather than along the 350 degree path in
    // the opposite direction. This matches WebKit's implementation but not
    // the current W3C spec. Fixing the spec to match this approach is discussed
    // at:
    // http://lists.w3.org/Archives/Public/www-style/2013May/0131.html
    var scale1 = 1;
    if (product < 0) {
      product = -product;
      scale1 = -1.0;
    }

    var epsilon = 1e-5;
    if (Math.abs(product - 1.0) < epsilon) {
      for (var i = 0; i < 4; ++i) interpolated[i] = this[i];
      return interpolated;
    }

    var denom = Math.sqrt(1 - product * product);
    var theta = Math.acos(product);
    var w = Math.sin(delta * theta) * (1 / denom);

    scale1 *= Math.cos(delta * theta) - product * w;
    var scale2 = w;

    return this.combine(v4q, scale1, scale2);
  }

});

module.exports = Vector4;
