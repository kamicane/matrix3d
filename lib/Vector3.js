/*
Vector3
This is mostly from three.js / transform_util source code.
*/'use strict';

var prime = require('prime');

var Vector3 = prime({

  constructor: function Vector3(x, y, z) {
    this[0] = x || 0;
    this[1] = y || 0;
    this[2] = z || 0;
  },

  clone: function() {
    return new Vector3(this[0], this[1], this[2]);
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

  // three.js
  length: function() {
    return Math.sqrt(this[0] * this[0] + this[1] * this[1] + this[2] * this[2]);
  },

  // three.js
  normalize: function() {
    var length = this.length();
    if (length === 0) return new Vector3();
    return new Vector3(this[0] / length, this[1] / length, this[2] / length);
  },

  // three.js
  dot: function(v3) {
    return this[0] * v3[0] + this[1] * v3[1] + this[2] * v3[2];
  },

  // three.js
  cross: function(v3) {
    var x = this[0], y = this[1], z = this[2];

    return new Vector3(
      y * v3[2] - z * v3[1],
      z * v3[0] - x * v3[2],
      x * v3[1] - y * v3[0]
    );
  },

  // trasform_util
  combine: function(v3, scale1, scale2) {
    var result = new Vector3;
    for (var i = 0; i < 3; i++) result[i] = this[i] * scale1 + v3[i] * scale2;
    return result;
  },

});

module.exports = Vector3;
