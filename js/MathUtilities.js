const TWO_PI = Math.PI * 2;
const PI_2 = Math.PI * 0.5;
const DEG2RAD = Math.PI / 180;
const RAD2DEG = 180 / Math.PI;

function degToRad(degrees)
{
	return degrees * DEG2RAD;
}

function radToDeg(radians)
{
	return radians * RAD2DEG;
}

function clamp(value, min, max)
{
	return Math.max(min, Math.min(max, value));
}

// Linear mapping from range <a1, a2> to range <b1, b2>
function mapLinear(x, a1, a2, b1, b2)
{
	return b1 + (x - a1) * (b2 - b1) / (a2 - a1);
}

// https://www.gamedev.net/tutorials/programming/general-and-gameplay-programming/inverse-lerp-a-super-useful-yet-often-overlooked-function-r5230/
function inverseLerp(x, y, value)
{
	return (value - x) / (y - x);
}

// https://en.wikipedia.org/wiki/Linear_interpolation
function lerp(x, y, t)
{
	return (1 - t) * x + t * y;
}

// http://en.wikipedia.org/wiki/Smoothstep
function smoothstep(x, min, max)
{
	if (x <= min) return 0;
	if (x >= max) return 1;

	x = (x - min) / (max - min);

	return x * x * (3 - 2 * x);
}

// Random integer from <low, high> interval
function randInt(low, high)
{
	return low + Math.floor(Math.random() * (high - low + 1));
}

// Random float from <low, high> interval
function randFloat(low, high)
{
	return low + Math.random() * (high - low);
}

// Random float from <-range/2, range/2> interval
function randFloatSpread(range)
{
	return range * (0.5 - Math.random());
}


// VECTOR 3 //////////////////////////////////////////////////////////////////////////////////////

class Vector3
{
	constructor(x = 0, y = 0, z = 0)
	{
		this.x = x;
		this.y = y;
		this.z = z;
	}

	set(x, y, z)
	{
		this.x = x;
		this.y = y;
		this.z = z;

		return this;
	}

	clone()
	{
		return new this.constructor(this.x, this.y, this.z);
	}

	copy(v)
	{
		this.x = v.x;
		this.y = v.y;
		this.z = v.z;

		return this;
	}

	add(v)
	{
		this.x += v.x;
		this.y += v.y;
		this.z += v.z;

		return this;
	}

	addScalar(s)
	{
		this.x += s;
		this.y += s;
		this.z += s;

		return this;
	}

	addVectors(a, b)
	{
		this.x = a.x + b.x;
		this.y = a.y + b.y;
		this.z = a.z + b.z;

		return this;
	}

	addScaledVector(v, s)
	{
		this.x += v.x * s;
		this.y += v.y * s;
		this.z += v.z * s;

		return this;
	}

	sub(v)
	{
		this.x -= v.x;
		this.y -= v.y;
		this.z -= v.z;

		return this;
	}

	subScalar(s)
	{
		this.x -= s;
		this.y -= s;
		this.z -= s;

		return this;
	}

	subVectors(a, b)
	{
		this.x = a.x - b.x;
		this.y = a.y - b.y;
		this.z = a.z - b.z;

		return this;
	}

	multiply(v)
	{
		this.x *= v.x;
		this.y *= v.y;
		this.z *= v.z;

		return this;
	}

	multiplyScalar(scalar)
	{
		this.x *= scalar;
		this.y *= scalar;
		this.z *= scalar;

		return this;
	}

	multiplyVectors(a, b)
	{
		this.x = a.x * b.x;
		this.y = a.y * b.y;
		this.z = a.z * b.z;

		return this;
	}

	applyAxisAngle(axis, angle)
	{
		return this.applyQuaternion(_quaternion.setFromAxisAngle(axis, angle));
	}

	applyMatrix4(m)
	{
		const x = this.x, y = this.y, z = this.z;
		const e = m.elements;

		const w = 1 / (e[3] * x + e[7] * y + e[11] * z + e[15]);

		this.x = (e[0] * x + e[4] * y + e[8] * z + e[12]) * w;
		this.y = (e[1] * x + e[5] * y + e[9] * z + e[13]) * w;
		this.z = (e[2] * x + e[6] * y + e[10] * z + e[14]) * w;

		return this;
	}

	applyQuaternion(q)
	{
		const x = this.x, y = this.y, z = this.z;
		const qx = q._x, qy = q._y, qz = q._z, qw = q._w;

		// calculate quat * vector

		const ix = qw * x + qy * z - qz * y;
		const iy = qw * y + qz * x - qx * z;
		const iz = qw * z + qx * y - qy * x;
		const iw = - qx * x - qy * y - qz * z;

		// calculate result * inverse quat

		this.x = ix * qw + iw * - qx + iy * - qz - iz * - qy;
		this.y = iy * qw + iw * - qy + iz * - qx - ix * - qz;
		this.z = iz * qw + iw * - qz + ix * - qy - iy * - qx;

		return this;
	}

	transformDirection(m)
	{
		// input: Matrix4 affine matrix
		// vector interpreted as a direction
		const x = this.x, y = this.y, z = this.z;
		const e = m.elements;

		this.x = e[0] * x + e[4] * y + e[8] * z;
		this.y = e[1] * x + e[5] * y + e[9] * z;
		this.z = e[2] * x + e[6] * y + e[10] * z;

		return this.normalize();
	}

	divide(v)
	{
		this.x /= v.x;
		this.y /= v.y;
		this.z /= v.z;

		return this;
	}

	divideScalar(scalar)
	{
		return this.multiplyScalar(1 / scalar);
	}

	min(v)
	{
		this.x = Math.min(this.x, v.x);
		this.y = Math.min(this.y, v.y);
		this.z = Math.min(this.z, v.z);

		return this;
	}

	max(v)
	{
		this.x = Math.max(this.x, v.x);
		this.y = Math.max(this.y, v.y);
		this.z = Math.max(this.z, v.z);

		return this;
	}

	clamp(min, max)
	{
		// assumes min < max, componentwise

		this.x = Math.max(min.x, Math.min(max.x, this.x));
		this.y = Math.max(min.y, Math.min(max.y, this.y));
		this.z = Math.max(min.z, Math.min(max.z, this.z));

		return this;
	}

	negate()
	{
		this.x = - this.x;
		this.y = - this.y;
		this.z = - this.z;

		return this;
	}

	dot(v)
	{
		return this.x * v.x + this.y * v.y + this.z * v.z;
	}

	lengthSq()
	{
		return this.x * this.x + this.y * this.y + this.z * this.z;
	}

	length()
	{
		return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
	}

	normalize()
	{
		return this.divideScalar(this.length() || 1);
	}

	setLength(length)
	{
		return this.normalize().multiplyScalar(length);
	}

	lerp(v, alpha)
	{
		this.x += (v.x - this.x) * alpha;
		this.y += (v.y - this.y) * alpha;
		this.z += (v.z - this.z) * alpha;

		return this;
	}

	lerpVectors(v1, v2, alpha)
	{
		this.x = v1.x + (v2.x - v1.x) * alpha;
		this.y = v1.y + (v2.y - v1.y) * alpha;
		this.z = v1.z + (v2.z - v1.z) * alpha;

		return this;
	}

	cross(v)
	{
		return this.crossVectors(this, v);
	}

	crossVectors(a, b)
	{
		const ax = a.x, ay = a.y, az = a.z;
		const bx = b.x, by = b.y, bz = b.z;

		this.x = ay * bz - az * by;
		this.y = az * bx - ax * bz;
		this.z = ax * by - ay * bx;

		return this;
	}

	projectOnVector(v)
	{
		const denominator = v.lengthSq();

		if (denominator === 0) return this.set(0, 0, 0);

		const scalar = v.dot(this) / denominator;

		return this.copy(v).multiplyScalar(scalar);
	}

	projectOnPlane(planeNormal)
	{
		_vector.copy(this).projectOnVector(planeNormal);

		return this.sub(_vector);
	}

	reflect(normal)
	{
		// reflect incident vector off plane orthogonal to normal
		// normal is assumed to have unit length

		return this.sub(_vector.copy(normal).multiplyScalar(2 * this.dot(normal)));
	}

	angleTo(v)
	{
		const denominator = Math.sqrt(this.lengthSq() * v.lengthSq());

		if (denominator === 0) return Math.PI / 2;

		const theta = this.dot(v) / denominator;

		// clamp, to handle numerical problems

		return Math.acos(MathUtils.clamp(theta, - 1, 1));
	}

	distanceTo(v)
	{
		return Math.sqrt(this.distanceToSquared(v));
	}

	distanceToSquared(v)
	{
		const dx = this.x - v.x, dy = this.y - v.y, dz = this.z - v.z;

		return dx * dx + dy * dy + dz * dz;
	}

	setFromSpherical(s)
	{
		return this.setFromSphericalCoords(s.radius, s.phi, s.theta);
	}

	setFromSphericalCoords(radius, phi, theta)
	{
		const sinPhiRadius = Math.sin(phi) * radius;

		this.x = sinPhiRadius * Math.sin(theta);
		this.y = Math.cos(phi) * radius;
		this.z = sinPhiRadius * Math.cos(theta);

		return this;
	}

	setFromMatrixPosition(m)
	{
		const e = m.elements;

		this.x = e[12];
		this.y = e[13];
		this.z = e[14];

		return this;
	}

	setFromMatrixScale(m)
	{
		const sx = this.setFromMatrixColumn(m, 0).length();
		const sy = this.setFromMatrixColumn(m, 1).length();
		const sz = this.setFromMatrixColumn(m, 2).length();

		this.x = sx;
		this.y = sy;
		this.z = sz;

		return this;
	}

	setFromMatrixColumn(m, index)
	{
		return this.fromArray(m.elements, index * 4);
	}

	setFromMatrix3Column(m, index)
	{
		return this.fromArray(m.elements, index * 3);
	}

	fromArray(array, offset = 0)
	{
		this.x = array[offset];
		this.y = array[offset + 1];
		this.z = array[offset + 2];

		return this;
	}

	equals(v)
	{
		return ((v.x === this.x) && (v.y === this.y) && (v.z === this.z));
	}

	random()
	{
		this.x = Math.random();
		this.y = Math.random();
		this.z = Math.random();

		return this;
	}

	randomDirection()
	{
		// Derived from https://mathworld.wolfram.com/SpherePointPicking.html

		const u = (Math.random() - 0.5) * 2;
		const t = Math.random() * Math.PI * 2;
		const f = Math.sqrt(1 - u ** 2);

		this.x = f * Math.cos(t);
		this.y = f * Math.sin(t);
		this.z = u;

		return this;
	}

	*[Symbol.iterator]()
	{
		yield this.x;
		yield this.y;
		yield this.z;
	}

}

Vector3.prototype.isVector3 = true;

const _vector = /*@__PURE__*/ new Vector3();


// QUATERNION /////////////////////////////////////////////////////////////////////////////////////

class Quaternion
{
	constructor(x = 0, y = 0, z = 0, w = 1)
	{
		this._x = x;
		this._y = y;
		this._z = z;
		this._w = w;
	}

	get x()
	{
		return this._x;
	}

	set x(value)
	{
		this._x = value;
		this._onChangeCallback();
	}

	get y()
	{
		return this._y;
	}

	set y(value)
	{
		this._y = value;
		this._onChangeCallback();
	}

	get z()
	{
		return this._z;
	}

	set z(value)
	{
		this._z = value;
		this._onChangeCallback();
	}

	get w()
	{
		return this._w;
	}

	set w(value)
	{
		this._w = value;
		this._onChangeCallback();
	}

	set(x, y, z, w)
	{
		this._x = x;
		this._y = y;
		this._z = z;
		this._w = w;

		this._onChangeCallback();

		return this;
	}

	clone()
	{
		return new this.constructor(this._x, this._y, this._z, this._w);
	}

	copy(quaternion)
	{
		this._x = quaternion._x;
		this._y = quaternion._y;
		this._z = quaternion._z;
		this._w = quaternion._w;

		this._onChangeCallback();

		return this;
	}

	setFromEuler(euler, update)
	{
		const x = euler._x, y = euler._y, z = euler._z;

		// http://www.mathworks.com/matlabcentral/fileexchange/
		// 	20696-function-to-convert-between-dcm-euler-angles-quaternions-and-euler-vectors/
		//	content/SpinCalc.m

		const cos = Math.cos;
		const sin = Math.sin;

		const c1 = cos(x / 2);
		const c2 = cos(y / 2);
		const c3 = cos(z / 2);

		const s1 = sin(x / 2);
		const s2 = sin(y / 2);
		const s3 = sin(z / 2);

		this._x = s1 * c2 * c3 + c1 * s2 * s3;
		this._y = c1 * s2 * c3 - s1 * c2 * s3;
		this._z = c1 * c2 * s3 + s1 * s2 * c3;
		this._w = c1 * c2 * c3 - s1 * s2 * s3;

		if (update !== false) this._onChangeCallback();

		return this;
	}

	setFromAxisAngle(axis, angle)
	{
		// http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm

		// assumes axis is normalized

		const halfAngle = angle / 2, s = Math.sin(halfAngle);

		this._x = axis.x * s;
		this._y = axis.y * s;
		this._z = axis.z * s;
		this._w = Math.cos(halfAngle);

		this._onChangeCallback();

		return this;
	}

	setFromRotationMatrix(m)
	{
		// http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm

		// assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

		const te = m.elements,

		m11 = te[0], m12 = te[4], m13 = te[8],
		m21 = te[1], m22 = te[5], m23 = te[9],
		m31 = te[2], m32 = te[6], m33 = te[10],

		trace = m11 + m22 + m33;

		if (trace > 0)
		{
			const s = 0.5 / Math.sqrt(trace + 1.0);

			this._w = 0.25 / s;
			this._x = (m32 - m23) * s;
			this._y = (m13 - m31) * s;
			this._z = (m21 - m12) * s;

		} else if (m11 > m22 && m11 > m33)
		{
			const s = 2.0 * Math.sqrt(1.0 + m11 - m22 - m33);

			this._w = (m32 - m23) / s;
			this._x = 0.25 * s;
			this._y = (m12 + m21) / s;
			this._z = (m13 + m31) / s;

		} else if (m22 > m33)
		{
			const s = 2.0 * Math.sqrt(1.0 + m22 - m11 - m33);

			this._w = (m13 - m31) / s;
			this._x = (m12 + m21) / s;
			this._y = 0.25 * s;
			this._z = (m23 + m32) / s;

		} else
		{
			const s = 2.0 * Math.sqrt(1.0 + m33 - m11 - m22);

			this._w = (m21 - m12) / s;
			this._x = (m13 + m31) / s;
			this._y = (m23 + m32) / s;
			this._z = 0.25 * s;
		}

		this._onChangeCallback();

		return this;
	}

	setFromUnitVectors(vFrom, vTo)
	{
		// assumes direction vectors vFrom and vTo are normalized

		let r = vFrom.dot(vTo) + 1;

		if (r < Number.EPSILON)
		{
			// vFrom and vTo point in opposite directions

			r = 0;

			if (Math.abs(vFrom.x) > Math.abs(vFrom.z))
			{
				this._x = - vFrom.y;
				this._y = vFrom.x;
				this._z = 0;
				this._w = r;
			} else
			{
				this._x = 0;
				this._y = - vFrom.z;
				this._z = vFrom.y;
				this._w = r;
			}
		} else
		{
			// crossVectors( vFrom, vTo ); // inlined to avoid cyclic dependency on Vector3

			this._x = vFrom.y * vTo.z - vFrom.z * vTo.y;
			this._y = vFrom.z * vTo.x - vFrom.x * vTo.z;
			this._z = vFrom.x * vTo.y - vFrom.y * vTo.x;
			this._w = r;
		}

		return this.normalize();
	}

	angleTo(q)
	{
		return 2 * Math.acos(Math.abs(MathUtils.clamp(this.dot(q), - 1, 1)));
	}


	identity()
	{
		return this.set(0, 0, 0, 1);
	}

	invert()
	{
		// quaternion is assumed to have unit length
		return this.conjugate();
	}

	conjugate()
	{
		this._x *= - 1;
		this._y *= - 1;
		this._z *= - 1;

		this._onChangeCallback();

		return this;
	}

	dot(v)
	{
		return this._x * v._x + this._y * v._y + this._z * v._z + this._w * v._w;
	}

	lengthSq()
	{
		return this._x * this._x + this._y * this._y + this._z * this._z + this._w * this._w;
	}

	length()
	{
		return Math.sqrt(this._x * this._x + this._y * this._y + this._z * this._z + this._w * this._w);
	}

	normalize()
	{
		let l = this.length();

		if (l === 0)
		{

			this._x = 0;
			this._y = 0;
			this._z = 0;
			this._w = 1;

		} else
		{

			l = 1 / l;

			this._x = this._x * l;
			this._y = this._y * l;
			this._z = this._z * l;
			this._w = this._w * l;

		}

		this._onChangeCallback();

		return this;
	}

	multiply(q)
	{
		return this.multiplyQuaternions(this, q);
	}

	premultiply(q)
	{
		return this.multiplyQuaternions(q, this);
	}

	multiplyQuaternions(a, b)
	{
		// from http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/code/index.htm

		const qax = a._x, qay = a._y, qaz = a._z, qaw = a._w;
		const qbx = b._x, qby = b._y, qbz = b._z, qbw = b._w;

		this._x = qax * qbw + qaw * qbx + qay * qbz - qaz * qby;
		this._y = qay * qbw + qaw * qby + qaz * qbx - qax * qbz;
		this._z = qaz * qbw + qaw * qbz + qax * qby - qay * qbx;
		this._w = qaw * qbw - qax * qbx - qay * qby - qaz * qbz;

		this._onChangeCallback();

		return this;
	}

	equals(quaternion)
	{
		return (quaternion._x === this._x) && (quaternion._y === this._y) && (quaternion._z === this._z) && (quaternion._w === this._w);
	}

	_onChange(callback)
	{
		this._onChangeCallback = callback;

		return this;
	}

	_onChangeCallback() { }

}

Quaternion.prototype.isQuaternion = true;

const _quaternion = /*@__PURE__*/ new Quaternion();



// MATRIX 4 ////////////////////////////////////////////////////////////////////////////////////////

class Matrix4
{
	constructor()
	{
		this.elements = [
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1
		];
	}

	set(n11, n12, n13, n14, n21, n22, n23, n24, n31, n32, n33, n34, n41, n42, n43, n44)
	{
		const te = this.elements;

		te[0] = n11; te[4] = n12; te[8] = n13; te[12] = n14;
		te[1] = n21; te[5] = n22; te[9] = n23; te[13] = n24;
		te[2] = n31; te[6] = n32; te[10] = n33; te[14] = n34;
		te[3] = n41; te[7] = n42; te[11] = n43; te[15] = n44;

		return this;
	}

	identity()
	{
		this.set(
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1
		);

		return this;
	}

	copy(m)
	{
		const te = this.elements;
		const me = m.elements;

		te[0] = me[0]; te[1] = me[1]; te[2] = me[2]; te[3] = me[3];
		te[4] = me[4]; te[5] = me[5]; te[6] = me[6]; te[7] = me[7];
		te[8] = me[8]; te[9] = me[9]; te[10] = me[10]; te[11] = me[11];
		te[12] = me[12]; te[13] = me[13]; te[14] = me[14]; te[15] = me[15];

		return this;
	}

	copyPosition(m)
	{
		const te = this.elements, me = m.elements;

		te[12] = me[12];
		te[13] = me[13];
		te[14] = me[14];

		return this;
	}


	extractBasis(xAxis, yAxis, zAxis)
	{
		xAxis.setFromMatrixColumn(this, 0);
		yAxis.setFromMatrixColumn(this, 1);
		zAxis.setFromMatrixColumn(this, 2);

		return this;
	}

	makeBasis(xAxis, yAxis, zAxis)
	{
		this.set(
			xAxis.x, yAxis.x, zAxis.x, 0,
			xAxis.y, yAxis.y, zAxis.y, 0,
			xAxis.z, yAxis.z, zAxis.z, 0,
			0, 0, 0, 1
		);

		return this;
	}

	extractRotation(m)
	{
		// this method does not support reflection matrices

		const te = this.elements;
		const me = m.elements;

		const scaleX = 1 / _v1.setFromMatrixColumn(m, 0).length();
		const scaleY = 1 / _v1.setFromMatrixColumn(m, 1).length();
		const scaleZ = 1 / _v1.setFromMatrixColumn(m, 2).length();

		te[0] = me[0] * scaleX;
		te[1] = me[1] * scaleX;
		te[2] = me[2] * scaleX;
		te[3] = 0;

		te[4] = me[4] * scaleY;
		te[5] = me[5] * scaleY;
		te[6] = me[6] * scaleY;
		te[7] = 0;

		te[8] = me[8] * scaleZ;
		te[9] = me[9] * scaleZ;
		te[10] = me[10] * scaleZ;
		te[11] = 0;

		te[12] = 0;
		te[13] = 0;
		te[14] = 0;
		te[15] = 1;

		return this;
	}

	makeRotationFromEuler(euler)
	{
		const te = this.elements;

		const x = euler._x, y = euler._y, z = euler._z;
		const a = Math.cos(x), b = Math.sin(x);
		const c = Math.cos(y), d = Math.sin(y);
		const e = Math.cos(z), f = Math.sin(z);
		
		const ae = a * e, af = a * f, be = b * e, bf = b * f;

		te[0] = c * e;
		te[4] = - c * f;
		te[8] = d;

		te[1] = af + be * d;
		te[5] = ae - bf * d;
		te[9] = - b * c;

		te[2] = bf - ae * d;
		te[6] = be + af * d;
		te[10] = a * c;
		
		// bottom row
		te[3] = 0;
		te[7] = 0;
		te[11] = 0;

		// last column
		te[12] = 0;
		te[13] = 0;
		te[14] = 0;
		te[15] = 1;

		return this;
	}


	makeRotationFromQuaternion(q)
	{
		return this.compose(_zero, q, _one);
	}

	lookAt(eye, target, up)
	{
		const te = this.elements;

		_z.subVectors(eye, target);

		if (_z.lengthSq() === 0)
		{
			// eye and target are in the same position
			_z.z = 1;
		}

		_z.normalize();
		_x.crossVectors(up, _z);

		if (_x.lengthSq() === 0)
		{
			// up and z are parallel
			if (Math.abs(up.z) === 1)
			{

				_z.x += 0.0001;

			} else
			{

				_z.z += 0.0001;

			}

			_z.normalize();
			_x.crossVectors(up, _z);
		}

		_x.normalize();
		_y.crossVectors(_z, _x);

		te[0] = _x.x; te[4] = _y.x; te[8] = _z.x;
		te[1] = _x.y; te[5] = _y.y; te[9] = _z.y;
		te[2] = _x.z; te[6] = _y.z; te[10] = _z.z;

		return this;
	}

	multiply(m)
	{
		return this.multiplyMatrices(this, m);
	}

	premultiply(m)
	{
		return this.multiplyMatrices(m, this);
	}

	multiplyMatrices(a, b)
	{
		const ae = a.elements;
		const be = b.elements;
		const te = this.elements;

		const a11 = ae[0], a12 = ae[4], a13 = ae[8], a14 = ae[12];
		const a21 = ae[1], a22 = ae[5], a23 = ae[9], a24 = ae[13];
		const a31 = ae[2], a32 = ae[6], a33 = ae[10], a34 = ae[14];
		const a41 = ae[3], a42 = ae[7], a43 = ae[11], a44 = ae[15];

		const b11 = be[0], b12 = be[4], b13 = be[8], b14 = be[12];
		const b21 = be[1], b22 = be[5], b23 = be[9], b24 = be[13];
		const b31 = be[2], b32 = be[6], b33 = be[10], b34 = be[14];
		const b41 = be[3], b42 = be[7], b43 = be[11], b44 = be[15];

		te[0] = a11 * b11 + a12 * b21 + a13 * b31 + a14 * b41;
		te[4] = a11 * b12 + a12 * b22 + a13 * b32 + a14 * b42;
		te[8] = a11 * b13 + a12 * b23 + a13 * b33 + a14 * b43;
		te[12] = a11 * b14 + a12 * b24 + a13 * b34 + a14 * b44;

		te[1] = a21 * b11 + a22 * b21 + a23 * b31 + a24 * b41;
		te[5] = a21 * b12 + a22 * b22 + a23 * b32 + a24 * b42;
		te[9] = a21 * b13 + a22 * b23 + a23 * b33 + a24 * b43;
		te[13] = a21 * b14 + a22 * b24 + a23 * b34 + a24 * b44;

		te[2] = a31 * b11 + a32 * b21 + a33 * b31 + a34 * b41;
		te[6] = a31 * b12 + a32 * b22 + a33 * b32 + a34 * b42;
		te[10] = a31 * b13 + a32 * b23 + a33 * b33 + a34 * b43;
		te[14] = a31 * b14 + a32 * b24 + a33 * b34 + a34 * b44;

		te[3] = a41 * b11 + a42 * b21 + a43 * b31 + a44 * b41;
		te[7] = a41 * b12 + a42 * b22 + a43 * b32 + a44 * b42;
		te[11] = a41 * b13 + a42 * b23 + a43 * b33 + a44 * b43;
		te[15] = a41 * b14 + a42 * b24 + a43 * b34 + a44 * b44;

		return this;
	}

	multiplyScalar(s)
	{
		const te = this.elements;

		te[0] *= s; te[4] *= s; te[8] *= s; te[12] *= s;
		te[1] *= s; te[5] *= s; te[9] *= s; te[13] *= s;
		te[2] *= s; te[6] *= s; te[10] *= s; te[14] *= s;
		te[3] *= s; te[7] *= s; te[11] *= s; te[15] *= s;

		return this;
	}

	determinant()
	{
		const te = this.elements;

		const n11 = te[0], n12 = te[4], n13 = te[8], n14 = te[12];
		const n21 = te[1], n22 = te[5], n23 = te[9], n24 = te[13];
		const n31 = te[2], n32 = te[6], n33 = te[10], n34 = te[14];
		const n41 = te[3], n42 = te[7], n43 = te[11], n44 = te[15];

		//TODO: make this more efficient
		//( based on http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm )

		return (
			n41 * (
				+ n14 * n23 * n32
				- n13 * n24 * n32
				- n14 * n22 * n33
				+ n12 * n24 * n33
				+ n13 * n22 * n34
				- n12 * n23 * n34
			) +
			n42 * (
				+ n11 * n23 * n34
				- n11 * n24 * n33
				+ n14 * n21 * n33
				- n13 * n21 * n34
				+ n13 * n24 * n31
				- n14 * n23 * n31
			) +
			n43 * (
				+ n11 * n24 * n32
				- n11 * n22 * n34
				- n14 * n21 * n32
				+ n12 * n21 * n34
				+ n14 * n22 * n31
				- n12 * n24 * n31
			) +
			n44 * (
				- n13 * n22 * n31
				- n11 * n23 * n32
				+ n11 * n22 * n33
				+ n13 * n21 * n32
				- n12 * n21 * n33
				+ n12 * n23 * n31
			)
		);
	}

	setPosition(x, y, z)
	{
		const te = this.elements;

		if (x.isVector3)
		{
			te[12] = x.x;
			te[13] = x.y;
			te[14] = x.z;
		} else
		{
			te[12] = x;
			te[13] = y;
			te[14] = z;
		}

		return this;
	}

	invert()
	{
		// based on http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
		const te = this.elements,

		n11 = te[0], n21 = te[1], n31 = te[2], n41 = te[3],
		n12 = te[4], n22 = te[5], n32 = te[6], n42 = te[7],
		n13 = te[8], n23 = te[9], n33 = te[10], n43 = te[11],
		n14 = te[12], n24 = te[13], n34 = te[14], n44 = te[15],

		t11 = n23 * n34 * n42 - n24 * n33 * n42 + n24 * n32 * n43 - n22 * n34 * n43 - n23 * n32 * n44 + n22 * n33 * n44,
		t12 = n14 * n33 * n42 - n13 * n34 * n42 - n14 * n32 * n43 + n12 * n34 * n43 + n13 * n32 * n44 - n12 * n33 * n44,
		t13 = n13 * n24 * n42 - n14 * n23 * n42 + n14 * n22 * n43 - n12 * n24 * n43 - n13 * n22 * n44 + n12 * n23 * n44,
		t14 = n14 * n23 * n32 - n13 * n24 * n32 - n14 * n22 * n33 + n12 * n24 * n33 + n13 * n22 * n34 - n12 * n23 * n34;

		const det = n11 * t11 + n21 * t12 + n31 * t13 + n41 * t14;

		if (det === 0) return this.set(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

		const detInv = 1 / det;

		te[0] = t11 * detInv;
		te[1] = (n24 * n33 * n41 - n23 * n34 * n41 - n24 * n31 * n43 + n21 * n34 * n43 + n23 * n31 * n44 - n21 * n33 * n44) * detInv;
		te[2] = (n22 * n34 * n41 - n24 * n32 * n41 + n24 * n31 * n42 - n21 * n34 * n42 - n22 * n31 * n44 + n21 * n32 * n44) * detInv;
		te[3] = (n23 * n32 * n41 - n22 * n33 * n41 - n23 * n31 * n42 + n21 * n33 * n42 + n22 * n31 * n43 - n21 * n32 * n43) * detInv;

		te[4] = t12 * detInv;
		te[5] = (n13 * n34 * n41 - n14 * n33 * n41 + n14 * n31 * n43 - n11 * n34 * n43 - n13 * n31 * n44 + n11 * n33 * n44) * detInv;
		te[6] = (n14 * n32 * n41 - n12 * n34 * n41 - n14 * n31 * n42 + n11 * n34 * n42 + n12 * n31 * n44 - n11 * n32 * n44) * detInv;
		te[7] = (n12 * n33 * n41 - n13 * n32 * n41 + n13 * n31 * n42 - n11 * n33 * n42 - n12 * n31 * n43 + n11 * n32 * n43) * detInv;

		te[8] = t13 * detInv;
		te[9] = (n14 * n23 * n41 - n13 * n24 * n41 - n14 * n21 * n43 + n11 * n24 * n43 + n13 * n21 * n44 - n11 * n23 * n44) * detInv;
		te[10] = (n12 * n24 * n41 - n14 * n22 * n41 + n14 * n21 * n42 - n11 * n24 * n42 - n12 * n21 * n44 + n11 * n22 * n44) * detInv;
		te[11] = (n13 * n22 * n41 - n12 * n23 * n41 - n13 * n21 * n42 + n11 * n23 * n42 + n12 * n21 * n43 - n11 * n22 * n43) * detInv;

		te[12] = t14 * detInv;
		te[13] = (n13 * n24 * n31 - n14 * n23 * n31 + n14 * n21 * n33 - n11 * n24 * n33 - n13 * n21 * n34 + n11 * n23 * n34) * detInv;
		te[14] = (n14 * n22 * n31 - n12 * n24 * n31 - n14 * n21 * n32 + n11 * n24 * n32 + n12 * n21 * n34 - n11 * n22 * n34) * detInv;
		te[15] = (n12 * n23 * n31 - n13 * n22 * n31 + n13 * n21 * n32 - n11 * n23 * n32 - n12 * n21 * n33 + n11 * n22 * n33) * detInv;

		return this;
	}

	scale(v)
	{
		const te = this.elements;
		const x = v.x, y = v.y, z = v.z;

		te[0] *= x; te[4] *= y; te[8] *= z;
		te[1] *= x; te[5] *= y; te[9] *= z;
		te[2] *= x; te[6] *= y; te[10] *= z;
		te[3] *= x; te[7] *= y; te[11] *= z;

		return this;
	}

	makeTranslation(x, y, z)
	{
		this.set(
			1, 0, 0, x,
			0, 1, 0, y,
			0, 0, 1, z,
			0, 0, 0, 1
		);

		return this;
	}

	makeRotationX(theta)
	{
		const c = Math.cos(theta), s = Math.sin(theta);

		this.set(
			1, 0, 0, 0,
			0, c, - s, 0,
			0, s, c, 0,
			0, 0, 0, 1
		);

		return this;
	}

	makeRotationY(theta)
	{
		const c = Math.cos(theta), s = Math.sin(theta);

		this.set(
			c, 0, s, 0,
			0, 1, 0, 0,
			- s, 0, c, 0,
			0, 0, 0, 1
		);

		return this;
	}

	makeRotationZ(theta)
	{
		const c = Math.cos(theta), s = Math.sin(theta);

		this.set(
			c, - s, 0, 0,
			s, c, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1
		);

		return this;
	}

	makeRotationAxis(axis, angle)
	{
		// Based on http://www.gamedev.net/reference/articles/article1199.asp

		const c = Math.cos(angle);
		const s = Math.sin(angle);
		const t = 1 - c;
		const x = axis.x, y = axis.y, z = axis.z;
		const tx = t * x, ty = t * y;

		this.set(
			tx * x + c,     tx * y - s * z, tx * z + s * y, 0,
			tx * y + s * z, ty * y + c,     ty * z - s * x, 0,
			tx * z - s * y, ty * z + s * x,  t * z * z + c, 0,
			0,              0,               0, 		1
		);

		return this;
	}

	makeScale(x, y, z)
	{
		this.set(
			x, 0, 0, 0,
			0, y, 0, 0,
			0, 0, z, 0,
			0, 0, 0, 1
		);

		return this;
	}

	makeShear(xy, xz, yx, yz, zx, zy)
	{
		this.set(
			1, yx, zx, 0,
			xy, 1, zy, 0,
			xz, yz, 1, 0,
			0, 0, 0, 1
		);

		return this;
	}

	compose(position, quaternion, scale)
	{
		const te = this.elements;

		const x = quaternion._x, y = quaternion._y, z = quaternion._z, w = quaternion._w;
		const x2 = x + x, y2 = y + y, z2 = z + z;
		const xx = x * x2, xy = x * y2, xz = x * z2;
		const yy = y * y2, yz = y * z2, zz = z * z2;
		const wx = w * x2, wy = w * y2, wz = w * z2;

		const sx = scale.x, sy = scale.y, sz = scale.z;

		te[0] = (1 - (yy + zz)) * sx;
		te[1] = (xy + wz) * sx;
		te[2] = (xz - wy) * sx;
		te[3] = 0;

		te[4] = (xy - wz) * sy;
		te[5] = (1 - (xx + zz)) * sy;
		te[6] = (yz + wx) * sy;
		te[7] = 0;

		te[8] = (xz + wy) * sz;
		te[9] = (yz - wx) * sz;
		te[10] = (1 - (xx + yy)) * sz;
		te[11] = 0;

		te[12] = position.x;
		te[13] = position.y;
		te[14] = position.z;
		te[15] = 1;

		return this;
	}

	decompose(position, quaternion, scale)
	{
		const te = this.elements;

		let sx = _v1.set(te[0], te[1], te[2]).length();
		const sy = _v1.set(te[4], te[5], te[6]).length();
		const sz = _v1.set(te[8], te[9], te[10]).length();

		// if determine is negative, we need to invert one scale
		const det = this.determinant();
		if (det < 0) sx = - sx;

		position.x = te[12];
		position.y = te[13];
		position.z = te[14];

		// scale the rotation part
		_m1.copy(this);

		const invSX = 1 / sx;
		const invSY = 1 / sy;
		const invSZ = 1 / sz;

		_m1.elements[0] *= invSX;
		_m1.elements[1] *= invSX;
		_m1.elements[2] *= invSX;

		_m1.elements[4] *= invSY;
		_m1.elements[5] *= invSY;
		_m1.elements[6] *= invSY;

		_m1.elements[8] *= invSZ;
		_m1.elements[9] *= invSZ;
		_m1.elements[10] *= invSZ;

		quaternion.setFromRotationMatrix(_m1);

		scale.x = sx;
		scale.y = sy;
		scale.z = sz;

		return this;
	}

	equals(matrix)
	{
		const te = this.elements;
		const me = matrix.elements;

		for (let i = 0; i < 16; i++)
		{
			if (te[i] !== me[i]) return false;
		}

		return true;
	}

}

Matrix4.prototype.isMatrix4 = true;

const _v1 = new Vector3();
const _m1 = new Matrix4();
const _zero = new Vector3(0, 0, 0);
const _one = new Vector3(1, 1, 1);
const _x = new Vector3();
const _y = new Vector3();
const _z = new Vector3();


// EULER //////////////////////////////////////////////////////////////////////////////////////////

class Euler
{
	constructor(x = 0, y = 0, z = 0, order = Euler.DefaultOrder)
	{
		this._x = x;
		this._y = y;
		this._z = z;
		this._order = order;
	}

	get x()
	{
		return this._x;
	}

	set x(value)
	{
		this._x = value;
		this._onChangeCallback();
	}

	get y()
	{
		return this._y;
	}

	set y(value)
	{
		this._y = value;
		this._onChangeCallback();
	}

	get z()
	{
		return this._z;
	}

	set z(value)
	{
		this._z = value;
		this._onChangeCallback();
	}

	set(x, y, z, order = this._order)
	{
		this._x = x;
		this._y = y;
		this._z = z;
		this._order = order;

		this._onChangeCallback();

		return this;
	}

	clone()
	{
		return new this.constructor(this._x, this._y, this._z, this._order);
	}

	copy(euler)
	{
		this._x = euler._x;
		this._y = euler._y;
		this._z = euler._z;
		this._order = euler._order;

		this._onChangeCallback();

		return this;
	}

	setFromRotationMatrix(m, order = this._order, update = true)
	{
		// assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

		const te = m.elements;
		const m11 = te[0], m12 = te[4], m13 = te[8];
		const m21 = te[1], m22 = te[5], m23 = te[9];
		const m31 = te[2], m32 = te[6], m33 = te[10];

		this._y = Math.asin(clamp(m13, - 1, 1));

		if (Math.abs(m13) < 0.9999999)
		{
			this._x = Math.atan2(- m23, m33);
			this._z = Math.atan2(- m12, m11);
		} else
		{
			this._x = Math.atan2(m32, m22);
			this._z = 0;
		}

		this._order = order;

		if (update === true) this._onChangeCallback();

		return this;
	}

	setFromQuaternion(q, order, update)
	{
		_matrix.makeRotationFromQuaternion(q);

		return this.setFromRotationMatrix(_matrix, order, update);
	}

	setFromVector3(v, order = this._order)
	{
		return this.set(v.x, v.y, v.z, order);
	}

	equals(euler)
	{
		return (euler._x === this._x) && (euler._y === this._y) && (euler._z === this._z) && (euler._order === this._order);
	}

	_onChange(callback)
	{
		this._onChangeCallback = callback;

		return this;
	}

	_onChangeCallback() { }
}

Euler.prototype.isEuler = true;

Euler.DefaultOrder = 'XYZ';
Euler.RotationOrders = ['XYZ', 'YZX', 'ZXY', 'XZY', 'YXZ', 'ZYX'];

const _matrix = /*@__PURE__*/ new Matrix4();



// EVENT DISPATCHER ///////////////////////////////////////////////////////////////////////////////

class EventDispatcher
{
	addEventListener(type, listener)
	{
		if (this._listeners === undefined) this._listeners = {};

		const listeners = this._listeners;

		if (listeners[type] === undefined)
		{
			listeners[type] = [];
		}

		if (listeners[type].indexOf(listener) === - 1)
		{
			listeners[type].push(listener);
		}
	}

	hasEventListener(type, listener)
	{
		if (this._listeners === undefined) return false;

		const listeners = this._listeners;

		return listeners[type] !== undefined && listeners[type].indexOf(listener) !== - 1;
	}

	removeEventListener(type, listener)
	{
		if (this._listeners === undefined) return;

		const listeners = this._listeners;
		const listenerArray = listeners[type];

		if (listenerArray !== undefined)
		{
			const index = listenerArray.indexOf(listener);

			if (index !== - 1)
			{
				listenerArray.splice(index, 1);
			}
		}
	}

	dispatchEvent(event)
	{
		if (this._listeners === undefined) return;

		const listeners = this._listeners;
		const listenerArray = listeners[event.type];

		if (listenerArray !== undefined)
		{
			event.target = this;

			// Make a copy, in case listeners are removed while iterating.
			const array = listenerArray.slice(0);

			for (let i = 0, l = array.length; i < l; i++)
			{
				array[i].call(this, event);
			}

			event.target = null;
		}
	}
}



// OBJECT 3D ///////////////////////////////////////////////////////////////////////////////////////

let _object3DId = 0;

class Object3D extends EventDispatcher
{
	constructor()
	{
		super();

		Object.defineProperty(this, 'id', { value: _object3DId++ });

		///this.uuid = MathUtils.generateUUID();

		this.name = '';
		this.type = 'Object3D';

		this.parent = null;
		this.children = [];

		this.up = Object3D.DefaultUp.clone();

		const position = new Vector3();
		const rotation = new Euler();
		const quaternion = new Quaternion();
		const scale = new Vector3(1, 1, 1);

		function onRotationChange()
		{
			quaternion.setFromEuler(rotation, false);
		}

		function onQuaternionChange()
		{
			rotation.setFromQuaternion(quaternion, undefined, false);
		}

		rotation._onChange(onRotationChange);
		quaternion._onChange(onQuaternionChange);

		Object.defineProperties(this, {
			position: {
				configurable: true,
				enumerable: true,
				value: position
			},
			rotation: {
				configurable: true,
				enumerable: true,
				value: rotation
			},
			quaternion: {
				configurable: true,
				enumerable: true,
				value: quaternion
			},
			scale: {
				configurable: true,
				enumerable: true,
				value: scale
			}
		});

		this.matrix = new Matrix4();
		this.matrixWorld = new Matrix4();

		this.matrixAutoUpdate = Object3D.DefaultMatrixAutoUpdate;
		this.matrixWorldNeedsUpdate = false;

		this.visible = true;
	}

	onBeforeRender( /* renderer, scene, camera, geometry, material, group */) { }

	onAfterRender( /* renderer, scene, camera, geometry, material, group */) { }

	applyMatrix4(matrix)
	{
		if (this.matrixAutoUpdate) this.updateMatrix();

		this.matrix.premultiply(matrix);

		this.matrix.decompose(this.position, this.quaternion, this.scale);
	}

	applyQuaternion(q)
	{
		this.quaternion.premultiply(q);

		return this;
	}

	setRotationFromAxisAngle(axis, angle)
	{
		// assumes axis is normalized

		this.quaternion.setFromAxisAngle(axis, angle);
	}

	setRotationFromEuler(euler)
	{
		this.quaternion.setFromEuler(euler, true);
	}

	setRotationFromMatrix(m)
	{
		// assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

		this.quaternion.setFromRotationMatrix(m);
	}

	setRotationFromQuaternion(q)
	{
		// assumes q is normalized

		this.quaternion.copy(q);
	}

	rotateOnAxis(axis, angle)
	{
		// rotate object on axis in object space
		// axis is assumed to be normalized

		_q1.setFromAxisAngle(axis, angle);

		this.quaternion.multiply(_q1);

		return this;
	}

	rotateOnWorldAxis(axis, angle)
	{
		// rotate object on axis in world space
		// axis is assumed to be normalized
		// method assumes no rotated parent

		_q1.setFromAxisAngle(axis, angle);

		this.quaternion.premultiply(_q1);

		return this;
	}

	rotateX(angle)
	{
		return this.rotateOnAxis(_xAxis, angle);
	}

	rotateY(angle)
	{
		return this.rotateOnAxis(_yAxis, angle);
	}

	rotateZ(angle)
	{
		return this.rotateOnAxis(_zAxis, angle);
	}

	translateOnAxis(axis, distance)
	{
		// translate object by distance along axis in object space
		// axis is assumed to be normalized

		_v1.copy(axis).applyQuaternion(this.quaternion);

		this.position.add(_v1.multiplyScalar(distance));

		return this;
	}

	translateX(distance)
	{
		return this.translateOnAxis(_xAxis, distance);
	}

	translateY(distance)
	{
		return this.translateOnAxis(_yAxis, distance);
	}

	translateZ(distance)
	{
		return this.translateOnAxis(_zAxis, distance);
	}

	localToWorld(vector)
	{
		return vector.applyMatrix4(this.matrixWorld);
	}

	worldToLocal(vector)
	{
		return vector.applyMatrix4(_m1.copy(this.matrixWorld).invert());
	}

	lookAt(x, y, z)
	{
		// This method does not support objects having non-uniformly-scaled parent(s)

		if (x.isVector3)
		{
			_target.copy(x);
		} else
		{
			_target.set(x, y, z);
		}

		const parent = this.parent;

		this.updateWorldMatrix(true, false);

		_position.setFromMatrixPosition(this.matrixWorld);

		if (this.isCamera || this.isLight)
		{
			_m1.lookAt(_position, _target, this.up);
		} else
		{
			_m1.lookAt(_target, _position, this.up);
		}

		this.quaternion.setFromRotationMatrix(_m1);

		if (parent)
		{
			_m1.extractRotation(parent.matrixWorld);
			_q1.setFromRotationMatrix(_m1);
			this.quaternion.premultiply(_q1.invert());
		}
	}

	add(object)
	{
		if (arguments.length > 1)
		{
			for (let i = 0; i < arguments.length; i++)
			{
				this.add(arguments[i]);
			}

			return this;
		}

		if (object === this)
		{
			console.error('Object3D.add: object can\'t be added as a child of itself.', object);
			return this;
		}

		if (object && object.isObject3D)
		{
			if (object.parent !== null)
			{
				object.parent.remove(object);
			}

			object.parent = this;
			this.children.push(object);

			object.dispatchEvent(_addedEvent);
		} else
		{
			console.error('Object3D.add: object not an instance of Object3D.', object);
		}

		return this;
	}

	remove(object)
	{
		if (arguments.length > 1)
		{
			for (let i = 0; i < arguments.length; i++)
			{
				this.remove(arguments[i]);
			}

			return this;
		}

		const index = this.children.indexOf(object);

		if (index !== - 1)
		{
			object.parent = null;
			this.children.splice(index, 1);

			object.dispatchEvent(_removedEvent);
		}

		return this;
	}

	getWorldPosition(target)
	{
		this.updateWorldMatrix(true, false);
		return target.setFromMatrixPosition(this.matrixWorld);
	}

	getWorldQuaternion(target)
	{
		this.updateWorldMatrix(true, false);
		this.matrixWorld.decompose(_position, target, _scale);
		return target;
	}

	getWorldScale(target)
	{
		this.updateWorldMatrix(true, false);
		this.matrixWorld.decompose(_position, _quaternion, target);
		return target;
	}

	getWorldDirection(target)
	{
		this.updateWorldMatrix(true, false);
		const e = this.matrixWorld.elements;
		return target.set(e[8], e[9], e[10]).normalize();
	}
	
	raycast( /* raycaster, intersects */) { }

	traverse(callback)
	{
		callback(this);

		const children = this.children;

		for (let i = 0, l = children.length; i < l; i++)
		{
			children[i].traverse(callback);
		}
	}

	updateMatrix()
	{
		this.matrix.compose(this.position, this.quaternion, this.scale);

		this.matrixWorldNeedsUpdate = true;
	}

	updateMatrixWorld(force)
	{
		if (this.matrixAutoUpdate) this.updateMatrix();

		if (this.matrixWorldNeedsUpdate || force)
		{
			if (this.parent === null)
			{
				this.matrixWorld.copy(this.matrix);
			} else
			{
				this.matrixWorld.multiplyMatrices(this.parent.matrixWorld, this.matrix);
			}

			this.matrixWorldNeedsUpdate = false;

			force = true;
		}

		// update children

		const children = this.children;

		for (let i = 0, l = children.length; i < l; i++)
		{
			children[i].updateMatrixWorld(force);
		}
	}

	updateWorldMatrix(updateParents, updateChildren)
	{
		const parent = this.parent;

		if (updateParents === true && parent !== null)
		{
			parent.updateWorldMatrix(true, false);
		}

		if (this.matrixAutoUpdate) this.updateMatrix();

		if (this.parent === null)
		{
			this.matrixWorld.copy(this.matrix);
		} else
		{
			this.matrixWorld.multiplyMatrices(this.parent.matrixWorld, this.matrix);
		}

		// update children

		if (updateChildren === true)
		{
			const children = this.children;

			for (let i = 0, l = children.length; i < l; i++)
			{
				children[i].updateWorldMatrix(false, true);
			}
		}
	}

}

Object3D.DefaultUp = new Vector3(0, 1, 0);
Object3D.DefaultMatrixAutoUpdate = true;

Object3D.prototype.isObject3D = true;


const _q1 = new Quaternion();

const _target = new Vector3();
const _position = new Vector3();
const _scale = new Vector3();

const _xAxis = new Vector3(1, 0, 0);
const _yAxis = new Vector3(0, 1, 0);
const _zAxis = new Vector3(0, 0, 1);

const _addedEvent = { type: 'added' };
const _removedEvent = { type: 'removed' };