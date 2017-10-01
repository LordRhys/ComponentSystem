#ifndef _MATRIX_H
#define _MATRIX_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <cmath>
#include <cfloat>

#include "MathUtils.h"
#include "Vector3.h"
#include "Vector4.h"
#include "Quaternion.h"
//#include "StringUtils.h"
//#include "Logger.h"

	class Matrix4
	{
	private:
		float value[16];
	public:
		Matrix4(void)
		{
			makeZero();
		}

		Matrix4(bool isZero)
		{
			if(isZero)
			{
				makeZero();
			}
			else
			{
				makeIdentity();
			}
		}

		Matrix4(float a, float b, float c, float d, float e, float f, float g, float h, float i, float j, float k, float l, float m, float n, float o, float p)
		{
			value[0] = a;
			value[1] = b;
			value[2] = c;
			value[3] = d;
			value[4] = e;
			value[5] = f;
			value[6] = g;
			value[7] = h;
			value[8] = i;
			value[9] = j;
			value[10] = k;
			value[11] = l;
			value[12] = m;
			value[13] = n;
			value[14] = o;
			value[15] = p;
		}

		~Matrix4(void)
		{
		}

		void makeZero()
		{
			value[0] = (float)0;
			value[1] = (float)0;
			value[2] = (float)0;
			value[3] = (float)0;
			value[4] = (float)0;
			value[5] = (float)0;
			value[6] = (float)0;
			value[7] = (float)0;
			value[8] = (float)0;
			value[9] = (float)0;
			value[10] = (float)0;
			value[11] = (float)0;
			value[12] = (float)0;
			value[13] = (float)0;
			value[14] = (float)0;
			value[15] = (float)0;
		}

		void makeIdentity()
		{
			value[0] = (float)1;
			value[1] = (float)0;
			value[2] = (float)0;
			value[3] = (float)0;
			value[4] = (float)0;
			value[5] = (float)1;
			value[6] = (float)0;
			value[7] = (float)0;
			value[8] = (float)0;
			value[9] = (float)0;
			value[10] = (float)1;
			value[11] = (float)0;
			value[12] = (float)0;
			value[13] = (float)0;
			value[14] = (float)0;
			value[15] = (float)1;
		}

		inline float determinant()
		{
			float a0 = value[0] * value[5] - value[1] * value[4];
			float a1 = value[0] * value[6] - value[2] * value[4];
			float a2 = value[0] * value[7] - value[3] * value[4];
			float a3 = value[1] * value[6] - value[2] * value[5];
			float a4 = value[1] * value[7] - value[3] * value[5];
			float a5 = value[2] * value[7] - value[3] * value[6];
			float b0 = value[8] * value[13] - value[9] * value[12];
			float b1 = value[8] * value[14] - value[10] * value[12];
			float b2 = value[8] * value[15] - value[11] * value[12];
			float b3 = value[9] * value[14] - value[10] * value[13];
			float b4 = value[9] * value[15] - value[11] * value[13];
			float b5 = value[10] * value[15] - value[11] * value[14];
			return a0 * b5 - a1 * b4 + a2 * b3 + a3 * b2 - a4 * b1 + a5* b0;
		}

		inline Matrix4 transpose()
		{
			Matrix4 result;
			result[0] = value[0];
			result[1] = value[4];
			result[2] = value[8];
			result[3] = value[12];
			result[4] = value[1];
			result[5] = value[5];
			result[6] = value[9];
			result[7] = value[13];
			result[8] = value[2];
			result[9] = value[6];
			result[10] = value[10];
			result[11] = value[14];
			result[12] = value[3];
			result[13] = value[7];
			result[14] = value[11];
			result[15] = value[15];
			return Matrix4(result[0],result[1],result[2],result[3],
				result[4],result[5],result[6],result[7],
				result[8],result[9],result[10],result[11],
				result[12],result[13],result[14],result[15]);
		}

		inline Matrix4 adjoint()
		{
			float a0 = value[ 0]*value[ 5] - value[ 1]*value[ 4];
			float a1 = value[ 0]*value[ 6] - value[ 2]*value[ 4];
			float a2 = value[ 0]*value[ 7] - value[ 3]*value[ 4];
			float a3 = value[ 1]*value[ 6] - value[ 2]*value[ 5];
			float a4 = value[ 1]*value[ 7] - value[ 3]*value[ 5];
			float a5 = value[ 2]*value[ 7] - value[ 3]*value[ 6];
			float b0 = value[ 8]*value[13] - value[ 9]*value[12];
			float b1 = value[ 8]*value[14] - value[10]*value[12];
			float b2 = value[ 8]*value[15] - value[11]*value[12];
			float b3 = value[ 9]*value[14] - value[10]*value[13];
			float b4 = value[ 9]*value[15] - value[11]*value[13];
			float b5 = value[10]*value[15] - value[11]*value[14];

			Matrix4 result;
		
			result[0] =	+value[ 5]*b5 - value[ 6]*b4 + value[ 7]*b3;
			result[1] =	-value[ 1]*b5 + value[ 2]*b4 - value[ 3]*b3;
			result[2] =	+value[13]*a5 - value[14]*a4 + value[15]*a3;
			result[3] =	-value[ 9]*a5 + value[10]*a4 - value[11]*a3;
			result[4] =	-value[ 4]*b5 + value[ 6]*b2 - value[ 7]*b1;
			result[5] =	+value[ 0]*b5 - value[ 2]*b2 + value[ 3]*b1;
			result[6] =	-value[12]*a5 + value[14]*a2 - value[15]*a1;
			result[7] =	+value[ 8]*a5 - value[10]*a2 + value[11]*a1;
			result[8] =	+value[ 4]*b4 - value[ 5]*b2 + value[ 7]*b0;
			result[9] =	-value[ 0]*b4 + value[ 1]*b2 - value[ 3]*b0;
			result[10] = +value[12]*a4 - value[13]*a2 + value[15]*a0;
			result[11] = -value[ 8]*a4 + value[ 9]*a2 - value[11]*a0;
			result[12] = -value[ 4]*b3 + value[ 5]*b1 - value[ 6]*b0;
			result[13] = +value[ 0]*b3 - value[ 1]*b1 + value[ 2]*b0;
			result[14] = -value[12]*a3 + value[13]*a1 - value[14]*a0;
			result[15] = +value[ 8]*a3 - value[ 9]*a1 + value[10]*a0;
		
			return Matrix4(result[0],result[1],result[2],result[3],
				result[4],result[5],result[6],result[7],
				result[8],result[9],result[10],result[11],
				result[12],result[13],result[14],result[15]);
		}

		inline Matrix4 inverse()
		{
			float determinant = this->determinant();

			Matrix4 inverse;

			if (std::fabs(determinant) >= FLT_EPSILON)
			{
				Matrix4 adjoint = this->adjoint();
				inverse = adjoint*(1/determinant);
				return inverse;
			}

			inverse.makeZero();
			return Matrix4(inverse[0],inverse[1],inverse[2],inverse[3],
				inverse[4],inverse[5],inverse[6],inverse[7],
				inverse[8],inverse[9],inverse[10],inverse[11],
				inverse[12],inverse[13],inverse[14],inverse[15]);
		}

		inline Matrix4 fastInverse()
		{
			Matrix3 outRot;

			outRot[0] = value[0];
			outRot[1] = value[4];
			outRot[2] = value[8];

			outRot[3] = value[1];
			outRot[4] = value[5];
			outRot[5] = value[9];

			outRot[6] = value[2];
			outRot[7] = value[6];
			outRot[8] = value[10];

			Vector3 outTrans;

			outTrans[0] = value[3];
			outTrans[1] = value[7];
			outTrans[2] = value[11];

			outTrans = outRot*outTrans;

			outTrans = -outTrans;

			Matrix4 result;

			result[0] = outRot[0];
			result[1] = outRot[1];
			result[2] = outRot[2];
			result[3] = outTrans[0];

			result[4] = outRot[3];
			result[5] = outRot[4];
			result[6] = outRot[5];
			result[7] = outTrans[1];

			result[8] = outRot[6];
			result[9] = outRot[7];
			result[10] = outRot[8];
			result[11] = outTrans[2];

			result[12] = 0;
			result[13] = 0;
			result[14] = 0;
			result[15] = 1;

			return Matrix4(result[0],result[1],result[2],result[3],
				result[4],result[5],result[6],result[7],
				result[8],result[9],result[10],result[11],
				result[12],result[13],result[14],result[15]);

		}

		void makeOrthoProjection(float const &left, float const &right, float const &bottom, float const &top, float const &zNear, float const &zFar)
		{
			value[0] = (float)(2) / (right - left);
			value[5] = (float)(2) / (top - bottom);
			value[10] = (float)(-2) / (zFar - zNear);

			value[3] = - ((right + left) / (right - left));
			value[7] = - ((top + bottom) / (top - bottom));
			value[11] = - ((zFar + zNear) / (zFar - zNear));

			value[15] = 1;

		}

		void makePerspectiveProjection(float const &fov, float const &aspect, float const &zNear, float const &zFar)
		{
			float range = std::tan(MathUtils::toRadians(fov / (float)(2))) * zNear;	
			float left = -range * aspect;
			float right = range * aspect;
			float bottom = -range;
			float top = range;

			value[0]  = (2.0f * zNear) / (right - left);
			value[1]  = 0;
			value[2]  = 0;
			value[3]  = 0;

			value[4]  = 0;
			value[5]  = (2.0f * zNear) / (top - bottom);
			value[6]  = 0;
			value[7]  = 0;

			value[8]  = 0;
			value[9]  = 0;
			value[10] = -(zFar + zNear) / (zFar - zNear);
			value[11] = -(2.0f * zFar * zNear) / (zFar - zNear);

			value[12] = 0;
			value[13] = 0;
			value[14] = -(1.0f);
			value[15] = 0;
		}

		void makeOrthoView(const Vector3& normal, const Vector3& origin, const Vector3& direction)
		{
			float dotND = normal.dot(direction);
			float dotNO = normal.dot(origin);
			value[0] = direction[0]*normal[0] - dotND;
			value[1] = direction[0]*normal[1];
			value[2] = direction[0]*normal[2];
			value[3] = -dotNO * direction[0];
			value[4] = direction[1]*normal[0];
			value[5] = direction[1]*normal[1] - dotND;
			value[6] = direction[1]*normal[2];
			value[7] = -dotNO * direction[1];
			value[8] = direction[2]*normal[0];
			value[9] = direction[2]*normal[1];
			value[10] = direction[2]*normal[2] - dotND;
			value[11] = -dotNO * direction[2];
			value[12] = 0.0f;
			value[13] = 0.0f;
			value[14] = 0.0f;
			value[15] = -dotND;
		}

		void makePerspectiveView(const Vector3& normal, const Vector3& origin, const Vector3& eye)
		{
			float dotND = normal.dot(eye - origin);
			value[0] = dotND - eye[0]*normal[0];
			value[1] = -eye[0]*normal[1];
			value[2] = -eye[0]*normal[2];
			value[3] = -(value[0]*eye[0] + value[1]*eye[1] + value[2]*eye[2]);
			value[4] = -eye[1]*normal[0];
			value[5] = dotND - eye[1]*normal[1];
			value[6] = -eye[1]*normal[2];
			value[7] = -(value[4]*eye[0] + value[5]*eye[1] + value[6]*eye[2]);
			value[8] = -eye[2]*normal[0];
			value[9] = -eye[2]*normal[1];
			value[10] = dotND- eye[2]*normal[2];
			value[11] = -(value[8]*eye[0] + value[9]*eye[1] + value[10]*eye[2]);
			value[12] = -normal[0];
			value[13] = -normal[1];
			value[14] = -normal[2];
			value[15] = normal.dot(eye);
		}

		void makeReflection(const Vector3& normal,const Vector3& origin)
		{

			float twoDotNO = ((float)2)*(normal.dot(origin));
			value[ 0] = (float)1 - ((float)2)*normal[0]*normal[0];
			value[ 1] = -((float)2)*normal[0]*normal[1];
			value[ 2] = -((float)2)*normal[0]*normal[2];
			value[ 3] = twoDotNO*normal[0];
			value[ 4] = -((float)2)*normal[1]*normal[0];
			value[ 5] = (float)1 - ((float)2)*normal[1]*normal[1];
			value[ 6] = -((float)2)*normal[1]*normal[2];
			value[ 7] = twoDotNO*normal[1];
			value[ 8] = -((float)2)*normal[2]*normal[0];
			value[ 9] = -((float)2)*normal[2]*normal[1];
			value[10] = (float)1 - ((float)2)*normal[2]*normal[2];
			value[11] = twoDotNO*normal[2];
			value[12] = (float)0;
			value[13] = (float)0;
			value[14] = (float)0;
			value[15] = (float)1;
		}

		inline Matrix3 toMatrix3()
		{
			return Matrix3(value[0],value[1],value[2],value[4],value[5],value[6],value[8],value[9],value[10]);
		}

		static inline Matrix4 toMatrix4(Quaternion q, Vector3 v, float scale)
		{
			Matrix3 rotation;

			q.toRotationMatrix(rotation);

			rotation *= scale;

			return Matrix4(rotation[0],rotation[1],rotation[2],v[0],
				rotation[3],rotation[4],rotation[5],v[1],
				rotation[6],rotation[7],rotation[8],v[2],
				0,0,0,1);
		}

		static inline Matrix4 toMatrix4(Quaternion q, Vector3 v, Matrix3 scale)
		{
			Matrix3 rotation;

			q.toRotationMatrix(rotation);

			rotation *= scale;

			return Matrix4(rotation[0],rotation[1],rotation[2],v[0],
				rotation[3],rotation[4],rotation[5],v[1],
				rotation[6],rotation[7],rotation[8],v[2],
				0,0,0,1);
		}

		static inline Matrix4 toMatrix4(Quaternion q, Vector3 v)
		{
			Matrix3 rotation;

			q.toRotationMatrix(rotation);

			return Matrix4(rotation[0],rotation[1],rotation[2],v[0],
				rotation[3],rotation[4],rotation[5],v[1],
				rotation[6],rotation[7],rotation[8],v[2],
				0,0,0,1);
		}

		static inline Matrix4 toMatrix4(Vector3 right, Vector3 up, Vector3 direction, Vector3 v)
		{
			return Matrix4(right[0],up[1],direction[2],v[0],
				right[3],up[4],direction[5],v[1],
				right[6],up[7],direction[8],v[2],
				0,0,0,1);
		}

		inline void operator +=(const Matrix4 &v)
		{
			(*this) = (*this) + v;
		}

		inline void operator -=(const Matrix4 &v)
		{
			(*this) = (*this) - v;
		}

		inline void operator *=(const Matrix4 &v)
		{
			(*this) = (*this) * v;
		}

		inline Matrix4 operator +(const Matrix4 &v)
		{
			Matrix4 result;
			result[0] = value[0] + v[0];
			result[1] = value[1] + v[1];
			result[2] = value[2] + v[2];
			result[3] = value[3] + v[3];
			result[4] = value[4] + v[4];
			result[5] = value[5] + v[5];
			result[6] = value[6] + v[6];
			result[7] = value[7] + v[7];
			result[8] = value[8] + v[8];
			result[9] = value[9] + v[9];
			result[10] = value[10] + v[10];
			result[11] = value[11] + v[11];
			result[12] = value[12] + v[12];
			result[13] = value[13] + v[13];
			result[14] = value[14] + v[14];
			result[15] = value[15] + v[15];
			return Matrix4(result[0],result[1],result[2],result[3],
				result[4],result[5],result[6],result[7],
				result[8],result[9],result[10],result[11],
				result[12],result[13],result[14],result[15]);
		}

		inline Matrix4 operator -(const Matrix4 &v)
		{
			Matrix4 result;
			result[0] = value[0] - v[0];
			result[1] = value[1] - v[1];
			result[2] = value[2] - v[2];
			result[3] = value[3] - v[3];
			result[4] = value[4] - v[4];
			result[5] = value[5] - v[5];
			result[6] = value[6] - v[6];
			result[7] = value[7] - v[7];
			result[8] = value[8] - v[8];
			result[9] = value[9] - v[9];
			result[10] = value[10] - v[10];
			result[11] = value[11] - v[11];
			result[12] = value[12] - v[12];
			result[13] = value[13] - v[13];
			result[14] = value[14] - v[14];
			result[15] = value[15] - v[15];
			return Matrix4(result[0],result[1],result[2],result[3],
				result[4],result[5],result[6],result[7],
				result[8],result[9],result[10],result[11],
				result[12],result[13],result[14],result[15]);
		}

		inline Matrix4 operator *(const Matrix4 &v)
		{
			Matrix4 result;

			result[0] = (value[0] * v[0]) + (value[1] * v[4]) + (value[2] * v[8]) + (value[3] * v[12]);
			result[4] = (value[0] * v[1]) + (value[1] * v[5]) + (value[2] * v[9]) + (value[3] * v[13]);
			result[8] = (value[0] * v[2]) + (value[1] * v[6]) + (value[2] * v[10]) + (value[3] * v[14]);
			result[12] = (value[0] * v[3]) + (value[1] * v[7]) + (value[2] * v[11]) + (value[3] * v[15]);
			result[1] = (value[4] * v[0]) + (value[5] * v[4]) + (value[6] * v[8]) + (value[7] * v[12]);
			result[5] = (value[4] * v[1]) + (value[5] * v[5]) + (value[6] * v[9]) + (value[7] * v[13]);
			result[9] = (value[4] * v[2]) + (value[5] * v[6]) + (value[6] * v[10]) + (value[7] * v[14]);
			result[13] = (value[4] * v[3]) + (value[5] * v[7]) + (value[6] * v[11]) + (value[7] * v[15]);
			result[2] = (value[8] * v[0]) + (value[9] * v[4]) + (value[10] * v[8]) + (value[11] * v[12]);
			result[6] = (value[8] * v[1]) + (value[9] * v[5]) + (value[10] * v[9]) + (value[11] * v[13]);
			result[10] = (value[8] * v[2]) + (value[9] * v[6]) + (value[10] * v[10]) + (value[11] * v[14]);
			result[14] = (value[8] * v[3]) + (value[9] * v[7]) + (value[10] * v[11]) + (value[11] * v[15]);
			result[3] = (value[12] * v[0]) + (value[13] * v[4]) + (value[14] * v[8]) + (value[15] * v[12]);
			result[7] = (value[12] * v[1]) + (value[13] * v[5]) + (value[14] * v[9]) + (value[15] * v[13]);
			result[11] = (value[12] * v[2]) + (value[13] * v[6]) + (value[14] * v[10]) + (value[15] * v[14]);
			result[15] = (value[12] * v[3]) + (value[13] * v[7]) + (value[14] * v[11]) + (value[15] * v[15]);

			return Matrix4(result[0],result[1],result[2],result[3],
				result[4],result[5],result[6],result[7],
				result[8],result[9],result[10],result[11],
				result[12],result[13],result[14],result[15]).transpose();
		}

		inline Matrix4 operator +(const float s)
		{
			Matrix4 result;
			result[0] = value[0] + s;
			result[1] = value[1] + s;
			result[2] = value[2] + s;
			result[3] = value[3] + s;
			result[4] = value[4] + s;
			result[5] = value[5] + s;
			result[6] = value[6] + s;
			result[7] = value[7] + s;
			result[8] = value[8] + s;
			result[9] = value[9] + s;
			result[10] = value[10] + s;
			result[11] = value[11] + s;
			result[12] = value[12] + s;
			result[13] = value[13] + s;
			result[14] = value[14] + s;
			result[15] = value[15] + s;
			return Matrix4(result[0],result[1],result[2],result[3],
				result[4],result[5],result[6],result[7],
				result[8],result[9],result[10],result[11],
				result[12],result[13],result[14],result[15]);
		}

		inline Matrix4 operator -(const float s)
		{
			Matrix4 result;
			result[0] = value[0] - s;
			result[1] = value[1] - s;
			result[2] = value[2] - s;
			result[3] = value[3] - s;
			result[4] = value[4] - s;
			result[5] = value[5] - s;
			result[6] = value[6] - s;
			result[7] = value[7] - s;
			result[8] = value[8] - s;
			result[9] = value[9] - s;
			result[10] = value[10] - s;
			result[11] = value[11] - s;
			result[12] = value[12] - s;
			result[13] = value[13] - s;
			result[14] = value[14] - s;
			result[15] = value[15] - s;
			return Matrix4(result[0],result[1],result[2],result[3],
				result[4],result[5],result[6],result[7],
				result[8],result[9],result[10],result[11],
				result[12],result[13],result[14],result[15]);
		}

		inline Matrix4 operator *(const float s)
		{
			Matrix4 result;
			result[0] = value[0] * s;
			result[1] = value[1] * s;
			result[2] = value[2] * s;
			result[3] = value[3] * s;
			result[4] = value[4] * s;
			result[5] = value[5] * s;
			result[6] = value[6] * s;
			result[7] = value[7] * s;
			result[8] = value[8] * s;
			result[9] = value[9] * s;
			result[10] = value[10] * s;
			result[11] = value[11] * s;
			result[12] = value[12] * s;
			result[13] = value[13] * s;
			result[14] = value[14] * s;
			result[15] = value[15] * s;
			return Matrix4(result[0],result[1],result[2],result[3],
				result[4],result[5],result[6],result[7],
				result[8],result[9],result[10],result[11],
				result[12],result[13],result[14],result[15]);
		}

		inline Matrix4 operator /(const float s)
		{
			Matrix4 result;
			result[0] = value[0] / s;
			result[1] = value[1] / s;
			result[2] = value[2] / s;
			result[3] = value[3] / s;
			result[4] = value[4] / s;
			result[5] = value[5] / s;
			result[6] = value[6] / s;
			result[7] = value[7] / s;
			result[8] = value[8] / s;
			result[9] = value[9] / s;
			result[10] = value[10] / s;
			result[11] = value[11] / s;
			result[12] = value[12] / s;
			result[13] = value[13] / s;
			result[14] = value[14] / s;
			result[15] = value[15] / s;
			return Matrix4(result[0],result[1],result[2],result[3],
				result[4],result[5],result[6],result[7],
				result[8],result[9],result[10],result[11],
				result[12],result[13],result[14],result[15]);
		}

		inline Vector4 operator* (const Vector4& v)
		{
			return Vector4(
				v[0]*value[0]+v[1]*value[1]+v[2]*value[2]+v[3]*value[3],
				v[0]*value[4]+v[1]*value[5]+v[2]*value[6]+v[3]*value[7],
				v[0]*value[8]+v[1]*value[9]+v[2]*value[10]+v[3]*value[11],
				v[0]*value[12]+v[1]*value[13]+v[2]*value[14]+v[3]*value[15]);
		}

		inline const float* getArray() const
		{
			return value;
		}

		inline float& operator [](unsigned int i)
		{
			return value[i];
		}

		inline float operator [](unsigned int i) const
		{
			return value[i];
		}

		/*void toString()
		{
			Logger::getInstance()->write(StringUtils::format("Matrix4[ %f, %f, %f, %f ]",value[0],value[1],value[2],value[3]));
			Logger::getInstance()->write(StringUtils::format("Matrix4[ %f, %f, %f, %f ]",value[4],value[5],value[6],value[7]));
			Logger::getInstance()->write(StringUtils::format("Matrix4[ %f, %f, %f, %f ]",value[8],value[9],value[10],value[11]));
			Logger::getInstance()->write(StringUtils::format("Matrix4[ %f, %f, %f, %f ]",value[12],value[13],value[14],value[15]));
		}*/
	};


#endif