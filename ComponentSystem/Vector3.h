#ifndef _VECTOR3_H
#define _VECTOR3_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <cmath>
#include <cfloat>

//#include "StringUtils.h"
//#include "Logger.h"

	class Vector3
	{
	private:
		float value[3];
	public:

		Vector3(void)
		{
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
		}

		Vector3(float _x, float _y, float _z)
		{
			value[0] = _x;
			value[1] = _y;
			value[2] = _z;
		}

		~Vector3()
		{
		}

		inline float angleBetween(Vector3 &v)
		{
			return acos(this->dot(v) / this->magnitude() / v.magnitude());
		}

		inline Vector3 projection(Vector3 &v)
		{
			return v*((this->dot(v))/(std::pow(v.magnitude(),2)));
		}

		inline float magnitude()
		{
			return std::sqrt((value[0]*value[0])+(value[1]*value[1])+(value[2]*value[2]));
		}

		inline Vector3 normalize()
		{
			return (*this) / this->magnitude();
		}

		inline float dot(const Vector3 &v) const
		{
			return value[0]*v[0]+value[1]*v[1]+value[2]*v[2];
		}

		inline Vector3 cross(const Vector3 &v)
		{
			Vector3 result;
			result[0] = value[1]*v[2] - value[2]*v[1]; 
			result[1] = value[2]*v[0] - value[0]*v[2];
			result[2] = value[0]*v[1] - value[1]*v[0];
			return Vector3(result[0],result[1],result[2]);
		}

		static inline Vector3 cross(const Vector3 &v0, const Vector3 &v1)
		{
			Vector3 result;
			result[0] = v0[1]*v1[2] - v0[2]*v1[1]; 
			result[1] = v0[2]*v1[0] - v0[0]*v1[2];
			result[2] = v0[0]*v1[1] - v0[1]*v1[0];
			return Vector3(result[0],result[1],result[2]);
		}

		inline Vector3 lerp(Vector3 &v, float step)
		{
			return ((*this) * (1.0f - step)) + (v * step);
		}

		inline Vector3 catMullRom(Vector3 &v0, Vector3 &v1, Vector3 &v2, float step)
		{
			float step2 = std::pow(step,2);
			float step3 = std::pow(step,3);

			Vector3 result;

			result[0] = 0.5f * ( ( 2.0f * v0[0] ) +
			( -value[0] + v1[0] ) * step +
			( 2.0f * value[0] - 5.0f * v0[0] + 4 * v1[0] - v2[0] ) * step2 +
			( -value[0] + 3.0f * v0[0] - 3.0f * v1[0] + v2[0] ) * step3 );

			result[1] = 0.5f * ( ( 2.0f * v0[1] ) +
			( -value[1] + v1[1] ) * step +
			( 2.0f * value[1] - 5.0f * v0[1] + 4 * v1[1] - v2[1] ) * step2 +
			( -value[1] + 3.0f * v0[1] - 3.0f * v1[1] + v2[1] ) * step3 );

			result[2] = 0.5f * ( ( 2.0f * v0[2] ) +
			( -value[2] + v1[2] ) * step +
			( 2.0f * value[1] - 5.0f * v0[2] + 4 * v1[2] - v2[2] ) * step2 +
			( -value[2] + 3.0f * v0[2] - 3.0f * v1[2] + v2[2] ) * step3 );

			return Vector3(result[0],result[1],result[2]);
		}

		void Vector3::Orthonormalize(Vector3& u, Vector3& v, Vector3& w)
		{

			u.normalize();

			float dot0 = u.dot(v); 
			v -= u*dot0;
			v.normalize();

			float dot1 = v.dot(w);
			dot0 = u.dot(w);
			w -= u*dot0 + v*dot1;
			w.normalize();
		}

		inline bool isZero()
		{
			if(value[0] == (float)0 && value[1] == (float)0 && value[2] == (float)0)
			{
				return true;
			}

			return false;
		}

		inline float getX()
		{
			return value[0];
		}

		inline float getY()
		{
			return value[1];
		}

		inline float getZ()
		{
			return value[2];
		}

		inline void setX(float _x)
		{
			value[0] = _x;
		}

		inline void setY(float _y)
		{
			value[1] = _y;
		}

		inline void setZ(float _z)
		{
			value[2] = _z;
		}

		inline void operator =(const Vector3 &v)
		{
			value[0] = v[0];
			value[1] = v[1];
			value[2] = v[2];
		}

		inline bool operator ==(const Vector3 &v)
		{
			if(value[0] == v[0] && value[1] == v[1] && value[2] == v[2]) { return true; }
			else { return false; }
		}

		inline bool operator !=(const Vector3 &v)
		{
			if(value[0] != v[0] || value[1] != v[1] || value[2] != v[2]) { return true; }
			else { return false; }
		}

		inline void operator +=(const float s)
		{
			value[0] = value[0] + s;
			value[1] = value[1] + s;
			value[2] = value[2] + s;
		}

		inline void operator -=(const float s)
		{
			value[0] = value[0] - s;
			value[1] = value[1] - s;
			value[2] = value[2] - s;
		}

		inline void operator *=(const float s)
		{
			value[0] = value[0] * s;
			value[1] = value[1] * s;
			value[2] = value[2] * s;
		}

		inline void operator /=(const float s)
		{
			value[0] = value[0] / s;
			value[1] = value[1] / s;
			value[2] = value[2] / s;
		}

		inline void operator +=(const Vector3 &v)
		{
			value[0] += v[0];
			value[1] += v[1];
			value[2] += v[2];
		}

		inline void operator -=(const Vector3 &v)
		{
			value[0] -= v[0];
			value[1] -= v[1];
			value[2] -= v[2];
		}

		inline Vector3 operator +(const Vector3 &v)
		{
			Vector3 result;
			result[0] = value[0] + v[0];
			result[1] = value[1] + v[1];
			result[2] = value[2] + v[2];
			return Vector3(result[0],result[1],result[2]);
		}

		inline Vector3 operator -(const Vector3 &v) const
		{
			Vector3 result;
			result[0] = value[0] - v[0];
			result[1] = value[1] - v[1];
			result[2] = value[2] - v[2];
			return Vector3(result[0],result[1],result[2]);
		}

		inline Vector3 operator -() const
		{
			Vector3 result;
			result[0] = -value[0];
			result[1] = -value[1];
			result[2] = -value[2];
			return Vector3(result[0],result[1],result[2]);
		}

		inline Vector3 operator +(const float s)
		{
			Vector3 result;
			result[0] = value[0] + s;
			result[1] = value[1] + s;
			result[2] = value[2] + s;
			return Vector3(result[0],result[1],result[2]);
		}

		inline Vector3 operator -(const float s)
		{
			Vector3 result;
			result[0] = value[0] - s;
			result[1] = value[1] - s;
			result[2] = value[2] - s;
			return Vector3(result[0],result[1],result[2]);
		}

		inline Vector3 operator *(const float s)
		{
			Vector3 result;
			result[0] = value[0] * s;
			result[1] = value[1] * s;
			result[2] = value[2] * s;
			return Vector3(result[0],result[1],result[2]);
		}

		inline Vector3 operator /(const float s)
		{
			Vector3 result;
			result[0] = value[0] / s;
			result[1] = value[1] / s;
			result[2] = value[2] / s;
			return Vector3(result[0],result[1],result[2]);
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
			Logger::getInstance()->write(StringUtils::format("Vector3[ %f, %f, %f ]",value[0],value[1],value[2]));
		}*/
	};



#endif

