#ifndef _VECTOR2_H
#define _VECTOR2_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <cmath>
#include <cfloat>

//#include "StringUtils.h"
//#include "Logger.h"

	class Vector2
	{
	private:
		float value[2];
	public:
		Vector2(void)
		{
			value[0] = 0;
			value[1] = 0;
		}

		Vector2(float _x, float _y)
		{
			value[0] = _x;
			value[1] = _y;
		}

		~Vector2()
		{
		}

		inline float angleBetween(Vector2 &v)
		{
			return acos(this->dot(v) / this->magnitude() / v.magnitude());
		}

		inline Vector2 projection(Vector2 &v)
		{
			return v*((this->dot(v))/(std::pow(v.magnitude(),2)));
		}

		inline float magnitude()
		{
			return std::sqrt((value[0]*value[0])+(value[1]*value[1]));
		}

		inline Vector2 normalize()
		{
			return (*this) / this->magnitude();
		}

		inline float dot(const Vector2 &v) const
		{
			return value[0]*v[0]+value[1]*v[1];
		}

		inline Vector2 lerp(Vector2 &v, float step)
		{
			return ((*this) * (1.0f - step)) + (v * step);
		}

		inline Vector2 catMullRom(Vector2 &v0, Vector2 &v1, Vector2 &v2, float step)
		{
			float step2 = std::pow(step,2);
			float step3 = std::pow(step,3);

			Vector2 result;

			result[0] = 0.5f * ( ( 2.0f * v0[0] ) +
			( -value[0] + v1[0] ) * step +
			( 2.0f * value[0] - 5.0f * v0[0] + 4 * v1[0] - v2[0] ) * step2 +
			( -value[0] + 3.0f * v0[0] - 3.0f * v1[0] + v2[0] ) *  step3);

			result[1] = 0.5f * ( ( 2.0f * v0[1] ) +
			( -value[1] + v1[1] ) * step +
			( 2.0f * value[1] - 5.0f * v0[1] + 4 * v1[1] - v2[1] ) * step2 +
			( -value[1] + 3.0f * v0[1] - 3.0f * v1[1] + v2[1] ) * step3 );

			return Vector2(result[0],result[1]);
		}

		void Vector2::Orthonormalize (Vector2& u, Vector2& v, Vector2& w)
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
			if(value[0] == (float)0 && value[1] == (float)0)
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

		inline void setX(float _x)
		{
			value[2] = _x;
		}

		inline void setY(float _y)
		{
			value[2] = _y;
		}


		inline void operator =(const Vector2 &v)
		{
			value[0] = v[0];
			value[1] = v[1];
		}

		inline bool operator ==(const Vector2 &v)
		{
			if(value[0] == v[0] && value[1] == v[1]) { return true; }
			else { return false; }
		}

		inline bool operator !=(const Vector2 &v)
		{
			if(value[0] != v[0] || value[1] != v[1]) { return true; }
			else { return false; }
		}

		inline void operator +=(const float s)
		{
			value[0] = value[0] + s;
			value[1] = value[1] + s;
		}

		inline void operator -=(const float s)
		{
			value[0] = value[0] - s;
			value[1] = value[1] - s;
		}

		inline void operator *=(const float s)
		{
			value[0] = value[0] * s;
			value[1] = value[1] * s;
		}

		inline void operator /=(const float s)
		{
			value[0] = value[0] / s;
			value[1] = value[1] / s;
		}

		inline void operator +=(const Vector2 &v)
		{
			value[0] += v[0];
			value[1] += v[1];
		}

		inline void operator -=(const Vector2 &v)
		{
			value[0] -= v[0];
			value[1] -= v[1];
		}

		inline Vector2 operator +(const Vector2 &v)
		{
			Vector2 result;
			result[0] = value[0] + v[0];
			result[1] = value[1] + v[1];
			return Vector2(result[0],result[1]);
		}

		inline Vector2 operator -(const Vector2 &v) const
		{
			Vector2 result;
			result[0] = value[0] - v[0];
			result[1] = value[1] - v[1];
			return Vector2(result[0],result[1]);
		}

		inline Vector2 operator -() const
		{
			Vector2 result;
			result[0] = -value[0];
			result[1] = -value[1];
			return Vector2(result[0],result[1]);
		}

		inline Vector2 operator +(const float s)
		{
			Vector2 result;
			result[0] = value[0] + s;
			result[1] = value[1] + s;
			return Vector2(result[0],result[1]);
		}

		inline Vector2 operator -(const float s)
		{
			Vector2 result;
			result[0] = value[0] - s;
			result[1] = value[1] - s;
			return Vector2(result[0],result[1]);
		}

		inline Vector2 operator *(const float s)
		{
			Vector2 result;
			result[0] = value[0] * s;
			result[1] = value[1] * s;
			return Vector2(result[0],result[1]);
		}

		inline Vector2 operator /(const float s)
		{
			Vector2 result;
			result[0] = value[0] / s;
			result[1] = value[1] / s;
			return Vector2(result[0],result[1]);
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
			Logger::getInstance()->write(StringUtils::format("Vector2[ %f, %f ]",value[0],value[1]));
		}*/
	};


#endif