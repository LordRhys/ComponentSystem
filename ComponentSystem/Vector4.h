#ifndef _Vector4_H
#define _Vector4_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <cmath>
#include <cfloat>

//#include "StringUtils.h"
//#include "Logger.h"

	class Vector4
	{
	private:
		float value[4];
	public:

		Vector4(void)
		{
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			value[3] = 1;
		}

		Vector4(float _x, float _y, float _z)
		{
			value[0] = _x;
			value[1] = _y;
			value[2] = _z;
			value[3] = 1;
		}

		Vector4(float _x, float _y, float _z, float _w)
		{
			value[0] = _x;
			value[1] = _y;
			value[2] = _z;
			value[3] = _w;
		}

		Vector4(Vector3 v)
		{
			value[0] = v[0];
			value[1] = v[1];
			value[2] = v[2];
			value[3] = (float)1;
		}

		~Vector4()
		{
		}

		inline Vector3 toVector3()
		{
			return Vector3(value[0],value[1],value[1]);
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

		inline float getW()
		{
			return value[3];
		}

		inline void setX(float _x)
		{
			value[2] = _x;
		}

		inline void setY(float _y)
		{
			value[2] = _y;
		}

		inline void setZ(float _z)
		{
			value[2] = _z;
		}

		inline void setW(float _w)
		{
			value[3] = _w;
		}

		inline void operator =(const Vector4 &v)
		{
			value[0] = v[0];
			value[1] = v[1];
			value[2] = v[2];
			value[3] = v[3];
		}

		inline bool operator ==(const Vector4 &v)
		{
			if(value[0] == v[0] && value[1] == v[1] && value[2] == v[2] && value[3] == v[3]) { return true; }
			else { return false; }
		}

		inline bool operator !=(const Vector4 &v)
		{
			if(value[0] != v[0] || value[1] != v[1] || value[2] != v[2] || value[3] != v[3]) { return true; }
			else { return false; }
		}

		inline void operator +=(const float s)
		{
			value[0] = value[0] + s;
			value[1] = value[1] + s;
			value[2] = value[2] + s;
			value[3] = value[3] + s;
		}

		inline void operator -=(const float s)
		{
			value[0] = value[0] - s;
			value[1] = value[1] - s;
			value[2] = value[2] - s;
			value[3] = value[3] - s;
		}

		inline void operator *=(const float s)
		{
			value[0] = value[0] * s;
			value[1] = value[1] * s;
			value[2] = value[2] * s;
			value[3] = value[3] * s;
		}

		inline void operator /=(const float s)
		{
			value[0] = value[0] / s;
			value[1] = value[1] / s;
			value[2] = value[2] / s;
			value[3] = value[3] / s;
		}

		inline void operator +=(const Vector4 &v)
		{
			value[0] += v[0];
			value[1] += v[1];
			value[2] += v[2];
			value[3] += v[3];
		}

		inline void operator -=(const Vector4 &v)
		{
			value[0] -= v[0];
			value[1] -= v[1];
			value[2] -= v[2];
			value[3] -= v[3];
		}

		inline Vector4 operator +(const Vector4 &v)
		{
			Vector4 result;
			result[0] = value[0] + v[0];
			result[1] = value[1] + v[1];
			result[2] = value[2] + v[2];
			result[3] = value[3] + v[3];
			return Vector4(result[0],result[1],result[2],result[3]);
		}

		inline Vector4 operator -(const Vector4 &v)
		{
			Vector4 result;
			result[0] = value[0] - v[0];
			result[1] = value[1] - v[1];
			result[2] = value[2] - v[2];
			result[3] = value[3] - v[3];
			return Vector4(result[0],result[1],result[2],result[3]);
		}

		inline Vector4 operator +(const float s)
		{
			Vector4 result;
			result[0] = value[0] + s;
			result[1] = value[1] + s;
			result[2] = value[2] + s;
			result[3] = value[3] + s;
			return Vector4(result[0],result[1],result[2],result[3]);
		}

		inline Vector4 operator -(const float s)
		{
			Vector4 result;
			result[0] = value[0] - s;
			result[1] = value[1] - s;
			result[2] = value[2] - s;
			result[3] = value[3] - s;
			return Vector4(result[0],result[1],result[2],result[3]);
		}

		inline Vector4 operator *(const float s)
		{
			Vector4 result;
			result[0] = value[0] * s;
			result[1] = value[1] * s;
			result[2] = value[2] * s;
			result[3] = value[3] * s;
			return Vector4(result[0],result[1],result[2],result[3]);
		}

		inline Vector4 operator /(const float s)
		{
			Vector4 result;
			result[0] = value[0] / s;
			result[1] = value[1] / s;
			result[2] = value[2] / s;
			result[3] = value[3] / s;
			return Vector4(result[0],result[1],result[2],result[3]);
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
			Logger::getInstance()->write(StringUtils::format("Vector4[ %f, %f, %f, %f ]",value[0],value[1],value[2],value[3]));
		}*/
	};


#endif