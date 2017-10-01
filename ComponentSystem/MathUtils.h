#ifndef _MATHUTILS_H
#define _MATHUTILS_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <cmath>

//#include "Random.h"
#include "Vector3.h"
#include "Vector2.h"

//#include "StringUtils.h"
//#include "Logger.h"
//
//#include "Colour.h"

#ifdef PARTICLESENGINE_DLL 
#define PARTICLESENGINE_API __declspec( dllexport )
#else
#define PARTICLESENGINE_API __declspec( dllimport )
#endif

namespace ParticlesEngine
{
	class PARTICLESENGINE_API MathUtils
	{
	public:
		static float toDegrees(float radians)
		{
			return (float)(radians*(180.0f/M_PI));
			return 0.0f;
		}

		static float toRadians(float degrees)
		{
			return (float)(degrees*(M_PI/180.0f));
			return 0.0f;
		}

		static float round(float r) 
		{
			return (r > 0.0f) ? floor(r + 0.5f) : ceil(r - 0.5f);
		}

		template <typename T>
		static T clamp(T min, T value, T max)
		{
			if(value < min)
			{
				return min;
			}
			else if(value > max)
			{
				return max;
			}

			return value;
		}

		template <typename T>
		static T minimum(T a, T b)
		{
			return (a < b ? a : b);
		}

		template <typename T>
		static T maximum(T a, T b)
		{
			return (a > b ? a : b);
		}

		static float lerp(float start, float end, float t)
		{
			return ((float)1-t)*start + t*end;
		}

		static Vector2 lerp(Vector2 start, Vector2 end, float dt)
		{
			float x,y;
			x = MathUtils::lerp(start[0],end[0],dt);
			y = MathUtils::lerp(start[1],end[1],dt);

			return Vector2(x,y);
		}

		static Vector3 lerp(Vector3 start, Vector3 end, float dt)
		{
			float x,y,z;
			x = MathUtils::lerp(start[0],end[0],dt);
			y = MathUtils::lerp(start[1],end[1],dt);
			z = MathUtils::lerp(start[1],end[1],dt);

			return Vector3(x,y,z);
		}

		/*static Colour lerp(Colour start, Colour end, float t)
		{
			Colour colour = Colour();
			colour.setR(((float)1-t)*start.getR() + t*end.getR());
			colour.setG(((float)1-t)*start.getG() + t*end.getG());
			colour.setB(((float)1-t)*start.getB() + t*end.getB());
			colour.setA(((float)1-t)*start.getA() + t*end.getA());

			return colour;
		}*/
		static Vector3 catMull(Vector3 point0, Vector3 point1, Vector3 point2, Vector3 point3, float delta)
		{
			Vector3 newPoint = (point1 * 2) + 
								(((-point0)+point2)*delta)+
								(((point0*2) - (point1*5)+(point2*4) - point3)*delta*delta) +
								((-point0 +(point1*3) - (point2*3) + point3)*delta*delta*delta);

			return newPoint;
		}
		/*static Vector3 randomPointOnALine(Vector3 start, Vector3 end)
		{
			float t = RandomUtils::randomf(0.0f,1.0f);

			return start.lerp(end,t);
		}

		static Vector3 randomPointFromAList(std::vector<Vector3> _list)
		{
			unsigned int count = _list.size();

			unsigned int i = RandomUtils::randomi(0,count-1);

			return _list[i];
		}

		static Vector3 randomPointInAPlane(Vector3 position, Vector3 normal, float radius)
		{
			 Vector3 random = Vector3();

			 do
			 {
			  random[0] = (float)2.0 * (float)RandomUtils::randomf((float)0.0,(float)1.0) - (float)1.0;
			  random[1] = (float)2.0 * (float)RandomUtils::randomf((float)0.0,(float)1.0) - (float)1.0;
			  random[2] = (float)2.0 * (float)RandomUtils::randomf((float)0.0,(float)1.0) - (float)1.0;
			  Vector3::cross(random, normal);
			 } while (random.isZero());

			 random.normalize();
			 random *= radius * (float)std::sqrt(RandomUtils::randomf((float)0.0,(float)1.0));
			 random += position;

			 return Vector3(random);
		}

		static Vector3 randomPointInACircleXZ(Vector3 center, float radius)
		{
			float randomRadius = radius*std::sqrt(RandomUtils::randomf((float)0.0,(float)1.0));
			float randomAngle = (float)2*(float)M_PI*RandomUtils::randomf((float)0.0,(float)1.0);
			float x = randomRadius * (float)std::cos(randomAngle);
			float z = randomRadius * (float)std::sin(randomAngle);

			return Vector3(center[0]+x,center[1],center[2]+z);
		}

		static Vector3 randomPointInACircleXY(Vector3 center, float radius)
		{
			float randomRadius = radius*std::sqrt(RandomUtils::randomf((float)0.0,(float)1.0));
			float randomAngle = (float)2*(float)M_PI*RandomUtils::randomf((float)0.0,(float)1.0);
			float x = randomRadius * (float)std::cos(randomAngle);
			float y = randomRadius * (float)std::sin(randomAngle);

			return Vector3(center[0]+x,center[1]+y,center[2]);
		}

		static Vector3 randomPointInACircleYZ(Vector3 center, float radius)
		{
			float randomRadius = radius*std::sqrt(RandomUtils::randomf((float)0.0,(float)1.0));
			float randomAngle = (float)2*(float)M_PI*RandomUtils::randomf((float)0.0,(float)1.0);
			float y = randomRadius * (float)std::cos(randomAngle);
			float z = randomRadius * (float)std::sin(randomAngle);

			return Vector3(center[0],center[1]+y,center[2]+z);
		}*/
	};
}

#endif