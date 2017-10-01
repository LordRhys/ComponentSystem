#ifndef _MATRIX3_H
#define _MATRIX3_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <cmath>
#include <cfloat>

#include "Vector3.h"
//#include "StringUtils.h"
//#include "Logger.h"

  class Matrix3
  {
  private:
    float value[9];
  public:

    Matrix3(void)
    {
      makeZero();
    }

    Matrix3(bool isZero)
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

    Matrix3(const Matrix3& m)
    {
      value[0] = m[0];
      value[1] = m[1];
      value[2] = m[2];
      value[3] = m[3];
      value[4] = m[4];
      value[5] = m[5];
      value[6] = m[6];
      value[7] = m[7];
      value[8] = m[8];
    }

    Matrix3(float a, float b, float c, float d, float e, float f, float g, float h, float i)
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
    }

    Matrix3(const float _value[9])
    {
      value[0] = _value[0];
      value[1] = _value[1];
      value[2] = _value[2];
      value[3] = _value[3];
      value[4] = _value[4];
      value[5] = _value[5];
      value[6] = _value[6];
      value[7] = _value[7];
      value[8] = _value[8];
    }

    Matrix3(const Vector3& u, const Vector3& v,const Vector3& w)
    {
      value[0] = u[0];
      value[1] = u[1];
      value[2] = u[2];
      value[3] = v[0];
      value[4] = v[1];
      value[5] = v[2];
      value[6] = w[0];
      value[7] = w[1];
      value[8] = w[2];
    }

    Matrix3(const Vector3& axis, float angle)
    {
      makeRotation(axis, angle);
    }

    Matrix3& operator= (const Matrix3& m)
    {
      value[0] = m[0];
      value[1] = m[1];
      value[2] = m[2];
      value[3] = m[3];
      value[4] = m[4];
      value[5] = m[5];
      value[6] = m[6];
      value[7] = m[7];
      value[8] = m[8];
      return *this;
    }

    Matrix3& makeZero()
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
      return *this;
    }

    Matrix3& makeIdentity()
    {
      value[0] = (float)1;
      value[1] = (float)0;
      value[2] = (float)0;
      value[3] = (float)0;
      value[4] = (float)1;
      value[5] = (float)0;
      value[6] = (float)0;
      value[7] = (float)0;
      value[8] = (float)1;
      return *this;
    }

    float getValue(int i)
    {
      return value[i];
    }

    Matrix3& makeRotation(const Vector3& axis,float angle)
    {
      float cs = std::cos(angle);
      float sn = std::sin(angle);
      float oneMinusCos = ((float)1) - cs;
      float x2 = axis[0]*axis[0];
      float y2 = axis[1]*axis[1];
      float z2 = axis[2]*axis[2];
      float xym = axis[0]*axis[1]*oneMinusCos;
      float xzm = axis[0]*axis[2]*oneMinusCos;
      float yzm = axis[1]*axis[2]*oneMinusCos;
      float xSin = axis[0]*sn;
      float ySin = axis[1]*sn;
      float zSin = axis[2]*sn;
    
      value[0] = x2*oneMinusCos + cs;
      value[1] = xym - zSin;
      value[2] = xzm + ySin;
      value[3] = xym + zSin;
      value[4] = y2*oneMinusCos + cs;
      value[5] = yzm - xSin;
      value[6] = xzm - ySin;
      value[7] = yzm + xSin;
      value[8] = z2*oneMinusCos + cs;
      return *this;
    }

    Matrix3 operator+ (const Matrix3& m)
    {
      Matrix3 result;
    
      result[0] = value[0] + m[0];
      result[1] = value[1] + m[1];
      result[2] = value[2] + m[2];
      result[3] = value[3] + m[3];
      result[4] = value[4] + m[4];
      result[5] = value[5] + m[5];
      result[6] = value[6] + m[6];
      result[7] = value[7] + m[7];
      result[8] = value[8] + m[8];
    
      return Matrix3(result[0],result[1],result[2],result[3],result[4],result[5],result[6],result[7],result[8]);
    }

    Matrix3 operator- (const Matrix3& m)
    {
      Matrix3 result;
    
      result[0] = value[0] - m[0];
      result[1] = value[1] - m[1];
      result[2] = value[2] - m[2];
      result[3] = value[3] - m[3];
      result[4] = value[4] - m[4];
      result[5] = value[5] - m[5];
      result[6] = value[6] - m[6];
      result[7] = value[7] - m[7];
      result[8] = value[8] - m[8];
    
      return Matrix3(result[0],result[1],result[2],result[3],result[4],result[5],result[6],result[7],result[8]);
    }

    Matrix3 operator* (float scalar)
    {
      Matrix3 result;
    
      result[1] = value[0] * scalar;
      result[2] = value[1] * scalar;
      result[3] = value[2] * scalar;
      result[4] = value[3] * scalar;
      result[5] = value[4] * scalar;
      result[6] = value[5] * scalar;
      result[7] = value[6] * scalar;
      result[8] = value[7] * scalar;
      result[9] = value[8] * scalar;
    
      return Matrix3(result[0],result[1],result[2],result[3],result[4],result[5],result[6],result[7],result[8]);
    }

    Matrix3 operator/ (float scalar) const
    {
      Matrix3 result;
    
      result[1] = value[0] / scalar;
      result[2] = value[1] / scalar;
      result[3] = value[2] / scalar;
      result[4] = value[3] / scalar;
      result[5] = value[4] / scalar;
      result[6] = value[5] / scalar;
      result[7] = value[6] / scalar;
      result[8] = value[7] / scalar;
      result[9] = value[8] / scalar;
    
      return Matrix3(result[0],result[1],result[2],result[3],result[4],result[5],result[6],result[7],result[8]);
    }

    Matrix3 operator- () const
    {
      Matrix3 result;
    
      result[1] = -value[0];
      result[2] = -value[1];
      result[3] = -value[2];
      result[4] = -value[3];
      result[5] = -value[4];
      result[6] = -value[5];
      result[7] = -value[6];
      result[8] = -value[7];
      result[9] = -value[8];
    
      return Matrix3(result[0],result[1],result[2],result[3],result[4],result[5],result[6],result[7],result[8]);
    }

    void operator+= (Matrix3& m)
    {
      value[0] += m.getValue(0);
      value[1] += m.getValue(1);
      value[2] += m.getValue(2);
      value[3] += m.getValue(3);
      value[4] += m.getValue(4);
      value[5] += m.getValue(5);
      value[6] += m.getValue(6);
      value[7] += m.getValue(7);
      value[8] += m.getValue(8);
    }

    void operator-= (Matrix3& m)
    {
      value[0] -= m.getValue(0);
      value[1] -= m.getValue(1);
      value[2] -= m.getValue(2);
      value[3] -= m.getValue(3);
      value[4] -= m.getValue(4);
      value[5] -= m.getValue(5);
      value[6] -= m.getValue(6);
      value[7] -= m.getValue(7);
      value[8] -= m.getValue(8);
    }

    void operator *= (float scalar)
    {
      value[0] *= scalar;
      value[1] *= scalar;
      value[2] *= scalar;
      value[3] *= scalar;
      value[4] *= scalar;
      value[5] *= scalar;
      value[6] *= scalar;
      value[7] *= scalar;
      value[8] *= scalar;
    }

    void operator/= (float scalar)
    {
      value[0] /= scalar;
      value[1] /= scalar;
      value[2] /= scalar;
      value[3] /= scalar;
      value[4] /= scalar;
      value[5] /= scalar;
      value[6] /= scalar;
      value[7] /= scalar;
      value[8] /= scalar;
    }


    Matrix3 transpose()
    {
      return Matrix3(value[0],
        value[3],
        value[6],
        value[1],
        value[4],
        value[7],
        value[2],
        value[5],
        value[8]);

    }

    void operator *=(Matrix3& m)
    {
      value[0] = value[0]*m.getValue(0) + value[1]*m.getValue(3) + value[2]*m.getValue(6);
      value[1] = value[0]*m.getValue(1) + value[1]*m.getValue(4) + value[2]*m.getValue(7);
      value[2] = value[0]*m.getValue(2) + value[1]*m.getValue(5) + value[2]*m.getValue(8);
      value[3] = value[3]*m.getValue(0) + value[4]*m.getValue(3) + value[5]*m.getValue(6);
      value[4] = value[3]*m.getValue(1) + value[4]*m.getValue(4) + value[5]*m.getValue(7);
      value[5] = value[3]*m.getValue(2) + value[4]*m.getValue(5) + value[5]*m.getValue(8);
      value[6] = value[6]*m.getValue(0) + value[7]*m.getValue(3) + value[8]*m.getValue(6);
      value[7] = value[6]*m.getValue(1) + value[7]*m.getValue(4) + value[8]*m.getValue(7);
      value[8] = value[6]*m.getValue(2) + value[7]*m.getValue(5) + value[8]*m.getValue(8);
    }

    Matrix3 operator* (Matrix3& m)
    {
      return Matrix3(
      value[0]*m.getValue(0) + value[1]*m.getValue(3) + value[2]*m.getValue(6),
      value[0]*m.getValue(1) + value[1]*m.getValue(4) + value[2]*m.getValue(7),
      value[0]*m.getValue(2) + value[1]*m.getValue(5) + value[2]*m.getValue(8),
      value[3]*m.getValue(0) + value[4]*m.getValue(3) + value[5]*m.getValue(6),
      value[3]*m.getValue(1) + value[4]*m.getValue(4) + value[5]*m.getValue(7),
      value[3]*m.getValue(2) + value[4]*m.getValue(5) + value[5]*m.getValue(8),
      value[6]*m.getValue(0) + value[7]*m.getValue(3) + value[8]*m.getValue(6),
      value[6]*m.getValue(1) + value[7]*m.getValue(4) + value[8]*m.getValue(7),
      value[6]*m.getValue(2) + value[7]*m.getValue(5) + value[8]*m.getValue(8));
    }

    Matrix3 inverse() 
    {
      Matrix3 inverse;

      float determinant = this->determinant();

      if (std::fabs(determinant) > FLT_EPSILON)
      {
        Matrix3 adjoinfloat = this->adjoint();
        inverse = adjoinfloat * (1/determinant);

        return inverse;
      }

      return inverse.makeZero();
    }

    Matrix3 adjoint()
    {
      return Matrix3(value[4]*value[8] - value[5]*value[7],
        value[2]*value[7] - value[1]*value[8],
        value[1]*value[5] - value[2]*value[4],
        value[5]*value[6] - value[3]*value[8],
        value[0]*value[8] - value[2]*value[6],
        value[2]*value[3] - value[0]*value[5],
        value[3]*value[7] - value[4]*value[6],
        value[1]*value[6] - value[0]*value[7],
        value[0]*value[4] - value[1]*value[3]);
    
    }

    float determinant() 
    {
      float co00 = value[4]*value[8] - value[5]*value[7];
      float co10 = value[5]*value[6] - value[3]*value[8];
      float co20 = value[3]*value[7] - value[4]*value[6];

      return value[0]*co00 + value[1]*co10 + value[2]*co20;
    }

    void ExtractAxisAngle(Vector3& axis, float& angle)
    {
      float trace = value[0] + value[4] + value[8];
      float cs = ((float)0.5)*(trace - (float)1);
      angle = std::acos(cs);

      if (angle > (float)0)
      {
        if (angle < M_PI)
        {
          axis[0] = value[7] - value[5];
          axis[1] = value[2] - value[6];
          axis[2] = value[3] - value[1];
          axis.normalize();
        }
        else
        {
          float halfInverse;
          if (value[0] >= value[4])
          {
            if (value[0] >= value[8])
            {
              axis[0] = ((float)0.5)*std::sqrt((float)1
                + value[0] - value[4] - value[8]);
              halfInverse = ((float)0.5)/axis[0];
              axis[1] = halfInverse*value[1];
              axis[2] = halfInverse*value[2];
            }
            else
            {
              axis[2] = ((float)0.5)*std::sqrt((float)1
                + value[8] - value[0] - value[4]);
              halfInverse = ((float)0.5)/axis[2];
              axis[0] = halfInverse*value[2];
              axis[1] = halfInverse*value[5];
            }
          }
          else
          {
            if (value[4] >= value[8])
            {
              axis[1] = ((float)0.5)*std::sqrt((float)1
                + value[4] - value[0] - value[8]);
              halfInverse  = ((float)0.5)/axis[1];
              axis[0] = halfInverse*value[1];
              axis[2] = halfInverse*value[5];
            }
            else
            {
              axis[2] = ((float)0.5)*std::sqrt((float)1
                + value[8] - value[0] - value[4]);
              halfInverse = ((float)0.5)/axis[2];
              axis[0] = halfInverse*value[2];
              axis[1] = halfInverse*value[5];
            }
          }
        }
      }
      else
      {
        axis[0] = (float)1;
        axis[1] = (float)0;
        axis[2] = (float)0;
      }
    }

    void Orthonormalize()
    {
      float invLength = 1/std::sqrt(value[0]*value[0] +value[3]*value[3] + value[6]*value[6]);

      value[0] *= invLength;
      value[3] *= invLength;
      value[6] *= invLength;

      float dot0 = value[0]*value[1] + value[3]*value[4] + value[6]*value[7];

      value[1] -= dot0*value[0];
      value[4] -= dot0*value[3];
      value[7] -= dot0*value[6];

      invLength = 1/std::sqrt(value[1]*value[1] + value[4]*value[4] + value[7]*value[7]);

      value[1] *= invLength;
      value[4] *= invLength;
      value[7] *= invLength;

      float dot1 = value[1]*value[2] + value[4]*value[5] + value[7]*value[8];

      dot0 = value[0]*value[2] + value[3]*value[5] + value[6]*value[8];

      value[2] -= dot0*value[0] + dot1*value[1];
      value[5] -= dot0*value[3] + dot1*value[4];
      value[8] -= dot0*value[6] + dot1*value[7];

      invLength = 1/std::sqrt(value[2]*value[2] + value[5]*value[5] + value[8]*value[8]);

      value[2] *= invLength;
      value[5] *= invLength;
      value[8] *= invLength;
    }

    inline Vector3 operator* (const Vector3& v)
    {
      Vector3 result;
    
      result[0] = v[0]*value[0] + v[1]*value[1] + v[2]*value[2];
      result[1] = v[0]*value[3] + v[1]*value[4] + v[2]*value[5];
      result[2] = v[0]*value[6] + v[1]*value[7] + v[2]*value[8];

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
      Logger::getInstance()->write(StringUtils::format("Matrix3[ %f, %f, %f ]",value[0],value[1],value[2]));
      Logger::getInstance()->write(StringUtils::format("Matrix3[ %f, %f, %f ]",value[3],value[4],value[5]));
      Logger::getInstance()->write(StringUtils::format("Matrix3[ %f, %f, %f ]",value[6],value[7],value[8]));
    }*/

  };


#endif