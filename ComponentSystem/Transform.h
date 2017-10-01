#pragma once

#include "Component.h"

#include "Matrix4.h"
#include "Matrix3.h"
#include "Vector3.h"
#include "Quaternion.h"

class GameObject;
class Transform : public Component
{
public:
  Quaternion m_rotation;
  Vector3 m_up;
  Vector3 m_right;
  Vector3 m_direction;
  Vector3 m_position;

  Matrix3 m_scale;

  Matrix4 m_localMatrix;
  Matrix4 m_worldMatrix;

  Matrix4 m_mv;
  Matrix4 m_mvp;

  Transform(void);
  ~Transform(void);
};
