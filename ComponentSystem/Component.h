#pragma once

#include <string>

class GameObject;
class Component
{
public:
  std::string m_name;

  GameObject* m_gameObject;

  bool m_isEnabled;

	Component(void);
  Component(std::string l_name);
	virtual ~Component(void);

  virtual void Destroy();
  virtual void Awake();
  virtual void Start();
  virtual void FixedUpdate();
  virtual void Update();
  virtual void LateUpdate();

};