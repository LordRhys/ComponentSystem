#pragma once

#include <string>

#include <objbase.h>
#include <guiddef.h>

class GameObject;
class Component
{
public:
  std::string m_name;

  GUID m_guid;

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

struct ComponentComparer
{
  std::string n_name;

  ComponentComparer(std::string l_name) 
    : n_name(l_name)
  {
  }

  bool operator()(Component* l_object)
  {
    return (l_object->m_name == n_name ? true : false);
  }
};