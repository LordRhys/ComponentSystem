#pragma once

#include <string>
#include <algorithm>
#include <typeinfo>
#include <vector>

#include "Component.h"

class GameObject
{
public:
  std::string m_name;
  std::string m_tag;
  GameObject* m_parent;

  typedef std::vector<Component*> component_vector;
  typedef component_vector::iterator component_vector_itr;
  typedef component_vector::const_iterator component_vector_const_itr;

  typedef std::vector<GameObject*> gameObject_vector;
  typedef gameObject_vector::iterator gameObject_vector_itr;
  typedef gameObject_vector::const_iterator gameObject_vector_const_itr;

  component_vector m_components;
  gameObject_vector m_children;

  GameObject(void);
  GameObject(std::string l_name);
	~GameObject(void);

  void Create();
  void Destroy();

  void AddComponent(Component* l_component);
  Component* FindComponentByName(std::string l_name);

  template <typename T>
  T* FindComponentByType()
  {
    for (component_vector_itr itr = m_components.begin(); itr != m_components.end(); ++itr)
    {
      if (T* l_type = dynamic_cast<T*>(itr))
      {
        return l_type;
      }
    }

    return NULL;
  };

  void AddChild(GameObject* l_gameObject);

  GameObject* FindChildByName(std::string l_name);
  std::vector<GameObject*> FindAllChildrenByName(std::string l_name);

};