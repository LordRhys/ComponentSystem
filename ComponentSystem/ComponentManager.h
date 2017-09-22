#pragma once

#include "Component.h"
#include "Renderer.h"
#include "GameObject.h"

class ComponentManager
{
public:
  typedef std::vector<Component*> component_vector;
  typedef component_vector::iterator component_vector_itr;
  typedef component_vector::const_iterator component_vector_const_itr;

  typedef std::vector<Renderer*> render_component_vector;
  typedef render_component_vector::iterator render_component_vector_itr;
  typedef render_component_vector::const_iterator render_component_vector_const_itr;

  component_vector m_components;
  render_component_vector m_renderComponents;

  int m_componentsSize;
  int m_currentComponentsSize;
  int m_currentComponentCount;

  ComponentManager(void);
  ~ComponentManager(void);

  void Clear();
  void Create();
  void Destroy();

  Component* AddComponent(Component* l_component);
  bool RemoveComponent(component_vector_const_itr l_component_itr);
  void RemoveComponentWithGameObject(GameObject* l_objects);

  void FixedUpdate();
  void Update();
  void LateUpdate();

  void Render();

  void ResizeComponents();

  // Singleton Implementation
  static ComponentManager* s_instance;

  static inline ComponentManager* getInstance(void)
  {
    if (s_instance == NULL)
    {
      s_instance = new ComponentManager();
    }

    return s_instance;
  }

  static inline bool exists(void)
  {
    return s_instance != 0;
  }
};