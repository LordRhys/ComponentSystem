#include "ComponentManager.h"
#include "GameObject.h"

  ComponentManager::ComponentManager(void)
  {
    m_componentsSize = 2000;
    m_currentComponentCount = 0;
    ResizeComponents();
  }
  
  ComponentManager::~ComponentManager(void)
  {
  }


  void ComponentManager::Clear()
  {
    for (component_vector_itr itr = this->m_components.begin(); itr != m_components.end(); ++itr)
    {
      delete (*itr);
    }

    m_components.clear();
    m_renderComponents.clear();
  }

  void ComponentManager::Create()
  {

  }

  void ComponentManager::Destroy()
  {
    Clear();
  }

  Component* ComponentManager::AddComponent(Component* l_component)
  {
    m_currentComponentCount++;

    if (m_currentComponentCount == m_currentComponentsSize - 1)
    {
      ResizeComponents();
    }

    if (Renderer* rendererable = dynamic_cast<Renderer*>(l_component))
    {
      this->m_renderComponents.push_back(rendererable);
    }

    this->m_components.push_back(l_component);
    int l_position = m_components.size() - 1;
    return m_components[l_position];
  }

  bool ComponentManager::RemoveComponent(component_vector_const_itr l_component_itr)
  {
    this->m_components.erase(l_component_itr);
    return true;
  }

  void ComponentManager::RemoveComponentWithGameObject(GameObject* l_objects)
  {
    component_vector components_to_remove;
    for (component_vector_const_itr itr = this->m_components.begin(); itr != this->m_components.end(); ++itr)
    {
      if ((*itr)->m_gameObject == l_objects)
      {
        components_to_remove.push_back((*itr));
      }
    }

    for (component_vector_const_itr itr = components_to_remove.begin(); itr != components_to_remove.end(); ++itr)
    {
      RemoveComponent(itr);
    }
  }

  void ComponentManager::FixedUpdate()
  {
    for (component_vector_itr itr = this->m_components.begin(); itr != this->m_components.end(); ++itr)
    {
      if ((*itr)->m_isEnabled)
      {
        (*itr)->FixedUpdate();
      }
    }
  }

  void ComponentManager::Update()
  {
    for (component_vector_itr itr = this->m_components.begin(); itr != this->m_components.end(); ++itr)
    {
      if ((*itr)->m_isEnabled)
      {
        (*itr)->Update();
      }
    }
  }

  void ComponentManager::LateUpdate()
  {
    for (component_vector_itr itr = this->m_components.begin(); itr != this->m_components.end(); ++itr)
    {
      if ((*itr)->m_isEnabled)
      {
        (*itr)->LateUpdate();
      }
    }
  }

  void ComponentManager::Render()
  {
    for (render_component_vector_itr itr = this->m_renderComponents.begin(); itr != this->m_renderComponents.end(); ++itr)
    {
      if ((*itr)->m_isEnabled)
      {
        (*itr)->Render();
      }
    }
  }

  void ComponentManager::ResizeComponents()
  {
    m_components.resize(m_currentComponentsSize + m_componentsSize);
    m_currentComponentsSize += m_componentsSize;
  }
