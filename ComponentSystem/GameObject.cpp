#include "GameObject.h"

  GameObject::GameObject(void)
  {
    this->m_name = "GameObject";
    this->m_tag = "";
    this->m_transform = NULL;
    this->m_parent = NULL;
  }
  
  GameObject::GameObject(std::string l_name)
  {
    this->m_name = l_name;
    this->m_tag = "";
    this->m_transform = NULL;
    this->m_parent = NULL;
  }

  GameObject::~GameObject(void)
  {
  }


  void GameObject::Create()
  {
    Transform* l_transform = new Transform();
    l_transform->m_gameObject = this;
    
  }

  void GameObject::Destroy()
  {

  }

  void GameObject::AddComponent(Component* l_component)
  {

  }

  Component* GameObject::FindComponentByName(std::string l_name)
  {

  }

  void GameObject::AddChild(GameObject * l_gameObject)
  {
  }

  GameObject* GameObject::FindChildByName(std::string l_name)
  {

  }

  std::vector<GameObject*> GameObject::FindAllChildrenByName(std::string l_name)
  {

  }
