#include "Component.h"
#include "GameObject.h"

Component::Component(void)
{
  this->m_name = "Component";
  this->m_isEnabled = true;
  this->m_gameObject = NULL;
  this->Awake();
  HRESULT l_result = CoCreateGuid(&m_guid);
}

Component::Component(std::string l_name)
{
  this->m_name = l_name;
  this->m_isEnabled = true;
  this->m_gameObject = NULL;
  this->Awake();
  HRESULT l_result = CoCreateGuid(&m_guid);
}

Component::~Component(void)
{

}

void Component::Destroy()
{
}


void Component::Awake()
{

}

void Component::Start()
{

}

void Component::FixedUpdate()
{

}

void Component::Update()
{
}

void Component::LateUpdate()
{

}
