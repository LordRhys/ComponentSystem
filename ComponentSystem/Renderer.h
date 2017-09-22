#pragma once

#include "Component.h"

class Renderer : public Component
{
public:
  Renderer(void);
  ~Renderer(void);

  void Render();
};
