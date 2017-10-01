#pragma once

#include <vector>

template<typename _InIt, typename _Pr> inline
std::vector<_InIt> find_all(_InIt _First, _InIt _Last, _Pr _Pred)
{
  std::vector<_InIt> l_result;

  for (; _First != _Last; ++_First)
  {
    if (_Pred(*_First))
    {
      l_result.push_back(_First);
    }
  }

  return l_result;
}