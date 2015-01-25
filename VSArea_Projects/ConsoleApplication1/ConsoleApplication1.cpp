// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

int get_something(const int &num1)
{
  return (num1-1);
}


int _tmain(int argc, _TCHAR* argv[])
{
  auto rez = get_something(10);
  return 0;
}

