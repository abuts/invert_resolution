// testSP.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "boost/shared_ptr.hpp"
#include <vector>

class reporter
{
  int index;
public:
  reporter():index(0){};
  reporter(int i0):index(i0){};
  void setInd(int ic)
  {
    index=ic;
  }
  int report()
  {
    return index;
  }
};

class holder
{
  size_t nObj;
  std::vector<boost::shared_ptr<reporter> > objects;
public: 

  holder(int nobj)
  {
   objects.resize(nobj);
   nObj = nobj;
   for(size_t i=0;i<10;i++)
   {
     objects[i]=boost::shared_ptr<reporter>(new reporter(int(i)));
   }
  }
  boost::shared_ptr<reporter> getObject(size_t index)
  {
    if ((index==0) || (index==4))
    {
      throw "InvalidIndex";
    }
    else
    {
      return objects[index];
    }
  }
};

int _tmain(int argc, _TCHAR* argv[])
{

  holder hobj(10);

  std::vector<int> result(10,-1);
  for(size_t i=0;i<10;i++)
  {
    boost::shared_ptr<reporter> spObj;
    try
    {
      spObj = hobj.getObject(i);
    }catch(...)
    {

    }
    if (spObj)
      result[i]=spObj->report();
  }

  return 0;
}

