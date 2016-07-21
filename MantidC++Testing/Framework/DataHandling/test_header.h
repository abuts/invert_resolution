#ifndef H_TEST_MAIN_DATAHANDLING_
#define H_TEST_MAIN_DATAHANDLING_

#include "MantidKernel/System.h"
#include "MantidAPI/FrameworkManager.h"

#include <cxxtest/TestSuite.h>

class test_sqw: public CxxTest::TestSuite
{
public:
      void testTMain(void)
      {
        Mantid::API::FrameworkManager::Instance();
        TS_WARN( "Test suite invoked" );
      }
};
#endif