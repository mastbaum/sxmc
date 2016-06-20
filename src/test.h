#ifndef __TEST_H__
#define __TEST_H__

/**
 * \file test.h
 *
 * Convergence/goodness-of-fit tests performed on a likelihood space.
*/

#include <string>

#include <sxmc/interval.h>

class LikelihoodSpace;

/** Test results. */
typedef enum { TEST_OK, TEST_FAIL, TEST_UNKNOWN } TestResult;


/**
 * Base class for tests on the likelihood space.
*/
class Test {
  public:
    /** Constructor. */
    Test() {}
 
    /** Destructor. */
    virtual ~Test() {}
 
    /** Perform the test. */
    virtual TestResult operator()(const LikelihoodSpace* ls) = 0;

    /** Get the status as a string. */
    std::string status(const TestResult& result) {
      if (result == TEST_OK) {
        return "OK";
      }
      else if (result == TEST_FAIL) {
        return "FAIL";
      }
      return "UNKNOWN";
    }
};

#endif  // __TEST_H__

