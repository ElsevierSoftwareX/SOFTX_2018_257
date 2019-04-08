/* TendonMech: An Open Source High Performance Code to Compute the Mechanical behavior of Tendon Fascicles
 *
 */

#ifndef _MYTIMER_H_
#define _MYTIMER_H_ 1

#include <sys/time.h>

double my_gettime()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}

class MyTimer
{
private:
  double t0, t1;

public:
  void start()
  {
    t0 = my_gettime();
  }

  void stop(string msg)
  {
    t1 = my_gettime();
#if defined(_PRINT_MYTIMER_)
    cout << msg << ": time = " << t1-t0 << " s" << endl;
#endif
  }
};

#endif
