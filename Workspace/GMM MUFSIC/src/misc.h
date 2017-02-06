#ifndef _MISC_H_
#define _MISC_H_

#include <sys/time.h>

class cap_tim {
    public:
        void               tic     (void);
        void               toc     (void);
        unsigned long long get_tim (void);
    private:
        struct timeval     tval_bef, tval_aft;
        unsigned long long tim = 0;
};

#ifdef USING_UNITTEST
#define GLOG() "[   GLOG   ] "
#endif

#endif
