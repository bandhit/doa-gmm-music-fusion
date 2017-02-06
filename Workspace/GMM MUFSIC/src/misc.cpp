#include <misc.h>
#include <iostream>

void cap_tim::tic (void) {
    tim = 0;
    gettimeofday(&tval_bef, NULL);
}

void cap_tim::toc (void) {
    gettimeofday(&tval_aft, NULL);
    tim = ((tval_aft.tv_sec * 1000000 + tval_aft.tv_usec) - (tval_bef.tv_sec * 1000000 + tval_bef.tv_usec));
}

unsigned long long cap_tim::get_tim (void) {
    return tim;
}
