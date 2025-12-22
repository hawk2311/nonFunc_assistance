#include <stdio.h>
#include <iostream>
#include "CycleTimer.h"
#include "ClockTimer.h"

int main() {
    ClockTimer t;
    t.start();
    printf("Hello world!\n");
    std::cout<<t.stop() << std::endl;
    return 0;
}
