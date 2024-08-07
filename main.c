#include <stdio.h>

#include "orientation_processing.h"

int main() {
    long int timestamp = 1581789429668;
    orientation_t or;
    measurement_t m;
    m.acc[0] = -80000;
    m.acc[1] = 60000;
    m.acc[2] = 255000;
    m.temp = 1775;
    for (int i = 0; i < GAUSS_NEWTON_WINDOW_SIZE * 6; i++) {
        int retval = process(&m, 5374047, timestamp + i, &or);
        if (retval == PROCESS_RESULT_NEW_ORIENTATION) {
            printf("%lu: ax = %6.4f ay = %6.4f\n", or.timestamp, or.angles.ax, or.angles.ay);
        }
    }
    return 0;
}