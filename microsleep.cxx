#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>

void microsleep(int sec, int usec)
{
    struct timeval tv;

    tv.tv_sec = sec;
    tv.tv_usec = usec;

    int retval = select(1, NULL, NULL, NULL, &tv);
    /* Don't rely on the value of tv now! */
}

