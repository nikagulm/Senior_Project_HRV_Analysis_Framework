#ifndef ABLATIONVARIANT_H
#define ABLATIONVARIANT_H

//which stage of the published algorithm is replaced
enum class AblationVariant {
    FixedThreshold,  //threshold set in the learning phase and never updated
    NoSearchBack,    //missed beats are never looked for
    MovingAverage    //moving average instead of the band-pass filter
};

#endif