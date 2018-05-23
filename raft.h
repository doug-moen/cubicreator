/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#ifndef RAFT_H
#define RAFT_H

#include "sliceDataStorage.h"

void generateRaft(SliceDataStorage& storage, int distance);
void generateRaft(SliceDataStorage& storage, int distance, int dashStride, int dashSize);

#endif//RAFT_H
