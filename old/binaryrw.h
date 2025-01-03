// binary read/write template
#ifndef BINARYRW_H
#define BINARYRW_H
#include <iostream>

using namespace std;


template<typename T>
ostream& bwrite(ofstream& stream, const T& value){
    return stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template<typename T>
istream& bread(istream& stream, T& value){
    return stream.read(reinterpret_cast<char*>(&value), sizeof(T));
}


#endif