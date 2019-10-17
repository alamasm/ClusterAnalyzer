#include <iostream>
#ifndef control_include
#define control_include
#include "Control.h"
#endif

class Interface
{
    public:
    Interface(Control control);
    void start();
    private:
    Control control;
    ifstream input;
    ofstream out;
    int from_file = 0;
    void parse();
};