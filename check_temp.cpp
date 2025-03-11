#include <iostream>  //for cin and cout
#include <math.h>  // for pow
#include <stdlib.h>  //for div(q,n).rem(quot),abs(int n)
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include "tensor_type.h"

int main(){

Mat_1_int TEMP = {0,1,3,12312,4,5}; 

cout<<TEMP.size()<<endl;

for(int i=0;i<TEMP.size();i++){
cout<<i<<" : "<<TEMP[i]<<endl;
}

return 0;
}
