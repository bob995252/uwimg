#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

void change_value(int *f, int SizeOfArray){
    
    for (int i=0; i<SizeOfArray; i++){
        printf("f = %d\n", f[i]);
    }
    for (int i = 0; i < SizeOfArray; ++i) f[i] = i * i;
    for (int i = 0; i < SizeOfArray; ++i) printf("(in function) f%d = %d\n", i, f[i]);
}

int main(){
    /*
    int flag[3], SizeOfArray;
    SizeOfArray = sizeof(flag)/sizeof(*flag);
    for (int i = 0; i<SizeOfArray; ++i) flag[i] = i;
    change_value(flag, SizeOfArray);
    printf("===== After change =====\n");
    for (int i = 0; i<SizeOfArray; ++i) printf("(in main) f%d = %d\n", i, flag[i]);*/

    return 0;
}