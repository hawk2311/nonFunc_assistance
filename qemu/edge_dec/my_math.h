// #pragma once 

// int absolut(int x){
//     if(x<0){
//         return (-x);
//     }else{
//         return x;
//     }
// }

// //find the value which has the clostest power of two to the value 
// int find_close(int y){
//     int i = 0;
//     int dis = y; 
//     while(1){
//         int new_dis = absolut(y-(i*i)); 
//         if(new_dis<= dis){ 
//             dis = new_dis; 
//             i++;
//         }else{
//             return i-1;
//         }
//     }
    
// }


// double my_sqrt(int x){
//     int val = find_close(x);
//     int pow_val = val*val;

//     if(x > pow_val){
//         return val + ((double)(x-pow_val)/(2*val));
//     } else {
//         return val - ((double)(pow_val-x)/(2*val));
//     }
// }

#pragma once

int absolut(int x) {
    return x < 0 ? -x : x;
}

// Gibt floor(sqrt(x)) zurück - reicht für Sobel magnitude
int my_sqrt(int x) {
    if (x <= 0) return 0;
    
    int i = 1;
    while (i * i <= x) i++;
    return i - 1;
}