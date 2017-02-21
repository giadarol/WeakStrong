#include "cmpx.h"

cmpx makecmpx(double r, double i){
    cmpx res;
    res.r = r;
    res.i = i;
    return res;
}

cmpx cadd(cmpx a, cmpx b){
    cmpx res;
    res.r = a.r + b.r;
    res.i = a.i + b.i;
    return res;
}

cmpx csub(cmpx a, cmpx b){
    cmpx res;
    res.r = a.r - b.r;
    res.i = a.i - b.i;
    return res;
}

cmpx cmul(cmpx a, cmpx b){
    cmpx res;
    res.r = a.r*b.r-a.i*b.i;
    res.i = a.i*b.r+a.r*b.i;
    return res;
}

cmpx cdiv(cmpx a, cmpx b){
    double den;
    cmpx res;
    den=b.r*b.r+b.i*b.i;
    res.r = (a.r*b.r+a.i*b.i)/den;
    res.i = (a.i*b.r-a.r*b.i)/den;
    return res;
}

cmpx expc(cmpx a){
    double mod;
    cmpx res;
    mod = exp(a.r);
    res.r = mod*cos(a.i);
    res.i = mod*sin(a.i);
    return res;
}

cmpx cscale(cmpx a, double f){
    cmpx res;
    res.r = f*a.r;
    res.i = f*a.i;
    return res;
}
