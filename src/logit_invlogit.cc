# include "cmath"

/* Fonction logit */
double logit (double x) {
    double Result=log(x)-log(1-x);
    return Result;
}
/* Fonction invlogit */
double invlogit (double x) {
    double Result=1/(1+exp(-x));
    return Result;
}
