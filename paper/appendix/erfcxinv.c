// From https://stackoverflow.com/questions/78832303/accurate-computation-of-the-inverse-of-the-scaled-complementary-error-function-u

#include <math.h>
double erfcx (double);

// In 1B tests: maximum errpr = 3.6145 ulps @ 1.7439789057636135
double erfcxinv (double a)
{
    const double c = 0.5641895835477562869481; // 1/sqrt(pi)
    double e, p, q, r, s, t, u, v, w, x;

    // special cases
    if (isnan (a)) return a + a;
    if (a < 0) return INFINITY / INFINITY; // NaN

    // Halley iteration does not work well around unity: rational approximation
    if ((a > 19.0/32) && (a < 51.0/32)) {
        t = 1 - a;
        p =            -9.9577362854734097e-1;
        p = fma (p, t,  2.5714995288801855e+1);
        p = fma (p, t, -2.0689956608976641e+2);
        p = fma (p, t,  7.4307246888729856e+2);
        p = fma (p, t, -1.3434194015527373e+3);
        p = fma (p, t,  1.2287551470371382e+3);
        p = fma (p, t, -4.9782808803242131e+2);
        p = fma (p, t,  4.9533724768600912e+1);
        q =         t  -2.7622396967308187e+1;
        q = fma (q, t,  2.4716925574847096e+2);
        q = fma (q, t, -1.0276851232159131e+3);
        q = fma (q, t,  2.2643646204587053e+3);
        q = fma (q, t, -2.7339512610053207e+3);
        q = fma (q, t,  1.7120977945481272e+3);
        q = fma (q, t, -4.3537299985712372e+2);
        return fma (p / q, t, t);
    }

    // starting approximation for Halley iteration
    if (a > 2.75) {
        x = -sqrt (log (a / 2));
    } else if (a < 0.375) {
        x = c / a;
    } else {
        t = a - 1;
        x = fma (fma (-0.32374940733, t, 0.87922595958), t, -1.003428026) * t;
    }

    // early exit if starting approximation already sufficiently accurate 
    if ((a < 1.6455e-8) || (a > 0x1.0p51)) return x; 

    // Halley iterations
    for (int i = 0; i < 3; i++) {
        e = erfcx (x);
        q = fma (e, x, -c);
        r = 1 / q;
        s = 0.5 * r;
        t = (a - e) * s;    // -f(x) / f'(x)
        u = fma (e, s, x);  // f"(x) / 2f'(x)
        v = fma (t, u, 1);
        w = 1 / v;
        x = fma (t, w, x);
    }
    return x;
}

double erfcx (double x)
{
    double a, d, e, m, p, q, r, s, t;

    a = fmax (x, 0.0 - x); // NaN preserving absolute value computation

    /* Compute q = (a-4)/(a+4) accurately. [0,INF) -> [-1,1] */
    m = a - 4.0;
    p = a + 4.0;
    r = 1.0 / p;
    q = m * r;
    t = fma (q + 1.0, -4.0, a); 
    e = fma (q, -a, t); 
    q = fma (r, e, q); 

    /* Approximate (1+2*a)*exp(a*a)*erfc(a) as p(q)+1 for q in [-1,1] */
    s = q * q;
    p =             0x1.edcad78fc8044p-31;  //  8.9820305531190140e-10
    t =             0x1.b1548f14735d1p-30;  //  1.5764464777959401e-09
    p = fma (p, s, -0x1.a1ad2e6c4a7a8p-27); // -1.2155985739342269e-08
    t = fma (t, s, -0x1.1985b48f08574p-26); // -1.6386753783877791e-08
    p = fma (p, s,  0x1.c6a8093ac4f83p-24); //  1.0585794011876720e-07
    t = fma (t, s,  0x1.31c2b2b44b731p-24); //  7.1190423171700940e-08
    p = fma (p, s, -0x1.b87373facb29fp-21); // -8.2040389712752056e-07
    t = fma (t, s,  0x1.3fef1358803b7p-22); //  2.9796165315625938e-07
    p = fma (p, s,  0x1.7eec072bb0be3p-18); //  5.7059822144459833e-06
    t = fma (t, s, -0x1.78a680a741c4ap-17); // -1.1225056665965572e-05
    p = fma (p, s, -0x1.9951f39295cf4p-16); // -2.4397380523258482e-05
    t = fma (t, s,  0x1.3be1255ce180bp-13); //  1.5062307184282616e-04
    p = fma (p, s, -0x1.a1df71176b791p-13); // -1.9925728768782324e-04
    t = fma (t, s, -0x1.8d4aaa0099bc8p-11); // -7.5777369791018515e-04
    p = fma (p, s,  0x1.49c673066c831p-8);  //  5.0319701025945277e-03
    t = fma (t, s, -0x1.0962386ea02b7p-6);  // -1.6197733983519948e-02
    p = fma (p, s,  0x1.3079edf465cc3p-5);  //  3.7167515521269866e-02
    t = fma (t, s, -0x1.0fb06dfedc4ccp-4);  // -6.6330365820039094e-02
    p = fma (p, s,  0x1.7fee004e266dfp-4);  //  9.3732834999538536e-02
    t = fma (t, s, -0x1.9ddb23c3e14d2p-4);  // -1.0103906603588378e-01
    p = fma (p, s,  0x1.16ecefcfa4865p-4);  //  6.8097054254651804e-02
    t = fma (t, s,  0x1.f7f5df66fc349p-7);  //  1.5379652102610957e-02
    p = fma (p, q, t);
    p = fma (p, q, -0x1.1df1ad154a27fp-3);  // -1.3962111684056208e-01
    p = fma (p, q,  0x1.dd2c8b74febf6p-3);  //  2.3299511862555250e-01

    /* Divide (1+p) by (1+2*a) ==> exp(a*a)*erfc(a) */
    d = a + 0.5;
    r = 1.0 / d;
    r = r * 0.5;
    q = fma (p, r, r); // q = (p+1)/(1+2*a)
    t = q + q;
    e = (p - q) + fma (t, -a, 1.0); // residual: (p+1)-q*(1+2*a)
    r = fma (e, r, q);

    /* Handle argument of infinity */
    if (a > 0x1.fffffffffffffp1023) r = 0.0;

    /* Handle negative arguments: erfcx(x) = 2*exp(x*x) - erfcx(|x|) */
    if (x < 0.0) {
        s = x * x;
        d = fma (x, x, -s);
        e = exp (s);
        r = e - r;
        r = fma (e, d + d, r); 
        r = r + e;
        if (e > 0x1.fffffffffffffp1023) r = e; // avoid creating NaN
    }
    return r;
}
