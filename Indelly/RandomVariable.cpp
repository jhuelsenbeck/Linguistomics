#include <cmath>
#include <iostream>
#include <limits>
#include "Msg.hpp"
#include "RandomVariable.hpp"



RandomVariable::RandomVariable(void) {

    uint32_t seed = (uint32_t)time(NULL);
    initialize(seed);
}

RandomVariable::RandomVariable(uint32_t seed) {

    initialize(seed);
}

double RandomVariable::beta(double a, double b) {

    return ( exp(lnGamma(a) + lnGamma(b) - lnGamma(a + b)) );
}

double RandomVariable::betaRv(double a, double b) {

    double z0 = rndGamma( a );
    double z1 = rndGamma( b );
    double sum = z0 + z1;
    double x = z0 / sum;
    return x;
}

double RandomVariable::betaPdf(double a, double b, double x) {

    double pdf;
    if ( x < 0.0 || 1.0 < x )
        pdf = 0.0;
    else
        pdf = pow(x, (a - 1.0)) * pow((1.0 - x), (b - 1.0)) / beta(a, b);
    return pdf;
}

double RandomVariable::lnBetaPdf(double a, double b, double x) {

    return ( (lnGamma(a + b) - lnGamma(a) - lnGamma(b)) + (a - 1.0) * log(x) + (b - 1.0) * log(1.0 - x) );
}

double RandomVariable::betaCdf(double a, double b, double x) {

    double cdf;
    if ( x <= 0.0 )
        cdf = 0.0;
    else if ( x <= 1.0 )
        cdf = incompleteBeta(a, b, x);
    else
        cdf = 1.0;
    return cdf;
}

double RandomVariable::betaQuantile(double alpha, double beta, double x) {

    int        i, nswitches;
    double    curPos, curFraction, increment;
    
    i = nswitches = 0;
    curPos = 0.5;
    bool stopIter = false;
    bool directionUp;
    increment = 0.25;
    curFraction = incompleteBeta(alpha, beta, curPos);
    if (curFraction > x)
    {
        directionUp = false;
    }
    else
    {
        directionUp = true;
    }
    
    while (!stopIter)
    {
        curFraction = incompleteBeta(alpha, beta, curPos);
        if (curFraction > x && directionUp == false)
        {
            /* continue going down */
            while (curPos - increment <= 0.0)
            {
                increment /= 2.0;
            }
            curPos -= increment;
        }
        else if (curFraction > x && directionUp == true)
        {
            /* switch directions, and go down */
            nswitches++;
            directionUp = false;
            while (curPos - increment <= 0.0)
            {
                increment /= 2.0;
            }
            increment /= 2.0;
            curPos -= increment;
        }
        else if (curFraction < x && directionUp == true)
        {
            /* continue going up */
            while (curPos + increment >= 1.0)
            {
                increment /= 2.0;
            }
            curPos += increment;
        }
        else if (curFraction < x && directionUp == false)
        {
            /* switch directions, and go up */
            nswitches++;
            directionUp = true;
            while (curPos + increment >= 1.0)
            {
                increment /= 2.0;
            }
            increment /= 2.0;
            curPos += increment;
        }
        else
        {
            stopIter = true;
        }
        if (i > 1000 || nswitches > 20)
            stopIter = true;
        i++;
    }
    
    return (curPos);
}

double RandomVariable::chiSquareQuantile(double prob, double v) {

    double         e = 0.5e-6, aa = 0.6931471805, p = prob,
                    a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0,
                    x = 0.0, b = 0.0;

    if (p < 0.000002 || p > 0.999998 || v <= 0.0)
        return (-1.0);
    double g = lnGamma(v/2.0);
    double xx = v/2.0;
    double c = xx - 1.0;
    double ch = pow((p*xx*exp(g+xx*aa)), 1.0/xx);
    if (v >= -1.24*log(p))
        goto l1;
    if (ch-e < 0)
        return (ch);
    goto l4;
    l1:
        if (v > 0.32)
            goto l3;
        ch = 0.4;
        a = log(1.0-p);
    l2:
        q = ch;
        p1 = 1.0+ch*(4.67+ch);
        p2 = ch*(6.73+ch*(6.66+ch));
        t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
        ch -= (1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
        if (fabs(q/ch-1.0)-0.01 <= 0.0)
            goto l4;
        else
            goto l2;
    l3:
        x = pointNormal (p);
        p1 = 0.222222/v;
        ch = v*pow((x*sqrt(p1)+1.0-p1), 3.0);
        if (ch > 2.2*v+6.0)
            ch = -2.0*(log(1.0-p)-c*log(0.5*ch)+g);
    l4:
        q = ch;
        p1 = 0.5*ch;
        if ( (t = incompleteGamma(p1, xx, g)) < 0.0 )
            {
            std::cerr << "Error in function \"IncompleteGamma" << std::endl;
            return (-1.0);
            }
        p2 = p-t;
        t = p2*exp(xx*aa+g+p1-c*log(ch));
        b = t/ch;
        a = 0.5*t-b*c;
        double s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
        double s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
        double s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
        double s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
        double s5 = (84.0+264.0*a+c*(175.0+606.0*a))/2520.0;
        double s6 = (120.0+c*(346.0+127.0*c))/5040.0;
        ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
        if (fabs(q/ch-1.0) > e)
            goto l4;
    return (ch);
}

void RandomVariable::dirichletRv(const std::vector<double> &a, std::vector<double> &z) {

    int n = (int)a.size();
    double sum = 0.0;
    for(int i=0; i<n; i++)
        {
        /* z[i] = rndGamma(a[i]) / 1.0; */
        z[i] = rndGamma(a[i]);
        sum += z[i];
        }
    for(int i=0; i<n; i++)
        z[i] /= sum;
}

void RandomVariable::discretizeGamma(std::vector<double>& catRate, double a, double b, int nCats, bool median) {

    double factor = a / b * nCats;

    if (median)
        {
        /* the median value for each category is used to represent all of the values
           in that category */
        double interval = 1.0 / (2.0 * nCats);
        for (int i=0; i<nCats; i++)
            catRate[i] = chiSquareQuantile((i * 2.0 + 1.0) * interval, 2.0 * a) / (2.0 * b);
        double t = 0.0;
        for (int i=0; i<nCats; i++)
            t += catRate[i];
        for (int i=0; i<nCats; i++)
            catRate[i] *= factor / t;
        }
    else
        {
        /* the mean value for each category is used to represent all of the values
           in that category */
        /* calculate the points in the gamma distribution */
        for (int i=0; i<nCats-1; i++)
            catRate[i] = chiSquareQuantile((i + 1.0) / nCats, 2.0 * a) / (2.0 * b);
        /* calculate the cumulative values */
        double lnGammaValue = lnGamma(a + 1.0);
        for (int i=0; i<nCats-1; i++)
            catRate[i] = incompleteGamma(catRate[i] * b, a + 1.0, lnGammaValue);
        catRate[nCats-1] = 1.0;
        /* calculate the relative values and rescale */
        for (int i=nCats-1; i>0; i--)
            {
            catRate[i] -= catRate[i-1];
            catRate[i] *= factor;
            }
        catRate[0] *= factor;
        }
}

double RandomVariable::exponentialRv(double lambda) {

    double u = uniformRv();
    return -(1.0/lambda) * log(u);
}

uint32_t RandomVariable::extractU32(void) {

    int i = index;
    if (index >= N)
        {
        twist();
        i = index;
        }

    uint32_t y = mt[i];
    index = i + 1;

    y ^= (mt[i] >> U);
    y ^= (y << S) & B;
    y ^= (y << T) & C;
    y ^= (y >> L);

    return y;
}

double RandomVariable::gammaRv(double a, double b) {

    return (rndGamma(a) / b);
}

double RandomVariable::incompleteBeta(double a, double b, double x) {
    
    double tol = 1.0E-07;
    
    double value;
    if ( x <= 0.0 )
    {
        value = 0.0;
        return value;
    }
    else if ( 1.0 <= x )
    {
        value = 1.0;
        return value;
    }
    
    /* change tail if necessary and determine S */
    double psq = a + b;
    
    double xx, cx, pp, qq;
    bool indx;
    if ( a < (a + b) * x )
    {
        xx = 1.0 - x;
        cx = x;
        pp = b;
        qq = a;
        indx = true;
    }
    else
    {
        xx = x;
        cx = 1.0 - x;
        pp = a;
        qq = b;
        indx = false;
    }
    
    double term = 1.0;
    int i = 1;
    value = 1.0;
    int ns = (int)(qq + cx * (a + b));
    
    /* use Soper's reduction formulas */
    double rx = xx / cx;
    
    double temp = qq - (double)i;
    if ( ns == 0 )
        rx = xx;
    
    int it = 0;
    int it_max = 1000;
    for (;;)
    {
        it++;
        if ( it_max < it )
        {
            //std::cerr << "Error in incompleteBeta: Maximum number of iterations exceeded!" << std::endl;
            return -1;
        }
        term = term * temp * rx / ( pp + ( double ) ( i ) );
        value = value + term;
        temp = fabs(term);
        if ( temp <= tol && temp <= tol * value )
            break;
        i++;
        ns--;
        if ( 0 <= ns )
        {
            temp = qq - (double)i;
            if ( ns == 0 )
                rx = xx;
        }
        else
        {
            temp = psq;
            psq = psq + 1.0;
        }
    }
    
    /* finish calculation */
    value = value * exp(pp * log(xx) + (qq - 1.0) * log(cx) - lnBeta(a, b) - log(pp));
    if ( indx )
        value = 1.0 - value;
    return value;
}

double RandomVariable::incompleteGamma (double x, double alpha, double LnGamma_alpha) {

    double            p = alpha, g = LnGamma_alpha,
                    accurate = 1e-8, overflow = 1e30,
                    rn = 0.0, a = 0.0, b = 0.0, an = 0.0,
                    gin, dif = 0.0, term = 0.0, pn[6];

    if (x == 0.0)
        return (0.0);
    if (x < 0 || p <= 0)
        return (-1.0);

    double factor = exp(p*log(x)-x-g);
    if (x > 1 && x >= p)
        goto l30;
    gin = 1.0;
    term = 1.0;
    rn = p;
    l20:
        rn++;
        term *= x/rn;
        gin += term;
        if (term > accurate)
            goto l20;
        gin *= factor/p;
        goto l50;
    l30:
        a = 1.0-p;
        b = a+x+1.0;
        term = 0.0;
        pn[0] = 1.0;
        pn[1] = x;
        pn[2] = x+1;
        pn[3] = x*b;
        gin = pn[2]/pn[3];
    l32:
        a++;
        b += 2.0;
        term++;
        an = a*term;
        for (int i=0; i<2; i++)
            pn[i+4] = b*pn[i+2]-an*pn[i];
        if (pn[5] == 0)
            goto l35;
        rn = pn[4]/pn[5];
        dif = fabs(gin-rn);
        if (dif>accurate)
            goto l34;
        if (dif<=accurate*rn)
            goto l42;
    l34:
        gin = rn;
    l35:
        for (int i=0; i<4; i++)
            pn[i] = pn[i+2];
        if (fabs(pn[4]) < overflow)
            goto l32;
        for (int i=0; i<4; i++)
            pn[i] /= overflow;
        goto l32;
    l42:
        gin = 1.0-factor*gin;
    l50:
        return (gin);
}

void RandomVariable::initialize(uint32_t seed) {

    mt[0] = seed;
    for (uint32_t i=1; i<N; i++)
        {
        mt[i] = (F * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i);
        }
    index = N;
    
    for (size_t i=0; i<10000; i++)
        extractU32();
}

/** Tests whether a double is finite */
bool RandomVariable::isFinite(double x) {
  
    return x > -std::numeric_limits<double>::infinity() && x < std::numeric_limits<double>::infinity();
}

bool RandomVariable::isNan(double x) {

    return x != x;
}

double RandomVariable::lnDirichletPdf(const std::vector<double> &a, const std::vector<double> &z) {

    int n = (int)a.size(); //!< we assume that a and z have the same size
    double alpha0 = 0.0;
    for (int i=0; i<n; i++)
        alpha0 += a[i];
    double lnP = lnGamma(alpha0);
    for (int i=0; i<n; i++)
        lnP -= lnGamma(a[i]);
    for (int i=0; i<n; i++)
        lnP += (a[i] - 1.0) * log(z[i]);
    return lnP;
}

double RandomVariable::lnExponentialPdf(double lambda, double x) {

    return (std::log(lambda) - lambda * x);
}

double RandomVariable::lnBeta(double a, double b)
{
    double corr, p, q;
    
    p = q = a;
    if (b < p) p = b;/* := min(a,b) */
    if (b > q) q = b;/* := max(a,b) */
    
    /* both arguments must be >= 0 */
    if (p < 0)
    {

        std::cout << "Cannot compute log-beta function for a = " << a << " and b = " << b;
    }
    else if (p == 0)
    {
        return -1;
        //return RbConstants::Double::inf;
    }
    else if (!isFinite(q))
    { /* q == +Inf */
        return -1;
        //return RbConstants::Double::neginf;
    }
    
    if (p >= 10) {
        /* p and q are big. */
        corr = lnGammacor(p) + lnGammacor(q) - lnGammacor(p + q);
        return log(q) * -0.5 + 0.918938533204672741780329736406 + corr
            + (p - 0.5) * log(p / (p + q)) + q * log1p(-p / (p + q));
    }
    else if (q >= 10) {
        /* p is small, but q is big. */
        corr = lnGammacor(q) - lnGammacor(p + q);
        return lnGamma(p) + corr + p - p * log(p + q) + (q - 0.5) * log1p(-p / (p + q));
    }
    else
        /* p and q are small: p <= q < 10. */
        return log(gamma(p) * (gamma(q) / gamma(p + q)));
    
}

double RandomVariable::log1p(double x) {

    if (x <= -1.0)
    {
        Msg::error("Invalid input argument for log1p: must be greater than -1.0");
    }

    if (fabs(x) > 1e-4)
    {
        // x is large enough that the obvious evaluation is OK
        return log(1.0 + x);
    }

    // Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
    // Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8

    return (-0.5*x + 1.0)*x;
}

double RandomVariable::lnGamma(double a) {

    double x = a;
    double f = 0.0;
    double z;
    if (x < 7)
        {
        f = 1.0;
        z = x - 1.0;
        while (++z < 7.0)
            f *= z;
        x = z;
        f = -log(f);
        }
    z = 1.0 / (x*x);
    return  (f + (x-0.5)*log(x) - x + 0.918938533204673 +
            (((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z +
            0.083333333333333)/x);
}

double RandomVariable::chebyshev_eval(double x, const double *a, const int n) {

    double b0, b1, b2, twox;
    int i;
    
    if (n < 1 || n > 1000)
    {
        Msg::error("Cannot compute chebyshev function");
    }
    
    if (x < -1.1 || x > 1.1)
    {
        Msg::error("Cannot compute chebyshev function");
    }
    
    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++)
    {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}

double RandomVariable::gamma(double x) {

    const static double gamcs[42] = {
        +.8571195590989331421920062399942e-2,
        +.4415381324841006757191315771652e-2,
        +.5685043681599363378632664588789e-1,
        -.4219835396418560501012500186624e-2,
        +.1326808181212460220584006796352e-2,
        -.1893024529798880432523947023886e-3,
        +.3606925327441245256578082217225e-4,
        -.6056761904460864218485548290365e-5,
        +.1055829546302283344731823509093e-5,
        -.1811967365542384048291855891166e-6,
        +.3117724964715322277790254593169e-7,
        -.5354219639019687140874081024347e-8,
        +.9193275519859588946887786825940e-9,
        -.1577941280288339761767423273953e-9,
        +.2707980622934954543266540433089e-10,
        -.4646818653825730144081661058933e-11,
        +.7973350192007419656460767175359e-12,
        -.1368078209830916025799499172309e-12,
        +.2347319486563800657233471771688e-13,
        -.4027432614949066932766570534699e-14,
        +.6910051747372100912138336975257e-15,
        -.1185584500221992907052387126192e-15,
        +.2034148542496373955201026051932e-16,
        -.3490054341717405849274012949108e-17,
        +.5987993856485305567135051066026e-18,
        -.1027378057872228074490069778431e-18,
        +.1762702816060529824942759660748e-19,
        -.3024320653735306260958772112042e-20,
        +.5188914660218397839717833550506e-21,
        -.8902770842456576692449251601066e-22,
        +.1527474068493342602274596891306e-22,
        -.2620731256187362900257328332799e-23,
        +.4496464047830538670331046570666e-24,
        -.7714712731336877911703901525333e-25,
        +.1323635453126044036486572714666e-25,
        -.2270999412942928816702313813333e-26,
        +.3896418998003991449320816639999e-27,
        -.6685198115125953327792127999999e-28,
        +.1146998663140024384347613866666e-28,
        -.1967938586345134677295103999999e-29,
        +.3376448816585338090334890666666e-30,
        -.5793070335782135784625493333333e-31
    };
    
    int i, n;
    double y;
    double sinpiy, value;
    
    /* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
     * (xmin, xmax) are non-trivial, see ./gammalims.c
     * xsml = exp(.01)*DBL_MIN
     * dxrel = sqrt(DBL_EPSILON) = 2 ^ -26
     */
# define ngam 22
# define xmin -170.5674972726612
# define xmax  171.61447887182298
# define xsml 2.2474362225598545e-308
# define dxrel 1.490116119384765696e-8
    
    if (isNan(x)) return x;
    
    /* If the argument is exactly zero or a negative integer
     * then return NaN. */
    if (x == 0 || (x < 0 && x == (long)x))
    {
        Msg::error("Cannot compute gamma function");
    }
    
    y = fabs(x);
    
    if (y <= 10) {
        
        /* Compute gamma(x) for -10 <= x <= 10
         * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
         * first of all. */
        
        n = int(x);
        if (x < 0) --n;
        y = x - n;/* n = floor(x)  ==>    y in [ 0, 1 ) */
        --n;
        value = chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
        if (n == 0)
            return value;/* x = 1.dddd = 1+y */
        
        if (n < 0) {
            /* compute gamma(x) for -10 <= x < 1 */
            
            /* exact 0 or "-n" checked already above */
            
            /* The answer is less than half precision */
            /* because x too near a negative integer. */
            if (x < -0.5 && fabs(x - (int)(x - 0.5) / x) < dxrel) {
                Msg::error("Cannot compute gamma function");
            }
            
            /* The argument is so close to 0 that the result would overflow. */
            if (y < xsml) {
                Msg::error("Cannot compute gamma function");
            }
            
            n = -n;
            
            for (i = 0; i < n; i++) {
                value /= (x + i);
            }
            return value;
        }
        else {
            /* gamma(x) for 2 <= x <= 10 */
            
            for (i = 1; i <= n; i++) {
                value *= (y + i);
            }
            return value;
        }
    }
    else {
        /* gamma(x) for     y = |x| > 10. */
        
        if (x > xmax) {            /* Overflow */
            Msg::error("Cannot compute gamma function");
        }
        
        if (x < xmin) {            /* Underflow */
            Msg::error("Cannot compute gamma function");
        }
        
        if (y <= 50 && y == (int)y) { /* compute (n - 1)! */
            value = 1.;
            for (i = 2; i < y; i++) value *= i;
        }
        else { /* normal case */
            value = exp((y - 0.5) * log(y) - y + 0.918938533204672741780329736406 + ((2*y == (int)2*y)? stirlerr(y) : lnGammacor(y)));
        }
        if (x > 0)
            return value;
        
        if (fabs((x - (int)(x - 0.5))/x) < dxrel){
            
            /* The answer is less than half precision because */
            /* the argument is too near a negative integer. */
            
            Msg::error("Cannot compute gamma function");
        }
        
        sinpiy = sin(3.141592653589793238462643383280 * y);
        if (sinpiy == 0) {        /* Negative integer arg - overflow */
            Msg::error("Cannot compute gamma function");
        }
        
        return -3.141592653589793238462643383280 / (y * sinpiy * value);
    }

}

double RandomVariable::lnGammacor(double x) {

    const static double algmcs[15] = {
        +.1666389480451863247205729650822e+0,
        -.1384948176067563840732986059135e-4,
        +.9810825646924729426157171547487e-8,
        -.1809129475572494194263306266719e-10,
        +.6221098041892605227126015543416e-13,
        -.3399615005417721944303330599666e-15,
        +.2683181998482698748957538846666e-17,
        -.2868042435334643284144622399999e-19,
        +.3962837061046434803679306666666e-21,
        -.6831888753985766870111999999999e-23,
        +.1429227355942498147573333333333e-24,
        -.3547598158101070547199999999999e-26,
        +.1025680058010470912000000000000e-27,
        -.3401102254316748799999999999999e-29,
        +.1276642195630062933333333333333e-30
    };
    
    double tmp;
    
    /* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
     *   xbig = 2 ^ 26.5
     *   xmax = DBL_MAX / 48 =  2^1020 / 3 */
#define nalgm 5
#define xbig  94906265.62425156
#undef  xmax
#define xmax  3.745194030963158e306
    
    if (x < 10)
    {
        Msg::error("Cannot compute log-gammacor function");
    }
    else if (x >= xmax) {
        Msg::error("Cannot compute log-gammacor function");
        /* allow to underflow below */
    }
    else if (x < xbig) {
        tmp = 10 / x;
        return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
    }
    return 1 / (x * 12);
}

double RandomVariable::lnGammaPdf(double a, double b, double x) {

    return a * log(b) - lnGamma(a) + (a - 1.0) * log(x) - x * b;
}

double RandomVariable::pointNormal(double prob) {

    double a0 = -0.322232431088;
    double a1 = -1.0;
    double a2 = -0.342242088547;
    double a3 = -0.0204231210245;
    double a4 = -0.453642210148e-4;
    double b0 = 0.0993484626060;
    double b1 = 0.588581570495;
    double b2 = 0.531103462366;
    double b3 = 0.103537752850;
    double b4 = 0.0038560700634;
    double p = prob;
    double p1 = ( p < 0.5 ? p : 1.0 - p);
    if (p1 < 1e-20)
        return (-9999.0);
    double y = sqrt( log(1.0/(p1*p1)) );
    double z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
    return ( p < 0.5 ? -z : z );
}

double RandomVariable::rndGamma(double s) {

    double r = 0.0;
    if (s <= 0.0)
        std::cout << "Gamma parameter less than zero" << std::endl;
    else if (s < 1.0)
        r = rndGamma1(s);
    else if (s > 1.0)
        r = rndGamma2(s);
    else
        r = -log(uniformRv());
    return (r);
}

double RandomVariable::rndGamma1(double s) {

    double            r, x = 0.0, small = 1e-37, w;
    static double   a, p, uf, ss = 10.0, d;
    
    if (s != ss)
        {
        a  = 1.0 - s;
        p  = a / (a + s * exp(-a));
        uf = p * pow(small / a, s);
        d  = a * log(a);
        ss = s;
        }
    for (;;)
        {
        r = uniformRv();
        if (r > p)
            {
            //x = a - log((1.0 - r) / (1.0 - p)), w = a * log(x) - d;  // The use of the comma operator here
            x = a - log((1.0 - r) / (1.0 - p));                        // and below, in the next conditional
            w = a * log(x) - d;                                        // statement, seems fishy.
            }
        else if (r>uf)
            {
            //x = a * pow(r / p, 1.0 / s), w = x;
            x = a * pow(r / p, 1.0 / s);
            w = x;
            }
        else
            return (0.0);
        r = uniformRv();
        if (1.0 - r <= w && r > 0.0)
        if (r*(w + 1.0) >= 1.0 || -log(r) <= w)
            continue;
        break;
        }
    
    return (x);
}

double RandomVariable::rndGamma2(double s) {

    double            r, d, f, g, x;
    static double    b, h, ss = 0.0;
    
    if (s != ss)
        {
        b  = s - 1.0;
        h  = sqrt(3.0 * s - 0.75);
        ss = s;
        }
    for (;;)
        {
        r = uniformRv();
        g = r - r * r;
        f = (r - 0.5) * h / sqrt(g);
        x = b + f;
        if (x <= 0.0)
            continue;
        r = uniformRv();
        d = 64.0 * r * r * g * g * g;
        if (d * x < x - 2.0 * f * f || log(d) < 2.0 * (b * log(x / b) - f))
            break;
        }
    
    return (x);
}

double RandomVariable::stirlerr(double n) {
    
#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */
    
    /*
     error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
     */
    const static double sferr_halves[31] = {
        0.0, /* n=0 - wrong, place holder only */
        0.1534264097200273452913848,  /* 0.5 */
        0.0810614667953272582196702,  /* 1.0 */
        0.0548141210519176538961390,  /* 1.5 */
        0.0413406959554092940938221,  /* 2.0 */
        0.03316287351993628748511048, /* 2.5 */
        0.02767792568499833914878929, /* 3.0 */
        0.02374616365629749597132920, /* 3.5 */
        0.02079067210376509311152277, /* 4.0 */
        0.01848845053267318523077934, /* 4.5 */
        0.01664469118982119216319487, /* 5.0 */
        0.01513497322191737887351255, /* 5.5 */
        0.01387612882307074799874573, /* 6.0 */
        0.01281046524292022692424986, /* 6.5 */
        0.01189670994589177009505572, /* 7.0 */
        0.01110455975820691732662991, /* 7.5 */
        0.010411265261972096497478567, /* 8.0 */
        0.009799416126158803298389475, /* 8.5 */
        0.009255462182712732917728637, /* 9.0 */
        0.008768700134139385462952823, /* 9.5 */
        0.008330563433362871256469318, /* 10.0 */
        0.007934114564314020547248100, /* 10.5 */
        0.007573675487951840794972024, /* 11.0 */
        0.007244554301320383179543912, /* 11.5 */
        0.006942840107209529865664152, /* 12.0 */
        0.006665247032707682442354394, /* 12.5 */
        0.006408994188004207068439631, /* 13.0 */
        0.006171712263039457647532867, /* 13.5 */
        0.005951370112758847735624416, /* 14.0 */
        0.005746216513010115682023589, /* 14.5 */
        0.005554733551962801371038690  /* 15.0 */
    };
    double nn;
    
    if (n <= 15.0) {
        nn = n + n;
        if (nn == (int)nn) return(sferr_halves[(int)nn]);
        return(lnGamma(n + 1.) - (n + 0.5)*log(n) + n - 0.918938533204672741780329736406);
    }
    
    nn = n*n;
    if (n>500) return((S0-S1/nn)/n);
    if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
    if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
    /* 15 < n <= 35 : */
    return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
}

void RandomVariable::twist(void) {

    for (uint32_t i=0; i<N; i++)
        {
        uint32_t x = (mt[i] & MASK_UPPER) + (mt[(i + 1) % N] & MASK_LOWER);
        uint32_t xA = x >> 1;

        if ( x & 0x1 )
            xA ^= A;

        mt[i] = mt[(i + M) % N] ^ xA;
        }
    index = 0;
}

double RandomVariable::uniformRv(void) {

    return (double)extractU32() / UINT32_MAX;
}
