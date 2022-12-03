#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <ostream>
#include <set>
#include <vector>

#ifdef PROFILE_FLAG
#  define PROFILE __attribute__((noinline))
#else
#  define PROFILE
#endif

using Scalar = std::int32_t;
using DoubleScalar = std::int64_t;

struct Number {//TODO: Cn
    std::vector<Scalar> primes;
    std::vector<Scalar> powersOfPrimes;
    std::vector<Scalar> divisors;
    std::vector<Scalar> coprimes;
    Scalar n;
    Scalar powerOfTwo = 0;
    Scalar phi = 0;
    Scalar lambda = 0;
    Scalar nskew = 0;
    bool squareFree = true;
};

std::vector<Number> numberCache = {};//TODO: non static?

Scalar getMaxPrime(const Number &number) {
    if (number.primes.empty()) {
        return 1;
    }
    return number.primes.back();
}

bool isPowerOfPrime(Number &number) {
    return number.primes.size() == 1;
}

PROFILE Scalar gcd(Scalar n, Scalar k) {
    while (k != 0) {
        const Scalar new_k = n % k;
        n = k;
        k = new_k;
    }
    return n;
}

PROFILE Scalar inverse(Scalar n, Scalar k) {
    Scalar i = 0;
    Scalar j = 1;
    while (k != 0) {
        const auto q = n / k;
        const auto new_k = n % k;
        n = k;
        k = new_k;
        const auto new_i = j;
        j = i - q * j;
        i = new_i;
    }
    if (n == 1) return i;
    return 0;
}

PROFILE Scalar gcd(const Number &number_n, const Number &number_k) {
    if (number_k.n == 0) {
        return number_n.n;
    }
    std::size_t ni = 0;//TODO: iterators?
    std::size_t ki = 0;
    const auto &n_primes = number_n.primes;
    const auto &k_primes = number_k.primes;
    const auto &n_powersOfPrimes = number_n.powersOfPrimes;
    const auto &k_powersOfPrimes = number_k.powersOfPrimes;
    const auto n_size = n_primes.size();
    const auto k_size = k_primes.size();
    Scalar result = 1;
    while (ni < n_size && ki < k_size) {
        const auto n_prime = n_primes[ni];
        const auto k_prime = k_primes[ki];
        if (n_prime < k_prime) {
            ++ni;
        } else if (n_prime > k_prime) {
            ++ki;
        } else {
            result *= std::min(n_powersOfPrimes[ni], k_powersOfPrimes[ki]);
            ++ki;
            ++ni;
        }
    }
    return result;
}

PROFILE bool isCoprime(const Number &number_n, const Number &number_k) {
    if (number_k.n == 0) {//TODO: dava to zmysel?
        return false;
    }
    const auto &n_primes = number_n.primes;
    const auto &k_primes = number_k.primes;
    auto ni = n_primes.begin();
    auto ki = k_primes.begin();
    const auto n_end = n_primes.end();
    const auto k_end = k_primes.end();
    while (ni != n_end && ki != k_end) {
        const auto n_prime = *ni;
        const auto k_prime = *ki;
        if (n_prime < k_prime) {
            ++ni;
        } else if (n_prime > k_prime) {
            ++ki;
        } else {
            return false;
        }
    }
    return true;
}

Scalar powerOfTwo(Scalar n) {
    return (n & (~(n - 1)));
}

PROFILE void factorize(Number& number, const std::vector<Scalar> &primes) {
    if (number.powerOfTwo != 0) {
        return;
    }
    auto n = number.n;
    number.powerOfTwo = powerOfTwo(n);
    Scalar divisorsCount = 1;
    Scalar powerOfPrime = 1;
    Scalar prime = n;
    if (n > 1) {
        for (const auto p: primes) {
            if (p * p > n) {
                break;
            }
            if (n % p == 0) {
                prime = p;
                powerOfPrime = p;
                n /= p;//TODO: optimize divisions?
                ++divisorsCount;

                while (n % p == 0) {
                    if (p != 2) {
                        number.squareFree = false;
                    }
                    powerOfPrime *= p;
                    n /= p;
                    ++divisorsCount;
                }
                break;
            }
        }
        if (n == prime) {
            powerOfPrime = n;
            n = 1;
            ++divisorsCount;
        }
        const auto &number_n_p = numberCache[n];
        const auto primesCount = number_n_p.primes.size() + 1;
        number.primes.reserve(primesCount);
        number.powersOfPrimes.reserve(primesCount);
        number.primes.push_back(prime);
        number.powersOfPrimes.push_back(powerOfPrime);
        if (n == 1) {
            number.divisors.reserve(divisorsCount);
            for (DoubleScalar power = 1; power <= powerOfPrime; power *= prime) {
                number.divisors.push_back(power);
            }
        } else {
            const auto &number_p_k = numberCache[powerOfPrime];
            number.primes.insert(number.primes.end(), number_n_p.primes.begin(), number_n_p.primes.end());
            number.powersOfPrimes.insert(number.powersOfPrimes.end(), number_n_p.powersOfPrimes.begin(), number_n_p.powersOfPrimes.end());
            number.squareFree = number.squareFree && number_n_p.squareFree;
            divisorsCount *= number_n_p.divisors.size();
            number.divisors.reserve(divisorsCount);
            for (const auto d: number_p_k.divisors) {
                for (const auto e: number_n_p.divisors) {
                    number.divisors.push_back(d * e);
                }
            }
            std::sort(number.divisors.begin(), number.divisors.end());//TODO: potrebujeme?
        }
    } else {
        number.divisors = {1};
    }
}

PROFILE void computePhi(Number &number) {
    if (number.phi != 0) {
        return;
    }//TODO: factorize
    number.phi = number.n;
    for (const auto p: number.primes) {
        number.phi -= number.phi / p;
    }
    if (number.primes.size() <= 1) {
        if (number.powerOfTwo > 4) {
            number.lambda = number.phi / 2;//TODO: >> 2??
        } else {
            number.lambda = number.phi;
        }
    } else {//TODO: optimize
        const auto p = number.primes.back();
        auto n_p = number.n / p;
        auto p_k = p;
        while (n_p % p == 0) {
            n_p /= p;
            p_k *= p;
        }
        auto &number_n_p = numberCache[n_p];//TODO: computePhi
        const auto phi_p_k = p_k - p_k / p;
        number.lambda = number_n_p.lambda * phi_p_k / gcd(numberCache[number_n_p.lambda], numberCache[phi_p_k]);//TODO: lcm
    }
}

PROFILE void prepareNumbers(Scalar N) {
    numberCache.resize(N + 1);
    std::vector<Scalar> primes;
    primes.reserve(1.25506 * 2 * std::sqrt(N) / log(N));
    for (Scalar i = 0; i <= N; ++i) {
        auto &number = numberCache[i];
        number.n = i;
        factorize(number, primes);
        if (isPowerOfPrime(number) && i == number.primes.back()) {//TODO: isPrime, computePhi to tiez zistuje
            if (i * DoubleScalar(i) <= N) {
                primes.push_back(i);
            }
        }
        computePhi(number);//TODO: lazy ak chceme pouzivat A, B
    }
}

void clear(Number &number);

PROFILE void computeCoprimes(Number &number) {
    if (!number.coprimes.empty()) {
        return;
    }
    const auto n = number.n;
    if (n == 0) {
        return;
    }
    computePhi(number);
    number.coprimes.reserve(number.phi);
    if (n == 1) {
        number.coprimes.push_back(0);
        return;
    }
    const auto p = getMaxPrime(number);
    if (n == p) {
        for (Scalar e = 1; e < n; ++e) {//TODO: spojit s p^k
            number.coprimes.push_back(e);
        }
    } else if (isPowerOfPrime(number)) {
        const auto n_p = n / p;
        auto &number_p = numberCache[n_p];
        computeCoprimes(number_p);
        for (Scalar m = 0; m < n; m += n_p) {
            for (const auto e: number_p.coprimes) {
                number.coprimes.push_back(m + e);
            }
        }
        clear(number_p);
    } else {
        const auto p_k = number.powersOfPrimes.back();
        const auto n_p_k = n / p_k;
        auto &number_n_p_k = numberCache[n_p_k];
        auto &number_p_k = numberCache[p_k];
        computeCoprimes(number_n_p_k);
        computeCoprimes(number_p_k);
        // l * n_p_k + a = e * n_p_k (mod p_k)
        const auto r = inverse(p_k, n_p_k % p_k);
        for (const auto a: number_n_p_k.coprimes) {
            const auto b = (n - a * r) % p_k;
            for (const auto e: number_p_k.coprimes) {
                auto l = e + b;
                if (l >= p_k) l -= p_k;
                number.coprimes.push_back(l * n_p_k + a);
            }
        }
        clear(number_n_p_k);
        clear(number_p_k);
    }
}

PROFILE Scalar compute_a(Scalar d, Scalar power_sum_e, Scalar rH, DoubleScalar power_sum, Scalar n_h) {
    return (((power_sum_e - d) / rH + n_h) * power_sum + d) % n_h;//TODO: optimize
}

PROFILE Scalar countCoprimeSolutions(const Number &number_a_b, const Number &number_n_b, const Number &number_gcd_nh_b) {//TODO: priamo cez rozklad
    // ax = b (%n_h)
    if (!isCoprime(number_n_b, number_a_b)) {
        return 0;
    }
    const auto gcd_nh_b_b = gcd(number_n_b, number_gcd_nh_b);
    const auto &number_gcd_nh_b_b = numberCache[gcd_nh_b_b];
    return gcd_nh_b_b * number_gcd_nh_b.phi / number_gcd_nh_b_b.phi;//TODO: can be precompted
}

PROFILE Scalar sEquals1(const Scalar d, const Number &number_n_div_d) {
    Scalar nskew = 0;
    std::vector<bool> visited_e(number_n_div_d.n, false);//TODO: global
    std::vector<Scalar> powers;//TODO: global
    powers.reserve(number_n_div_d.n);
    for (const auto n_h : number_n_div_d.divisors) {
        if (n_h < d) {//TODO: range from second / except last / skip based on predicate
            continue;
        }
        auto &number_n_h = numberCache[n_h];
        if (number_n_h.lambda % d != 0) {
            continue;
        }
        for (auto e = 2; e < n_h; ++e) {
            if (visited_e[e]) {
                continue;
            }
            if (!isCoprime(number_n_h, numberCache[e])) {
                continue;
            }
            DoubleScalar power_e = e;
            DoubleScalar power_sum_e = 1;
            powers.push_back(1);
            while (power_e != 1) {//TODO: until visited, power sum reuse?
                powers.push_back(power_e);
                power_e = (power_e * e) % number_n_h.n;//TODO: alias n_h, n_div_d
            }
            const Scalar order = powers.size();
            if (order % d == 0) {
                const auto step = order / d;
                if (!visited_e[powers[step]]) {
                    for (auto i = step; i < order; i += step) {
                        power_sum_e += powers[i];
                    }
                    power_sum_e = power_sum_e % number_n_h.n;
                    if (power_sum_e == 0) {
                        nskew += numberCache[d].phi * number_n_h.phi;
                    }
                }
            }
            for (const auto power: powers) {
                visited_e[power] = true;
            }
            powers.clear();
        }
        for (auto e = 2; e < n_h; ++e) {
            visited_e[e] = false;
        }
    }
    return nskew;
}

PROFILE Scalar sOtherThan1(const Scalar d, const Number &number_n_div_d) {
    Scalar nskew = 0;
    for (const auto gcd_b: number_n_div_d.divisors) {
        if (gcd_b == number_n_div_d.n) {
            continue;
        }
        if (gcd_b < d) {//TODO: nejaky lepsi check? vnutri? toto skoro nic nerobi aj tak sa to dost skoro zahodi
            continue;
        }
        auto &number_n_b = numberCache[number_n_div_d.n / gcd_b];
        computeCoprimes(number_n_b);
        std::vector<bool> visited_s(number_n_div_d.n, false);//TODO: global
        std::vector<bool> visited_e(number_n_div_d.n, false);//TODO: global
        std::vector<Scalar> powers;//TODO: global
        powers.reserve(number_n_div_d.n);
        for (const auto coprime_b: number_n_b.coprimes) {
            const auto b = gcd_b * coprime_b;
            const auto s = b + 1;
            if (visited_s[s]) {
                continue;
            }
            if (!isCoprime(number_n_div_d, numberCache[s])) {
                continue;
            }
            DoubleScalar power_s = s;
            DoubleScalar power_sum = 1;
            powers.push_back(1);
            while (power_s != 1) {
                power_sum += power_s;
                powers.push_back(power_s);
                power_s = (power_s * s) % number_n_div_d.n;
            }
            const Scalar rH = powers.size();
            auto &number_rH = numberCache[rH];
            computeCoprimes(number_rH);
            for (const auto coprime_rH: number_rH.coprimes) {
                visited_s[powers[coprime_rH]] = true;
            }
            if (number_rH.n != number_n_b.n) {
                clear(number_rH);
            }
            powers.clear();
            power_sum = power_sum % number_n_div_d.n;//TODO: treba toto modulovat? alias n_div_d
            const auto power_sum_b = power_sum / number_n_b.n;
            const auto min_r = d * rH;//TODO: better estimate?
            const auto &number_gcd_b = numberCache[gcd_b];
            for (const auto gcd_nh_b: number_gcd_b.divisors) {
                if (gcd_nh_b < d) {
                    continue;
                }
                const auto n_h = number_n_b.n * gcd_nh_b;
                if (n_h < min_r) {//TODO: necessary?
                    continue;
                }
                if (power_sum_b % gcd_nh_b == 0) {//TODO: gcd ratame o riadok nizsie
                    continue;
                }
                const auto &number_gcd_nh_b = numberCache[gcd_nh_b];
                const auto g_p = gcd(numberCache[power_sum_b], number_gcd_nh_b);
                const auto g_r = gcd_nh_b / g_p;
                if (g_r < d) {
                    continue;
                }
                const auto r = rH * g_r;
                auto &number_r = numberCache[r];
                if (number_r.lambda % d != 0) {
                    continue;
                }
                for (auto e = rH + 1; e < r; e += rH) {
                    if (visited_e[e]) {
                        continue;
                    }
                    if (!isCoprime(number_r, numberCache[e])) {
                        continue;
                    }
                    DoubleScalar power_e = e;
                    DoubleScalar power_sum_e = 1;
                    powers.push_back(1);
                    while (power_e != 1) {//TODO: until visited, power sum reuse?
                        powers.push_back(power_e);
                        power_e = (power_e * e) % number_r.n;//TODO: alias n_h, n_div_d
                    }
                    const Scalar order = powers.size();
                    if (order % d == 0) {
                        const auto step = order / d;
                        if (!visited_e[powers[step]]) {
                            for (auto i = step; i < order; i += step) {
                                power_sum_e += powers[i];
                            }
                            power_sum_e = power_sum_e % number_r.n;

                            if (power_sum_e % gcd_nh_b == 0) {
                                const auto a = compute_a(d, power_sum_e, rH, power_sum, n_h);
                                const auto &number_a_b = numberCache[a / gcd_nh_b];
                                nskew += numberCache[d].phi * number_rH.phi * countCoprimeSolutions(number_a_b, number_n_b, number_gcd_nh_b);
                            }
                        }
                    }
                    for (const auto power: powers) {
                        visited_e[power] = true;
                    }
                    powers.clear();
                }
                for (auto e = rH + 1; e < r; e += rH) {
                    visited_e[e] = false;
                }
            }
        }
        clear(number_n_b);
    }
    return nskew;
}

PROFILE void countSkewmorphisms(Number &number) {
    if (number.nskew != 0) {
        return;
    }
    const auto n = number.n;
    const auto phi = number.phi;
    auto nskew = phi;

    if (gcd(number, numberCache[phi]) > 1) {//TODO: toto viem nejako vyuzit a predratat?
        const auto maxPrime = getMaxPrime(number);
        const auto &number_n_div_maxPrime = numberCache[n / maxPrime];
        for (const auto d: number_n_div_maxPrime.divisors) {//TODO: iterate over n_d first
            if (d == 1) {
                continue;
            }
            const auto n_div_d = n / d;
            auto &number_n_div_d = numberCache[n_div_d];

            nskew += sEquals1(d, number_n_div_d);
            nskew += sOtherThan1(d, number_n_div_d);
        }
    }
    number.nskew = nskew;
}

void print(const Number &number) {
    const auto n = number.n;
    const auto phi = number.phi;
    const auto nskew = number.nskew;
    printf("Total number of skew morphisms of C_%d is %d\n", n, nskew);
    printf("Check sub-total of automorphisms of C_%d is %d\n", n, phi);
    auto perc = static_cast<Scalar>(std::round(float(100 * phi) / float(nskew)));
    printf("Automorphisms account for ~%d%% of all skew morphisms.\n\n", perc);
}

void fprint(const Number &number, std::ostream& stream) {
    const auto n = number.n;
    const auto phi = number.phi;
    const auto nskew = number.nskew;
    stream << "Total number of skew morphisms of C_" << n << " is " << nskew << "\n";
    stream << "Check sub-total of automorphisms of C_" << n << " is " << phi << "\n";
    auto perc = static_cast<Scalar>(std::round(float(100 * phi) / float(nskew)));
    stream << "Automorphisms account for ~" << perc << "% of all skew morphisms.\n\n";
}

void clear(Number &number) {
    if (!number.coprimes.empty()) {//TODO: neratame nieco viackrat? jasne, ze hej, mozme si predratat reference counter?
        number.coprimes = std::vector<Scalar>{};
    }
}

int main(int argc, char *argv[]) {
    Scalar N = 10000;
    if (argc > 1) {
        N = std::stoi(argv[1]);
    }
    Scalar A = 1;
    if (argc > 2) {
        A = std::stoi(argv[2]);
    }

    prepareNumbers(N);

    for (auto n = A; n <= N; ++n) {
        auto &number = numberCache[n];
        if (number.powerOfTwo <= 16 && number.squareFree) {
            countSkewmorphisms(number);
//            fprint(number, std::cerr);
        }
        if (number.nskew != 0) {
            print(number);
        }
    }

    return 0;
}
