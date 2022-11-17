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
    std::vector<Scalar> orders;// TODO: s v deliteli n
    std::vector<Scalar> order2OrderIndex;
    std::vector<Scalar> orderIndex2CoprimesBegin;
    std::vector<Scalar> divisors;
    std::vector<Scalar> coprimes;//TODO: vector vectorov
    std::vector<Scalar> inverses;
    std::vector<Scalar> powerSums;
    std::vector<Scalar> gcdPowerSums;
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

PROFILE Scalar gcd(const Number &number_n, const Number &number_k) {
    if (number_k.n == 0) {
        return number_n.n;
    }
    std::size_t ni = 0;
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
            for (Scalar power = 1; power <= powerOfPrime; power *= prime) {
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

std::set<int> toCountWithMultiples = {};//TODO

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
            toCountWithMultiples.insert(i);
        }
        computePhi(number);//TODO: lazy ak chceme pouzivat A, B
    }
    toCountWithMultiples.insert(1);
}

void computeOrders(Number &number);

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
    toCountWithMultiples.insert(number.n);
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
    } else {
        const auto p_k = number.powersOfPrimes.back();
        const auto n_p_k = n / p_k;
        auto &number_n_p_k = numberCache[n_p_k];
        auto &number_p_k = numberCache[p_k];
        computeCoprimes(number_n_p_k);
        computeOrders(number_p_k);//TODO
        // l * n_p_k + a = e * n_p_k (mod p_k)
        const auto r = number_p_k.inverses[n_p_k % p_k];//TODO: jediny usage inverzov
        for (const auto a: number_n_p_k.coprimes) {
            const auto b = (n - a * r) % p_k;
            for (const auto e: number_p_k.coprimes) {
                auto l = e + b;
                if (l >= p_k) l -= p_k;
                number.coprimes.push_back(l * n_p_k + a);
            }
        }
    }
}

PROFILE void computeOrders(Number &number) {
    if (!number.orders.empty()) {
        return;
    }
    const auto n = number.n;
    if (n == 0) {
        return;
    }
    computeCoprimes(number);
    if (n == 1) {
        number.orders = {1};
        number.powerSums = {0};
        number.gcdPowerSums = {0};
        return;
    }
    number.orders.resize(n, 0);
    number.orders[1] = 1;
    if (isPowerOfPrime(number)) {
        number.inverses.resize(n, 0);
        number.inverses[1] = 1;
    }
    number.powerSums.resize(n, 0);//TODO: len pre s co dava zmysel
    number.powerSums[1] = 1;
    number.gcdPowerSums.resize(n, 0);
    number.gcdPowerSums[1] = 1;
    std::vector<Scalar> powers; //TODO: global
    powers.reserve(number.lambda);
    for (const auto c: number.coprimes) {
        if (number.orders[c] != 0) {
            continue;
        }
        powers.clear();
        DoubleScalar power = 1;
        do {
            powers.push_back(power);
            power = (power * c) % n;//TODO: pre vsetky delitele potom zratat, a inverse je phi - 1
        } while (power != 1); //TODO: number.orders[power] == 0
        const Scalar order = powers.size() * number.orders[power];//TODO: power == 1
        const auto &orderNumber = numberCache[order];
        for (const auto o: orderNumber.divisors) {
            const auto divisor = order / o;
            auto &oNumber = numberCache[o];
            computeCoprimes(oNumber);
            DoubleScalar powerSum = 0;
            for (Scalar i = 0; i < order; i += divisor) {
                powerSum = powerSum + powers[i];//TODO: reuse powerSum suborbity
            }
            powerSum = powerSum % n;
            const auto gcdPowerSum = gcd(number, numberCache[powerSum]);//TODO: handle 0?
            for (const auto cop: oNumber.coprimes) {
                const auto div = cop * divisor;//TODO: rename
                const auto e = powers[div];
                number.orders[e] = o;
                number.powerSums[e] = powerSum;
                number.gcdPowerSums[e] = gcdPowerSum;
            }
        }
        if (isPowerOfPrime(number) && powers.size() == number.coprimes.size()) {
            for (std::size_t i = 1; i < powers.size(); ++i) {
                number.inverses[powers[i]] = powers[powers.size() - i];
            }
        }
    }
}

PROFILE void computeCoprimesSortedByOrder(Number &number) {
    if (!number.order2OrderIndex.empty()) {
        return;
    }
    computeOrders(number);
    number.order2OrderIndex.resize(number.lambda + 1, -1);
    const auto &possibleOrders = numberCache[number.lambda].divisors;
    number.orderIndex2CoprimesBegin.resize(possibleOrders.size() + 1);
    number.orderIndex2CoprimesBegin[0] = 0;
    number.orderIndex2CoprimesBegin[possibleOrders.size()] = number.coprimes.size();

    std::vector<std::size_t> knownOrderBegins;
    std::vector<std::size_t> nextKnownOrderBegins;
    knownOrderBegins.reserve(possibleOrders.size());
    nextKnownOrderBegins.reserve(possibleOrders.size());

    nextKnownOrderBegins.push_back(0);
    nextKnownOrderBegins.push_back(possibleOrders.size());

    while (knownOrderBegins.size() < possibleOrders.size()) {
        knownOrderBegins.swap(nextKnownOrderBegins);
        nextKnownOrderBegins.clear();

        auto begin = knownOrderBegins[0];
        nextKnownOrderBegins.push_back(begin);
        for (std::size_t i = 1; i < knownOrderBegins.size(); ++i) {
            const auto end = knownOrderBegins[i];
            const auto middle = (begin + end) / 2;
            if (middle != begin) {
                nextKnownOrderBegins.push_back(middle);
                const auto middleOrder = possibleOrders[middle];

                const auto &beginIterator = number.coprimes.begin() + number.orderIndex2CoprimesBegin[begin];
                const auto &endIterator = number.coprimes.begin() + number.orderIndex2CoprimesBegin[end];
                const auto &middleIterator = std::partition(beginIterator, endIterator,
                                                            [middleOrder, &number](const auto &item) {
                                                                return number.orders[item] < middleOrder;
                                                            });

                number.orderIndex2CoprimesBegin[middle] = middleIterator - number.coprimes.begin();//TODO: store iterators
            }
            nextKnownOrderBegins.push_back(end);
            begin = end;
        }
    }
    for (Scalar i = 0; i < possibleOrders.size(); ++i) {
        number.order2OrderIndex[possibleOrders[i]] = i;
    }
}

Scalar isPQ(const Number &number) {//TODO: Scalar -> int, aj isAB
    if (number.squareFree && number.powerOfTwo <= 2 && number.primes.size() == 2) {
        return (number.primes[0] - 1) * (number.primes[1] - 1);
    }
    return 0;
}

void countSkewmorphisms(Number &number);

Scalar isAB(const Number &number) {
    for (const auto divisor: number.divisors) {
        auto &a = numberCache[divisor];
        auto &b = numberCache[number.n / a.n];
        if (a.n == 1) {
            continue;
        }
        if (a.n > b.n) {
            break;
        }
        if (gcd(a, b) == 1 && gcd(a, numberCache[b.phi]) == 1 && gcd(b, numberCache[a.phi]) == 1) {//TODO: is coprime
            countSkewmorphisms(a);
            countSkewmorphisms(b);
            return a.nskew * b.nskew;
        }
    }
    return 0;
}

PROFILE Scalar computeHelpSumsOrder(Scalar s, const Number &number_n_h) {
    const auto n_h = number_n_h.n;
    const auto s_mod = s % n_h;
    const auto orderSumsCount = number_n_h.orders[s_mod];
    const auto h_value = number_n_h.gcdPowerSums[s_mod];
    if (h_value == n_h) {//TODO:
        return orderSumsCount;
    } else {
        return orderSumsCount * (n_h / h_value);
    }
}

PROFILE Scalar scitaj(Scalar d, Scalar e, const Number &number_r, Scalar rH, Scalar power_sum, Scalar n_h) {
    return (((number_r.powerSums[e] - d) / rH + n_h) * static_cast<DoubleScalar>(power_sum) + d) % n_h;
}

PROFILE Scalar countCoprimeSolutions(Scalar a, const Number &number_n_h, Scalar gcd_b) {//TODO: priamo cez rozklad
    // ax = b (%n_h)
    const auto gcd_a = gcd(number_n_h, numberCache[a]);
    if (gcd_b != gcd_a) {
        return 0;
    }
    auto &number = numberCache[number_n_h.n / gcd_a];
    const auto &number_ga = numberCache[gcd_a];
    const auto gaa = gcd(number, number_ga);
    const auto &number_gaa = numberCache[gaa];
    return gaa * number_ga.phi / number_gaa.phi;
}

PROFILE bool remainder1(Scalar a, Scalar b) {
    return a % b != 1;
}

PROFILE Scalar sEquals1(const Scalar d, const Number &number_n_div_d) {
    Scalar nskew = 0;
    for (const auto n_h : number_n_div_d.divisors) {
        if (n_h == 1) {//TODO: range from second / except last
            continue;
        }
        auto &number_n_h = numberCache[n_h];
        if (d > number_n_h.lambda) {
            continue;
        }
        auto &number_r = number_n_h;
        computeCoprimesSortedByOrder(number_r);
        const auto order_index = number_r.order2OrderIndex[d];//TODO:
        if (order_index == -1) {
            continue;
        }
        const auto begin = number_r.orderIndex2CoprimesBegin[order_index];//TODO: len pre d ktore davaju zmysel?
        const auto end = number_r.orderIndex2CoprimesBegin[order_index + 1];
        for (std::size_t i = begin; i < end; ++i) {
            const auto e = number_r.coprimes[i];
            if (number_r.powerSums[e] == 0) {//TODO: toto je po castiach rovnake
                nskew += number_n_h.phi;
            }
        }
    }
    return nskew;
}

PROFILE Scalar sOtherThan1(const Scalar d, const Number &number_n_div_d, const Scalar s) {
    Scalar nskew = 0;
    const auto power_sum = number_n_div_d.powerSums[s];
    const auto gcd_power_sum = number_n_div_d.gcdPowerSums[s];
    const Scalar rH = number_n_div_d.orders[s];//TODO: sorted by order?
    bool bComputed = false;
    Scalar b = s - 1;
    Scalar gcd_b;
    for (const auto n_h : number_n_div_d.divisors) {
        if (n_h == 1) {
            continue;
        }
        if (gcd_power_sum % n_h == 0) {
            continue;
        }
        const auto small_gcd_n_h = number_n_div_d.n / n_h;
        if (b % small_gcd_n_h != 0) {
            continue;
        }
        auto &number_n_h = numberCache[n_h];
        computeOrders(number_n_h);
        const auto r = computeHelpSumsOrder(s, number_n_h);
//        if (r < d * rH) {//TODO: toto nie je blbost?
//            continue;
//        }

//if (r > n) {
//throw "";
//}
        auto &number_r = numberCache[r];//TODO: nemoze byt vacsie ako N?
        if (d > number_r.lambda) {
            continue;
        }
        computeCoprimesSortedByOrder(number_r);
        const auto order_index = number_r.order2OrderIndex[d];
        if (order_index == -1) {
            continue;
        }
        const auto begin = number_r.orderIndex2CoprimesBegin[order_index];//TODO: len pre d ktore davaju zmysel?
        const auto end = number_r.orderIndex2CoprimesBegin[order_index + 1];
        for (std::size_t i = begin; i < end; ++i) {
            const auto e = number_r.coprimes[i];//TODO: iterate e % rH == 1
            if (remainder1(e, rH)) {
                continue;
            }
            if (!bComputed) {
                gcd_b = gcd(number_n_div_d, numberCache[b]);
                bComputed = true;
            }
            const auto big_vysledok = scitaj(d, e, number_r, rH, power_sum, n_h);
            nskew += countCoprimeSolutions(big_vysledok, number_n_h, gcd_b / small_gcd_n_h);//TODO: bez nutnosti delenia
        }
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
        if (const auto p_1_q_1 = isPQ(number)) {
            nskew += p_1_q_1;
        } else if (const auto a_b = isAB(number)) {
            nskew = a_b;
        } else {
            const auto maxPrime = getMaxPrime(number);
            const auto &number_n_div_maxPrime = numberCache[n / maxPrime];
            for (const auto d: number_n_div_maxPrime.divisors) {
                if (d == 1) {
                    continue;
                }
                const auto n_div_d = n / d;
                auto &number_n_div_d = numberCache[n_div_d];
                computeCoprimesSortedByOrder(number_n_div_d);

                nskew += sEquals1(d, number_n_div_d);
                for (const auto s: number_n_div_d.coprimes) {
                    if (s == 1) {
                        continue;
                    }
                    nskew += sOtherThan1(d, number_n_div_d, s);
                }
            }
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
    if (!number.coprimes.empty()) {
        number.orders = std::vector<Scalar>{};
        number.coprimes = std::vector<Scalar>{};
        number.order2OrderIndex = std::vector<Scalar>{};
        number.orderIndex2CoprimesBegin = std::vector<Scalar>{};//TODO: neratame nieco viackrat? jasne, ze hej, mozme si predratat reference counter?
        number.inverses = std::vector<Scalar>{};
        number.powerSums = std::vector<Scalar>{};
        number.gcdPowerSums = std::vector<Scalar>{};
        toCountWithMultiples.erase(number.n);
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

    while (!toCountWithMultiples.empty()) {
        const auto p = *toCountWithMultiples.rbegin();
        const auto lowerMultiplier = (A + p - 1) / p;
        const auto upperMultiplier = N / p;
        for (auto multiplier = upperMultiplier; multiplier >= lowerMultiplier; --multiplier) {
            const auto n = multiplier * p;

            auto &number = numberCache[n];
            if (number.powerOfTwo <= 16 && number.squareFree) {
                countSkewmorphisms(number);
//                fprint(number, std::cerr);
            }

            clear(number);
        }
//        for (auto multiplier = 1; multiplier < lowerMultiplier; ++multiplier) {
//            const auto n = multiplier * p;
//            auto &number = numberCache[n];
//            clear(number);
//        }
        toCountWithMultiples.erase(p);
    }
    for (auto n = A; n <= N; ++n) {
        const auto &number = numberCache[n];
        if (number.nskew != 0) {
            print(number);
        }
    }

    return 0;
}
