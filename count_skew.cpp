#include <iostream>
#include <cstdio>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <set>
#include <boost/multiprecision/cpp_int.hpp>

//#define PROFILE __attribute__((noinline))
#define PROFILE

using Scalar = std::int32_t;
using DoubleScalar = std::int64_t;
using QuadScalar = boost::multiprecision::int128_t;

constexpr Scalar N = 100000;
constexpr Scalar N_1 = N + 1;
//constexpr Scalar N_1_2 = N_1 * N_1;

std::vector<Scalar> primes = {};

//std::vector<Scalar> gcdCache = std::vector<Scalar>(N_1_2, 0);
//std::array<Scalar, N_1_2> gcdCache = {};
//Scalar gcdCache[N_1_2] = {};
Scalar orderSumsCount = 0;

void computePrimes() {
    for (Scalar i = 2; i <= N; ++i) {
        bool isPrime = true;
        for (const auto p: primes) {
            if (i % p == 0) {
                isPrime = false;
                break;
            }
        }
        if (isPrime) {
            primes.push_back(i);
        }
    }
}


PROFILE Scalar gcdCompute(Scalar n, Scalar k) {
    while (k != 0) {
        const Scalar new_k = n % k;
        n = k;
        k = new_k;
    }
    return n;
}

Scalar gcd(Scalar n, Scalar k) {
    return gcdCompute(n, k);
//    const auto index = n * N_1 + k;
//    if (gcdCache[index] == 0) {
//        if (k == 0) {
//            gcdCache[index] = n;
//        } else {
//            gcdCache[index] = gcd(k, n % k);
//        }
//    }
//    return gcdCache[index];
}


struct Number {//TODO: Cn
    std::vector<Scalar> oddFactors;//TODO: gcd ratat z faktorizacie
    std::vector<Scalar> orders;// TODO: s v deliteli n
    std::vector<Scalar> order2OrderIndex;
    std::vector<Scalar> orderIndex2CoprimesBegin;
    std::vector<Scalar> divisors;
    std::vector<Scalar> coprimes;
    std::vector<Scalar> inverses;
    std::vector<Scalar> powerSums;
    Scalar n;
    Scalar powerOfTwo = 0;
    Scalar phi = 0;
    Scalar nskew = 0;
    bool squareFree = true;
};

void computeOrders(Number &number);

Scalar powerOfTwo(Scalar n) {
    return (n & (~(n - 1)));
}

//std::vector<Number> numberCache = std::vector<Number>(N_1);
//std::array<Number, N_1> numberCache = {};
Number numberCache[N_1] = {};

PROFILE Number &getNumber(Scalar n) {//TODO: const
    if (n > N) {
        throw "";
    }
    auto &number = numberCache[n];
    computeOrders(number);
    return number;
}

void factorize(Number& number) {
    if (number.powerOfTwo != 0) {
        return;
    }
    auto n = number.n;
    number.powerOfTwo = powerOfTwo(n);
    Scalar divisorsCount = 1;
    auto powerOfTwo = number.powerOfTwo;
    while (powerOfTwo >>= 1)
    {
        divisorsCount++;
    }
    if (n != 0) {
        n /= number.powerOfTwo;
        for (const auto p: primes) {
            Scalar power = 1;
            if (n % p == 0) {
                number.oddFactors.push_back(p);
                n /= p;
                ++power;
            }
            while (n % p == 0) {
                number.squareFree = false;
                n /= p;
                ++power;
            }
            divisorsCount *= power;
            if (n == 1) {
                break;
            }
        }
    }
    number.divisors.reserve(divisorsCount);
    for (Scalar i = 1; i <= number.n; ++i) {
        if (number.n % i == 0) {
            number.divisors.push_back(i);
        }
    }
}

void computePhi(Number &number) {
    if (number.phi != 0) {
        return;
    }
    factorize(number);
    number.phi = number.n;
    if (number.powerOfTwo > 1) {
        number.phi -= number.phi / 2;
    }
    for (const auto p: number.oddFactors) {
        number.phi -= number.phi / p;
    }
}

Scalar getMaxPrime(const Number &number) {
    if (number.oddFactors.empty()) {
        if (number.powerOfTwo == 1) {
            return 1;
        }
        return 2;
    }
    return number.oddFactors.back();
}

std::set<int> toCountWithMultiples = {};

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
        return;
    }
    toCountWithMultiples.insert(number.n);
//    if (++number.toCountWithMultiples > 1) {
//        printf("toCountWithMultiples %d %d\n", number.n, number.toCountWithMultiples);
//    }
    for (Scalar e = 1; e <= n; ++e) {//TODO: computeDivisors, computeCoprimes, computeOrders
        if (number.powerOfTwo > 1) {
            if (e % 2 == 0) {//TOOO: isCoprime, teoreticky even faster cez faktorizaciu oboch? takisto gcd faster?
                continue;
            }
        }
        bool coprime = true;
        for (const auto p : number.oddFactors) {
            if (e % p == 0) {
                coprime = false;
                break;
            }
        }
        if (coprime) {
            number.coprimes.push_back(e);
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
        return;
    }
    number.orders.resize(n, 0);
    number.orders[1] = 1;
    number.inverses.resize(n, 0);
    number.inverses[1] = 1;
    number.powerSums.resize(n, 0);//TODO: len pre s co dava zmysel
    number.powerSums[1] = 1;//TODO: coprimes 1 je 1
    std::vector<Scalar> powers; //TODO: global
    powers.reserve(number.phi);
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
        const Scalar order = powers.size() * number.orders[power];
        const auto &orderNumber = numberCache[order];
        for (const auto o: orderNumber.divisors) {
            const auto divisor = order / o;
            auto &oNumber = numberCache[o];
            computeCoprimes(oNumber);
            Scalar powerSum = 0;
            for (Scalar i = 0; i < order; i += divisor) {
                powerSum = (powerSum + powers[i]) % n;
            }
            for (const auto cop: oNumber.coprimes) {
                const auto div = cop * divisor;//TODO: rename
                const auto e = powers[div];
                number.orders[e] = o;
                number.powerSums[e] = powerSum;
            }
        }
        for (std::size_t i = 1; i < powers.size(); ++i) {
            number.inverses[powers[i]] = powers[powers.size() - i];
        }
    }
}

PROFILE void reorganizeCoprimes(Number &number) {
    if (!number.order2OrderIndex.empty()) {
        return;
    }
    computeOrders(number);
    number.order2OrderIndex.resize(number.phi + 1, -1);
    const auto& number_phi = numberCache[number.phi];
    number.orderIndex2CoprimesBegin.reserve(number_phi.divisors.size() + 1);//TODO: nemusia tam byt vsetky
    number.orderIndex2CoprimesBegin.push_back(0);
    std::vector<Scalar> coprimes; //TODO: global
    coprimes.reserve(number.phi);
    for (const auto o: number_phi.divisors) {
        number.order2OrderIndex[o] = number.orderIndex2CoprimesBegin.size() - 1;
        for (const auto coprime: number.coprimes) {
            if (number.orders[coprime] == o) {//TODO: N^2
                coprimes.push_back(coprime);
            }
        }
        number.orderIndex2CoprimesBegin.push_back(coprimes.size());
    }
    coprimes.swap(number.coprimes);
}

Scalar isPQ(const Number &number) {//TODO: Scalar -> int, aj isAB
    if (number.powerOfTwo == 1 && number.oddFactors.size() == 2) {
        return (number.oddFactors[0] - 1) * (number.oddFactors[1] - 1);
    } else if (number.powerOfTwo == 2 && number.oddFactors.size() == 1) {
        return number.oddFactors[0] - 1;
    }
    return 0;
}

Scalar count(const Number &number);

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
        if (gcd(a.n, b.n) == 1 && gcd(a.n, b.phi) == 1 && gcd(a.phi, b.n) == 1) {//TODO: is coprime
            if (a.nskew == 0) {
                a.nskew = count(a);//TODO:
            }
            if (b.nskew == 0) {
                b.nskew = count(b);
            }
            return a.nskew * b.nskew;
        }
    }
    return 0;
}

DoubleScalar pow(QuadScalar x, unsigned int y, DoubleScalar p)
{
    QuadScalar res = 1;     // Initialize result

    x = x % p; // Update x if it is more than or
    // equal to p

    if (x == 0) return 0; // In case x is divisible by p;

    while (y > 0) {
        // If y is odd, multiply x with result
        if (y & 1)
            res = (res * x) % p;

        // y must be even now
        y = y >> 1; // y = y/2
        x = (x * x) % p;
    }
    return res.convert_to<DoubleScalar>();
}

PROFILE bool computeHelpSums(const Number &number, Scalar possible_d, Scalar max_d) {//TODO: niektore s nepotrebujeme predratat, napr 1
    const auto n = number.n;
    bool atLeastOne = possible_d != 1;
    if (!atLeastOne) {
        if (n % 2 == 0) {
            atLeastOne = true;
        }
        for (const auto &p: number.oddFactors) {//TODO: all factors
            if (p > max_d) {
                break;
            }
            atLeastOne = true;
        }
    }
    return atLeastOne;
}

PROFILE Scalar computeHelpSumsOrder(Scalar s, const Number &number_n_h) {
//    orderSumsCount = order(n / gcd_n_h, s);//TODO
    const auto n_h = number_n_h.n;
    const auto s_mod = s % n_h;
    orderSumsCount = number_n_h.orders[s_mod];
    const auto h_value = number_n_h.powerSums[s_mod];
    if (h_value == 0) {
        return orderSumsCount;
    } else {
        return orderSumsCount * (n_h / gcd(n_h, h_value));
    }
}

PROFILE Scalar scitaj(Scalar d, Scalar e, Scalar n_div_d, Scalar r, Scalar s, Scalar n, const Number &number_n_h) {
    Scalar mocnina = 1;
    Scalar vysledok = 0;
    const auto h_value = number_n_h.powerSums[s % number_n_h.n];
    for (Scalar i = 0; i < d; i++) {
        const auto increment = s == 1 ? mocnina : ((pow(s, mocnina % orderSumsCount, static_cast<DoubleScalar>(n) * (s - 1)) - 1) / (s - 1) +
                    ((h_value * static_cast<DoubleScalar>(mocnina / orderSumsCount)))) % n_div_d;//TODO
        vysledok += increment;
        mocnina = (mocnina * e) % r;
    }
    return vysledok;
}

PROFILE bool overScitane(DoubleScalar big_vysledok, Scalar s, Scalar n_div_d, Scalar small_small_h) {
    return ((s - 1) - big_vysledok * small_small_h) % n_div_d == 0;//TODO: gcd(s-1, n_div_d)
}

PROFILE Scalar countCoprimeSolutions(DoubleScalar big_vysledok, Scalar n_div_d, const Number &number_n_h, Scalar b, Scalar gcd_b) {
    const auto a = big_vysledok % n_div_d; // TODO: reuse
    // ax = b (%n_div_d)
    const auto gcd_a = gcd(n_div_d, a);
    if (gcd_b % gcd_a != 0) {
        return 0;
    }
    const auto& number = getNumber(n_div_d / gcd_a);//TODO:
    const auto x = b == 0 ? 0 : (number.inverses[a / gcd_a] * b / gcd_a) % number.n;
    Scalar nskew = 0;
    for (Scalar sol = x; sol < number_n_h.n; sol += number.n) {
        if (number_n_h.orders[sol] != 0) {
            ++nskew;
        }
    }
    return nskew;
}

PROFILE bool divisible3(Scalar a, Scalar b) {
    return a % b != 0;
}

std::vector<Scalar> possible_ds = {};

PROFILE Scalar countForE(DoubleScalar big_vysledok, Scalar s, Scalar n_div_d, const Number &number_n_h) {
    Scalar nskew = 0;
    for (const auto small_small_h: number_n_h.coprimes) {
        if (overScitane(big_vysledok, s, n_div_d, small_small_h)) {
            nskew++;//TODO: neda sa toto rychlejsie?
        }
    }
    return nskew;
}

Scalar count(const Number &number) {
    const auto n = number.n;
    const auto phi = number.phi;
    auto nskew = phi;

    if (gcd(n, phi) > 1) {
        if (const auto p_1_q_1 = isPQ(number)) {
            nskew += p_1_q_1;
        } else if (const auto a_b = isAB(number)) {
            nskew = a_b;
        } else {
            const auto maxPrime = getMaxPrime(number);
            for (Scalar s = 1; s < (n + 1) / 2; ++s) {
                const auto gcd_n_s = gcd(n, s);

                auto possible_d = gcd_n_s;

                if (possible_d % maxPrime == 0) {
                    continue;
                }

                if (possible_d % 2 == 0) {
                    possible_d *= powerOfTwo(n) / powerOfTwo(possible_d);
                }

                const auto max_d = static_cast<Scalar>(std::ceil(static_cast<float>(n) / static_cast<float>(s))) - 1;
                if (possible_d > max_d) {
                    continue;
                }

                if (!computeHelpSums(number, possible_d, max_d)) {//TODO: strasne vela checkov tu je, ale my uz iterujeme inak, tak asi az tak netreba
                    continue;
                }

                const auto &number_n_div_possible_d = numberCache[n / maxPrime / possible_d];
                possible_ds.clear();
                for (const auto small_d: number_n_div_possible_d.divisors) {
                    const auto d = small_d * possible_d;
                    if (d == 1) {
                        continue;
                    }
                    if (d > max_d) {
                        break;
                    }
                    possible_ds.push_back(d);
                }
                for (const auto d: possible_ds) {
                    const auto n_div_d = n / d;
//                    const Scalar rH = order(n_div_d, s);//TODO:
                    const auto &number_n_div_d = getNumber(n_div_d);
                    const Scalar rH = number_n_div_d.orders[s];
                    bool bComputed = false;
                    Scalar b;
                    Scalar gcd_b;
                    for (const auto small_gcd_n_h : number_n_div_d.divisors) {
                        if (small_gcd_n_h == number_n_div_d.n) {
                            break;
                        }
                        const auto gcd_n_h = small_gcd_n_h * d;
                        const auto &number_n_h = getNumber(n / gcd_n_h);
                        const auto r = computeHelpSumsOrder(s, number_n_h);
                        auto &number_r = numberCache[r];
                        if (d > number_r.phi) {
                            continue;
                        }
                        reorganizeCoprimes(number_r);
                        const auto order_index = number_r.order2OrderIndex[d];
                        if (order_index == -1) {
                            continue;
                        }
                        const auto begin = number_r.orderIndex2CoprimesBegin[order_index];
                        const auto end = number_r.orderIndex2CoprimesBegin[order_index + 1];
                        for (std::size_t i = begin; i < end; ++i) {
                            const auto e = number_r.coprimes[i];
                            if (divisible3(e - 1, rH)) {
                                continue;
                            }
                            if (!bComputed) {
                                b = (s - 1) % n_div_d;
                                gcd_b = gcd(n_div_d, b);
                                bComputed = true;
                            }
                            const auto big_vysledok = scitaj(d, e, n_div_d, r, s, n, number_n_h) * static_cast<DoubleScalar>(small_gcd_n_h);
                            nskew += countCoprimeSolutions(big_vysledok, n_div_d, number_n_h, b, gcd_b);
                        }
                    }
                }
            }
        }
    }
    return nskew;
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
        toCountWithMultiples.erase(number.n);
    }
}

int main() {
    computePrimes();
    possible_ds.reserve(N);
    for (Scalar i = 0; i <= N; ++i) {
        auto &number = numberCache[i];
        number.n = i;
        computePhi(number);//TODO: lazy ak chceme pouzivat A, B
    }
    const Scalar A = 1;
    const Scalar B = N;

    toCountWithMultiples.insert(1);
    for (const auto p : primes) {
        toCountWithMultiples.insert(p);
    }

    while (!toCountWithMultiples.empty()) {
        const auto p = *toCountWithMultiples.rbegin();
        const auto lowerMultiplier = (A + p - 1) / p;
        const auto upperMultiplier = B / p;
        for (auto multiplier = upperMultiplier; multiplier >= lowerMultiplier; --multiplier) {
            const auto n = multiplier * p;

            auto &number = numberCache[n];
            if (number.powerOfTwo <= 16 && number.squareFree && number.nskew == 0) {
                number.nskew = count(number);//TODO: tento assignment je cudny
                fprint(number, std::cerr);
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
    for (auto n = A; n <= B; ++n) {
        const auto &number = numberCache[n];
        if (number.nskew != 0) {
            print(number);
        }
    }

    return 0;
}
