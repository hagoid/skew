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

struct Number {//TODO: Cn
    std::vector<Scalar> primes;//TODO: gcd ratat z faktorizacie
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

//std::vector<Number> numberCache = std::vector<Number>(N_1);
//std::array<Number, N_1> numberCache = {};
Number numberCache[N_1] = {};

void computePrimes() {
    for (Scalar n = 2; n <= N; ++n) {
        bool isPrime = true;
        for (const auto p: primes) {
            if (p * p > n) {
                break;
            }
            if (n % p == 0) {
                isPrime = false;
                break;
            }
        }
        if (isPrime) {
            primes.push_back(n);
        }
    }
}


PROFILE Scalar gcdCompute(Scalar n, Scalar k) {
    if (n < N_1) {
        const auto &number = numberCache[n];
        if (number.squareFree && k < N_1) {
            if (k == 0) {
                return n;
            }
            const auto &kumber = numberCache[k];
            std::size_t ni = 0;
            std::size_t ki = 0;
            const auto &nrimes = number.primes;
            const auto &krimes = kumber.primes;
            const auto nize = nrimes.size();
            const auto kize = krimes.size();
            Scalar result = std::min(number.powerOfTwo, kumber.powerOfTwo);
            if (result > 1) {
                ++ni;
                ++ki;
            }
            while (ni < nize && ki < kize) {
                if (nrimes[ni] < krimes[ki]) {
                    ++ni;
                } else if (nrimes[ni] > krimes[ki]) {
                    ++ki;
                } else {
                    result *= nrimes[ni];
                    ++ki;
                    ++ni;
                }
            }
            return result;
        }
        const auto &divisors = number.divisors;
        const auto size = divisors.size();
        for (std::size_t i = 0; i < size; ++i) {
            const auto divisor = divisors[size - i - 1];
            if (k % divisor == 0) {
                return divisor;
            }
        }
    }
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

void computeOrders(Number &number);

Scalar powerOfTwo(Scalar n) {
    return (n & (~(n - 1)));
}

PROFILE Number &getNumber(Scalar n) {//TODO: const
    if (n > N) {
        throw "";
    }
    auto &number = numberCache[n];
    computeOrders(number);
    return number;
}

PROFILE void factorize(Number& number) {
    if (number.powerOfTwo != 0) {
        return;
    }
    auto n = number.n;
    number.powerOfTwo = powerOfTwo(n);
    Scalar divisorsCount = 1;
    if (n > 1) {
        for (const auto p: primes) {
            Scalar power = 1;
            if (p * p > n) {
                number.primes.push_back(n);
                n = 1;
                ++power;
            }
            if (n % p == 0) {
                number.primes.push_back(p);
                n /= p;
                ++power;
            }
            while (n % p == 0) {
                if (p != 2) {
                    number.squareFree = false;
                }
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
    for (const auto p: number.primes) {
        number.phi -= number.phi / p;
    }
}

Scalar getMaxPrime(const Number &number) {
    if (number.primes.empty()) {
        return 1;
    }
    return number.primes.back();
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
        bool coprime = true;
        for (const auto p : number.primes) {//TOOO: isCoprime, teoreticky even faster cez faktorizaciu oboch? takisto gcd faster?
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
    if (number.squareFree && number.powerOfTwo <= 2 && number.primes.size() == 2) {
        return (number.primes[0] - 1) * (number.primes[1] - 1);
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

PROFILE Scalar scitaj(Scalar d, Scalar e, const Number &number_r, Scalar rH, Scalar power_sum, Scalar n, const Number &number_n_h) {
    return (((number_r.powerSums[e] - d) / rH + n) * static_cast<DoubleScalar>(power_sum) + d) % number_n_h.n;
}

PROFILE Scalar countCoprimeSolutions(Scalar a, Scalar n_div_d, const Number &number_n_h, Scalar b, Scalar gcd_b) {
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
            const auto &number_n_div_maxPrime = numberCache[n / maxPrime];
            for (const auto d: number_n_div_maxPrime.divisors) {
                if (d == 1 || d == n) {
                    continue;
                }
                const auto n_div_d = n / d;
                auto &number_n_div_d = numberCache[n_div_d];
                reorganizeCoprimes(number_n_div_d);//TODO: rename
                for (const auto s: number_n_div_d.coprimes) {
                    const auto power_sum = number_n_div_d.powerSums[s];
                    const auto gcd_power_sum = gcd(n_div_d, power_sum);//TODO: gcd s powersumom vyuzivame
                    const Scalar rH = number_n_div_d.orders[s];//TODO: reorganized?
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
                        const auto &number_n_h = getNumber(n_h);
                        const auto r = computeHelpSumsOrder(s, number_n_h);
                        if (r < d * rH) {
                            continue;
                        }
                        auto &number_r = numberCache[r];
                        if (d > number_r.phi) {//TODO: lambda?
                            continue;
                        }
                        reorganizeCoprimes(number_r);
                        const auto order_index = number_r.order2OrderIndex[d];
                        if (order_index == -1) {
                            continue;
                        }
                        const auto begin = number_r.orderIndex2CoprimesBegin[order_index];//TODO: len pre d ktore davaju zmysel?
                        const auto end = number_r.orderIndex2CoprimesBegin[order_index + 1];
                        const auto small_gcd_n_h = n_div_d / n_h;
                        for (std::size_t i = begin; i < end; ++i) {
                            const auto e = number_r.coprimes[i];
                            if (divisible3(e - 1, rH)) {
                                continue;
                            }
                            if (!bComputed) {
                                gcd_b = gcd(n_div_d, b);//TODO: number
                                bComputed = true;
                            }
                            const auto big_vysledok = scitaj(d, e, number_r, rH, power_sum, n, number_n_h) * small_gcd_n_h;
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
