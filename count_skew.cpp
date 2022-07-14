#include <iostream>
#include <cstdio>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <set>

//#define PROFILE __attribute__((noinline))
#define PROFILE

constexpr int N = 10000;
constexpr int N_1 = N + 1;
//constexpr int N_1_2 = N_1 * N_1;

std::vector<int> primes = {};

//std::vector<int> gcdCache = std::vector<int>(N_1_2, 0);
//std::vector<int> helpSums = std::vector<int>(N_1, 0);
//std::array<int, N_1_2> gcdCache = {};
//std::array<int, N_1> helpSums = {};
//int gcdCache[N_1_2] = {};
int helpSums[N_1] = {};
int helpSumsCount = 0;
int orderSumsCount = 0;

void computePrimes() {
    for (int i = 2; i < N; ++i) {
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


PROFILE int gcdCompute(int n, int k) {
    while (k != 0) {
        const int new_k = n % k;
        n = k;
        k = new_k;
    }
    return n;
}

int gcdPrecompute(int n, int k) {
    return gcdCompute(n, k);
//    const auto index = n * N_1 + k;
//    if (gcdCache[index] == 0) {
//        if (k == 0) {
//            gcdCache[index] = n;
//        } else {
//            gcdCache[index] = gcdPrecompute(k, n % k);
//        }
//    }
//    return gcdCache[index];
}

int gcd(int n, int k) {
//    return gcdCache[n * N_1 + k];
    return gcdPrecompute(n, k);
}

struct CoprimesByOrder {
    std::vector<int> coprimes;
    int order;
};

struct Number {
    std::vector<int> oddFactors;//TODO: gcd ratat z faktorizacie
    std::vector<int> orders;// TODO: s v deliteli n
    std::vector<int> orders_index;
    std::vector<CoprimesByOrder> coprimesByOrder;
    std::vector<int> divisors;
    std::vector<int> coprimes;
    int n;
    int powerOfTwo;
    int phi;
    int nskew;
    bool squareFree;
};

int powerOfTwo(int n) {
    return (n & (~(n - 1)));
}

//std::vector<Number> numberCache = std::vector<Number>(N_1);
//std::array<Number, N_1> numberCache = {};
Number numberCache[N_1] = {};

Number &getNumber(int n) {
    return numberCache[n];
}

Number factorize(int n) {
    Number number = {.n = n, .powerOfTwo = powerOfTwo(n), .phi = 0, .nskew = 0, .squareFree = true};
    if (n != 0) {
        n /= number.powerOfTwo;
        for (const auto p: primes) {
            if (n % p == 0) {
                number.oddFactors.push_back(p);
                n /= p;
            }
            while (n % p == 0) {
                number.squareFree = false;
                n /= p;
            }
            if (n == 1) {
                break;
            }
        }
    }
    return number;
}

void computePhi(Number &number) {
    number.phi = number.n;
    if (number.powerOfTwo > 1) {
        number.phi -= number.phi / 2;
    }
    for (const auto p: number.oddFactors) {
        number.phi -= number.phi / p;
    }
}


int findMaxPrime(const Number &number) {
    if (number.oddFactors.empty()) {
        if (number.powerOfTwo == 1) {
            return 1;
        }
        return 2;
    }
    return number.oddFactors.back();
}

void computeOrders(Number &number) {
    const auto n = number.n;
    if (n == 0) {
        return;
    }
    computePhi(number);
    number.coprimes.reserve(number.phi);
    if (n == 1) {
        return;
    }
    for (int e = 1; e <= n; ++e) {
        const auto d = gcdPrecompute(n, e);
        if (d == e) {
            number.divisors.push_back(e);
        }
        if (d == 1) {
            number.coprimes.push_back(e);
        }
    }
    number.divisors.shrink_to_fit();//TODO: count?
    number.orders.resize(n, 0);
    number.orders[1] = 1;
    std::set<int> orders;
    for (const auto c: number.coprimes) {
        if (number.orders[c] != 0) {
            continue;
        }
        std::vector<int> powers;
        powers.reserve(number.phi);
        int power = 1;
        do {
            powers.push_back(power);
            power = (power * c) % n;
        } while (power != 1); //TODO: number.orders[power] == 0
        const int order = powers.size() * number.orders[power];
        const auto &orderNumber = getNumber(order);
        for (const auto divisor: orderNumber.divisors) {
            const auto o = order / divisor;
            const auto &oNumber = getNumber(o);
            orders.insert(o);
            for (const auto cop: oNumber.coprimes) {
                const auto div = cop * divisor;//TODO: rename
                if (div >= powers.size()) {
                    break;
                }
                number.orders[powers[div]] = o;
            }
        }
    }
    number.orders_index.resize(number.phi + 1, -1);
    for (const auto o: orders) {
        number.orders_index[o] = number.coprimesByOrder.size();
        number.coprimesByOrder.emplace_back();
        auto &coprimesByOrder = number.coprimesByOrder.back();
        coprimesByOrder.order = o;
        for (const auto coprime: number.coprimes) {
            if (number.orders[coprime] == o) {
                coprimesByOrder.coprimes.push_back(coprime);
            }
        }
    }
}

int isPQ(const Number &number) {
    if (number.powerOfTwo == 1 && number.oddFactors.size() == 2) {
        return (number.oddFactors[0] - 1) * (number.oddFactors[1] - 1);
    } else if (number.powerOfTwo == 2 && number.oddFactors.size() == 1) {
        return number.oddFactors[0] - 1;
    }
    return 0;
}

int isAB(const Number &number) {
    for (const auto divisor: number.divisors) {
        const auto &a = getNumber(divisor);
        const auto &b = getNumber(number.n / a.n);
        if (a.n == 1) {
            continue;
        }
        if (a.n > b.n) {
            break;
        }
        if (gcdPrecompute(a.n, b.n) == 1 && gcdPrecompute(a.n, b.phi) == 1 && gcdPrecompute(a.phi, b.n) == 1) {
            return a.nskew * b.nskew;
        }
    }
    return 0;
}

PROFILE bool computeHelpSums(int s, const Number &number, int possible_d,
                             int max_d) {//TODO: niektore s nepotrebujeme predratat, napr 1
    const auto n = number.n;
    const auto &number_n_d = getNumber(n / possible_d);
    helpSumsCount = number_n_d.orders[s];
    bool atLeastOne = possible_d != 1;
    if (!atLeastOne) {
        int maxHelpSumsCount = 0;
        if (n % 2 == 0) {
            const auto &number_n_p = getNumber(n / 2);
            maxHelpSumsCount = number_n_p.orders[s];
            atLeastOne = true;
        }
        for (const auto &p: number.oddFactors) {
            if (maxHelpSumsCount == helpSumsCount || p > max_d) {
                break;
            }
            const auto &number_n_p = getNumber(n / p);
            maxHelpSumsCount = std::max(maxHelpSumsCount, number_n_p.orders[s]);
            atLeastOne = true;
        }
        helpSumsCount = maxHelpSumsCount;
    }
    if (atLeastOne) {
        helpSums[0] = 0;
        for (int i = 1; i <= helpSumsCount; ++i) {
            helpSums[i] = (helpSums[i - 1] * s + 1) % n;//TODO: 1 na konci nie na zaciatku
        }
    }
    return atLeastOne;
}

PROFILE int computeHelpSumsOrder(int s, const Number &number_n_h) {//TODO: toto treba vobec robit?
//    orderSumsCount = order(n / gcd_n_h, s);//TODO
    orderSumsCount = number_n_h.orders[s % number_n_h.n];
    const auto h_value = helpSums[orderSumsCount] % number_n_h.n;
    if (h_value == 0) {
        return orderSumsCount;
    } else {
        return orderSumsCount * number_n_h.n / gcd(number_n_h.n, h_value);
    }
}

PROFILE int scitaj(int d, int e, int n_div_d, int r) {
    int mocnina = 1;
    int vysledok = 0;
    for (int i = 0; i < d; i++) {
        vysledok = (vysledok + helpSums[mocnina % orderSumsCount] +
                    (helpSums[orderSumsCount] * (mocnina / orderSumsCount))) % n_div_d;//TODO
        mocnina = (mocnina * e) % r;
    }
    return vysledok;
}

PROFILE bool overScitane(int big_vysledok, int s, int n_div_d, int small_small_h) {
    return ((s - 1) - big_vysledok * small_small_h) % n_div_d == 0;//TODO: gcd(s-1, n_div_d)
}

PROFILE bool divisible3(int a, int b) {
    return a % b != 0;
}

int count(const Number &number) {
    const auto n = number.n;
    const int phi = number.phi;
    auto nskew = phi;

    if (gcd(n, phi) > 1) {
        if (const auto p_1_q_1 = isPQ(number)) {
            nskew += p_1_q_1;
        } else if (const auto a_b = isAB(number)) {
            nskew = a_b;
        } else {
            const auto maxPrime = findMaxPrime(number);
            for (int s = 1; s < (n + 1) / 2; ++s) {
                const auto gcd_n_s = gcd(n, s);

                auto possible_d = gcd_n_s;

                if (possible_d % maxPrime == 0) {
                    continue;
                }

                if (possible_d % 2 == 0) {
                    possible_d *= powerOfTwo(n) / powerOfTwo(possible_d);
                }

                const auto max_d = static_cast<int>(std::ceil(static_cast<float>(n) / static_cast<float>(s))) - 1;
                if (possible_d > max_d) {
                    continue;
                }

                if (!computeHelpSums(s, number, possible_d,
                                     max_d)) {//TODO: strasne vela checkov tu je, ale my uz iterujeme inak, tak asi az tak netreba
                    continue;
                }

                const auto &number_n_div_possible_d = getNumber(n / maxPrime / possible_d);
                std::vector<int> possible_ds = {};
                for (const auto &small_d: number_n_div_possible_d.divisors) {
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
//                    const int rH = order(n_div_d, s);//TODO:
                    const auto &number_n_div_d = getNumber(n_div_d);
                    const int rH = number_n_div_d.orders[s];
//                    for (const auto& small_gcd_n_h : number_n_div_d.divisors) {
                    const auto size = number_n_div_d.divisors.size() - 1;
                    for (size_t iii = 0; iii < size; ++iii) {
                        const auto small_gcd_n_h = number_n_div_d.divisors[iii];
                        const auto gcd_n_h = small_gcd_n_h * d;
                        const auto &number_n_h = getNumber(n / gcd_n_h);
                        const auto r = computeHelpSumsOrder(s, number_n_h);
                        const auto &number_r = getNumber(r);
                        if (d > number_r.phi) {
                            continue;
                        }
                        const auto order_index = number_r.orders_index[d];
                        if (order_index == -1) {
                            continue;
                        }
                        const auto &coprimesByOrder = number_r.coprimesByOrder[order_index];
                        for (const auto e: coprimesByOrder.coprimes) {
                            if (divisible3(e - 1, rH)) {
                                continue;
                            }
                            const auto big_vysledok = scitaj(d, e, n_div_d, r) * small_gcd_n_h;
                            for (const auto small_small_h: number_n_h.coprimes) {
                                if (overScitane(big_vysledok, s, n_div_d, small_small_h)) {
                                    nskew++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    printf("Total number of skew morphisms of C_%d is %d\n", n, nskew);
    printf("Check sub-total of automorphisms of C_%d is %d\n", n, phi);
    int perc = static_cast<int>(std::round(float(100 * phi) / float(nskew)));
    printf("Automorphisms account for ~%d%% of all skew morphisms.\n\n", perc);
    return nskew;
}

int main() {
    computePrimes();
    for (int i = 0; i <= N; ++i) {
        numberCache[i] = factorize(i);
        computeOrders(numberCache[i]);
    }
    for (int n = 1; n <= N; n++) {
        auto &number = getNumber(n);
        if (number.powerOfTwo > 16 || !number.squareFree) continue;

        number.nskew = count(number);
    }

    return 0;
}
