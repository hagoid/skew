#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <ostream>
#include <set>
#include <sstream>
#include <string_view>
#include <unordered_set>
#include <vector>

#ifdef PROFILE_FLAG
#  define PROFILE __attribute__((noinline))
#else
#  define PROFILE
#endif

using Scalar = std::int32_t;
using DoubleScalar = std::int64_t;
using Index = Scalar;

using Orbit = std::vector<Scalar>;
using Orbits = std::vector<Orbit>;
using Function = std::vector<Scalar>;

struct OrbitPlace {
    Index orbitIndex = -1;
    Index indexOnOrbit = 0;
};
using OrbitPlaces = std::vector<OrbitPlace>;

struct Permutation {
    Orbits orbits;
    OrbitPlaces places;
};

struct SkewMorphism {
    Permutation permutation;
    Function pi;
    std::vector<Index> free_x;
    Scalar n = 0;
    Scalar d = 0;
    Scalar h = 0;
    Scalar r = 0;
    Scalar s = 0;
    Scalar c = 0;
    Scalar max_orbits = 0;
    bool inverseOrbit = false;
    bool powerOfInverseOrbit = false;
};

struct CompactSkewMorphism {
    Orbit orbit1;
    Function pi;
};

using Class = std::vector<Index>;
using Classes = std::vector<Class>;

PROFILE const Orbit& getOrbit1(const SkewMorphism &skewMorphism) {
    return skewMorphism.permutation.orbits[0];
}

PROFILE void computeMaxOrbits(SkewMorphism &skewMorphism) {//TODO: premenovat este to aj free_x rata a shrinkuje a inverseOrbituje
    auto &orbits = skewMorphism.permutation.orbits;
    if (orbits.empty()) {
        skewMorphism.max_orbits = skewMorphism.permutation.places.size();
        return;
    }

    skewMorphism.max_orbits = std::count_if(orbits.begin(), orbits.end(),
                                            [r = skewMorphism.r](const auto &orbit) { return orbit.size() == r; });//TODO: od zaciatku idu

    const auto &orbit1 = getOrbit1(skewMorphism);
    std::set<Index> set;
    for (std::size_t i = 1; i < orbit1.size(); ++i) {//TODO: preco ideme od 1?
        const auto index = orbit1[i] % skewMorphism.d;
        set.insert(index);
    }
    std::copy(set.begin(), set.end(), std::back_inserter(skewMorphism.free_x));

    const auto inverse1 = skewMorphism.n - 1;
    for (const auto o: orbit1) {
        if (o == inverse1) {
            skewMorphism.inverseOrbit = true;
        }
    }

    skewMorphism.permutation.orbits.shrink_to_fit();
    for (auto &orbit: skewMorphism.permutation.orbits) {
        orbit.shrink_to_fit();
    }
    skewMorphism.permutation.places.shrink_to_fit();
    skewMorphism.pi.shrink_to_fit();
    skewMorphism.free_x.shrink_to_fit();
}

PROFILE Scalar apply(const SkewMorphism &skewMorphism, Scalar n, Scalar a) {
    const auto &place = skewMorphism.permutation.places[a];
    if (place.orbitIndex == -1) {
        return a;
    }
    const auto &orbit = skewMorphism.permutation.orbits[place.orbitIndex];
    return orbit[(place.indexOnOrbit + n) % Scalar(orbit.size())];
}

PROFILE Scalar apply2(const SkewMorphism &skewMorphism, Scalar n, Scalar a) {//TODO: rename
    const auto &place = skewMorphism.permutation.places[a];
    if (place.orbitIndex == -1) {
        return a;
    }
    const auto &orbit = skewMorphism.permutation.orbits[place.orbitIndex];
    auto indexOnOrbit = place.indexOnOrbit + n;
    const Scalar orbitSize = orbit.size();
    if (indexOnOrbit >= orbitSize) indexOnOrbit -= orbitSize;
    return orbit[indexOnOrbit];
}

PROFILE Scalar apply1(const SkewMorphism &skewMorphism, Scalar a) {
    const auto &place = skewMorphism.permutation.places[a];
    if (place.orbitIndex == -1) {
        return a;
    }
    const auto &orbit = skewMorphism.permutation.orbits[place.orbitIndex];
    auto indexOnOrbit = place.indexOnOrbit + 1;
    const Scalar orbitSize = orbit.size();
    if (indexOnOrbit >= orbitSize) indexOnOrbit -= orbitSize;
    return orbit[indexOnOrbit];
}

struct SkewMorphisms {
    std::vector<SkewMorphism> automorphisms;
    std::vector<SkewMorphism> properCosetPreserving;
    std::vector<SkewMorphism> properNonPreserving;
    Classes autoClasses;
    Classes preservingClasses;
    Classes nonPreservingClasses;
};

struct Number {//TODO: Cn
    std::vector<Scalar> primes;
    std::vector<Scalar> powersOfPrimes;
    std::vector<Scalar> divisors;
    std::vector<Scalar> coprimes;
    SkewMorphisms skewMorphisms;
    Scalar n = 0;
    Scalar powerOfTwo = 0;
    Scalar phi = 0;
    Scalar lambda = 0;
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
//        clear(number_p);
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
            if (b < 0) throw "";
            for (const auto e: number_p_k.coprimes) {
                auto l = e + b;
                if (l >= p_k) l -= p_k;
                number.coprimes.push_back(l * n_p_k + a);
            }
        }
//        clear(number_n_p_k);
//        clear(number_p_k);
    }
}

PROFILE Scalar compute_a(Scalar d, Scalar power_sum_e, Scalar rH, DoubleScalar power_sum, Scalar n_h) {
    return (((power_sum_e - d) / rH + n_h) * power_sum + d) % n_h;//TODO: optimize
}

PROFILE std::string to_short_string(const Orbit& orbit1, const Function& pi) {
    std::string result;
    result += "(";
    bool first = true;
    for (const auto e: orbit1) {
        result += (first ? "" : ", ") + std::to_string(e);
        first = false;
    }
    result += ")";
    result += "[";
    first = true;
    for (const auto e: pi) {
        result += (first ? "" : ", ") + std::to_string(e);
        first = false;
    }
    result += "]";
    return result;
}

PROFILE std::string to_short_string(const SkewMorphism &skewMorphism) {
    if (skewMorphism.permutation.orbits.empty()) {
        return "(1)[0]";
    }
    return to_short_string(getOrbit1(skewMorphism), skewMorphism.pi);
}

PROFILE bool computeFunction(Scalar n, const Function &pi, const Orbit &orbit1, Function &function) {
    Scalar oldO = 1;
    for (std::size_t i = 1; i < orbit1.size(); ++i) {
        const auto o = orbit1[i];
        if (o != 0 && oldO != 0) {
            function[oldO] = -o;
        }
        oldO = o;
    }
    if (oldO != 0) {
        function[oldO] = -1;
    }

    for (std::size_t i1 = 0; i1 < n; i1 += pi.size()) {
        for (std::size_t i2 = 1; i2 <= pi.size(); ++i2) {
            const auto i = i1 + i2;
            if (i >= n) {
                return true;
            }
            auto value = function[i - 1] + orbit1[pi[i2 - 1]];
            if (value >= n) value -= n;

            if (function[i] < 0 && function[i] != -value) {
                for (const auto o: orbit1) {
                    function[o] = 0;
                }
                return false;
            }
            function[i] = value;
        }
    }
    return true;
}

PROFILE Function computeFunction(Scalar n, const Function &pi, const Orbit &orbit1) {
    Function function(n, 0);
    if (computeFunction(n, pi, orbit1, function)) {
        return function;
    }
    return {};
}

PROFILE void computeOrbits(SkewMorphism &skewMorphism, const Function &function) {//TODO: zjednotit
    const auto n = skewMorphism.n;
    skewMorphism.permutation.orbits.reserve(n);
    skewMorphism.permutation.places.resize(n);
    for (Scalar i = 0; i < n; ++i) {
        auto &orbitPlace = skewMorphism.permutation.places[i];
        if (orbitPlace.orbitIndex != -1) {
            continue;
        }
        auto c = function[i];
        if (c == i) {
            orbitPlace.indexOnOrbit = i;
            continue;
        }
        const auto orbitIndex = skewMorphism.permutation.orbits.size();
        skewMorphism.permutation.orbits.emplace_back();
        auto &orbit = skewMorphism.permutation.orbits.back();
        orbit.reserve(skewMorphism.r);
        orbitPlace.orbitIndex = orbitIndex;
        orbitPlace.indexOnOrbit = orbit.size();
        orbit.push_back(i);
        do {
            auto &orbitPlaceC = skewMorphism.permutation.places[c];
            orbitPlaceC.orbitIndex = orbitIndex;
            orbitPlaceC.indexOnOrbit = orbit.size();
            orbit.push_back(c);
            c = function[c];
        } while (c != i);
    }

    auto &orbits = skewMorphism.permutation.orbits;
    if (!orbits.empty()) {
        std::sort(orbits.begin() + 1, orbits.end(), [](const auto &lhs, const auto &rhs) {
            if (lhs.size() == rhs.size()) {
                return lhs[0] < rhs[0];
            }
            return lhs.size() > rhs.size();
        });
        for (Index i = 0; i < orbits.size(); ++i) {
            for (const auto e: orbits[i]) {
                skewMorphism.permutation.places[e].orbitIndex = i;
            }
        }
    }
}

PROFILE void finish(SkewMorphism &skewMorphism, const Orbit& orbit1) {
    const auto n = skewMorphism.n;
    const auto function = computeFunction(n, skewMorphism.pi, orbit1);
    computeOrbits(skewMorphism, function);
    computeMaxOrbits(skewMorphism);
}

PROFILE void addSkewMorphism(SkewMorphism skewMorphism, SkewMorphisms &skewMorphisms, bool newClass) {
    const auto add = [](Classes &classes, std::vector<SkewMorphism> &skews, SkewMorphism skew, bool newClass) {
        if (newClass) {
            classes.emplace_back();
        }
        classes.back().push_back(skews.size());
        skews.emplace_back(std::move(skew));
    };
    if (skewMorphism.d == 1) {
        add(skewMorphisms.autoClasses, skewMorphisms.automorphisms, std::move(skewMorphism), newClass);
    } else if (skewMorphism.h % skewMorphism.d == 0) {
        add(skewMorphisms.preservingClasses, skewMorphisms.properCosetPreserving, std::move(skewMorphism), newClass);
    } else {
        add(skewMorphisms.nonPreservingClasses, skewMorphisms.properNonPreserving, std::move(skewMorphism), newClass);
    }
}

template <typename T>
PROFILE std::size_t hash_vector(const std::vector<T> &vector) {
    std::string_view bytestring(reinterpret_cast<const char *>(vector.data()), vector.size() * sizeof(T));
    return std::hash<std::string_view>{}(bytestring);
}

template <typename T, std::size_t N>
PROFILE std::size_t hash_array(const std::array<T, N> &array) {
    std::string_view bytestring(reinterpret_cast<const char *>(array.data()), N * sizeof(T));
    return std::hash<std::string_view>{}(bytestring);
}

struct HashCompact {
    PROFILE std::size_t operator()(const CompactSkewMorphism &compactSkew) const {
        std::array<std::size_t, 2> sub_hash = {hash_vector(compactSkew.orbit1), hash_vector(compactSkew.pi)};
        return hash_array(sub_hash);
    }
};

struct EqualCompact {
    PROFILE bool operator()(const CompactSkewMorphism &lhs, const CompactSkewMorphism &rhs) const {
        return lhs.orbit1 == rhs.orbit1 && lhs.pi == rhs.pi;
    }
};

using SkewSet = std::unordered_set<CompactSkewMorphism, HashCompact, EqualCompact>;

void addSkewClassByRepresentant(SkewSet &foundSkews,
                                std::vector<std::vector<std::vector<std::size_t>>> &roots,
                                const SkewMorphism &ro,
                                const SkewMorphism &phi);

std::size_t getPreservingSkewCount(const SkewMorphisms &skewMorphisms);//TODO: remove

PROFILE void sEquals1(const Scalar d, const Number &number_n_div_d, SkewMorphisms &skewMorphisms) {
//    SkewSet foundSkews;
//    std::vector<std::vector<std::vector<std::size_t>>> roots;
    std::vector<bool> visited_e(number_n_div_d.n, false);//TODO: global
    std::vector<Scalar> powers;//TODO: global
    powers.reserve(number_n_div_d.n);
    auto &number_d = numberCache[d];
    computeCoprimes(number_d);
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
                        std::vector<Scalar> orbit_d1;
                        const auto n = d * number_n_div_d.n;
                        const auto gcd_h = n / number_n_h.n;
                        for (Scalar i = 0; i < n; i += gcd_h) {
                            orbit_d1.push_back(i);
                        }
                        computeCoprimes(number_n_h);
                        bool newClass = true;
                        for (const auto coprime_h: number_n_h.coprimes) {
                            Orbit orbit1;
                            std::transform(orbit_d1.begin(), orbit_d1.end(), std::back_inserter(orbit1),
                                           [coprime_h, n](const auto d_h) {
                                               return (1 + d_h * coprime_h) % n;
                                           });
                            for (const auto coprime_e: number_d.coprimes) {
                                SkewMorphism skew;
                                skew.n = n;
                                skew.d = d;
                                skew.h = orbit1[1] - 1;
                                skew.s = 1;
                                skew.c = 2;
                                skew.r = n_h;
                                for (Scalar id = 0; id < d; ++id) {
                                    skew.pi.push_back(powers[((coprime_e * id) % d) * step]);
                                }
                                finish(skew, orbit1);
                                addSkewMorphism(std::move(skew), skewMorphisms, newClass);
                                newClass = false;
                            }
                        }
//                        clear(number_n_h);
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
}

PROFILE void sOtherThan1(const Scalar d, const Number &number_n_div_d, SkewMorphisms &skewMorphisms) {
    SkewSet foundSkews;
    std::vector<std::vector<std::vector<std::size_t>>> roots;
    auto &number_d = numberCache[d];
    computeCoprimes(number_d);
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
        std::vector<Scalar> powers_s;//TODO: global
        std::vector<Scalar> powers_e;//TODO: global
        powers_s.reserve(number_n_div_d.n);
        powers_e.reserve(number_n_div_d.n);
//        const auto number_n_b_coprimes = std::move(number_n_b.coprimes);
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
            powers_s.clear();
            powers_s.push_back(1);
            while (power_s != 1) {
                power_sum += power_s;
                powers_s.push_back(power_s);
                power_s = (power_s * s) % number_n_div_d.n;
            }
            const Scalar rH = powers_s.size();
            auto &number_rH = numberCache[rH];
            computeCoprimes(number_rH);
            for (const auto coprime_rH: number_rH.coprimes) {
                visited_s[powers_s[coprime_rH]] = true;
            }
//            if (number_rH.n != number_n_b.n) {
//                clear(number_rH);
//            }
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
                    powers_e.push_back(1);
                    while (power_e != 1) {//TODO: until visited, power sum reuse?
                        powers_e.push_back(power_e);
                        power_e = (power_e * e) % number_r.n;//TODO: alias n_h, n_div_d
                    }
                    const Scalar order = powers_e.size();
                    if (order % d == 0) {
                        const auto step = order / d;
                        if (!visited_e[powers_e[step]]) {
                            for (auto i = step; i < order; i += step) {
                                power_sum_e += powers_e[i];
                            }
                            power_sum_e = power_sum_e % number_r.n;

                            if (power_sum_e % gcd_nh_b == 0) {
                                const auto a = compute_a(d, power_sum_e, rH, power_sum, n_h);
                                const auto &number_a_b = numberCache[a / gcd_nh_b];
                                if (isCoprime(number_n_b, number_a_b)) {
                                    const auto gcd_nh_b_b = gcd(number_n_b, number_gcd_nh_b);
                                    const auto &number_gcd_nh_b_b = numberCache[gcd_nh_b_b];
                                    const auto increment = numberCache[d].phi * number_rH.phi * gcd_nh_b_b * number_gcd_nh_b.phi / number_gcd_nh_b_b.phi;//TODO: can be precompted;

                                    auto increment2 = 0;
                                    const auto n = d * number_n_div_d.n;
                                    auto &number_n_h = numberCache[n_h];
                                    computeCoprimes(number_n_h);
                                    computeCoprimes(number_rH);
                                    const auto gcd_h = n / n_h;
                                    for (const auto coprime_s: number_rH.coprimes) {
                                        std::vector<Scalar> orbit_d1;
                                        Scalar sum = 0;
                                        for (Scalar i = 0; i < g_r; ++i) {
                                            for (Scalar is = 0; is < rH; ++is) {
                                                orbit_d1.push_back(sum);
                                                sum += gcd_h * powers_s[(coprime_s * is) % rH];  // TODO: iterate over coprimes and multiply by it instead of applying automorphism
                                            }
                                        }
                                        const auto x = (inverse(number_n_b.n, number_a_b.n) * ((powers_s[coprime_s] - 1) / gcd_b)) % number_n_b.n;
                                        for (const auto coprime_h: number_n_h.coprimes) {
                                            if ((coprime_h - x) % number_n_b.n != 0) {
                                                continue;
                                            }
                                            Orbit orbit1;
                                            std::transform(orbit_d1.begin(), orbit_d1.end(), std::back_inserter(orbit1),
                                                           [coprime_h, n](const auto d_h) {
                                                               return (1 + d_h * coprime_h) % n;
                                                           });
                                            for (const auto coprime_e: number_d.coprimes) {
                                                SkewMorphism skew;
                                                skew.n = n;
                                                skew.d = d;
                                                skew.h = orbit1[1] - 1;
                                                skew.s = powers_s[coprime_s];
                                                skew.c = 2;
                                                skew.r = r;
                                                for (Scalar id = 0; id < d; ++id) {
                                                    skew.pi.push_back(powers_e[((coprime_e * id) % d) * step]);
                                                }
                                                finish(skew, orbit1);

                                                SkewMorphism ro;
                                                ro.n = skew.r;
                                                ro.d = 1;
                                                ro.h = skew.pi[1] - 1;
                                                ro.s = skew.pi[1];
                                                ro.c = 1;
                                                ro.r = skew.d;
                                                ro.pi = {1};
                                                finish(ro, skew.pi);

                                                addSkewClassByRepresentant(foundSkews, roots, ro, skew);
                                                ++increment2;
                                            }
                                        }
                                    }
                                    if (increment != increment2) {
                                        throw "";
                                    }
//                                    clear(number_n_h);
                                }
                            }
                        }
                    }
                    for (const auto power: powers_e) {
                        visited_e[power] = true;
                    }
                    powers_e.clear();
                }
                for (auto e = rH + 1; e < r; e += rH) {
                    visited_e[e] = false;
                }
            }
        }
//        clear(number_n_b);
    }
}

PROFILE void computeAutomorphisms(Number &number) {
    SkewSet foundSkews;
    std::vector<std::vector<std::vector<std::size_t>>> roots;
    computeCoprimes(number);
    const auto n = number.n;
    const auto one = 1 % n;
    std::vector<bool> visited_s(number.n, false);//TODO: global
    std::vector<Scalar> powers_s;//TODO: global
    powers_s.reserve(number.n);
    for (const auto s: number.coprimes) {
        if (visited_s[s]) {
            continue;
        }
        DoubleScalar power_s = s;
        powers_s.clear();
        powers_s.push_back(one);
        while (power_s != one) {
            powers_s.push_back(power_s);
            power_s = (power_s * s) % number.n;
        }
        const Scalar rH = powers_s.size();
        auto &number_rH = numberCache[rH];
        computeCoprimes(number_rH);
        for (const auto coprime_rH: number_rH.coprimes) {
            visited_s[powers_s[coprime_rH]] = true;
        }
        bool newClass = true;

        for (const auto coprime_s: number_rH.coprimes) {
            Orbit orbit1;
            orbit1.reserve(rH);
            for (Scalar i = 0; i < rH; ++i) {
                orbit1.push_back(powers_s[(coprime_s * i) % rH]);
            }

            SkewMorphism automorphism;
            automorphism.n = n;//TODO: unify filling these
            automorphism.d = 1;
            automorphism.h = powers_s[coprime_s] - one;
            automorphism.s = powers_s[coprime_s];
            automorphism.c = automorphism.h == 0 ? 0 : 1;
            automorphism.r = orbit1.size();
            automorphism.pi.resize(1, 1 % automorphism.r);

            finish(automorphism, orbit1);

            addSkewMorphism(std::move(automorphism), number.skewMorphisms, newClass);
            newClass = false;
        }
    }
//    clear(number);
}

std::string to_string(const Orbit &orbit) {
    std::string result;
    bool first = true;
    for (const auto e: orbit) {
        result += (first ? "" : ", ") + std::to_string(e);
        first = false;
    }
    return result;
}

std::string to_string(const SkewMorphism &skewMorphism) {
    if (skewMorphism.permutation.orbits.empty()) {
        return "Id(sn)";
    }
    std::string result;
    for (const auto &orbit: skewMorphism.permutation.orbits) {
        result += "(" + to_string(orbit) + ")";
    }
    return result;
}

std::string to_json(const Orbit &orbit) {
    std::string result;
    bool first = true;
    result += "[";
    for (const auto e: orbit) {
        result += (first ? "" : ",") + std::to_string(e);
        first = false;
    }
    result += "]";
    return result;
}

std::string to_json(const SkewMorphism &skewMorphism) {
    std::string result;
    bool first = true;
    result += "[";
    for (const auto &orbit: skewMorphism.permutation.orbits) {
        result += (first ? "" : ",") + to_json(orbit);
        first = false;
    }
    result += "]";
    return result;
}

PROFILE bool isSkewMorphism(const SkewMorphism &skewMorphism, const Function& function) {
    const auto n = skewMorphism.n;
    bool good = true;
    Function c_j(n, 0);
    std::vector<bool> visited_j(skewMorphism.d, false);
    for (std::size_t ij = 0; ij < skewMorphism.d; ++ij) {
        const auto j = skewMorphism.pi[ij];
        for (Scalar i = 1; i < n; ++i) {
            c_j[i] = apply(skewMorphism, skewMorphism.r - j, i);//TODO: skewMorphism.r je order?
        }
        std::size_t ijj = 0;
        for (; ijj < skewMorphism.d; ++ijj) {
            if (visited_j[ijj]) {
                continue;
            }
            const auto e = c_j[ijj];
            good = true;
            for (std::size_t x = 1; x < n; ++x) {
                if ((c_j[(function[x] + ijj) % n] + (n - x)) % n != e) {
                    good = false;
                    break;
                }
            }
            if (good) {
                visited_j[ijj] = true;
                break;
            }
        }
        if (!good) {
            break;
        }
        for (std::size_t i = ijj + skewMorphism.d; i < n; i += skewMorphism.d) {
            const auto e = c_j[i];
            for (std::size_t x = 1; x < n; ++x) {
                if ((c_j[(function[x] + i) % n] + (n - x)) % n != e) {
                    good = false;
                    break;
                }
            }
            if (!good) {
                break;
            }
        }
        if (!good) {
            break;
        }
    }
    for (const auto v: visited_j) {
        if (!v) {
            good = false;
            break;
        }
    }
    return good;
}

PROFILE void periodicallyFillOrbit(std::size_t p, std::size_t i, const SkewMorphism &psi, Orbit &t) {
    const std::size_t orbitSize = psi.r;
    const auto e = t[i];

    const auto &place = psi.permutation.places[e];
    const auto &orbit = psi.permutation.orbits[place.orbitIndex];

    std::size_t index = place.indexOnOrbit;
    for (std::size_t j = p + i; j < t.size(); j += p) {
        ++index;
        if (index == orbitSize) index = 0;
        t[j] = orbit[index];
    }
}

PROFILE bool isPermutation(Function &function) {//TODO: spojit s pocitanim funkcie
    std::vector<bool> visited(function.size(), false);
    for (const auto f: function) {
        if (visited[f]) {
            return false;
        }
        visited[f] = true;
    }
    return true;
}

PROFILE bool computeOrbit1(const Function &function, Orbit &orbit1) {
    Scalar c = 1;
    std::size_t index = 0;
    do {
        orbit1[index] = c;
        c = function[c];
        ++index;
    } while (c != 1 && index < orbit1.size());
    return c == 1 && index == orbit1.size();
}

PROFILE bool compareOrbits(const Orbit &sparseOrbit, const Orbit &orbit2) {
    for (std::size_t i = 0; i < sparseOrbit.size(); ++i) {
        if (sparseOrbit[i] == 0) {
            continue;
        }
        if (sparseOrbit[i] != orbit2[i]) {
            return false;
        }
    }
    return true;
}

PROFILE bool checkPsiCycles(const Scalar p, const SkewMorphism &psi, const Orbit &orbit1) {
    for (std::size_t i = p; i < orbit1.size(); ++i) {
        if (orbit1[i] != apply1(psi, orbit1[i - p])) {
            return false;
        }
    }
    return true;
}

PROFILE bool checkFirstPMod(const SkewMorphism &ro, const Orbit &orbit1) {
    for (std::size_t i = 1; i < ro.d; ++i) {
        if (orbit1[i] % ro.r != ro.pi[i]) {
            return false;
        }
    }
    return true;
}

PROFILE Scalar powerToSkew(Scalar p, const SkewMorphism &ro) {
    const auto &orbit1_ro = getOrbit1(ro);
    for (const auto prime: numberCache[p].primes) {
        if (p == prime) {
            break;
        }
        Scalar period = 1;
        while (orbit1_ro[period] % prime != 1) {
            ++period;
        }
        if (ro.r % period != 0) {
            continue;
        }
        bool ok = true;
        for (std::size_t i = period + 1; i < orbit1_ro.size(); ++i) {
            if ((orbit1_ro[i] - orbit1_ro[i - period]) % prime != 0) {
                ok = false;
                break;
            }
        }
        if (!ok) {
            continue;
        }
        const auto modulo = ro.pi[0] % period;//TODO: copy-paste, zjednotit
        bool all = true;
        for (std::size_t j = 1; j < ro.pi.size(); ++j) {
            if (ro.pi[j] % period != modulo) {
                all = false;
                break;
            }
        }
        if (all) {
            return prime;
        }
    }
    return p;
}

PROFILE SkewMorphism &getSkewByIndex(SkewMorphisms &skewMorphisms, std::size_t index) {
    if (index < skewMorphisms.automorphisms.size()) {
        return skewMorphisms.automorphisms[index];
    }
    index -= skewMorphisms.automorphisms.size();
    if (index < skewMorphisms.properCosetPreserving.size()) {
        return skewMorphisms.properCosetPreserving[index];
    }
    index -= skewMorphisms.properCosetPreserving.size();
    return skewMorphisms.properNonPreserving[index];
}

PROFILE std::size_t getSkewCount(const SkewMorphisms &skewMorphisms) {
    return skewMorphisms.automorphisms.size() +
           skewMorphisms.properCosetPreserving.size() +
           skewMorphisms.properNonPreserving.size();
}

PROFILE std::size_t getClassesCount(const SkewMorphisms &skewMorphisms) {
    return skewMorphisms.autoClasses.size() +
           skewMorphisms.preservingClasses.size() +
           skewMorphisms.nonPreservingClasses.size();
}

PROFILE const SkewMorphism &getPreservingSkewByIndex(const SkewMorphisms &skewMorphisms, std::size_t index) {  // TODO: iterators
    if (index < skewMorphisms.automorphisms.size()) {
        return skewMorphisms.automorphisms[index];
    }
    index -= skewMorphisms.automorphisms.size();
    return skewMorphisms.properCosetPreserving[index];
}

PROFILE std::size_t getPreservingSkewCount(const SkewMorphisms &skewMorphisms) {
    return skewMorphisms.automorphisms.size() + skewMorphisms.properCosetPreserving.size();
}

PROFILE std::size_t getPreservingClassIndex(const SkewMorphisms &skewMorphisms, std::size_t index) {
    if (index < skewMorphisms.autoClasses.size()) {
        return skewMorphisms.autoClasses[index][0];
    }
    index -= skewMorphisms.autoClasses.size();
    return skewMorphisms.automorphisms.size() + skewMorphisms.preservingClasses[index][0];
}

PROFILE std::size_t getPreservingClassesCount(const SkewMorphisms &skewMorphisms) {
    return skewMorphisms.autoClasses.size() + skewMorphisms.preservingClasses.size();
}

PROFILE const SkewMorphism &getProperSkewByIndex(const SkewMorphisms &skewMorphisms, std::size_t index) {
    if (index < skewMorphisms.properCosetPreserving.size()) {
        return skewMorphisms.properCosetPreserving[index];
    }
    index -= skewMorphisms.properCosetPreserving.size();
    return skewMorphisms.properNonPreserving[index];
}

PROFILE std::size_t getProperSkewCount(const SkewMorphisms &skewMorphisms) {
    return skewMorphisms.properCosetPreserving.size() + skewMorphisms.properNonPreserving.size();
}

PROFILE std::size_t getProperClassIndex(const SkewMorphisms &skewMorphisms, std::size_t index) {
    if (index < skewMorphisms.preservingClasses.size()) {
        return skewMorphisms.preservingClasses[index][0];
    }
    index -= skewMorphisms.preservingClasses.size();
    return skewMorphisms.properCosetPreserving.size() + skewMorphisms.nonPreservingClasses[index][0];
}

PROFILE std::size_t getProperClassesCount(const SkewMorphisms &skewMorphisms) {
    return skewMorphisms.preservingClasses.size() + skewMorphisms.nonPreservingClasses.size();
}

PROFILE bool isPowerOf(const SkewMorphism &base, Scalar exponent, const SkewMorphism &power) {
    for (const auto &orbit: power.permutation.orbits) {
        Scalar oldO = orbit[0];
        for (std::size_t i = 1; i < orbit.size(); ++i) {
            const auto o = orbit[i];
            if (apply(base, exponent, oldO) != o) {
                return false;
            }
            oldO = o;
        }
        if (apply(base, exponent, oldO) != orbit[0]) {
            return false;
        }
    }
    for (const auto &place: power.permutation.places) {//TODO: fixed points as part of Permutation
        if (place.orbitIndex == -1 && apply(base, exponent, place.indexOnOrbit) != place.indexOnOrbit) {
            return false;
        }
    }
    return true;
}

PROFILE Function inverse(const Function &function, Scalar n) {
    Function result(n, -1);
    for (Scalar x = 0; x < function.size(); ++x) {
        result[function[x]] = x;
    }
    return result;
}

PROFILE void initializeSplitIndex(Orbit &t, std::vector<Index> &splitIndex, const std::vector<std::vector<Scalar>> &values, const std::vector<Index> &free_x_index, Scalar exponent, const SkewMorphism &power) {
    for (std::size_t i = 0; i <= values.size(); ++i) {
        splitIndex[i] = 0;
    }
    t[0] = 1;
    periodicallyFillOrbit(exponent, 0, power, t);
    for (std::size_t i = 0; i < free_x_index.size(); ++i) {
        const auto index = free_x_index[i];
        const auto valueIndex = splitIndex[i];
        const auto value = values[i][valueIndex];
        t[index] = value;
        //TODO: pozri size_t
        periodicallyFillOrbit(exponent, index, power, t);
    }
}

PROFILE bool isValidSplitIndex(std::vector<Index> &splitIndex, const std::vector<std::vector<Scalar>> &values) {
    return splitIndex[values.size()] == 0;
}

PROFILE void incrementSplitIndex(Orbit &t, std::vector<Index> &splitIndex, const std::vector<std::vector<Scalar>> &values, const std::vector<Index> &free_x_index, Scalar exponent, const SkewMorphism &power) {
    ++splitIndex[0];
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (splitIndex[i] == 0 || splitIndex[i] == values[i].size()) {
            if (splitIndex[i] == values[i].size()) {
                splitIndex[i] = -static_cast<Index>(values[i].size());
            }
            ++splitIndex[i + 1];
        } else {
            const auto index = free_x_index[i];
            const auto valueIndex = splitIndex[i] >= 0 ? splitIndex[i] :  -(splitIndex[i] + 1);
            const auto value = values[i][valueIndex];
            t[index] = value;
            periodicallyFillOrbit(exponent, index, power, t);
            break;
        }
    }
}

PROFILE void clearOrbit(Orbit &orbit) {
    for (auto &o: orbit) {
        o = 0;
    }
}

PROFILE void addSkewClassByRepresentant(SkewSet &foundSkews,  // TODO: pridaj do SkewMorphisms
                                        std::vector<std::vector<std::vector<std::size_t>>> &roots,
                                        const SkewMorphism &ro,
                                        const SkewMorphism &phi) {
    const auto n = phi.n;
    const auto r = phi.r;
    const auto d = phi.d;

    auto &number = numberCache[n];
    auto &number_r = numberCache[r];

    CompactSkewMorphism compactSkewMorphism = {getOrbit1(phi), phi.pi};  // TODO: do not copy

    if (foundSkews.find(compactSkewMorphism) != foundSkews.end()) {
        return;
    }

//    if (phi.h % phi.d != 0) {
        for (std::size_t j = 0; j < roots.size(); ++j) {
            auto &other = getSkewByIndex(number.skewMorphisms, j);
            if (phi.r % other.r != 0) {
                continue;
            }
            const auto e = phi.r / other.r;
            if (isPowerOf(phi, e, other)) {
                if (phi.inverseOrbit) {
                    other.powerOfInverseOrbit = true;
                }
                roots[j][e].push_back(roots.size());
            }
        }
        roots.emplace_back(n);
//    }

    foundSkews.insert(compactSkewMorphism);
    addSkewMorphism(phi, number.skewMorphisms, true);

    CompactSkewMorphism compactSkewMorphism2;
    auto &orbit1 = compactSkewMorphism2.orbit1;
    orbit1.resize(r, 1);
    auto &orbit1_ro = compactSkewMorphism2.pi;
    orbit1_ro.resize(d, 1);

    computeCoprimes(number);
    computeCoprimes(number_r);
    for (const auto coprime_r: number_r.coprimes) {
        const auto coprime_r_inverse = (r + inverse(r, coprime_r)) % r;
        for (const auto coprime_n: number.coprimes) {
            const auto coprime_n_inverse = (n + inverse(n, coprime_n)) % n;
            const auto coprime_d = coprime_n % d;
            SkewMorphism phi2;
            phi2.n = n;
            phi2.d = d;
            phi2.r = r;
            const auto &coprime_n_place = phi.permutation.places[coprime_n];
            const auto &coprime_n_orbit = phi.permutation.orbits[coprime_n_place.orbitIndex];
            std::size_t index = coprime_n_place.indexOnOrbit;
            for (std::size_t i = 1; i < r; ++i) {
                index += coprime_r;
                if (index >= r) index -= r;
                orbit1[i] = (coprime_n_orbit[index] * coprime_n_inverse) % n;
            }
            const auto &coprime_r_place = ro.permutation.places[coprime_r];
            const auto &coprime_r_orbit = ro.permutation.orbits[coprime_r_place.orbitIndex];
            index = coprime_r_place.indexOnOrbit;
            for (std::size_t i = 1; i < d; ++i) {
                index += coprime_d;
                if (index >= d) index -= d;
                orbit1_ro[i] = (coprime_r_orbit[index] * coprime_r_inverse) % r;
            }

            if (foundSkews.find(compactSkewMorphism2) != foundSkews.end()) {
                continue;
            }

            phi2.h = orbit1[1] - 1;
            phi2.pi = orbit1_ro;

            finish(phi2, orbit1);

            phi2.s = apply1(phi2, d) / d;
            phi2.c = ro.c + 1;

//            if (phi2.h % phi2.d != 0) {
                for (std::size_t j = 0; j < roots.size(); ++j) {
                    auto &other = getSkewByIndex(number.skewMorphisms, j);
                    if (phi2.r % other.r != 0) {
                        continue;
                    }
                    const auto e = phi2.r / other.r;
                    if (isPowerOf(phi2, e, other)) {
                        if (phi2.inverseOrbit) {
                            other.powerOfInverseOrbit = true;
                        }
                        roots[j][e].push_back(roots.size());
                    }
                }
                roots.emplace_back(n);
//            }

            foundSkews.insert(compactSkewMorphism2);
            addSkewMorphism(std::move(phi2), number.skewMorphisms, false);
        }
    }
}

//TODO: zoradit podla radu
//TODO: kanonicke poradie autov, aj cosetov

//TODO: nemozu nejake nasobenia scitania pretiect?
PROFILE void computeProperNotPreserving(Number &number) {
    SkewSet foundSkews;
    const auto n = number.n;
    const auto number_nlambda = numberCache[n * number.lambda];

    const auto maxPrime = getMaxPrime(number);
    auto n_div_maxPrime = n / maxPrime;

    std::vector<std::vector<std::vector<std::size_t>>> roots(getPreservingSkewCount(number.skewMorphisms));
    for (std::size_t i = 0; i < roots.size(); ++i) {
        auto &r = roots[i];
        r.resize(n);
        auto &power = getSkewByIndex(number.skewMorphisms, i);
        for (std::size_t j = 0; j < roots.size(); ++j) {
            const auto &root = getSkewByIndex(number.skewMorphisms, j);
            if (root.r % power.r != 0) {
                continue;
            }
            const auto exponent = root.r / power.r;
            if (isPowerOf(root, exponent, power)) {
                if (root.inverseOrbit) {
                    power.powerOfInverseOrbit = true;
                }
                r[exponent].push_back(j);
            }
        }
    }

    Function function(n, 0);
    std::vector<std::vector<OrbitPlace>> moduloOrbits(n);
    std::vector<Index> splitIndex(n + 1, 0);

    CompactSkewMorphism compactSkewMorphism;

    for (const auto m: number_nlambda.divisors) {
        if (m >= n) {
            continue;
        }
        const auto &number_m = numberCache[m];
        if (isCoprime(number, number_m)) {
            continue;
        }
        const auto r = m;

        Orbit t(r, 0);
        t[0] = 1;
        Orbit &orbit1 = compactSkewMorphism.orbit1;
        orbit1.clear();
        orbit1.resize(r, 0);
//        for (std::size_t ro_class_index = 0; ro_class_index < getProperClassesCount(number_m.skewMorphisms); ++ro_class_index) {
//            const auto ro_index = getProperClassIndex(number_m.skewMorphisms, ro_class_index);

        for (std::size_t ro_index = 0; ro_index < getProperSkewCount(number_m.skewMorphisms); ++ro_index) {
            const auto &ro = getProperSkewByIndex(number_m.skewMorphisms, ro_index);

            const auto d = ro.r;

            if (n % (d * maxPrime) != 0) {
                continue;
            }

            const auto n_div_d = n / d;

            if (isCoprime(numberCache[n_div_d], number_m)) {
                continue;
            }

            const auto p = ro.d;
            const auto ord_psi = m / p;

            const auto &orbit1_ro = getOrbit1(ro);

            compactSkewMorphism.pi = orbit1_ro;  // TODO: do not copy

            for (std::size_t psi_class_index = 0; psi_class_index < getPreservingClassesCount(number.skewMorphisms); ++psi_class_index) {
                const auto psi_index = getPreservingClassIndex(number.skewMorphisms, psi_class_index);

//            for (std::size_t psi_index = 0; psi_index < getPreservingSkewCount(number.skewMorphisms); ++psi_index) {
                const auto &psi = getPreservingSkewByIndex(number.skewMorphisms, psi_index);

                if (psi.h == 0) {
                    continue;
                }
                if (psi.r != ord_psi) {
                    continue;
                }
                if (psi.h % d != 0) {
                    continue;
                }
                if (psi.max_orbits < p) {
                    continue;
                }

                const auto exponent = powerToSkew(p, ro);
                const auto p_exponent = p / exponent;

                const auto &possiblePowerIndices = roots[psi_index][p_exponent];

                std::vector<Index> free_x_index;
                std::set<Index> set;
                for (const auto index: ro.free_x) {
                    set.insert(index % exponent);
                }
                std::copy(set.begin(), set.end(), std::back_inserter(free_x_index));

                clearOrbit(t);
                const auto inverse_ro_pi = inverse(ro.pi, d);
                for (const auto power_index: possiblePowerIndices) {
                    const auto power = getSkewByIndex(number.skewMorphisms, power_index);  // TODO: do not copy

                    const auto &power_orbits = power.permutation.orbits;

                    for (std::size_t i = 0; i < exponent; ++i) {
                        moduloOrbits[i].clear();
                    }
                    for (Index i = 0; i < power.max_orbits; ++i) {
                        const auto &orbit = power_orbits[i];
                        const auto modulo = orbit[0] % d;
                        auto index = inverse_ro_pi[modulo];
                        if (index == -1) {
                            continue;
                        }
                        bool all = true;
                        for (std::size_t j = 1; j < orbit.size(); ++j) {
                            index += exponent;
                            if (index >= p) index -= p;
                            if (orbit[j] % d != ro.pi[index]) {
                                all = false;
                                break;
                            }
                        }
                        if (!all) {
                            continue;
                        }
                        index = inverse_ro_pi[modulo];
                        const auto smallestIndex = index % exponent;
                        const auto offset = index < exponent ? 0 : p_exponent - index / exponent;
                        moduloOrbits[smallestIndex].push_back(OrbitPlace{.orbitIndex = i, .indexOnOrbit = offset});
                    }
                    bool allModulos = true;
                    for (Index i = 0; i < exponent; ++i) {
                        if (moduloOrbits[i].empty()) {
                            allModulos = false;
                            break;
                        }
                    }
                    if (!allModulos) {
                        continue;
                    }

                    std::vector<std::vector<Scalar>> free_x_values;
                    free_x_values.reserve(free_x_index.size());
                    for (const auto index: free_x_index) {
                        free_x_values.emplace_back();
//                        free_x_values.back().reserve(n_div_d);
                        for (const auto moduloOrbitPlace: moduloOrbits[index]) {
                            const auto &orbit = power_orbits[moduloOrbitPlace.orbitIndex];
                            for (std::size_t i = moduloOrbitPlace.indexOnOrbit; i < orbit.size(); i += p_exponent) {
                                free_x_values.back().push_back(orbit[i]);
                            }
                        }
                    }

                    for (initializeSplitIndex(t, splitIndex, free_x_values, free_x_index, exponent, power); isValidSplitIndex(splitIndex, free_x_values); incrementSplitIndex(t, splitIndex, free_x_values, free_x_index, exponent, power)) {
                        const Function &pi = orbit1_ro;
                        if (!computeFunction(n, pi, t, function)) {
                            continue;
                        }
                        if (!isPermutation(function)) {
                            continue;
                        }

                        if (!computeOrbit1(function, orbit1)) {
                            continue;
                        }

                        if (!checkFirstPMod(ro, orbit1)) {
                            continue;
                        }
                        //TODO: ktore checky treba robit?
                        if (!compareOrbits(t, orbit1)) {
                            continue;
                        }
                        if (!checkPsiCycles(exponent, power, orbit1)) {
                            continue;
                        }

                        if (foundSkews.find(compactSkewMorphism) != foundSkews.end()) {
                            continue;
                        }
                        SkewMorphism phi;
                        phi.n = n;
                        phi.d = d;
                        phi.h = t[1] - 1;
                        phi.r = r;
                        computeOrbits(phi, function);
                        phi.pi = pi;

                        if (!isSkewMorphism(phi, function)) {
                            continue;
                        }

                        computeMaxOrbits(phi);
                        phi.s = apply1(phi, d) / d;
                        phi.c = ro.c + 1;

                        addSkewClassByRepresentant(foundSkews, roots, ro, phi);
                    }
                }
            }
        }
    }
    //printf("\r");
//    printf("%d %ld\n", n, foundSkews.size());
}

PROFILE void countSkewmorphisms(Number &number) {
    if (getSkewCount(number.skewMorphisms) != 0) {
        return;
    }
    const auto n = number.n;
    const auto phi = number.phi;
    computeAutomorphisms(number);

    if (gcd(number, numberCache[phi]) > 1) {//TODO: toto viem nejako vyuzit a predratat?
        const auto maxPrime = getMaxPrime(number);
        const auto &number_n_div_maxPrime = numberCache[n / maxPrime];
        for (const auto d: number_n_div_maxPrime.divisors) {//TODO: iterate over n_d first
            if (d == 1) {
                continue;
            }
            const auto n_div_d = n / d;
            auto &number_n_div_d = numberCache[n_div_d];

            sEquals1(d, number_n_div_d, number.skewMorphisms);
            sOtherThan1(d, number_n_div_d, number.skewMorphisms);
        }
        if (number.powerOfTwo > 16 || !number.squareFree) {
            computeProperNotPreserving(number);
        }
    }
}

void print(const Number &number) {
    const auto n = number.n;
    const auto phi = number.phi;
    const auto nskew = getSkewCount(number.skewMorphisms);
    const auto nproper = getProperSkewCount(number.skewMorphisms);
    const auto npreserving = getPreservingSkewCount(number.skewMorphisms) - phi;
    const auto nother = nproper - npreserving;
    const auto nclasses = getProperClassesCount(number.skewMorphisms);
    const auto nclasses_a = number.skewMorphisms.autoClasses.size();
    const auto nclasses_c = number.skewMorphisms.preservingClasses.size();
    const auto nclasses_o = number.skewMorphisms.nonPreservingClasses.size();
//    printf("Total number of skew morphisms of C_%d is %d\n", n, nskew);
//    printf("Check sub-total of automorphisms of C_%d is %d\n", n, phi);
//    auto perc = static_cast<Scalar>(std::round(float(100 * phi) / float(nskew)));
//    printf("Automorphisms account for ~%d%% of all skew morphisms.\n\n", perc);
//    printf("%d    %ld + %d (%ld + %ld)   %ld %ld %ld\n", n, nproper, phi, npreserving, nother, nclasses_a, nclasses_c, nclasses_o);
    printf("%d    %ld + %d (%ld + %ld)\n", n, nproper, phi, npreserving, nother);
}

void fprint(const Number &number, std::ostream& stream) {
    const auto n = number.n;
    const auto phi = number.phi;
    const auto nskew = getSkewCount(number.skewMorphisms);
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

PROFILE SkewMorphism from_short_string(const std::string &string, Scalar n) {
    SkewMorphism result;
    std::istringstream stream;
    stream.str(string);
    char c = '\0';
    Orbit orbit1;
    while (c != ')') {
        stream.get(c);
        Scalar i = 0;
        stream >> i;
        orbit1.push_back(i);
        stream.get(c);
    }
    Function pi;
    while (c != ']') {
        stream.get(c);
        Scalar i = 0;
        stream >> i;
        pi.push_back(i);
        stream.get(c);
    }
    result.d = pi.size();
    result.r = orbit1.size();
    result.n = n;
    if (orbit1.size() == 1) {
        result.h = 0;
    } else {
        result.h = orbit1[1] - 1;
    }
    result.pi = std::move(pi);
    finish(result, orbit1);
    result.s = apply1(result, result.d % result.n) / result.d;
//    result.c = ???;
    return result;
}

SkewMorphism quotient(const SkewMorphism &skew) {
    SkewMorphism ro;
    ro.n = skew.r;
    ro.r = skew.d;

    ro.d = 0;
    Scalar a = 1 % skew.d;
    const auto pi1 = skew.pi[a];
    do {
        ++ro.d;
        ro.pi.push_back(a);
        a = apply1(skew, a) % skew.d;
    } while (skew.pi[a] != pi1);

    finish(ro, skew.pi);

    ro.h = apply1(ro, 1 % ro.n) - 1 % ro.n;
    ro.s = apply1(ro, ro.d % ro.n) / ro.d;
    ro.c = skew.c - 1;

    return ro;
}

bool readSkewMorphisms(Scalar n, SkewMorphisms &skewMorphisms) {
    std::string line;
    std::ifstream ifile;

    SkewSet foundSkews;
    std::vector<std::vector<std::vector<std::size_t>>> roots;

    ifile = std::ifstream{"../skew/" + std::to_string(n) + "_auto.txt"};
    if (!ifile.is_open()) {
        return false;
    }
    while (std::getline(ifile, line)) {
        auto skew = from_short_string(line, n);
        skew.c = skew.h == 0 ? 0 : 1;
        addSkewMorphism(skew, skewMorphisms, true);
        computeCoprimes(numberCache[skew.r]);
        Orbit orbit12(skew.r, 1);
        for (const auto coprime: numberCache[skew.r].coprimes) {
            if (coprime <= 1) {
                continue;
            }
            SkewMorphism skew2;
            skew2.n = skew.n;
            skew2.r = skew.r;
            skew2.d = skew.d;
            skew2.pi = skew.pi;
            const auto &orbit1 = getOrbit1(skew);
            for (std::size_t i = 1; i < skew.r; ++i) {
                orbit12[i] = orbit1[(i * coprime) % skew.r];
            }

            skew2.h = orbit12[1] - 1;
            skew2.s = orbit12[1];
            skew2.c = skew.c;

            finish(skew2, orbit12);

            addSkewMorphism(std::move(skew2), skewMorphisms, false);
        }
    }
    ifile = std::ifstream{"../skew/" + std::to_string(n) + "_coset.txt"};
    while (std::getline(ifile, line)) {
        auto skew = from_short_string(line, n);
        skew.c = 2;
        addSkewClassByRepresentant(foundSkews, roots, quotient(skew), skew);
    }
    ifile = std::ifstream{"../skew/" + std::to_string(n) + "_other.txt"};
    while (std::getline(ifile, line)) {
        auto skew = from_short_string(line, n);
        for (SkewMorphism ro = skew; ro.h != 0; ro = quotient(ro)) {
            ++skew.c;
        }
        addSkewClassByRepresentant(foundSkews, roots, quotient(skew), skew);
    }
    return true;
}

void printMato(std::ofstream &output, const SkewMorphism &skew) {
    output << "\nSkew morphism " << to_string(skew);
    output << "\n of order " << std::to_string(skew.r);
    output << "\n with kernel of order " << std::to_string(skew.n / skew.d);
    output << "\n and smallest kernel generator " << std::to_string(skew.d);
    if (skew.d == 1) {
        output << "\n and power function values [ 1 ]" << std::endl;
    } else {
        output << "\n and power function values [ " << to_string(skew.pi) << " ]" << std::endl;
    }
}

void printHtml(std::ofstream &output, const SkewMorphism &skew, Index index, const std::string &additionalClass, Scalar classSize) {
    output << "\n<li class='skew " << additionalClass <<" n-" << skew.n << " d-" << skew.d << " h-" << skew.h << " r-" << skew.r << " s-" << skew.s << " c-" << skew.c;
    if (skew.inverseOrbit) {
        output << " inv-1";
    }
    if (skew.powerOfInverseOrbit) {
        output << " pow-inv-1";
    }
    output << "' id='n-" << skew.n <<"-id-" << index << "' data-repr='" << to_json(skew) << "'>";
    output << "\nSkew morphsim class of size " << classSize << " <a class='link' href='#n-" << skew.n <<"-id-" << index << "'>#</a>";
    output << "\n<ul>";
    output << "\n<li class='repr'></li>";
    output << "\n<li>of order " << std::to_string(skew.r) << "</li>";
    output << "\n<li>with kernel of order " << std::to_string(skew.n / skew.d) << "</li>";
    output << "\n<li>and smallest kernel generator " << std::to_string(skew.d) << "</li>";
    if (skew.d == 1) {
        output << "\n<li>and power function values [ 1 ]</li>" << std::endl;
    } else {
        output << "\n<li>and power function values [ " << to_string(skew.pi) << " ]</li>";
    }
    output << "\n<li>and with periodicity " << std::to_string(quotient(skew).d) << "</li>";

    output << "\n</ul>" << std::endl;
    output << "\n</li>" << std::endl;
}

std::string factorizationString(Scalar n) {
    if (n == 1) {
        return "1";
    }
    std::string result;
    const auto &number = numberCache[n];
    for (std::size_t i = 0; i < number.primes.size(); ++i) {
        const Scalar prime = number.primes[i];
        const Scalar exponent = numberCache[number.powersOfPrimes[i]].divisors.size() - 1;
        if (i > 0) {
            result += "x";
        }
        result += std::to_string(prime);
        if (exponent > 1) {
            result += "^" + std::to_string(exponent);
        }
    }
    return result;
}

void printHtml(std::ofstream &output, const Number &number) {
    const auto n = number.n;
    auto &skewMorphisms = number.skewMorphisms;
    output << "\n<div class='n' id='n-" << n << "'>";
    output << "\nn = " << n << " = " << factorizationString(n);
    output << "\n<ul>";
    output << "\n<li>Number of selected skew morphism classes of C_" << n << " is <span class='count'>0</span></li>";
    output << "\n<li>Total number of skew morphisms of C_" << n << " is " << getSkewCount(skewMorphisms) << "</li>";
    output << "\n<li>Total number of skew morphism classes of C_" << n << " is " << getClassesCount(skewMorphisms) << "</li>";
    output << "\n<li>Total number of automorphisms of C_" << n << " is " << number.phi << "</li>";
    output << "\n</ul>";

    Index index = 0;
    output << "\n<ol>";
    for (const auto &c: number.skewMorphisms.autoClasses) {
        const auto &skew = number.skewMorphisms.automorphisms[c[0]];
        printHtml(output, skew, ++index, "auto", c.size());
    }
    for (const auto &c: number.skewMorphisms.preservingClasses) {
        const auto &skew = number.skewMorphisms.properCosetPreserving[c[0]];
        printHtml(output, skew, ++index, "coset", c.size());
    }
    for (const auto &c: number.skewMorphisms.nonPreservingClasses) {
        const auto &skew = number.skewMorphisms.properNonPreserving[c[0]];
        printHtml(output, skew, ++index, "other", c.size());
    }
    output << "\n</ol>";

    output << "\n</div>" << std::endl;
}

int main(int argc, char *argv[]) {
    Scalar N = 100;
    if (argc > 1) {
        N = std::stoi(argv[1]);
    }
    Scalar A = 1;
    if (argc > 2) {
        A = std::stoi(argv[2]);
    }

    //TODO: compute coprimes ale len do N
    prepareNumbers(N * N);

    std::vector<Scalar> computeVector;
    std::set<Scalar> computeSet;
    for (Scalar n = A; n <= N; ++n) {
        computeVector.push_back(n);
        computeSet.insert(n);
    }
    for (std::size_t i = 0; i < computeVector.size(); ++i) {
        const auto a = computeVector[i];
        for (const auto d: numberCache[a * numberCache[a].lambda].divisors) {
            if (d >= a) {
                continue;
            }
            if (isCoprime(numberCache[a], numberCache[d])) {
                continue;
            }
            if (computeSet.find(d) == computeSet.end()) {
                computeVector.push_back(d);
                computeSet.insert(d);
            }
        }
    }

//    for (auto n = A; n <= N; ++n) {
    for (auto n: computeSet) {
        auto &number = numberCache[n];
//        if (number.powerOfTwo <= 16 && number.squareFree) {
            if (!readSkewMorphisms(n, number.skewMorphisms)) {
                countSkewmorphisms(number);
            }
//            fprint(number, std::cerr);
//        }
        if (getSkewCount(number.skewMorphisms) != number.phi) {
            print(number);
        }
        std::ofstream file;
        file = std::ofstream{"../skew/" + std::to_string(number.n) + "_auto.txt"};
        for (const auto &c: number.skewMorphisms.autoClasses) {
            const auto &skew = number.skewMorphisms.automorphisms[c[0]];
            file << to_short_string(skew) << std::endl;
        }
        file = std::ofstream{"../skew/" + std::to_string(number.n) + "_coset.txt"};
        for (const auto &c: number.skewMorphisms.preservingClasses) {
            const auto &skew = number.skewMorphisms.properCosetPreserving[c[0]];
            file << to_short_string(skew) << std::endl;
        }
        file = std::ofstream{"../skew/" + std::to_string(number.n) + "_other.txt"};
        for (const auto &c: number.skewMorphisms.nonPreservingClasses) {
            const auto &skew = number.skewMorphisms.properNonPreserving[c[0]];
            file << to_short_string(skew) << std::endl;
        }
    }

//    std::ofstream output("skew" + std::to_string(N) + "output.txt");
//    for (auto n = N; n <= N; ++n) {
//        auto &number = numberCache[n];
//        output << "\nn = " << n << std::endl;
//
//        auto &skewMorphisms = number.skewMorphisms;
//
//        std::sort(skewMorphisms.automorphisms.begin(), skewMorphisms.automorphisms.end(), [](const auto& lhs, const auto &rhs) { return lhs.h < rhs.h; });
//
//        for (const auto &skew: skewMorphisms.automorphisms) {
//            printMato(output, skew);
//        }
//
//        std::sort(skewMorphisms.properCosetPreserving.begin(), skewMorphisms.properCosetPreserving.end(), [](const auto& lhs, const auto &rhs) {return lhs.d < rhs.d; });
//
//        for (const auto &skew: skewMorphisms.properCosetPreserving) {
//            printMato(output, skew);
//        }
//
//        std::sort(skewMorphisms.properNonPreserving.begin(), skewMorphisms.properNonPreserving.end(), [](const auto& lhs, const auto &rhs) {return lhs.d < rhs.d; });
//
//        for (const auto &skew: skewMorphisms.properNonPreserving) {
//            printMato(output, skew);
//        }
//
//        output << "\nTotal number of skew morphisms of C_" << n << " is " << getSkewCount(skewMorphisms);
//        output << "\nCheck sub-total of automorphisms of C_" << n << " is " << number.phi << std::endl;
//
//        output << "\n.............................................................................." << std::endl;
//    }

    for (auto n = A; n <= N; ++n) {
        std::ofstream output("../html/skew/" + std::to_string(n) + ".html");
        const auto prev = n > 1 ? n - 1 : n;
        output <<
               "<html>\n"
               "  <head>\n"
               "    <script src='../jquery-3.6.4.min.js'></script>\n"
               "    <script src='../script.js'></script>\n"
               "    <link rel='stylesheet' href='../stylesheet.css'>\n"
               "  </head>\n"
               "  <body>\n"
               "    <form action='javascript:void(0);'>\n"
               "      <a id='prev' href='" << prev  << ".html'>" << prev << "</a>\n"
               "      <div>"
               "        <input type='checkbox' id='checkbox-auto' name='auto'><label for='checkbox-auto'>automorphisms</label>\n"
               "        <input type='checkbox' id='checkbox-coset' name='coset'><label for='checkbox-coset'>proper coset preserving</label>\n"
               "        <input type='checkbox' id='checkbox-other' name='other' checked><label for='checkbox-other'>proper coset non-preserving</label>\n"
               "        <br />\n"
               "        <input type='checkbox' id='checkbox-inv' name='inv'><label for='checkbox-inv'>power of -1 on orbit 1</label>\n"
               "        <br />\n"
               "        <input type='number' id='modulo' name='mod' value='" << n << "'> <label for='modulo'>modulo</label>\n"
               "      </div>\n"
               "      <a id='next' href='" << n + 1 << ".html'>" << n + 1 << "</a>\n"
               "    </form>\n"
                ;
        auto &number = numberCache[n];

        printHtml(output, number);
        output <<
               "  </body>\n"
               "</html>\n"
                ;
    }


    return 0;
}
