#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <ostream>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>
#include <numeric>
#include <CLI11.hpp>

#ifdef PROFILE_FLAG
#  define PROFILE __attribute__((noinline))
#else
#  define PROFILE
#endif

using Scalar = std::int32_t;
using DoubleScalar = std::int64_t;
using Index = Scalar;

using OrbitContainer = std::vector<Scalar>;
using OrbitsContainer = std::vector<OrbitContainer>;
using Function = std::vector<Scalar>;

struct OrbitPlace {
    Index orbitIndex = -1;
    Index indexOnOrbit = 0;
};
using OrbitPlaces = std::vector<OrbitPlace>;

struct Permutation {
    OrbitsContainer orbits;
    OrbitPlaces places;
};

struct CompactSkewMorphism {
    OrbitContainer orbit1;
    OrbitContainer pi;
};

class OrbitWrapper;

namespace old {
class OrbitWrapper {
public:
    friend class ::OrbitWrapper;

    class Iterator {
    public:
        Iterator(const OrbitContainer &orbit, Index index);
        Iterator(const OrbitContainer &orbit, Index index, std::size_t size);

        const Scalar &operator*() const;
        Iterator &operator++();

        friend bool operator!=(const Iterator &lhs, const Iterator &rhs);

    private:
        const OrbitContainer &orbit;
        std::size_t index;
        std::size_t size;
    };

    OrbitWrapper(const OrbitsContainer& orbits, const OrbitPlace &place);

    std::size_t size() const;
    Iterator begin() const noexcept;
    Iterator begin(std::size_t size) const noexcept;
    Iterator end() const noexcept;

    bool isSingleElement() const;

private:
    const OrbitsContainer &orbits;
    const OrbitPlace &place;
};

OrbitWrapper::Iterator::Iterator(const OrbitContainer &orbit, Index index) : Iterator(orbit, index, orbit.size()) {}

OrbitWrapper::Iterator::Iterator(const OrbitContainer &orbit, Index index, std::size_t size) : orbit(orbit), index(index), size(size) {}

const Scalar &OrbitWrapper::Iterator::operator*() const {
    return orbit[index];
}

OrbitWrapper::Iterator &OrbitWrapper::Iterator::operator++() {
    ++index;
    if (index == size) {
        index = 0;
    }
    return *this;
}

bool operator!=(const OrbitWrapper::Iterator &lhs, const OrbitWrapper::Iterator &rhs) {
    return lhs.index != rhs.index;
}

OrbitWrapper::OrbitWrapper(const OrbitsContainer &orbits, const OrbitPlace &place) : orbits(orbits), place(place) {}

std::size_t OrbitWrapper::size() const {
    return orbits[place.orbitIndex].size();
}

OrbitWrapper::Iterator OrbitWrapper::begin() const noexcept {
    return {orbits[place.orbitIndex], place.indexOnOrbit};
}

OrbitWrapper::Iterator OrbitWrapper::begin(std::size_t size) const noexcept {
    return {orbits[place.orbitIndex], place.indexOnOrbit, size};
}

OrbitWrapper::Iterator OrbitWrapper::end() const noexcept {
    return begin();
}

bool OrbitWrapper::isSingleElement() const {
    return place.orbitIndex == -1;
}
}

class OrbitWrapper {
public:
    friend class Iterator;
    class Iterator {
    public:
        Iterator(const OrbitWrapper &orbitWrapper, Index index, bool first);

        Scalar operator*() const;
        Iterator &operator++();

        friend bool operator!=(const Iterator &lhs, const Iterator &rhs);

    private:
        std::reference_wrapper<const OrbitWrapper> orbitWrapper;
        Index index;
        bool first;
    };

    OrbitWrapper(const OrbitContainer& orbit, Index indexOnOrbit, Scalar multiplier, Scalar exponent, Scalar n);
    OrbitWrapper(Scalar element, Scalar n);
    OrbitWrapper(const old::OrbitWrapper &orbitWrapper, Scalar n);

    OrbitWrapper &operator*=(Scalar m);
    OrbitWrapper &operator^=(Scalar e);
    Scalar operator[](std::size_t index) const;

    friend bool operator!=(const OrbitWrapper &lhs, const OrbitWrapper &rhs);
    friend bool operator<(const OrbitWrapper &lhs, const OrbitWrapper &rhs);

    explicit operator OrbitContainer() const;

    std::size_t size() const;
    Iterator begin() const noexcept;
    Iterator end() const noexcept;

    bool isSingleElement() const;

private:
    const OrbitContainer *orbit{nullptr};
    Index indexOnOrbit;
    Scalar multiplier{1};
    Scalar exponent{1};
    Scalar n;
};

PROFILE OrbitWrapper::Iterator::Iterator(const OrbitWrapper &orbitWrapper, Index index, bool first) : orbitWrapper(orbitWrapper), index(index), first(first) {}

PROFILE Scalar OrbitWrapper::Iterator::operator*() const {
    DoubleScalar value;
    if (orbitWrapper.get().orbit == nullptr) {
        value = index;
    } else {
        value = (*orbitWrapper.get().orbit)[index];
    }
    return (value * orbitWrapper.get().multiplier) % orbitWrapper.get().n;
}

PROFILE OrbitWrapper::Iterator &OrbitWrapper::Iterator::operator++() {
    if (orbitWrapper.get().orbit != nullptr) {
        index += orbitWrapper.get().exponent;
        if (index >= orbitWrapper.get().orbit->size()) {
            index -= orbitWrapper.get().orbit->size();
        }
    }
    first = false;
    return *this;
}

PROFILE bool operator!=(const OrbitWrapper::Iterator &lhs, const OrbitWrapper::Iterator &rhs) {
    return lhs.index != rhs.index || lhs.first != rhs.first || &lhs.orbitWrapper.get() != &rhs.orbitWrapper.get();
}

PROFILE OrbitWrapper::OrbitWrapper(const OrbitContainer &orbit, Index indexOnOrbit, Scalar multiplier, Scalar exponent, Scalar n) : orbit(&orbit), indexOnOrbit(indexOnOrbit), multiplier(multiplier), exponent(exponent), n(n) {
    exponent %= size();
}

PROFILE OrbitWrapper::OrbitWrapper(Scalar element, Scalar n) : indexOnOrbit(element), n(n) {}

PROFILE OrbitWrapper &OrbitWrapper::operator*=(Scalar m) {
    multiplier = (multiplier * static_cast<DoubleScalar>(m)) % n;
    return *this;
}

PROFILE OrbitWrapper &OrbitWrapper::operator^=(Scalar e) {
    exponent = (exponent * static_cast<DoubleScalar>(e)) % size();
    return *this;
}

PROFILE Scalar OrbitWrapper::operator[](std::size_t index) const {
    if (orbit == nullptr) {
        return *begin();
    }
    return *Iterator{*this, static_cast<Index>((indexOnOrbit + static_cast<DoubleScalar>(index) * exponent) % size()), true};
}

PROFILE bool operator!=(const OrbitWrapper &lhs, const OrbitWrapper &rhs) {
    if (lhs.size() != rhs.size()) {
        return true;
    }
    auto literator = lhs.begin();
    auto riterator = rhs.begin();
    while (literator != lhs.end()) {
        if (*literator != *riterator) {
            return true;
        }
        ++literator;
        ++riterator;
    }
    return false;
}

PROFILE bool operator<(const OrbitWrapper &lhs, const OrbitWrapper &rhs) {
    auto literator = lhs.begin();
    auto riterator = rhs.begin();
    while (literator != lhs.end() || riterator != rhs.end()) {
        if (*literator < *riterator) {
            return true;
        }
        if (*literator > *riterator) {
            return false;
        }
        ++literator;
        ++riterator;
    }
    if (literator != lhs.end()) {
        return false;
    }
    if (riterator != rhs.end()) {
        return true;
    }
    return false;
}

PROFILE OrbitWrapper::operator OrbitContainer() const {  // TODO: bez tohto a bez templatov
    OrbitContainer result;
    result.reserve(size());
    std::transform(begin(), end(), std::back_inserter(result), [](const auto &a) {return a;});
    return result;
}

PROFILE std::size_t OrbitWrapper::size() const {
    if (orbit == nullptr) {
        return 1;
    }
    return orbit->size();
}

PROFILE OrbitWrapper::Iterator OrbitWrapper::begin() const noexcept {
    return {*this, indexOnOrbit, true};
}

PROFILE OrbitWrapper::Iterator OrbitWrapper::end() const noexcept {
    return {*this, indexOnOrbit, false};
}

PROFILE bool OrbitWrapper::isSingleElement() const {
    return orbit == nullptr;
}

class SkewMorphism;

class WeakClassRepresentant {
public:
    WeakClassRepresentant(Permutation permutation, const SkewMorphism *quotient);

    const OrbitsContainer &orbits() const;
    const OrbitContainer &orbit1() const;
    Scalar orbit1(std::size_t index) const;
    old::OrbitWrapper orbitOf(Scalar element) const;
    OrbitWrapper quotientOrbitOf(Scalar element) const;

    Scalar max_orbits() const;

    Scalar n() const;
    Scalar r() const;
    Scalar d() const;
    Scalar ordAut() const;
    Scalar c() const;
private:
    Permutation _permutation;
    const SkewMorphism *quotient;
    Scalar _ordAut;
};

PROFILE const OrbitsContainer &WeakClassRepresentant::orbits() const {
    return _permutation.orbits;
}

PROFILE const OrbitContainer & WeakClassRepresentant::orbit1() const {
    return _permutation.orbits[0];
}

PROFILE Scalar WeakClassRepresentant::orbit1(std::size_t index) const {
    return orbit1()[index];
}

PROFILE old::OrbitWrapper WeakClassRepresentant::orbitOf(Scalar element) const {
    return {_permutation.orbits, _permutation.places[element]};
}

PROFILE OrbitWrapper::OrbitWrapper(const old::OrbitWrapper &orbitWrapper, Scalar n) : OrbitWrapper(orbitWrapper.place.indexOnOrbit, n)
{
    if (!orbitWrapper.isSingleElement()) {
        orbit = &orbitWrapper.orbits[orbitWrapper.place.orbitIndex];
    }
}

PROFILE Scalar WeakClassRepresentant::n() const {
    return _permutation.places.size();
}

PROFILE Scalar WeakClassRepresentant::r() const {
    return _permutation.orbits[0].size();
}

class SkewMorphism {
public:
    SkewMorphism(const WeakClassRepresentant &representant, Scalar a, Scalar e);

    OrbitWrapper orbit1() const;
    Scalar orbit1(std::size_t index) const;
    OrbitWrapper orbitOf(Scalar element) const;

    OrbitWrapper pi() const;
    Scalar pi(std::size_t index) const;

    Scalar n() const;
    Scalar d() const;
    Scalar ordAut() const;
    Scalar r() const;
    Scalar h() const;
    Scalar s() const;
    Scalar c() const;

    Scalar phi_1() const;

    bool inverseOrbit() const;

    std::unordered_map<Scalar, std::vector<Index>> roots;
    std::unordered_map<Scalar, Index> preservingSubgroups;
    bool powerOfInverseOrbit = false;
    std::uint64_t hash;  // TODO: checkni ci nie su rovnake

private:
    const WeakClassRepresentant &representant;
    Scalar a;
    Scalar aInv;
    Scalar e;
    Scalar eInv;
    Scalar _phi_1;
};

PROFILE WeakClassRepresentant::WeakClassRepresentant(Permutation permutation, const SkewMorphism *quotient) : _permutation(std::move(permutation)), quotient(quotient){
    if (quotient == nullptr || quotient->r() == 1) {
        _ordAut = r();
    } else {
        _ordAut = quotient->ordAut();
    }
}

PROFILE Scalar WeakClassRepresentant::d() const {
    if (quotient == nullptr) {
        return 1;
    }
    return quotient->orbit1().size();
}

Scalar WeakClassRepresentant::ordAut() const {
    return _ordAut;
}

PROFILE OrbitWrapper WeakClassRepresentant::quotientOrbitOf(Scalar element) const {
    if (quotient == nullptr) {
        return {0, 1};
    }
    return quotient->orbitOf(element);
}

PROFILE OrbitWrapper SkewMorphism::orbit1() const {
    return orbitOf(1);
}

PROFILE Scalar SkewMorphism::orbit1(std::size_t index) const {
    return orbit1()[index];
}

PROFILE OrbitWrapper SkewMorphism::orbitOf(Scalar element) const {
    OrbitWrapper result{representant.orbitOf((static_cast<DoubleScalar>(element) * a) % n()), n()};
    result *= aInv;
    result ^= e;
    return result;
}

PROFILE OrbitWrapper SkewMorphism::pi() const {
    auto result = representant.quotientOrbitOf(e % r());
    result *= eInv;
    result ^= a;
    return result;
}

PROFILE Scalar SkewMorphism::pi(std::size_t index) const {
    return pi()[index];
}

PROFILE Scalar SkewMorphism::n() const {
    return representant.n();
}

PROFILE Scalar SkewMorphism::d() const {
    return representant.d();
}

PROFILE Scalar SkewMorphism::ordAut() const {
    return representant.ordAut();
}

PROFILE Scalar SkewMorphism::r() const {
    return representant.r();
}

PROFILE Scalar SkewMorphism::h() const {
    const Scalar one = 1 % n();
    const Scalar one_pi = 1 % r();

    return orbit1(one_pi) - one;
}

PROFILE Scalar SkewMorphism::s() const {
    DoubleScalar _s = 0;
    for (const auto e: this->pi()) {
        _s += orbit1(e);
    }

    return (_s % this->n()) / d();
}

std::uint64_t hash(const CompactSkewMorphism &skew, Scalar n);

struct HashCompact {
    PROFILE std::size_t operator()(const CompactSkewMorphism &compactSkew) const {
        return hash(compactSkew, n);
    }

    Scalar n{0};
};

struct EqualCompact {
    PROFILE bool operator()(const CompactSkewMorphism &lhs, const CompactSkewMorphism &rhs) const {
        return lhs.orbit1 == rhs.orbit1 && lhs.pi == rhs.pi;
    }
};

//using SkewSet = std::unordered_set<CompactSkewMorphism, HashCompact, EqualCompact>;
using SkewIndexMap = std::unordered_map<CompactSkewMorphism, Index, HashCompact, EqualCompact>;

struct Class {
    Scalar size() const { return end - begin; }
    std::unique_ptr<WeakClassRepresentant> weakClassRepresentant;
    Index begin;
    Index end;
    Index representant;
};
using Classes = std::vector<Class>;

template<typename SkewType>
PROFILE bool isIdentity(const SkewType &skewMorphism) {
    return skewMorphism.r() == 1;
}

struct SkewMorphisms {
    std::vector<std::unique_ptr<SkewMorphism>> skews;
    std::vector<Classes> classes;
    SkewIndexMap skewIndexMap;
    Scalar automorphismsEnd;
    Scalar preservingEnd;
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

PROFILE void computePreservingSubgroups(SkewMorphism &skewMorphism) {
    const auto &number = numberCache[skewMorphism.n()];
    for (const auto d: number.divisors) {
        const auto &number_d = numberCache[skewMorphism.n() / d];
        const auto &orbit = skewMorphism.orbitOf(d % skewMorphism.n());
        const auto one = 1 % number_d.n;
        CompactSkewMorphism compact;
        if (orbit.isSingleElement()) {
            compact.orbit1 = {one};
            compact.pi = {0};
        } else {
            auto orbitIterator = ++orbit.begin();
            const auto orbitEnd = orbit.end();
            if (*orbitIterator % d != 0) {
                continue;
            }
            compact.orbit1.reserve(orbit.size());
            compact.orbit1.push_back(one);
            while (orbitIterator != orbitEnd) {
                compact.orbit1.push_back(*orbitIterator / d);
                ++orbitIterator;
            }
            const auto one_pi = 1 % orbit.size();
            compact.pi.push_back(one_pi);
            for (std::size_t i = d % skewMorphism.d(); skewMorphism.pi(i) % orbit.size() != one_pi; i = (i + d) % skewMorphism.d()) {
                compact.pi.push_back(skewMorphism.pi(i) % orbit.size());
            }
        }
        const auto index = number_d.skewMorphisms.skewIndexMap.find(compact)->second;
        skewMorphism.preservingSubgroups.insert({d, index});
    }
}

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
        number.skewMorphisms.skewIndexMap = SkewIndexMap{0, HashCompact{.n = i}};
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

template <typename IterableOrbit>
PROFILE std::string to_short_string(const IterableOrbit& orbit1, const IterableOrbit& pi) {
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
    return to_short_string(skewMorphism.orbit1(), skewMorphism.pi());
}

PROFILE void cleanup(const OrbitContainer &orbit1, Function &function) {
    for (const auto o: orbit1) {
        function[o] = 0;
    }
}

PROFILE void initialize(const OrbitContainer &orbit1, Function &function) {
    Scalar oldO = orbit1[0];
    for (std::size_t i = 1; i < orbit1.size(); ++i) {
        const auto o = orbit1[i];
        function[oldO] = -o;
        oldO = o;
    }
    function[oldO] = -orbit1[0];
    function[0] = 0;
}

PROFILE bool computeFunction(Scalar n, const OrbitContainer &pi, const OrbitContainer &orbit1, Function &function) {
    initialize(orbit1, function);

    Scalar value = 0;
    auto functionIterator = ++function.begin();
    auto piIterator = pi.begin();
    const auto fend = function.end();
    const auto pend = pi.end();
    while (functionIterator != fend) {
        value += orbit1[*piIterator];
        if (value >= n) value -= n;

        if (*functionIterator < 0 && *functionIterator != -value) {
            cleanup(orbit1, function);
            return false;
        }
        *functionIterator = value;

        ++functionIterator;
        ++piIterator;
        if (!(piIterator != pend)) {//TODO
            piIterator = pi.begin();
        }
    }
    return true;
}

PROFILE Function computeFunction(Scalar n, const OrbitContainer &pi, const OrbitContainer &orbit1) {
    Function function(n, 0);
    if (computeFunction(n, pi, orbit1, function)) {
        return function;
    }
    return {};
}

PROFILE Permutation computePermutation(const OrbitContainer &orbit1, const Function &function) {//TODO: zjednotit
    Permutation permutation;
    const auto n = static_cast<Scalar>(function.size());
    const auto r = static_cast<Scalar>(orbit1.size());
    permutation.orbits.reserve(n);
    permutation.places.resize(n);
    const auto one = 1 % n;
    for (Index i = 0; i < orbit1.size(); ++i) {
        auto &orbitPlaceC = permutation.places[orbit1[i]];
        orbitPlaceC.orbitIndex = 0;
        orbitPlaceC.indexOnOrbit = i;
    }
    permutation.orbits.emplace_back(orbit1);
    for (Scalar i = 0; i < n; ++i) {
        auto &orbitPlace = permutation.places[i];
        if (orbitPlace.orbitIndex != -1) {
            continue;
        }
        auto c = function[i];
        if (c == i) {
            orbitPlace.indexOnOrbit = i;
            continue;
        }
        const auto orbitIndex = permutation.orbits.size();
        permutation.orbits.emplace_back();
        auto &orbit = permutation.orbits.back();
        orbit.reserve(r);
        orbitPlace.orbitIndex = orbitIndex;
        orbitPlace.indexOnOrbit = orbit.size();
        orbit.push_back(i);
        do {
            auto &orbitPlaceC = permutation.places[c];
            orbitPlaceC.orbitIndex = orbitIndex;
            orbitPlaceC.indexOnOrbit = orbit.size();
            orbit.push_back(c);
            c = function[c];
        } while (c != i);
    }

    auto &orbits = permutation.orbits;
    std::sort(orbits.begin() + 1, orbits.end(), [](const auto &lhs, const auto &rhs) {
        if (lhs.size() == rhs.size()) {
            return lhs[0] < rhs[0];
        }
        return lhs.size() > rhs.size();
    });
    for (Index i = 0; i < orbits.size(); ++i) {
        for (const auto e: orbits[i]) {
            permutation.places[e].orbitIndex = i;
        }
    }

    permutation.orbits.shrink_to_fit();
    for (auto &orbit: permutation.orbits) {
        orbit.shrink_to_fit();
    }
    permutation.places.shrink_to_fit();

    return permutation;
}

PROFILE Permutation computePermutation(Scalar n, const OrbitContainer &pi, const OrbitContainer &orbit1) {
    return computePermutation(orbit1, computeFunction(n, pi, orbit1));
}

PROFILE CompactSkewMorphism toCompact(const SkewMorphism &skewMorphism) {
    return {OrbitContainer{skewMorphism.orbit1()}, OrbitContainer{skewMorphism.pi()}};
}

SkewMorphism &getSkewByIndex(SkewMorphisms &skewMorphisms, std::size_t index);
const SkewMorphism &quotient(const SkewMorphism &skew);
const SkewMorphism *quotientPtr(const CompactSkewMorphism &compact);

PROFILE SkewMorphism &powerOf(const SkewMorphism &skewMorphism, const std::pair<Scalar, Index> &sub, SkewMorphisms &skewMorphisms) {
    const auto e = sub.first;
    const auto r = skewMorphism.r() / e;
    const auto &orbit1 = skewMorphism.orbit1();
    const auto &roOther = getSkewByIndex(numberCache[r].skewMorphisms, sub.second);
    CompactSkewMorphism compactOther;
    compactOther.orbit1.reserve(r);
    for (std::size_t i = 0; i < skewMorphism.r(); i += e) {
        compactOther.orbit1.push_back(orbit1[i]);
    }
    compactOther.pi = OrbitContainer{roOther.orbit1()};
    return getSkewByIndex(skewMorphisms, skewMorphisms.skewIndexMap.find(compactOther)->second);
}

PROFILE SkewMorphism &mod(const SkewMorphism &skewMorphism, Scalar n) {
    auto &skewMorphisms = numberCache[n].skewMorphisms;
    if (n == 1) {
        return getSkewByIndex(skewMorphisms, 0);
    }
    const auto orbit1 = OrbitContainer{skewMorphism.orbit1()};
    const auto pi = OrbitContainer{skewMorphism.pi()};
    CompactSkewMorphism compactOther;
    compactOther.orbit1.push_back(1);
    if (orbit1.size() > 1) {
        for (std::size_t i = 1; i < orbit1.size() && orbit1[i] % n != 1; ++i) {
            compactOther.orbit1.push_back(orbit1[i] % n);
        }
    }
    const Scalar r = compactOther.orbit1.size();
    if (r > 1) {
        compactOther.pi.push_back(1);
        for (std::size_t i = 1; i < pi.size() && pi[i] % r != 1; ++i) {
            compactOther.pi.push_back(pi[i] % r);
        }
    } else {
        compactOther.pi.push_back(0);
    }
    compactOther.orbit1.shrink_to_fit();
    compactOther.pi.shrink_to_fit();
    return getSkewByIndex(skewMorphisms, skewMorphisms.skewIndexMap.find(compactOther)->second);
}

PROFILE std::uint64_t cyrb53(const std::vector<Scalar> &vector, std::uint32_t seed = 0) {
    std::uint32_t h1 = 0xdeadbeef ^ seed;
    std::uint32_t h2 = 0x41c6ce57 ^ seed;
    for (const auto element: vector) {
        h1 = (h1 ^ element) * 2654435761;
        h2 = (h2 ^ element) * 1597334677;
    }

    h1 = ((h1 ^ (h1 >> 16)) * 2246822507) ^ ((h2 ^ (h2 >> 13)) * 3266489909);
    h2 = ((h2 ^ (h2 >> 16)) * 2246822507) ^ ((h1 ^ (h1 >> 13)) * 3266489909);

    return 4294967296 * (2097151 & h2) + (h1 >> 0);
}

PROFILE std::vector<Scalar> toVector(const CompactSkewMorphism &skew, Scalar n) {
    std::vector<Scalar> vector = skew.orbit1;
    vector.push_back(n);
    vector.insert(vector.end(), skew.pi.begin(), skew.pi.end());
    return vector;
}

PROFILE std::vector<Scalar> toVector(const SkewMorphism &skew) {
    return toVector(toCompact(skew), skew.n());
}

PROFILE std::uint64_t hash(const SkewMorphism &skew) {
    return cyrb53(toVector(skew));
}

PROFILE std::uint64_t hash(const CompactSkewMorphism &skew, Scalar n) {
    return cyrb53(toVector(skew, n));
}

PROFILE void computeHashForClass(const Class &c, SkewMorphisms &skewMorphisms) {
    auto &skews = skewMorphisms.skews;
    const auto least = std::min_element(skews.begin() + c.begin, skews.begin() + c.end, [](const auto &lhs, const auto &rhs) {
        return toVector(*lhs) < toVector(*rhs);
    });
    const auto h = hash(**least);
    for (std::size_t i = c.begin; i < c.end; ++i) {
        skews[i]->hash = h;
    }
}

PROFILE void computeRoots(SkewMorphisms &skewMorphisms, Index index) {
    auto &skewMorphism = getSkewByIndex(skewMorphisms, index);
    computePreservingSubgroups(skewMorphism);
    const auto &ro = quotient(skewMorphism);
    for (const auto sub: ro.preservingSubgroups) {
        const auto e = sub.first;
        if (skewMorphism.r() % e != 0) {
            continue;
        }
        auto &other = powerOf(skewMorphism, sub, skewMorphisms);
        if (skewMorphism.powerOfInverseOrbit) {
            other.powerOfInverseOrbit = true;
        }
        other.roots[e].push_back(index);  // TODO: aktualne nevyuzivam
    }
}

PROFILE const WeakClassRepresentant &addClassRepresentant(const CompactSkewMorphism &compact, Permutation permutation, SkewMorphisms &skewMorphisms) {
    auto representant = std::make_unique<WeakClassRepresentant>(std::move(permutation), quotientPtr(compact));

    const auto c = representant->c();
    while (skewMorphisms.classes.size() <= c) {
        skewMorphisms.classes.emplace_back();
    }
    auto &classes = skewMorphisms.classes[c];
    const Index index = skewMorphisms.skews.size();
    classes.emplace_back(Class{std::move(representant), index, index, index});

    return *classes.back().weakClassRepresentant;
}

PROFILE void addSkewMorphism(SkewMorphism skewMorphism, SkewMorphisms &skewMorphisms) {
    const auto c = skewMorphism.c();
    while (skewMorphisms.classes.size() <= c) {
        skewMorphisms.classes.emplace_back();
    }
    auto &skews = skewMorphisms.skews;
    auto &classes = skewMorphisms.classes[c];
    auto &skewIndexMap = skewMorphisms.skewIndexMap;
    Index index = skews.size();

    skewIndexMap.insert({toCompact(skewMorphism), index});
    skews.emplace_back(std::make_unique<SkewMorphism>(std::move(skewMorphism)));

    computeRoots(skewMorphisms, index);

    classes.back().end = skews.size();
}

void addSkewClassByRepresentant(const CompactSkewMorphism &compact, Permutation permutation, Scalar n);

void addSkewClassByRepresentant(const CompactSkewMorphism &compact, Scalar n) {
    addSkewClassByRepresentant(compact, computePermutation(n, compact.pi, compact.orbit1), n);
}

std::size_t getPreservingSkewCount(const SkewMorphisms &skewMorphisms);//TODO: remove

PROFILE SkewMorphism::SkewMorphism(const WeakClassRepresentant &representant, Scalar a, Scalar e)
        : representant(representant)
        , a(a)
        , e(e) {
    aInv = (n() + inverse(n(), a)) % n();
    eInv = (r() + inverse(r(), e)) % r();
    powerOfInverseOrbit = inverseOrbit();

    _phi_1 = orbit1(1);
}

PROFILE Scalar SkewMorphism::c() const {
    return representant.c();
}

PROFILE Scalar SkewMorphism::phi_1() const {
    return _phi_1;
}

PROFILE Scalar WeakClassRepresentant::c() const {
    if (r() == 1) {
        return 0;
    }
    return quotient->c() + 1;
}

PROFILE Scalar WeakClassRepresentant::max_orbits() const {
    if (isIdentity(*this)) {
        return n();
    }
    Scalar result = 0;
    for (const auto &orbit: orbits()) {
        if (orbit.size() != r()) {
            return result;
        }
        ++result;
    }
    return result;
}

PROFILE bool SkewMorphism::inverseOrbit() const {
    if (isIdentity(*this)) {
        return true;
    }
    const auto inverse1 = n() - 1;
    for (const auto o: orbit1()) {
        if (o == inverse1) {
            return true;
        }
    }
    return false;
}

PROFILE void sEquals1(const Scalar d, const Number &number_n_div_d, SkewMorphisms &skewMorphisms) {
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
                        CompactSkewMorphism compact;
                        compact.orbit1.reserve(n_h);
                        compact.pi.reserve(d);
                        const auto n = d * number_n_div_d.n;
                        const auto gcd_h = n / number_n_h.n;
                        for (Scalar i = 0; i < n; i += gcd_h) {
                            compact.orbit1.push_back(i + 1);
                        }
                        for (Scalar id = 0; id < d; ++id) {
                            compact.pi.push_back(powers[id * step]);
                        }
                        addSkewClassByRepresentant(compact, n);
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
    auto &number_d = numberCache[d];
    computeCoprimes(number_d);
    std::vector<bool> visited_s(number_n_div_d.n, false);//TODO: global
    for (const auto &autoClass: number_n_div_d.skewMorphisms.classes[1]) {
        for (std::size_t autoIndex = autoClass.begin; autoIndex < autoClass.end; ++autoIndex) {
            const auto &powers_s = number_n_div_d.skewMorphisms.skews[autoIndex]->orbit1();
            const auto s = powers_s[1];
            if (visited_s[s]) {
                continue;
            }
            const auto b = s - 1;
            const auto gcd_b = gcd(number_n_div_d, numberCache[b]);
            const auto n_b = number_n_div_d.n / gcd_b;
            const auto coprime_b = b / gcd_b;
            if (gcd_b == number_n_div_d.n) {
                continue;
            }
            if (gcd_b < d) {//TODO: nejaky lepsi check? vnutri? toto skoro nic nerobi aj tak sa to dost skoro zahodi
                continue;
            }
            auto &number_n_b = numberCache[n_b];
            computeCoprimes(number_n_b);
            std::vector<bool> visited_e(number_n_div_d.n, false);//TODO: global
            std::vector<Scalar> powers_e;//TODO: global
            powers_e.reserve(number_n_div_d.n);
            auto power_sum = std::accumulate(powers_s.begin(), powers_s.end(), DoubleScalar{0});
            const Scalar rH = powers_s.size();
            auto &number_rH = numberCache[rH];
            computeCoprimes(number_rH);
            for (const auto coprime_rH: number_rH.coprimes) {
                visited_s[powers_s[coprime_rH]] = true;
            }
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
                                            CompactSkewMorphism compact;
                                            compact.orbit1.reserve(orbit_d1.size());
                                            compact.pi.reserve(d);
                                            std::transform(orbit_d1.begin(), orbit_d1.end(), std::back_inserter(compact.orbit1),
                                                           [coprime_h, n](const auto d_h) {
                                                               return (1 + d_h * coprime_h) % n;
                                                           });
                                            for (const auto coprime_e: number_d.coprimes) {
                                                compact.pi.clear();
                                                for (Scalar id = 0; id < d; ++id) {
                                                    compact.pi.push_back(powers_e[((coprime_e * id) % d) * step]);
                                                }

                                                addSkewClassByRepresentant(compact, n);
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
    computeCoprimes(number);
    const auto n = number.n;
    const auto one = 1 % n;
    std::vector<bool> visited_s(number.n, false);//TODO: global
    CompactSkewMorphism compact;
    std::vector<Scalar> powers_s;
    powers_s.reserve(number.n);
    compact.pi = {1};
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
        const Scalar r = powers_s.size();
        auto &number_r = numberCache[r];
        for (const auto rH: number_r.divisors) {
            const auto g = r / rH;
            if (visited_s[powers_s[g % r]]) {
                continue;
            }
            auto &number_rH = numberCache[rH];

            computeCoprimes(number_rH);
            for (const auto coprime_rH: number_rH.coprimes) {
                visited_s[powers_s[g * coprime_rH]] = true;
            }

            compact.orbit1.clear();
            for (Index i = 0; i < r; i += g) {
                compact.orbit1.push_back(powers_s[i]);
            }
            compact.pi[0] = 1 % rH;

            addSkewClassByRepresentant(compact, n);
        }
    }
//    clear(number);
}

template <typename IterableOrbit>
std::string to_string(const IterableOrbit &orbit) {
    std::string result;
    bool first = true;
    for (const auto e: orbit) {
        result += (first ? "" : ", ") + std::to_string(e);
        first = false;
    }
    return result;
}

std::string to_string(const SkewMorphism &skewMorphism) {
    if (isIdentity(skewMorphism)) {
        return "Id(sn)";
    }
    std::string result;
    auto orbits = computePermutation(skewMorphism.n(), OrbitContainer{skewMorphism.pi()}, OrbitContainer{skewMorphism.orbit1()}).orbits;
    for (const auto &orbit: orbits) {
        result += "(" + to_string(orbit) + ")";
    }
    return result;
}

template <typename IterableOrbit>
std::string to_json(const IterableOrbit &orbit) {
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
    auto orbits = computePermutation(skewMorphism.n(), OrbitContainer{skewMorphism.pi()}, OrbitContainer{skewMorphism.orbit1()}).orbits;
    for (const auto &orbit: orbits) {
        result += (first ? "" : ",") + to_json(orbit);
        first = false;
    }
    result += "]";
    return result;
}

PROFILE bool isSkewMorphism(const CompactSkewMorphism &compact, const OrbitsContainer &orbits, const Function& function) {
    const auto n = static_cast<Scalar>(function.size());
    bool good = true;
    Function c_j(n, 0);
    std::vector<bool> visited_j(compact.pi.size(), false);
    for (std::size_t ij = 0; ij < compact.pi.size(); ++ij) {
        const auto j = compact.pi[ij];
        for (const auto &orbit: orbits) {
            for (std::size_t oi = 0; oi < orbit.size(); ++oi) {
                c_j[orbit[oi]] = orbit[(oi + compact.orbit1.size() - j) % orbit.size()];  // TODO: compact.orbit1.size() je order?
            }
        }
        for (Scalar aa = 1; aa < n; ++aa) {//TODO: potrebujem?
            if (c_j[aa] == 0) {
                c_j[aa] = aa;
            }
        }
        std::size_t ijj = 0;
        for (; ijj < compact.pi.size(); ++ijj) {
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
        for (std::size_t i = ijj + compact.pi.size(); i < n; i += compact.pi.size()) {
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

PROFILE void periodicallyFillOrbit(std::size_t p, std::size_t i, const WeakClassRepresentant &power, OrbitContainer &t) {
    const std::size_t orbitSize = power.r();
    const auto e = t[i];

    const auto &orbit = power.orbitOf(e);

    auto orbitIterator = orbit.begin();
    for (std::size_t j = p + i; j < t.size(); j += p) {
        ++orbitIterator;
        t[j] = *orbitIterator;
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

PROFILE bool computeOrbit1(const Function &function, OrbitContainer &orbit1) {
    Scalar c = 1;
    std::size_t index = 0;
    do {
        orbit1[index] = c;
        c = function[c];
        ++index;
    } while (c != 1 && index < orbit1.size());
    return c == 1 && index == orbit1.size();
}

PROFILE OrbitContainer computeOrbit1(const Function &function) {
    OrbitContainer orbit1;
    Scalar c = 1;
    std::size_t index = 0;
    do {
        orbit1.push_back(c);
        c = function[c];
        ++index;
    } while (c != 1);
    return orbit1;
}

PROFILE bool compareOrbits(const OrbitContainer &sparseOrbit, const OrbitContainer &orbit2) {
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

PROFILE bool checkPowerCycles(const Scalar p, const WeakClassRepresentant &power, const OrbitContainer &orbit1) {
    for (std::size_t a = 0; a < p; ++a) {
        const auto &orbit = power.orbitOf(orbit1[a]);
        auto orbitIterator = orbit.begin();
        for (std::size_t i = a + p; i < orbit1.size(); i += p) {
            ++orbitIterator;
            if (orbit1[i] != *orbitIterator) {
                return false;
            }
        }
    }
    return true;
}

PROFILE bool checkFirstPMod(const SkewMorphism &ro, const OrbitContainer &orbit1) {
    for (std::size_t i = 1; i < ro.d(); ++i) {
        if (orbit1[i] % ro.r() != ro.pi(i)) {
            return false;
        }
    }
    return true;
}


PROFILE SkewMorphism &getSkewByIndex(SkewMorphisms &skewMorphisms, std::size_t index) {
    return *skewMorphisms.skews[index];
}

PROFILE std::size_t getSkewCount(const SkewMorphisms &skewMorphisms) {
    return skewMorphisms.skews.size();
}

PROFILE const Class &getClass(const SkewMorphisms &skewMorphisms, std::size_t index) {
    std::size_t i = 0;
    while (index >= skewMorphisms.classes[i].size()) {
        index -= skewMorphisms.classes[i].size();
        ++i;
    }
    return skewMorphisms.classes[i][index];
}

PROFILE std::size_t getClassesCount(const SkewMorphisms &skewMorphisms) {
    std::size_t result = 0;
    for (const auto &classes: skewMorphisms.classes) {
        result += classes.size();
    }
    return result;
}

PROFILE const SkewMorphism &getPreservingSkewByIndex(const SkewMorphisms &skewMorphisms, std::size_t index) {  // TODO: iterators
    return *skewMorphisms.skews[index];
}

PROFILE std::size_t getPreservingSkewCount(const SkewMorphisms &skewMorphisms) {
    return skewMorphisms.preservingEnd;
}

PROFILE std::size_t getPreservingClassIndex(const SkewMorphisms &skewMorphisms, std::size_t index) {
    if (index == 0) {
        return skewMorphisms.classes[0][0].begin;
    }
    index -= 1;
    if (index < skewMorphisms.classes[1].size()) {
        return skewMorphisms.classes[1][index].begin;
    }
    index -= skewMorphisms.classes[1].size();
    return skewMorphisms.classes[2][index].begin;
}

PROFILE std::size_t getPreservingClassesCount(const SkewMorphisms &skewMorphisms) {
    return 1 + skewMorphisms.classes[1].size() + skewMorphisms.classes[2].size();
}

PROFILE const SkewMorphism &getProperSkewByIndex(const SkewMorphisms &skewMorphisms, std::size_t index) {
    return *skewMorphisms.skews[skewMorphisms.automorphismsEnd + index];
}

PROFILE std::size_t getProperSkewCount(const SkewMorphisms &skewMorphisms) {
    return skewMorphisms.skews.size() - skewMorphisms.automorphismsEnd;
}

PROFILE std::size_t getProperClassIndex(const SkewMorphisms &skewMorphisms, std::size_t index) {
    std::size_t i = 2;
    while (index >= skewMorphisms.classes[i].size()) {
        index -= skewMorphisms.classes[i].size();
        ++i;
    }
    return skewMorphisms.classes[i][index].begin;
}

PROFILE std::size_t getProperClassesCount(const SkewMorphisms &skewMorphisms) {
    std::size_t result = 0;
    for (std::size_t i = 2; i < skewMorphisms.classes.size(); ++i) {
        result += skewMorphisms.classes[i].size();
    }
    return result;
}

PROFILE Function positionOnOrbit(const OrbitWrapper &orbit, Scalar n) {
    Function result(n, -1);
    for (Scalar x = 0; x < orbit.size(); ++x) {
        result[orbit[x]] = x;
    }
    return result;
}

PROFILE void initializeSplitIndex(OrbitContainer &t, std::vector<Index> &splitIndex, const std::vector<std::vector<Scalar>> &values, const std::vector<Index> &free_x_index, Scalar exponent, const WeakClassRepresentant &power) {
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

PROFILE void incrementSplitIndex(OrbitContainer &t, std::vector<Index> &splitIndex, const std::vector<std::vector<Scalar>> &values, const std::vector<Index> &free_x_index, Scalar exponent, const WeakClassRepresentant &power) {
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

PROFILE void clearOrbit(OrbitContainer &orbit) {
    for (auto &o: orbit) {
        o = 0;
    }
}

PROFILE void addSkewClassByRepresentant(const CompactSkewMorphism &compact, Permutation permutation, Scalar n) {
    const Scalar r = compact.orbit1.size();
    const Scalar d = compact.pi.size();

    auto &number = numberCache[n];
    auto &number_r = numberCache[r];

    const auto &skewIndexMap = number.skewMorphisms.skewIndexMap;
    if (skewIndexMap.find(compact) != skewIndexMap.end()) {
        return;
    }

    const auto &phiRepresentant = addClassRepresentant(compact, std::move(permutation), number.skewMorphisms);

    const Scalar one = 1 % n;
    const Scalar one_pi = 1 % r;

    addSkewMorphism(SkewMorphism(phiRepresentant, one, one_pi), number.skewMorphisms);

    computeCoprimes(number_r);
    std::vector<Scalar> coprimes_n;
    if (phiRepresentant.d() == 1) {
        coprimes_n = {one};
    } else if (phiRepresentant.c() == 2 && number.skewMorphisms.skews.back()->s() == 1) {
        computeCoprimes(numberCache[d]);
        coprimes_n = numberCache[d].coprimes;
        for (auto &c: coprimes_n) {
            while (!isCoprime(number, numberCache[c])) {
                c += d;
            }
        }
    } else {
        computeCoprimes(number);
        coprimes_n = number.coprimes;
    }
    for (const auto coprime_r: number_r.coprimes) {
        for (const auto coprime_n: coprimes_n) {
            auto phi2 = SkewMorphism(phiRepresentant, coprime_n, coprime_r);

            if (skewIndexMap.find({OrbitContainer{phi2.orbit1()}, OrbitContainer{phi2.pi()}}) != skewIndexMap.end()) {
                continue;
            }

            addSkewMorphism(std::move(phi2), number.skewMorphisms);
        }
    }
    auto &c = number.skewMorphisms.classes[phiRepresentant.c()].back();
    computeHashForClass(c, number.skewMorphisms);
}

PROFILE bool quotientEquals(const WeakClassRepresentant &skew, const SkewMorphism &quotient) {
    if (skew.quotientOrbitOf(1) != quotient.orbit1()) {
        return false;
    }
    const auto &orbit1 = skew.orbit1();
    const auto d = skew.d();
    for (std::size_t i = 0; i < quotient.d(); ++i) {
        if (orbit1[i] % d != quotient.pi(i)) {
            return false;
        }
    }
    return true;
}

void computeOdd(Number &number, Scalar c);
void computeEven(Number &number, Scalar c);

//TODO: zoradit podla radu
//TODO: kanonicke poradie autov, aj cosetov

//TODO: nemozu nejake nasobenia scitania pretiect?
PROFILE void computeProperNotPreserving(Number &number) {
    const auto n = number.n;
    const auto &number_nlambda = numberCache[n * number.lambda];

    Scalar maxC = 3;
    for (const auto m: number_nlambda.divisors) {
        if (m >= n) {
            continue;
        }
        const auto &number_m = numberCache[m];
        const Scalar m_maxC = number_m.skewMorphisms.classes.size();
        maxC = std::max(maxC, m_maxC);
    }
    for (Scalar c = 3; c <= maxC; ++c) {
        if (c % 2 == 1) {
            computeOdd(number, c);
        } else {
            computeEven(number, c);
        }
    }
}

std::set<std::pair<Scalar, Scalar>> complexities;

PROFILE void computeOdd(Number &number, Scalar c) {
    const auto n = number.n;
    const auto &number_nlambda = numberCache[n * number.lambda];

    const auto maxPrime = getMaxPrime(number);
    auto n_div_maxPrime = n / maxPrime;

    Function function(n, 0);

    CompactSkewMorphism compactSkewMorphism;
    const auto &skewIndexMap = number.skewMorphisms.skewIndexMap;

    for (const auto m: number_nlambda.divisors) {
        if (m >= n) {
            continue;
        }
        const auto &number_m = numberCache[m];
        if (isCoprime(number, number_m)) {
            continue;
        }
        const auto r = m;

        OrbitContainer t(r, 0);
        OrbitContainer &orbit1 = compactSkewMorphism.orbit1;
        orbit1.clear();
        orbit1.resize(r, 0);

        for (std::size_t ro_index = 0; ro_index < getProperSkewCount(number_m.skewMorphisms); ++ro_index) {
            const auto &ro = getProperSkewByIndex(number_m.skewMorphisms, ro_index);

            if (ro.c() + 1 != c) {
                continue;
            }

            const auto d = ro.r();

            if (n % (d * maxPrime) != 0) {
                continue;
            }

            const auto n_div_d = n / d;

            if (isCoprime(numberCache[n_div_d], number_m)) {
                continue;
            }

            const auto p = ro.d();

            const auto &orbit1_ro = OrbitContainer{ro.orbit1()};

            compactSkewMorphism.pi = orbit1_ro;  // TODO: do not copy

            const auto exponent = ro.ordAut();
            const auto ord_power = m / exponent;

            clearOrbit(t);

            const auto &powerQuotient = getSkewByIndex(numberCache[ord_power].skewMorphisms, ro.preservingSubgroups.find(exponent)->second);

            for (std::size_t power_class_index = 0; power_class_index < getClassesCount(number.skewMorphisms); ++power_class_index) {
                const auto &powerClass = getClass(number.skewMorphisms, power_class_index);
                const auto &power = *powerClass.weakClassRepresentant;

                if (power.r() != ord_power) {
                    continue;
                }
//                if (power.orbit1(p_exponent) % d != 1) {
//                    continue;
//                }
                if (power.max_orbits() < exponent) {
                    continue;
                }

                if (!quotientEquals(power, powerQuotient)) {
                    continue;
                }

                for (Scalar a = ro.pi(1); a < n; a += d) {
                    const auto &orbit_a = power.orbitOf(a);
                    if (orbit_a.isSingleElement() || orbit_a.size() != ord_power) {
                        continue;
                    }
                    t[1] = a;
                    periodicallyFillOrbit(exponent, 1, power, t);
                    const auto &pi = orbit1_ro;
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
                    if (!checkPowerCycles(exponent, power, orbit1)) {
                        continue;
                    }

                    if (skewIndexMap.find(compactSkewMorphism) != skewIndexMap.end()) {
                        continue;
                    }

                    auto permutation = computePermutation(orbit1, function);

                    if (!isSkewMorphism(compactSkewMorphism, permutation.orbits, function)) {
                        continue;
                    }

                    complexities.emplace(c, power.c());
                    addSkewClassByRepresentant(compactSkewMorphism, std::move(permutation), n);
                }
            }
        }
    }
    //printf("\r");
//    printf("%d %ld\n", n, foundSkews.size());
}

PROFILE std::vector<std::pair<Scalar, OrbitContainer>>
bbb(const Scalar n, const Scalar e, const Scalar n_div_e, const WeakClassRepresentant &psi) {
    std::vector<std::pair<Scalar, OrbitContainer>> ts;
    const auto psi_r = psi.r();

    Function psiFunction = computeFunction(n_div_e, OrbitContainer{psi.quotientOrbitOf(1)}, psi.orbit1());
    Function psiPower(n_div_e);
    std::iota(std::begin(psiPower), std::end(psiPower), 0);


    for (Scalar ro_1_mod_psi_r = 0; ro_1_mod_psi_r < psi_r; ++ro_1_mod_psi_r) {
        for (Scalar k = 1; k < n_div_e; ++k) {
            auto b = k;
            OrbitContainer t;
            t.reserve(n);
            t.push_back(1);
//                    t.resize(1);
            for (std::size_t i = 1; i < n; ++i) {
                t.push_back(b * e + 1);
                b = (psiPower[b] + k);
                if (b >= n_div_e) b -= n_div_e;
                if (b == 0) {
                    break;
                }
            }
            const Scalar r = t.size();
            if (r == n) {
                continue;
            }
            ts.emplace_back(std::make_pair(ro_1_mod_psi_r, std::move(t)));
        }
        for (std::size_t i = 0; i < n_div_e; ++i) {
            psiPower[i] = psiFunction[psiPower[i]];
        }
    }
    std::sort(ts.begin(), ts.end(),[](const auto &lhs, const auto &rhs) { return lhs.second.size() < rhs.second.size(); });
    return ts;
}

PROFILE void computeEven(Number &number, Scalar c) {
    const auto n = number.n;
    const auto &number_nlambda = numberCache[n * number.lambda];

    const auto maxPrime = getMaxPrime(number);
    auto n_div_maxPrime = n / maxPrime;

    Function function(n, 0);

    CompactSkewMorphism compactSkewMorphism;
    OrbitContainer &orbit1 = compactSkewMorphism.orbit1;
    orbit1.reserve(n);
    const auto &skewIndexMap = number.skewMorphisms.skewIndexMap;

//    OrbitContainer t;
//    t.reserve(n);
//    t.push_back(1);

    for (const auto e: number.divisors) {
        if (e == 1 || e == n) {
            continue;
        }
        const auto n_div_e = n / e;
        const auto &number_n_div_e = numberCache[n_div_e];

        for (std::size_t psi_class_index = 0; psi_class_index < getClassesCount(number_n_div_e.skewMorphisms); ++psi_class_index) {
            const auto &psiClass = getClass(number_n_div_e.skewMorphisms, psi_class_index);
            const auto &psi = *psiClass.weakClassRepresentant;
            const auto psi_r = psi.r();

            std::vector<std::pair<Scalar, OrbitContainer>> ts = bbb(n, e, n_div_e, psi);

            for (const auto &pair: ts) {
                const auto ro_1_mod_psi_r = pair.first;
                const auto &t = pair.second;
                const Scalar r = t.size();
                auto &number_r = numberCache[r];
                const auto psi_1 = t[1];

                if (number_r.skewMorphisms.classes.size() < c) {
                    continue;
                }

                for (const auto &roClass: number_r.skewMorphisms.classes[c - 1]) {
                    for (auto ro_index = roClass.begin; ro_index < roClass.end; ++ro_index) {
                        const auto &ro = getSkewByIndex(number_r.skewMorphisms, ro_index);

                        if (ro.ordAut() != e) {
                            continue;
                        }

                        const auto d = ro.r();

                        if (n % (d * maxPrime) != 0) {
                            continue;
                        }

                        const auto n_div_d = n / d;

                        if (isCoprime(numberCache[n_div_d], number_r)) {
                            continue;
                        }

                        if (ro_1_mod_psi_r != ro.phi_1() % psi_r) {
                            continue;
                        }

                        if (!checkFirstPMod(ro, t)) {
                            continue;
                        }
                        const auto &pi = OrbitContainer{ro.orbit1()};
                        if (!computeFunction(n, pi, t, function)) {
                            continue;
                        }
                        if (!isPermutation(function)) {
                            continue;
                        }

                        orbit1.resize(r);
                        if (!computeOrbit1(function, orbit1)) {
                            continue;
                        }

                        //TODO: ktore checky treba robit?
                        if (!compareOrbits(t, orbit1)) {
                            continue;
                        }
//                        if (!checkPowerCycles(exponent, power, orbit1)) {
//                            continue;
//                        }

                        compactSkewMorphism.pi = pi;  // TODO: do not copy
                        if (skewIndexMap.find(compactSkewMorphism) != skewIndexMap.end()) {
                            continue;
                        }

                        auto permutation = computePermutation(orbit1, function);

                        if (!isSkewMorphism(compactSkewMorphism, permutation.orbits, function)) {
                            continue;
                        }
                        complexities.emplace(c, psi.c());

                        addSkewClassByRepresentant(compactSkewMorphism, std::move(permutation), n);
                    }
                }
            }
        }
    }
    //printf("\r");
//    printf("%d %ld\n", n, foundSkews.size());
}

struct MatoComparator {
    bool operator()(const SkewMorphism &lhs, const SkewMorphism &rhs) {
        if (lhs.n() != rhs.n()) {
            return lhs.n() < rhs.n();
        }
        if (std::min(lhs.c(), 3) != std::min(rhs.c(), 3)) {
            return lhs.c() < rhs.c();
        }
        if (lhs.d() != rhs.d()) {
            return lhs.d() < rhs.d();
        }
        if (lhs.c() < 3) {
            if (lhs.s() != rhs.s()) {
                return lhs.s() < rhs.s();
            }
            if (lhs.d() == 1) {
                return false;
            }
            if (lhs.h() != rhs.h()) {
                return lhs.h() < rhs.h();
            }
            if (lhs.pi(1) != rhs.pi(1)) {
                return lhs.pi(1) < rhs.pi(1);
            }
            return false;
        }
        if (lhs.r() != rhs.r()) {
            return lhs.r() < rhs.r();
        }
        auto &number_r = numberCache[lhs.r()];
        const auto &lq = quotient(lhs);
        const auto &rq = quotient(rhs);
//        auto lleastq = lq;
//        auto rleastq = rq;
//        Scalar lc = 1;
//        Scalar rc = 1;
//        auto &lnumber_r = numberCache[lq.r];
//        computeCoprimes(lnumber_r);
//        for (const auto coprime_r: lnumber_r.coprimes) {
//            const auto r = lq.r;
//            const auto d = lq.d;
//            const auto coprime_r_inverse = (r + inverse(r, coprime_r)) % r;
//            CompactSkewMorphism compact;
//            auto &orbit1 = compact.orbit1;
//            auto &orbit1_ro = compact.pi;
//            orbit1.resize(r, 1);
//            orbit1_ro.resize(d, 1);
//            const auto &orbit = getOrbit1(lq);
//            std::size_t index = 0;
//            for (std::size_t i = 1; i < r; ++i) {
//                index += coprime_r;
//                if (index >= r) index -= r;
//                orbit1[i] = orbit[index];
//            }
//            const auto &ro = quotient(lq);
//            const auto &coprime_r_place = ro.permutation.places[coprime_r];
//            const auto &coprime_r_orbit = ro.permutation.orbits[coprime_r_place.orbitIndex];
//            index = coprime_r_place.indexOnOrbit;
//            for (std::size_t i = 1; i < d; ++i) {
//                ++index;
//                if (index >= d) index -= d;
//                orbit1_ro[i] = (coprime_r_orbit[index] * coprime_r_inverse) % r;
//            }
//
//            const auto &lq2 = number_r.skewMorphisms.skews[number_r.skewMorphisms.skewIndexMap.find(compact)->second];
//            if (operator()(*lq2, lleastq)) {
//                lleastq = *lq2;
//                lc = coprime_r_inverse;
//            }
//        }
//        auto &rnumber_r = numberCache[rq.r];
//        computeCoprimes(rnumber_r);
//        for (const auto coprime_r: rnumber_r.coprimes) {
//            const auto r = rq.r;
//            const auto d = rq.d;
//            const auto coprime_r_inverse = (r + inverse(r, coprime_r)) % r;
//            CompactSkewMorphism compact;
//            auto &orbit1 = compact.orbit1;
//            auto &orbit1_ro = compact.pi;
//            orbit1.resize(r, 1);
//            orbit1_ro.resize(d, 1);
//            const auto &orbit = getOrbit1(rq);
//            std::size_t index = 0;
//            for (std::size_t i = 1; i < r; ++i) {
//                index += coprime_r;
//                if (index >= r) index -= r;
//                orbit1[i] = orbit[index];
//            }
//            const auto &ro = quotient(rq);
//            const auto &coprime_r_place = ro.permutation.places[coprime_r];
//            const auto &coprime_r_orbit = ro.permutation.orbits[coprime_r_place.orbitIndex];
//            index = coprime_r_place.indexOnOrbit;
//            for (std::size_t i = 1; i < d; ++i) {
//                ++index;
//                if (index >= d) index -= d;
//                orbit1_ro[i] = (coprime_r_orbit[index] * coprime_r_inverse) % r;
//            }
//
//            const auto &rq2 = number_r.skewMorphisms.skews[number_r.skewMorphisms.skewIndexMap.find(compact)->second];
//            if (operator()(*rq2, rleastq)) {
//                rleastq = *rq2;
//                rc = coprime_r_inverse;
//            }
//        }
//        if (operator()(lleastq, rleastq)) {
//            return true;
//        }
//        if (operator()(rleastq, lleastq)) {
//            return false;
//        }
        if (operator()(lq, rq)) {
            return true;
        }
        if (operator()(rq, lq)) {
            return false;
        }
        auto &skewMorphisms = numberCache[lhs.n()].skewMorphisms;
        const auto &lcp = powerOf(lhs, *lq.preservingSubgroups.find(lq.d()), skewMorphisms);
        const auto &rcp = powerOf(rhs, *rq.preservingSubgroups.find(rq.d()), skewMorphisms);
        if (operator()(lcp, rcp)) {
            return true;
        }
        if (operator()(rcp, lcp)) {
            return false;
        }

        return lhs.orbit1() < rhs.orbit1();
    }
};

PROFILE void countSkewmorphisms(Number &number) {
    if (getSkewCount(number.skewMorphisms) != 0) {
        return;
    }
    const auto n = number.n;
    const auto phi = number.phi;
    number.skewMorphisms.classes.resize(4);
    computeAutomorphisms(number);
    number.skewMorphisms.automorphismsEnd = number.skewMorphisms.skews.size();

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
        number.skewMorphisms.preservingEnd = number.skewMorphisms.skews.size();

        if (number.powerOfTwo > 16 || !number.squareFree) {
            computeProperNotPreserving(number);
        }
    } else {
        number.skewMorphisms.preservingEnd = number.skewMorphisms.skews.size();
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
    const auto nclasses_a = 1 + number.skewMorphisms.classes[1].size();
    const auto nclasses_c = number.skewMorphisms.classes[2].size();
    const auto nclasses_o = getProperClassesCount(number.skewMorphisms) - nclasses_c;
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

PROFILE CompactSkewMorphism fromString(const std::string &string, Scalar n) {
    CompactSkewMorphism result;
    std::istringstream stream;
    stream.str(string);
    char c = '\0';
    while (c != ')') {
        stream.get(c);
        Scalar i = 0;
        stream >> i;
        result.orbit1.push_back(i % n);
        stream.get(c);
    }
    while (c != ']') {
        stream.get(c);
        Scalar i = 0;
        stream >> i;
        result.pi.push_back(i % result.orbit1.size());
        stream.get(c);
    }
    return result;
}

PROFILE CompactSkewMorphism compactQuotient(const CompactSkewMorphism &compact) {
    const auto &orbit1 = compact.orbit1;
    OrbitContainer ro_pi;
    const auto one = 1 % compact.pi.size();
    ro_pi.push_back(one);
    if (orbit1.size() > 1) {
        for (std::size_t i = 1; orbit1[i] % compact.pi.size() != one; ++i) {
            ro_pi.push_back(orbit1[i] % compact.pi.size());
        }
    }
    ro_pi.shrink_to_fit();
    return {compact.pi, ro_pi};
}

PROFILE CompactSkewMorphism compactQuotient(const SkewMorphism &skew) {
    return compactQuotient(toCompact(skew));
}

PROFILE const SkewMorphism &quotient(const CompactSkewMorphism &compact) {
    const auto quotient = compactQuotient(compact);
    const auto &number_r = numberCache[compact.orbit1.size()];
    const auto index = number_r.skewMorphisms.skewIndexMap.find(quotient)->second;
    return *number_r.skewMorphisms.skews[index];
}

PROFILE const SkewMorphism &quotient(const SkewMorphism &skew) {
    return quotient(toCompact(skew));
}

PROFILE const SkewMorphism *quotientPtr(const CompactSkewMorphism &compact) {
    if (compact.orbit1[0] == 0) {
        return nullptr;
    }
    return &quotient(compact);
}

bool readSkewMorphisms(Scalar n, SkewMorphisms &skewMorphisms) {
    std::string line;
    std::ifstream ifile;

    if (skewMorphisms.classes.size() < 4) {
        skewMorphisms.classes.resize(4);
    }

    ifile = std::ifstream{"../skew/" + std::to_string(n) + "_auto.txt"};
    if (!ifile.is_open()) {
        return false;
    }
    while (std::getline(ifile, line)) {
        auto skew = fromString(line, n);
        addSkewClassByRepresentant(skew, n);
    }
    skewMorphisms.automorphismsEnd = skewMorphisms.skews.size();
    ifile = std::ifstream{"../skew/" + std::to_string(n) + "_coset.txt"};
    while (std::getline(ifile, line)) {
        auto skew = fromString(line, n);
        addSkewClassByRepresentant(skew, n);
    }
    skewMorphisms.preservingEnd = skewMorphisms.skews.size();

    ifile = std::ifstream{"../skew/" + std::to_string(n) + "_other.txt"};
    while (std::getline(ifile, line)) {
        auto skew = fromString(line, n);
        addSkewClassByRepresentant(skew, n);
    }
    return true;
}

std::string oldSplit(const std::string &line) {
    std::string result;
    std::stringstream ss(line);
    std::string w;
    std::string last;
    std::string d;
    while (getline(ss, w, ' ')) {
        if (last.size() + d.size() <= 80) {
            last += d;
        } else {
            result += "\n" + last;
            last = "";
            d = "";
        }
        if (last.size() + w.size() <= 80) {
            last += w;
            d = " ";
        } else {
            result += "\n" + last;
            last = w;
        }
    }
    if (!last.empty()) {
        result += "\n" + last;
    }
    return result;
}

void printMato(std::ofstream &output, const SkewMorphism &skew, bool old = false) {
    if (!old) {
        output << "\nSkew morphism " << to_string(skew);
        output << "\n of order " << std::to_string(skew.r());
    } else {
        std::string line = "Skew morphism ";
        if (isIdentity(skew)) {
            line += "Id(sn)";
        } else {
            auto orbits = computePermutation(skew.n(), OrbitContainer{skew.pi()}, OrbitContainer{skew.orbit1()}).orbits;
            std::sort(orbits.begin(), orbits.end(), [](const auto &lhs, const auto &rhs) {return lhs[0] < rhs[0];});
            for (const auto &orbit: orbits) {
                line += "(" + to_string(orbit) + ")";
            }
        }
        output << oldSplit(line);
    }
    output << "\n with kernel of order " << std::to_string(skew.n() / skew.d());
    if (!old) {
        output << "\n and smallest kernel generator " << std::to_string(skew.d());
        if (skew.d() == 1) {
            output << "\n and power function values [ 1 ]" << std::endl;
        } else {
            output << "\n and power function values [ " << to_string(skew.pi()) << " ]" << std::endl;
        }
        output << "\n  and with periodicity " << std::to_string(quotient(skew).d());
    } else {
        std::string line = " and power function values [ ";
        if (skew.d() == 1) {
            line += "1";
            for (int i = 1; i < skew.n(); ++i) {
                line += ", 1";
            }
        } else {
            line += to_string(skew.pi());
            for (int i = 1; i < skew.n() / skew.d(); ++i) {
                line += ", " + to_string(skew.pi());
            }
        }
        line += " ]";
        output << oldSplit(line) << std::endl;
    }
}

void printHtml(std::ofstream &output, const SkewMorphism &skew, const std::string &additionalClass, Scalar classSize) {
    const auto &q = quotient(skew);
    output << "\n<li class='skew " << additionalClass <<" n-" << skew.n() << " d-" << skew.d() << " h-" << skew.h() << " r-" << skew.r() << " s-" << skew.s() << " c-" << skew.c();
    if (skew.inverseOrbit()) {
        output << " inv-1";
    }
    if (skew.powerOfInverseOrbit) {
        output << " pow-inv-1";
    }
    auto &number = numberCache[skew.n()];
    for (const auto &sub: q.preservingSubgroups) {
        const auto &power = powerOf(skew, sub, number.skewMorphisms);
        output << " roots-" << power.hash;
    }
    output << "' id='" << skew.hash << "' data-repr='" << to_json(skew) << "' data-pi='" << to_json(skew.pi()) << "'><pre>";
    output << "\nSkew morphsim class of size " << classSize << " <a class='link' href='#" << skew.hash << "'>#</a>";
    output << "\n  <span class='repr'></span>";
    output << "\n  of order " << std::to_string(skew.r());
    output << " <a href='?auto=true&coset=true&other=true&inv=true#roots-" << skew.hash << "'>&radic;</a>";
    output << "\n  with kernel of order " << std::to_string(skew.n() / skew.d());
    output << "\n  and smallest kernel generator " << std::to_string(skew.d());
    if (skew.d() == 1) {
        output << "\n  and power function values [ 1 ]";
    } else {
        output << "\n  and power function values [ " << to_string(skew.pi()) << " ]";
    }
    output << " <a href='" << skew.r() << ".html?auto=true&coset=true&other=true&inv=true#" << q.hash << "'>#</a>";
    output << "\n  and with periodicity " << std::to_string(q.d());

    output << "\n\n</pre>" << std::endl;
    output << "\n</li>" << std::endl;
}

std::string factorizationHtml(Scalar n) {
    if (n == 1) {
        return "1";
    }
    std::string result;
    const auto &number = numberCache[n];
    for (std::size_t i = 0; i < number.primes.size(); ++i) {
        const Scalar prime = number.primes[i];
        const Scalar exponent = numberCache[number.powersOfPrimes[i]].divisors.size() - 1;
        if (i > 0) {
            result += "&centerdot;";
        }
        result += std::to_string(prime);
        if (exponent > 1) {
            result += "<sup>" + std::to_string(exponent) + "</sup>";
        }
    }
    return result;
}

void printHtml(std::ofstream &output, const Number &number) {
    const auto n = number.n;
    auto &skewMorphisms = number.skewMorphisms;
    output << "\n<div class='n' id='n-" << n << "'>";
    output << "\n<pre>";
    output << "\nNumber of selected skew morphism classes of C_" << n << " is <span class='count'>0</span>";
    output << "\nTotal number of skew morphisms of C_" << n << " is " << getSkewCount(skewMorphisms);
    output << "\nTotal number of skew morphism classes of C_" << n << " is " << getClassesCount(skewMorphisms);
    output << "\nTotal number of automorphisms of C_" << n << " is " << number.phi;
    output << "\n</pre>";

    output << "\n<ol>";
    for (std::size_t i = 0; i < 2; ++i) {
        for (const auto &c: number.skewMorphisms.classes[i]) {
            const auto &skew = *number.skewMorphisms.skews[c.representant];
            printHtml(output, skew, "auto", c.size());
        }
    }
    for (const auto &c: number.skewMorphisms.classes[2]) {
        const auto &skew = *number.skewMorphisms.skews[c.representant];
        printHtml(output, skew, "coset", c.size());
    }
    for (std::size_t i = 3; i < number.skewMorphisms.classes.size(); ++i) {
        for (const auto &c: number.skewMorphisms.classes[i]) {
            const auto &skew = *number.skewMorphisms.skews[c.representant];
            printHtml(output, skew, "other", c.size());
        }
    }
    output << "\n</ol>";

    output << "\n</div>" << std::endl;
}

int main(int argc, char *argv[]) {
    CLI::App app{"Compute list of skew morphisms"};

    Scalar N = 100;
    app.add_option("N", N, "Compute up to");
    Scalar A = 1;
    app.add_option("A", A, "Compute from");
    bool recompute = false;
    app.add_flag("--recompute", recompute, "Recompute instead of loading");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
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
            if (isCoprime(numberCache[a], numberCache[d]) && d != 1) {
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
            if (recompute || !readSkewMorphisms(n, number.skewMorphisms)) {
                countSkewmorphisms(number);
            }
//            fprint(number, std::cerr);
//        }
        if (getSkewCount(number.skewMorphisms) != number.phi) {
            print(number);
        }


        const auto cubeFree = [](const Number &number) -> bool {
            for (std::size_t i = 0; i < number.primes.size(); ++i) {
                if (number.primes[i] == 2) {
                    continue;
                }
                if (number.powersOfPrimes[i] > number.primes[i] * number.primes[i]) {
                    return false;
                }
                if (number.powersOfPrimes[i] > number.primes[i] && !isCoprime(number, numberCache[number.primes[i] - 1])) {
                    return false;
                }
            }
            return true;
        };

        if ((number.skewMorphisms.classes.size() == 4) != (number.powerOfTwo < 64 && cubeFree(number))) {
            printf("!!!! %d\n", n);
        }

        auto found = std::vector<bool>(number.skewMorphisms.skews.size(), false);
        for (std::size_t i = 0; i < number.skewMorphisms.skews.size(); ++i) {
            found[i] = number.skewMorphisms.skews[i]->c() < 2;
        }
        for (const auto &beta_: number.skewMorphisms.skews) {
            if (n == 1) { //!!!!
                continue;
            }
            const auto &beta = *beta_;
            const auto ord_beta = beta.r();
            for (Scalar m = 2; m <= n / ord_beta; ++m) {//!!!! m >= 2
                for (const auto &alpha_: numberCache[m * ord_beta].skewMorphisms.skews) {
                    const auto &alpha = * alpha_;
                    if (alpha.c() % 2 != 0 || alpha.ordAut() != m) {
                        continue;
                    }
                    for (Scalar h = 0; h < n; ++h) {
                        Function phi_function(n, 0);
                        for (std::size_t i = 1; i <= n; ++i) {
                            phi_function[i % n] = (phi_function[i - 1] + beta.orbitOf(h)[(alpha.orbit1(i - 1) - 1) / m]) % n;
                        }

                        bool a = true;
                        auto visited = std::vector<bool>(n, false);
                        Scalar t = 1;
                        for (Scalar x = 0; !visited[t] || x % alpha.pi().size() != 0; ++x) {
                            visited[t] = true;
                            if ((t % alpha.r()) != alpha.pi(x)) {
                                a = false;
                                break;
                            }
                            t = phi_function[t];
                        }

                        const auto &powerQuotient = getSkewByIndex(numberCache[ord_beta].skewMorphisms, alpha.preservingSubgroups.find(m)->second);
                        bool b = &quotient(beta) == &powerQuotient;

                        bool c = true;

                        for (Scalar x = 0; x < n; ++x) {
                            Scalar phi_m_x = x;
                            for (Scalar i = 0; i < m; ++i) {
                                phi_m_x = phi_function[phi_m_x];
                            }
                            if (phi_m_x != beta.orbitOf(x)[1]) {
                                c = false;
                                break;
                            }
                        }

                        if (!a) {
                            continue;
                        }
                        if (!b) {
                            continue;
                        }
                        if (!c) {
                            continue;
                        }

                        if (!isPermutation(phi_function)) {
                            std::cout << "perm" << std::endl;
                            std::cout << "n: " << n << ", m: " << m << ", h: " << h << std::endl;
                            std::cout << "alpha: " << to_json(alpha) << std::endl;
                            std::cout << "beta: " << to_json(beta) << std::endl;
                            std::cout << "phi_function: " << to_json(phi_function) << std::endl;
                            throw "";
                        }
                        const auto orbit1 = computeOrbit1(phi_function);
                        const auto permutation = computePermutation(orbit1, phi_function);

                        if (!isSkewMorphism({orbit1, OrbitContainer{alpha.orbit1()}}, permutation.orbits, phi_function)) {
                            std::cout << "skew" << std::endl;
                            std::cout << "n: " << n << ", m: " << m << ", h: " << h << std::endl;
                            std::cout << "alpha: " << to_json(alpha) << std::endl;
                            std::cout << "beta: " << to_json(beta) << std::endl;
                            std::cout << "phi_function: " << to_json(phi_function) << std::endl;
                            throw "";
                        }
                        const auto index = number.skewMorphisms.skewIndexMap[{orbit1, OrbitContainer{alpha.orbit1()}}];
                        if (found[index]) {
                            std::cout << "perm" << std::endl;
                            std::cout << "n: " << n << ", m: " << m << ", h: " << h << std::endl;
                            std::cout << "alpha: " << to_json(alpha) << std::endl;
                            std::cout << "beta: " << to_json(beta) << std::endl;
                            std::cout << "phi_function: " << to_json(phi_function) << std::endl;
                            throw "";
                        }
                        found[index] = true;
                    }
                }
            }
        }
        for (const Scalar m : number.divisors) {
            if (m == 1) {
                continue;
            }
            for (const auto &beta_: numberCache[n / m].skewMorphisms.skews) {
                const auto &beta = *beta_;
                for (Scalar l = 1; l < n; ++l) {
                    for (const auto &alpha_ : numberCache[l].skewMorphisms.skews) {
                        const auto &alpha = *alpha_;
                        if (alpha.c() % 2 != 1 || alpha.ordAut() != m) {
                            continue;
                        }
                        for (Scalar f = 0; f < n / m; ++f) {
                            OrbitContainer psi(l + 1, 0);
                            psi[0] = 1;
                            Scalar ord_psi = 0;
                            for (Scalar i = 1; i < l + 1; ++i) {
                                psi[i] = (psi[i - 1] + beta.orbitOf(f)[alpha.orbitOf(i - 1)[1]] * m) % n;
                                if (ord_psi == 0 && psi[i] == 1)
                                {
                                    ord_psi = i;
                                }
                            }

                            bool e = l == ord_psi;
                            if (!e) {
                                continue;
                            }

                            // OrbitContainer psi;
                            // psi.push_back(1);
                            // for (Scalar i = 1; i < n; ++i) {
                            //     Scalar t = (1 + (beta.orbitOf((psi.back() - 1) / m)[alpha.orbit1(1)] + f) * m) % n;
                            //     if (t == 1) {
                            //         break;
                            //     }
                            //     psi.push_back(t);
                            // }
                            Function phi_function(n, 0);
                            for (std::size_t i = 1; i <= n; ++i) {
                                phi_function[i % n] = (phi_function[i - 1] + psi[alpha.orbit1(i - 1)]) % n;
                            }

                            bool a = true;
                            for (Scalar x = 0; x < psi.size(); ++x) {
                                if ((psi[x] % alpha.r()) != alpha.pi(x)) {
                                    a = false;
                                    break;
                                }
                            }

                            bool b = true;
                            const auto &alpha_m = powerOf(alpha, *quotient(alpha).preservingSubgroups.find(m), numberCache[l].skewMorphisms);
                            if (alpha_m.preservingSubgroups.find(beta.r()) != alpha_m.preservingSubgroups.end()) {
                                const auto &powerQuotient = mod(alpha_m, beta.r());;
                                b = &quotient(beta) == &powerQuotient;
                            } else {
                                b = false;
                            }

                            bool c = true;
                            for (Scalar x = 0; x < n / m; ++x) {
                                if (phi_function[x * m] != beta.orbitOf(x)[1] * m) {
                                    c = false;
                                    break;
                                }
                            }

                            bool d = true;
                            for (Scalar x = 0; x < n; ++x) {
                                Scalar phi_m = m;
                                for (Scalar i = 0; i < alpha.orbit1(x); ++i) {
                                    phi_m = phi_function[phi_m];
                                }
                                if (phi_function[(x + m) % n] != (phi_function[x] + phi_m) % n) {
                                    d = false;
                                    break;
                                }
                            }

                            // bool d = true;
                            // Scalar old = 1;
                            // for (std::size_t index = 1; index < psi.size(); ++index) {
                            //     const auto &x = psi[index];
                            //     if (phi_function[old] != x) {
                            //         d = false;
                            //         break;
                            //     }
                            //     old = x;
                            // }
                            // d = d && (phi_function[old] == 1);

                            // bool e = l == psi.size();

                            if (!a) {
                                continue;
                            }
                            if (!b) {
                                if (c && d && e) {
                                    std::cout << "B" << std::endl;
                                }
                                continue;
                            }
                            if (!c) {
                                continue;
                            }
                            if (!d) {
                                continue;
                            }

                            if (!isPermutation(phi_function)) {
                                std::cout << "perm" << std::endl;
                                std::cout << "n: " << n << ", m: " << m << ", l: " << l << ", f: " << f << std::endl;
                                std::cout << "alpha: " << to_json(alpha) << std::endl;
                                std::cout << "beta: " << to_json(beta) << std::endl;
                                std::cout << "phi_function: " << to_json(phi_function) << std::endl;
                                throw "";
                            }
                            const auto orbit1 = computeOrbit1(phi_function);
                            const auto permutation = computePermutation(orbit1, phi_function);

                            if (!isSkewMorphism({orbit1, OrbitContainer{alpha.orbit1()}}, permutation.orbits, phi_function)) {
                                std::cout << "skew" << std::endl;
                                std::cout << "n: " << n << ", m: " << m << ", l: " << l << ", f: " << f << std::endl;
                                std::cout << "alpha: " << to_json(alpha) << std::endl;
                                std::cout << "beta: " << to_json(beta) << std::endl;
                                std::cout << "phi_function: " << to_json(phi_function) << std::endl;
                                throw "";
                            }
                            const auto index = number.skewMorphisms.skewIndexMap[{orbit1, OrbitContainer{alpha.orbit1()}}];
                            if (found[index]) {
                                std::cout << "single" << std::endl;
                                std::cout << "n: " << n << ", m: " << m << ", l: " << l << ", f: " << f << std::endl;
                                std::cout << "alpha: " << to_json(alpha) << std::endl;
                                std::cout << "beta: " << to_json(beta) << std::endl;
                                std::cout << "phi_function: " << to_json(phi_function) << std::endl;
                                throw "";
                            }
                            found[index] = true;
                        }
                    }
                }
            }
        }
        if (n != 1) {
            if (!std::ranges::all_of(found, std::identity())) {
                std::cout << "all" << std::endl;
                std::cout << "n: " << n << std::endl;
                throw "";
            }
        }
        
        std::ofstream file;
        file = std::ofstream{"../skew/" + std::to_string(number.n) + "_auto.txt"};
        for (std::size_t i = 0; i < 2; ++i) {
            for (const auto &c: number.skewMorphisms.classes[i]) {
                const auto &skew = *number.skewMorphisms.skews[c.representant];
                file << to_short_string(skew) << std::endl;
            }
        }
        file = std::ofstream{"../skew/" + std::to_string(number.n) + "_coset.txt"};
        for (const auto &c: number.skewMorphisms.classes[2]) {
            const auto &skew = *number.skewMorphisms.skews[c.representant];
            file << to_short_string(skew) << std::endl;
        }
        file = std::ofstream{"../skew/" + std::to_string(number.n) + "_other.txt"};
        for (std::size_t i = 3; i < number.skewMorphisms.classes.size(); ++i) {
            for (const auto &c: number.skewMorphisms.classes[i]) {
                const auto &skew = *number.skewMorphisms.skews[c.representant];
                file << to_short_string(skew) << std::endl;
            }
        }
    }

//    std::ofstream output("skew" + std::to_string(N) + "output.txt");
//    for (auto n = A; n <= N; ++n) {
//        auto &number = numberCache[n];
//        output << "\nn = " << n << std::endl;
//
//        std::vector<SkewMorphism> skews;
//        std::transform(number.skewMorphisms.skews.begin(), number.skewMorphisms.skews.end(), std::back_inserter(skews), [](const auto &skew) {return *skew;});
//
//        std::sort(skews.begin(), skews.end(), MatoComparator{});
//
//        for (const auto &skew: skews) {
//            printMato(output, skew, true);
//        }
//
//        output << "\nTotal number of skew morphisms of C_" << n << " is " << skews.size();
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
               "        <h3>n = " << n << " = " << factorizationHtml(n) << "</h3>"
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

    for (const auto &complexity: complexities) {
        printf("%d,%d\n", complexity.first, complexity.second);
    }

    return 0;
}
