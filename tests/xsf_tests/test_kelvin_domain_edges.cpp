#include "../testing_utils.h"
#include <xsf/kelvin.h>

TEST_CASE("ker subnormal-input", "[ker][xsf_tests]") {
    using test_case = std::tuple<double, double>;
    // klvna forms log(x/2). For subnormal x the halving is inexact, and at the
    // smallest subnormal it underflows to zero, so the logarithm is -inf and the
    // ker/kei terms that multiply it by an underflowed ber/bei give inf*0 = NaN.
    //
    // Odd multiples of 2^-1074 round onto the neighbouring even one, which is why
    // main returns the same value for 3, 4 and 5 units in the last place.
    //
    // Reference values computed with mpmath with 60 digits of precision:
    //
    // import mpmath as mp
    // mp.mp.dps = 60
    // u = 4.9406564584124654e-324
    // for k in [1, 2, 3, 5, 7, 9, 11, 13, 17, 25]:
    //     print(k * u, mp.nstr(mp.ker(0, mp.mpf(k * u)), 18))
    auto [x, ref] = GENERATE(
        test_case{4.9406564584124654e-324, 744.556003437039675},
        test_case{9.8813129168249309e-324, 743.862856256479729},
        test_case{1.4821969375237396e-323, 743.457391148371565},
        test_case{2.4703282292062327e-323, 742.946565524605574},
        test_case{3.4584595208887258e-323, 742.610093287984361},
        test_case{4.4465908125712189e-323, 742.358778859703455}, test_case{5.434722104253712e-323, 742.158108164241304},
        test_case{6.4228533959362051e-323, 741.991054079578138},
        test_case{8.3991159793011913e-323, 741.722790092983459},
        test_case{1.2351641146031164e-322, 741.337127612171474},
        test_case{9.9998886718268301e-321, 736.943172406632319},
        test_case{9.9999999999999694e-311, 713.917310343812578}, test_case{2.2250738585072014e-308, 708.512350047922519}
    );
    const double w = xsf::ker(x);
    const double rel_error = xsf::extended_relative_error(w, ref);

    CAPTURE(x, w, ref, rel_error);
    REQUIRE(rel_error <= 5e-16);
}

TEST_CASE("kei subnormal-input", "[kei][xsf_tests]") {
    // kei tends to -pi/4 as x tends to zero, and klvna already returns exactly that
    // at x == 0. The same log(x/2) underflow makes it NaN at the smallest subnormal.
    const double x = GENERATE(
        4.9406564584124654e-324, 9.8813129168249309e-324, 1.4821969375237396e-323, 2.4703282292062327e-323,
        3.4584595208887258e-323, 4.4465908125712189e-323, 9.9998886718268301e-321, 9.9999999999999694e-311,
        2.2250738585072014e-308
    );
    const double ref = -0.78539816339744831; // -pi/4, to 60 digits from mpmath
    const double w = xsf::kei(x);
    const double rel_error = xsf::extended_relative_error(w, ref);

    CAPTURE(x, w, ref, rel_error);
    REQUIRE(rel_error <= 5e-16);
}

TEST_CASE("kelvin subnormal-input finite", "[ker][kei][kerp][keip][xsf_tests]") {
    // kerp overflows in this range and -inf is the right answer; what must not
    // happen is NaN. Same for keip, which underflows to zero.
    const double x = GENERATE(
        4.9406564584124654e-324, 9.8813129168249309e-324, 1.4821969375237396e-323, 2.4703282292062327e-323,
        3.4584595208887258e-323
    );
    const double her = xsf::kerp(x);
    const double hei = xsf::keip(x);

    CAPTURE(x, her, hei);
    REQUIRE_FALSE(std::isnan(her));
    REQUIRE_FALSE(std::isnan(hei));
}

TEST_CASE("kelvin infinite-input", "[ker][kei][kerp][keip][xsf_tests]") {
    // ker and kei decay to zero at infinity exactly as k0 and k1 do, and
    // xsf::cephes::k0(inf) already returns 0.0. ber and bei oscillate with
    // unbounded amplitude, so they have no limit and stay NaN.
    const double inf = std::numeric_limits<double>::infinity();

    REQUIRE(xsf::ker(inf) == 0.0);
    REQUIRE(xsf::kei(inf) == 0.0);
    REQUIRE(xsf::kerp(inf) == 0.0);
    REQUIRE(xsf::keip(inf) == 0.0);

    REQUIRE(std::isnan(xsf::ber(inf)));
    REQUIRE(std::isnan(xsf::bei(inf)));
    REQUIRE(std::isnan(xsf::berp(inf)));
    REQUIRE(std::isnan(xsf::beip(inf)));
}

TEST_CASE("kelvin zero-input unchanged", "[ber][bei][ker][kei][xsf_tests]") {
    // The explicit x == 0 branch is untouched by the above; pinned so that it stays
    // that way.
    REQUIRE(xsf::ber(0.0) == 1.0);
    REQUIRE(xsf::bei(0.0) == 0.0);
    REQUIRE(std::isinf(xsf::ker(0.0)));
    REQUIRE(xsf::kei(0.0) == -0.25 * 3.141592653589793);
    REQUIRE(xsf::berp(0.0) == 0.0);
    REQUIRE(xsf::beip(0.0) == 0.0);
    REQUIRE(std::isinf(xsf::kerp(0.0)));
    REQUIRE(xsf::keip(0.0) == 0.0);
}
