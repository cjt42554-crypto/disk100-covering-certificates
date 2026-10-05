#!/usr/bin/env python3
"""Fresh finite Fraction guards for the new uniform actual-tau papers.

Written but NOT EXECUTED by the author. The root is the only intended
executor, once independent static clearance has been recorded.
This file imports no mathematical source, certificate, or old checker.
It evaluates a fixed finite list of exact rational expressions, not a
parameter grid, interval optimizer, trig function, or covering problem.
Polynomial-to-trigonometric remainder signs remain paper premises.
"""

from fractions import Fraction
import json
from pathlib import Path


Q = Fraction


def cosine_upper4(x):
    return 1 - x**2 / 2 + x**4 / 24


def cosine_lower6(x):
    return cosine_upper4(x) - x**6 / 720


def sine_upper5(x):
    return x - x**3 / 6 + x**5 / 120


def sine_lower3(x):
    return x - x**3 / 6


def fraction_pair(value):
    return {"numerator": value.numerator, "denominator": value.denominator}


def main():
    records = []

    def strict(name, margin, paper):
        records.append({"name": name, "relation": "margin > 0",
                        "margin": fraction_pair(margin),
                        "pass": margin > 0, "paper": paper})

    def equal(name, actual, expected, paper):
        records.append({"name": name, "relation": "actual == expected",
                        "actual": fraction_pair(actual),
                        "expected": fraction_pair(expected),
                        "difference": fraction_pair(actual - expected),
                        "pass": actual == expected, "paper": paper})

    r0 = Q(5661, 50000)
    rmax = Q(29, 250)
    d_pro = Q(9, 400)
    d_local = Q(223, 10000)

    # Mean-envelope Section 4: the .3018 arc-sum interface through .116.
    p = "q2_tau_mean_arc_envelope/PROOF.md Section 4, first corollary"
    beta = Q(1257, 2000) + Q(100, 99) * (rmax - r0)
    equal("arc3018_beta_endpoint", beta, Q(124999, 198000), p)
    strict("arc3018_beta_below_06314", Q(".6314") - beta, p)
    x = Q(".6314")
    equal("arc3018_cube_identity", x**3, Q(".251717687144"), p)
    strict("arc3018_cube_lower", x**3 - Q(".2517"), p)
    strict("arc3018_cube_upper", Q(".25172") - x**3, p)
    strict("arc3018_square_upper", Q(".399") - x**2, p)
    equal("arc3018_fifth_product_identity", Q(".25172") * Q(".399"),
          Q(".10043628"), p)
    strict("arc3018_fifth_upper", Q(".10044") - x**5, p)
    equal("arc3018_coarse_sine_identity",
          x - Q(".2517") / 6 + Q(".10044") / 120,
          Q(".590287"), p)
    strict("arc3018_P5_below_sine_upper", Q(".5903") - sine_upper5(x), p)
    base_beta = Q(".6285")
    strict("arc_common_fourth_power_upper", Q(".1561") - base_beta**4, p)
    strict("arc_common_P4_below_cosine_upper",
           Q(".809") - cosine_upper4(base_beta), p)
    denominator = 1 - rmax * Q(".809")
    equal("arc3018_tangent_denominator", denominator, Q(".906156"), p)
    equal("arc3018_tangent_numerator", rmax * Q(".5903"), Q(".0684748"), p)
    margin = denominator * Q(".07559") - rmax * Q(".5903")
    equal("arc3018_tangent_margin_identity", margin, Q(".00002153204"), p)
    strict("arc3018_tangent_margin", margin, p)
    tangent_lower = Q(".07545") + Q(".07545")**3 / 3
    equal("arc3018_tangent_lower_identity", tangent_lower,
          Q(".075593171467875"), p)
    strict("arc3018_tangent_lower_gap", tangent_lower - Q(".07559"), p)
    equal("arc3018_sum_identity", 4 * Q(".07545"), Q(".3018"), p)

    # Mean-envelope Section 4: retain total .294 through .11364.
    p = "q2_tau_mean_arc_envelope/PROOF.md Section 4, second corollary"
    r_arc = Q(".11364")
    beta_arc = base_beta + Q(100, 99) * (r_arc - r0)
    strict("arc0294_beta_below_062893", Q(".62893") - beta_arc, p)
    x = Q(".62893")
    equal("arc0294_cube_identity", x**3, Q(".248775113635957"), p)
    strict("arc0294_cube_lower", x**3 - Q(".24874"), p)
    strict("arc0294_cube_upper", Q(".24878") - x**3, p)
    strict("arc0294_square_upper", Q(".39556") - x**2, p)
    equal("arc0294_fifth_product_identity", Q(".24878") * Q(".39556"),
          Q(".0984074168"), p)
    strict("arc0294_fifth_upper", Q(".0986") - x**5, p)
    equal("arc0294_coarse_sine_identity",
          x - Q(".24874") / 6 + Q(".0986") / 120,
          Q(".588295"), p)
    strict("arc0294_P5_below_sine_upper", Q(".5883") - sine_upper5(x), p)
    tangent_lower = Q(".0735") + Q(".0735")**3 / 3
    equal("arc0294_tangent_lower_identity", tangent_lower,
          Q(".073632355125"), p)
    margin = tangent_lower - r_arc * (Q(".5883") + Q(".809") * tangent_lower)
    equal("arc0294_tangent_margin_identity", margin,
          Q(".000008570228348355"), p)
    strict("arc0294_tangent_margin", margin, p)
    equal("arc0294_sum_identity", 4 * Q(".0735"), Q(".294"), p)

    # New radius transfer arithmetic: the eight old slabs are premises,
    # and are neither read nor recalculated by this guard.
    p = "q2_radius_uniform_tau/PROOF.md Sections 4-5"
    equal("radius_interval_width", rmax - r0, Q(139, 50000), p)
    equal("radius_lipschitz_payment", Q(23, 20) * (rmax - r0),
          Q(3197, 1000000), p)
    residual = Q(1, 4000) - Q(23, 20) * (rmax - r0) + Q(1, 125) * Q(13, 30)
    equal("radius_slab_transfer_residual_identity", residual,
          Q(1559, 3000000), p)
    strict("radius_slab_transfer_residual", residual, p)
    strict("radius_signed_tangent_direction",
           Q(163, 1000) - Q(27, 20) * rmax, p)
    strict("radius_uniform_Hessian_rational_bound",
           Q(7, 4) - Q(14, 970) - Q(175, 102), p)
    coefficient = (Q(14, 1000) + rmax * Q(104, 25)) * Q(225, 169)
    equal("radius_remainder_coefficient_identity", coefficient,
          Q(111726, 169000), p)
    strict("radius_remainder_below_two_thirds", Q(2, 3) - coefficient, p)
    tau_cost = Q(213, 1000) + Q(2, 3) * d_pro
    equal("radius_tau_cost_identity", tau_cost, Q(57, 250), p)
    strict("radius_tau_cost_below_023", Q(23, 100) - tau_cost, p)

    # New local radius endpoint guards and its final squared-norm gap.
    p = "q2_local_radius_uniform_tau/RADIUS_UNIFORM_ACTUAL_TAU_PROOF.md Sections 4-7"
    local_rmax = Q(".1135")
    equal("local_peak_product_identity", local_rmax * Q(".948662"),
          Q(".107673137"), p)
    equal("local_peak_sine_lower_identity", sine_lower3(Q(".108")),
          Q(".107790048"), p)
    strict("local_peak_margin", sine_lower3(Q(".108"))
           - local_rmax * Q(".948662"), p)
    tip_cost = 2 * local_rmax * Q(".989")
    equal("local_tip_cost_identity", tip_cost, Q(".224503"), p)
    strict("local_tip_cost_below_0225", Q(".225") - tip_cost, p)
    numerator = Q(".114") * Q("8.03") * Q("1.221") * d_local
    equal("local_path_perturbation_numerator_identity", numerator,
          Q(".024925330386"), p)
    strict("local_path_perturbation_margin", Q(".061") * Q(".4096")
           - numerator, p)
    strict("local_fullpath_cost_reserve", Q(".30") - Q(".225") - Q(".061"), p)
    upper_norm = 1 + r0**2 - Q("1.62") * r0
    equal("local_mean_upper_identity", upper_norm, Q(".8294023684"), p)
    equal("local_sine_cubic_identity", sine_lower3(Q(".0735")),
          Q(".0734338224375"), p)
    strict("local_sine_cubic_above_007343", sine_lower3(Q(".0735"))
           - Q(".07343"), p)
    equal("local_root_radicand_upper_identity", local_rmax**2 - Q(".07343")**2,
          Q(".0074902851"), p)
    root_square = Q(".086548875")**2
    equal("local_root_upper_square_identity", root_square,
          Q(".007490707763765625"), p)
    strict("local_root_square_margin", root_square
           - (local_rmax**2 - Q(".07343")**2), p)
    equal("local_norm_floor_identity", Q(".997298875") - Q(".086548875"),
          Q(".91075"), p)
    lower_norm = Q(".91075")**2
    equal("local_norm_floor_square_identity", lower_norm, Q(".8294655625"), p)
    gap = lower_norm - upper_norm
    equal("local_final_gap_identity", gap, Q(".0000631941"), p)
    strict("local_final_gap", gap, p)

    # New relaxed-branch margin; the actual unit-anchor equations,
    # angle/lens domain and choice of normal are paper obligations.
    p = "q2_tau_geometry_thresholds/ANALYSIS.md Section 6"
    cosine_margin = cosine_lower6(Q(73, 50)) - Q(11, 100)
    equal("branch_cosine_P6_margin_identity", cosine_margin,
          Q(783848711, 11250000000000), p)
    strict("branch_cosine_P6_margin", cosine_margin, p)
    beta_n = Q(6, 7) - Q(27, 1040)
    equal("branch_normal_lower_identity", beta_n, Q(6051, 7280), p)
    k_n = Q(11, 100) * beta_n - Q(27, 1040)
    equal("branch_dot_lower_identity", k_n, Q(47661, 728000), p)
    strict("branch_u_domain", 1 - rmax**2 - Q(".993")**2, p)
    cross = 2 * Q(".993") * k_n
    equal("branch_cross_lower_identity", cross, Q(47327373, 364000000), p)
    strict("branch_cross_above_013", cross - Q(13, 100), p)
    extra = r0 * (Q(".13") - rmax)
    equal("branch_uniform_extra_margin_identity", extra, Q(".00158508"), p)
    strict("branch_extra_above_relaxed_cap", extra - Q(3, 2000), p)

    # Parameter transport constants: central geometry is an external
    # hypothesis. No support peak, global deficit, or changed target is
    # computed. q^2/8 is the paper completed-square residual formula.
    p = "q2_tau_mean_arc_envelope/PROOF.md Section 5 and q2_parameter_uniform_tau/PRO_PARAMETER_COROLLARY.md"
    a, b = Q(".0003"), Q(".0002")
    gross_loss = Q("4.352") * a + Q("4.63") * b
    equal("pro_parameter_gross_loss_identity", gross_loss, Q(".0022316"), p)
    bootstrap = Q(".025") - gross_loss
    equal("pro_parameter_bootstrap_identity", bootstrap, Q(".0227684"), p)
    strict("pro_parameter_bootstrap_above_central_fee", bootstrap - d_pro, p)
    q = 4 * a + b
    strict("pro_parameter_monotonic_transport", 1 - b - q / Q(".3"), p)
    fee = d_pro - Q(".6") * a - Q(".1725") * b
    equal("pro_parameter_wider_fee_identity", fee, Q(44571, 2000000), p)
    balanced_fee = d_pro - Q(".6") * b - Q(".1725") * b
    equal("pro_parameter_balanced_fee_identity", balanced_fee,
          Q(44691, 2000000), p)
    equal("pro_parameter_peak_displacement_identity", q / 4, Q(7, 20000), p)
    equal("pro_parameter_support_residual_identity", q**2 / 8,
          Q(49, 200000000), p)
    strict("pro_parameter_residual_positive", q**2 / 8, p)

    p = "q2_parameter_uniform_tau/PARAMETER_UNIFORM_LEMMA.md Sections 3-6"
    a, b = Q(".0004"), Q(".0002")
    gross_loss = Q("4.352") * a + Q("4.63") * b
    equal("local_parameter_gross_loss_identity", gross_loss, Q(".0026668"), p)
    strict("local_parameter_bootstrap_above_central_fee",
           Q(".025") - gross_loss - d_local, p)
    q = 4 * a + b
    strict("local_parameter_monotonic_transport", 1 - b - q / Q(".2"), p)
    local_fee = d_local - Q(".6") * a - Q(".1723") * b
    equal("local_parameter_broad_fee_identity", local_fee,
          Q(1101277, 50000000), p)
    narrow_fee = d_local - Q(".6") * Q(".00002") - Q(".1723") * Q(".00002")
    equal("local_parameter_narrow_fee_identity", narrow_fee,
          Q(11142277, 500000000), p)
    equal("local_parameter_peak_displacement_identity", q / 4, Q(9, 20000), p)
    equal("local_parameter_support_residual_identity", q**2 / 8,
          Q(81, 200000000), p)
    strict("local_parameter_residual_positive", q**2 / 8, p)

    all_pass = all(row["pass"] for row in records)
    result = {"status": "PASS" if all_pass else "FAIL",
              "scope": "Fresh finite rational constant guards only",
              "old_48_guard_replays": 0,
              "mathematical_source_imports": 0,
              "numerical_parameter_scans": 0,
              "global_covering_certified": False,
              "guards": records}
    output = Path(__file__).resolve().parent / "verification.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    for row in records:
        print(("PASS " if row["pass"] else "FAIL ") + row["name"])
    print("OVERALL " + result["status"] + " fresh finite constant guards")
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
