# SPDX-License-Identifier: MIT
"""Screen annual Bornholm CfD terms against fixed outside options."""

import argparse
from dataclasses import dataclass
from math import isfinite
from pathlib import Path

import numpy as np
import pandas as pd

REFERENCE_STRIKE = 75.0
REFERENCE_DE_SHARE = 0.7
TOLERANCE = 1e-8  # million EUR/year, for endpoint arithmetic


@dataclass(frozen=True)
class AgreementCase:
    """
    Incremental annual gains before the candidate's CfD, in million EUR.

    The reference already includes its own settlement. Its contract stays
    fixed throughout the candidate's sweep. Generation is in TWh/year.
    """

    generation: float
    capture_price: float
    dk: float
    de: float
    developer: float

    def __post_init__(self):
        values = (self.generation, self.capture_price, self.dk, self.de, self.developer)
        if not all(map(isfinite, values)) or self.generation <= 0:
            raise ValueError("Use finite annual accounts and positive generation.")

    @classmethod
    def from_anchor(cls, generation, capture_price, dk, de, developer):
        """Remove the candidate's €75/MWh, 70% DE settlement from observed gains."""
        support = generation * (REFERENCE_STRIKE - capture_price)
        return cls(
            generation,
            capture_price,
            dk + (1 - REFERENCE_DE_SHARE) * support,
            de + REFERENCE_DE_SHARE * support,
            developer - support,
        )

    def support(self, strike):
        return self.generation * (strike - self.capture_price)

    def strike(self, support):
        return self.capture_price + support / self.generation

    def settlement(self, strike, german_share):
        if not 0 <= german_share <= 1 or not isfinite(strike) or strike < 0:
            raise ValueError("Use a nonnegative strike and a German share in [0, 1].")
        total = self.support(strike)
        de_payment = german_share * total
        dk_payment = total - de_payment
        return dict(
            strike=strike,
            german_share=german_share,
            dk_support=dk_payment,
            de_support=de_payment,
            net_support=total,
            dk_gain=self.dk - dk_payment,
            de_gain=self.de - de_payment,
            developer_gain=self.developer + total,
        )


@dataclass(frozen=True)
class Participation:
    """Annual incremental welfare floors and positive net CfD payment ceilings."""

    dk_floor: float = 0.0
    de_floor: float = 0.0
    developer_floor: float = 0.0
    dk_cap: float = float("inf")
    de_cap: float = float("inf")

    def __post_init__(self):
        if not all(map(isfinite, (self.dk_floor, self.de_floor, self.developer_floor))):
            raise ValueError("Participation floors must be finite.")
        if not (self.dk_cap >= 0 and self.de_cap >= 0):
            raise ValueError("Funding ceilings must be nonnegative.")

    def payment_limits(self, case):
        return min(case.dk - self.dk_floor, self.dk_cap), min(
            case.de - self.de_floor, self.de_cap
        )

    def minimum_support(self, case):
        return max(self.developer_floor - case.developer, case.support(0))


def support_interval(case, requirements, *, fixed_german_share=None):
    """
    Return exact feasible net support bounds or None.

    A common share allocates support and repayments with the same sign.
    If a country needs a credit, feasibility therefore requires net repayments.
    """
    lower = requirements.minimum_support(case)
    dk_limit, de_limit = requirements.payment_limits(case)
    if fixed_german_share is None:
        if min(dk_limit, de_limit) < 0:
            upper = min(0, dk_limit) + min(0, de_limit)
        else:
            upper = dk_limit + de_limit
    else:
        if not 0 <= fixed_german_share <= 1:
            raise ValueError("German funding share must be in [0, 1].")
        upper = float("inf")
        for share, limit in (
            (1 - fixed_german_share, dk_limit),
            (fixed_german_share, de_limit),
        ):
            if share > 0:
                upper = min(upper, limit / share)
            elif limit < 0:
                return None
    return (lower, max(lower, upper)) if lower <= upper + TOLERANCE else None


def german_payment_interval(support, requirements, case):
    dk_limit, de_limit = requirements.payment_limits(case)
    return max(min(0, support), support - dk_limit), min(max(0, support), de_limit)


def balanced_agreement(case, requirements, *, dk_weight=1.0):
    """
    At minimum feasible support, maximise the smaller weighted country gain.

    Germany's weight is one. We compare gains above the participation floors.
    A DK weight of two seeks a 2:1 DK/DE split of those gains when attainable.
    """
    if not isfinite(dk_weight) or dk_weight <= 0:
        raise ValueError("Danish gain weight must be finite and positive.")
    interval = support_interval(case, requirements)
    if interval is None:
        return None
    support = interval[0]
    lower, upper = german_payment_interval(support, requirements, case)
    dk_room = case.dk - requirements.dk_floor
    de_room = case.de - requirements.de_floor
    target_de_payment = (dk_weight * de_room - dk_room + support) / (1 + dk_weight)
    de_payment = min(max(target_de_payment, lower), upper)
    share = min(1.0, max(0.0, de_payment / support)) if support != 0 else 0.5
    return case.settlement(max(0.0, case.strike(support)), share)


def summarise(case, requirements, *, dk_weight=1.0):
    free = support_interval(case, requirements)
    fixed = support_interval(case, requirements, fixed_german_share=REFERENCE_DE_SHARE)
    result = {
        "participation_surplus": (
            case.dk
            + case.de
            + case.developer
            - requirements.dk_floor
            - requirements.de_floor
            - requirements.developer_floor
        )
    }
    for label, interval in (("variable_share", free), ("fixed_70_percent_de", fixed)):
        result[f"{label}_feasible"] = interval is not None
        result[f"{label}_min_strike"] = (
            max(0, case.strike(interval[0])) if interval else None
        )
        result[f"{label}_max_strike"] = case.strike(interval[1]) if interval else None
    agreement = balanced_agreement(case, requirements, dk_weight=dk_weight)
    if agreement is not None:
        result.update(agreement)
    return result


def save_figure(cases, requirements, output, dk_weight, *, provisional):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        len(cases) // 2, 2, figsize=(10, 3.4 * (len(cases) // 2)), squeeze=False
    )
    strike_points = [REFERENCE_STRIKE]
    for _, case in cases:
        interval = support_interval(case, requirements)
        if interval is not None:
            strike_points.extend(case.strike(value) for value in interval)
    padding = max(0.25, (max(strike_points) - min(strike_points)) * 0.1)
    strike_limits = (max(0, min(strike_points) - padding), max(strike_points) + padding)
    for ax, (key, case) in zip(axes.flat, cases, strict=True):
        free = support_interval(case, requirements)
        if free is None:
            ax.text(
                0.5,
                0.5,
                "No feasible strike and funding share",
                transform=ax.transAxes,
                ha="center",
            )
        else:
            for start, end in ((free[0], min(0, free[1])), (max(0, free[0]), free[1])):
                if start > end or (start == end == 0):
                    continue
                supports = np.linspace(start, end, 201)
                lower, upper = np.array(
                    [
                        german_payment_interval(total, requirements, case)
                        for total in supports
                    ]
                ).T
                # A vertical line below represents the zero-payment case.
                divisor = np.where(supports == 0, np.nan, supports)
                lo = np.clip(np.minimum(lower / divisor, upper / divisor), 0, 1)
                hi = np.clip(np.maximum(lower / divisor, upper / divisor), 0, 1)
                strikes = case.strike(supports)
                ax.fill_between(
                    strikes, 100 * lo, 100 * hi, alpha=0.25, color="#007D78"
                )
                ax.plot(strikes, 100 * lo, strikes, 100 * hi, color="#007D78")
            if free[0] <= 0 <= free[1]:
                ax.vlines(case.capture_price, 0, 100, color="#007D78")
            agreement = balanced_agreement(case, requirements, dk_weight=dk_weight)
            ax.scatter(
                agreement["strike"],
                100 * agreement["german_share"],
                color="#293A52",
                zorder=3,
            )
        ax.axhline(70, linestyle="--", color="#777777", linewidth=1)
        ax.set(
            title=f"{key['design_label']} · {key['allocation']}",
            xlabel="Strike [EUR/MWh]",
            ylabel="German CfD share [%]",
            ylim=(-2, 102),
        )
        ax.set_xlim(*strike_limits)
    title = "Bornholm 2035 · Feasible annual CfD terms versus split radial"
    fig.suptitle(title + (" (provisional)" if provisional else ""))
    fig.text(
        0.5,
        0.02,
        "Shading shows feasible terms. The dot marks the selected split. The dashed line marks 70%.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    fig.savefig(output, dpi=180)
    plt.close(fig)


def assess_matched(root, output, requirements, *, dk_weight, allow_provisional=False):
    from scripts.analyse_bornholm_matched import compare

    runs, changes = compare(root)
    provisional = (
        not runs.status.eq("checks_passed").all() or changes.residual.abs().max() > 0.01
    )
    if provisional and not allow_provisional:
        raise ValueError(
            "Review the matched solves or use --allow-provisional for diagnostics."
        )
    run = runs.set_index("design")
    anchor = changes.query(
        "reference == 'central_split_radial' and strike == @REFERENCE_STRIKE"
    )
    cases, rows = [], []
    for row in anchor.itertuples():
        design = run.loc[row.design]
        case = AgreementCase.from_anchor(
            design.generation_twh, design.capture_price, row.DK, row.DE, row.developer
        )
        key = dict(
            design=row.design,
            reference=row.reference,
            design_label="Matched hybrid"
            if row.design == "central_hybrid_matched"
            else "Treaty hybrid",
            allocation="Equal cable costs"
            if row.pooled_cost_share
            else "Landing-country costs",
        )
        cases.append((key, case))
        rows.append(
            dict(
                **key,
                provisional=provisional,
                system_gain=row.system,
                others_gain=row.others,
                generation_twh=case.generation,
                capture_price=case.capture_price,
                dk_before_support=case.dk,
                de_before_support=case.de,
                developer_before_support=case.developer,
                dk_floor=requirements.dk_floor,
                de_floor=requirements.de_floor,
                developer_floor=requirements.developer_floor,
                dk_cap=requirements.dk_cap,
                de_cap=requirements.de_cap,
                dk_weight=dk_weight,
                **summarise(case, requirements, dk_weight=dk_weight),
            )
        )
    output.mkdir(parents=True, exist_ok=True)
    runs.to_csv(output / "matched_runs.csv", index=False)
    changes.to_csv(output / "welfare_changes.csv", index=False)
    pd.DataFrame(rows).to_csv(output / "agreement_screen.csv", index=False)
    save_figure(
        cases,
        requirements,
        output / "strike_and_funding_share.png",
        dk_weight,
        provisional=provisional,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("transfer_dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dk-floor", type=float, default=0.0)
    parser.add_argument("--de-floor", type=float, default=0.0)
    parser.add_argument("--developer-floor", type=float, default=0.0)
    parser.add_argument("--dk-cap", type=float, default=float("inf"))
    parser.add_argument("--de-cap", type=float, default=float("inf"))
    parser.add_argument("--dk-weight", type=float, default=1.0)
    parser.add_argument("--allow-provisional", action="store_true")
    args = parser.parse_args()
    assess_matched(
        args.transfer_dir,
        args.output,
        Participation(
            args.dk_floor, args.de_floor, args.developer_floor, args.dk_cap, args.de_cap
        ),
        dk_weight=args.dk_weight,
        allow_provisional=args.allow_provisional,
    )
