import csv
import math


YAZAR_CSV = "results/followup-mature-fulldepth-yazar-cell-drift.csv"
GAO_CSV = "results/followup-mature-fulldepth-gao-cell-drift.csv"
MOUSE_DONORS = 8
YAZAR_TOTAL_MATURE_DONORS = 932
NORMAL_MDE80_MULTIPLIER = 2.8


def read_scatter_anchor(filename):
    with open(filename, newline="", encoding="utf-8") as handle:
        matches = [
            row
            for row in csv.DictReader(handle)
            if row["metric"] == "scatter"
            and row["depth_mode"] == "full_depth"
            and row["model_spec"] == "stratified_global_primary"
            and row["group"] == "Stratified_CellType_Global"
            and row["is_pooled"] == "True"
        ]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one scatter anchor in {filename}, found {len(matches)}")
    return matches[0]


def main():
    yazar = read_scatter_anchor(YAZAR_CSV)
    gao = read_scatter_anchor(GAO_CSV)
    yazar_se = float(yazar["se_sd"])

    rounded_se_mouse_932 = 0.0049 * math.sqrt(YAZAR_TOTAL_MATURE_DONORS / MOUSE_DONORS)
    exact_se_mouse_932 = yazar_se * math.sqrt(YAZAR_TOTAL_MATURE_DONORS / MOUSE_DONORS)
    exact_se_mouse_925 = yazar_se * math.sqrt(int(yazar["n_donors"]) / MOUSE_DONORS)

    rounded_mde_932 = NORMAL_MDE80_MULTIPLIER * rounded_se_mouse_932
    headline_mde = NORMAL_MDE80_MULTIPLIER * exact_se_mouse_932
    exact_mde_925 = NORMAL_MDE80_MULTIPLIER * exact_se_mouse_925

    print(
        "Yazar anchor: "
        f"slope={yazar['slope_sd_per_decade']} SD/decade, "
        f"SE={yazar['se_sd']} SD/decade, p={yazar['pvalue']}, "
        f"n_donors={yazar['n_donors']}"
    )
    print(
        "Gao anchor: "
        f"slope={gao['slope_sd_per_decade']} SD/decade, "
        f"p={gao['pvalue']}, n_donors={gao['n_donors']}"
    )
    print("Mouse design: 4 donors at 18 months + 4 donors at 30 months = 8 donors")
    print(
        "Rounded SE, n=932: "
        f"SE_mouse = 0.0049 * sqrt(932/8) = {rounded_se_mouse_932:.10f}; "
        f"MDE80 = 2.8 * SE_mouse = {rounded_mde_932:.10f}"
    )
    print(
        "Exact SE, n=932 (HEADLINE): "
        f"SE_mouse = {yazar_se:.16f} * sqrt(932/8) = {exact_se_mouse_932:.10f}; "
        f"MDE80 = 2.8 * SE_mouse = {headline_mde:.10f}"
    )
    print(
        "Exact SE, n=925: "
        f"SE_mouse = {yazar_se:.16f} * sqrt(925/8) = {exact_se_mouse_925:.10f}; "
        f"MDE80 = 2.8 * SE_mouse = {exact_mde_925:.10f}"
    )

    gao_ratio = abs(float(gao["slope_sd_per_decade"])) / rounded_mde_932
    yazar_ratio = abs(float(yazar["slope_sd_per_decade"])) / rounded_mde_932
    print(
        f"Gao / audit MDE80 = {abs(float(gao['slope_sd_per_decade'])):.10f} "
        f"/ {rounded_mde_932:.10f} = {gao_ratio:.10f} ({gao_ratio:.1%})"
    )
    print(
        f"Yazar / audit MDE80 = {abs(float(yazar['slope_sd_per_decade'])):.10f} "
        f"/ {rounded_mde_932:.10f} = {yazar_ratio:.10f} ({yazar_ratio:.1%})"
    )

    t_multiplier = 2.4469 + 0.9057
    t_sensitivity_mde = t_multiplier * exact_se_mouse_932
    print(
        "SENSITIVITY ONLY, NOT THE TEX VALUE: "
        "t(0.975, df=6) + t(0.80, df=6) = 2.4469 + 0.9057 = "
        f"{t_multiplier:.4f}; MDE80 = {t_sensitivity_mde:.10f} "
        f"(approximately {t_sensitivity_mde:.4f})"
    )

    assert round(headline_mde, 3) == 0.148, (
        f"Headline MDE80 does not round to 0.148: {headline_mde:.10f}"
    )
    assert round(rounded_mde_932, 3) == 0.148
    assert round(exact_mde_925, 3) == 0.148
    print(f"Filed mouse scatter detection floor: {headline_mde:.3f} SD/decade")


if __name__ == "__main__":
    main()
