import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

METRICS = ["lddt", "bb_lddt", "tm_score", "rmsd"]


def compute_af3_metrics(preds, evals, name):
    metrics = {}

    top_model = None
    top_confidence = -1000
    for model_id in range(5):
        # Load confidence file
        confidence_file = (
            Path(preds) / f"seed-1_sample-{model_id}" / "summary_confidences.json"
        )
        with confidence_file.open("r") as f:
            confidence_data = json.load(f)
            confidence = confidence_data["ranking_score"]
            if confidence > top_confidence:
                top_model = model_id
                top_confidence = confidence

        # Load eval file
        eval_file = Path(evals) / f"{name}_model_{model_id}.json"
        with eval_file.open("r") as f:
            eval_data = json.load(f)
            for metric_name in METRICS:
                if metric_name in eval_data:
                    metrics.setdefault(metric_name, []).append(eval_data[metric_name])

            if "dockq" in eval_data and eval_data["dockq"] is not None:
                metrics.setdefault("dockq_>0.23", []).append(
                    np.mean(
                        [float(v > 0.23) for v in eval_data["dockq"] if v is not None]
                    )
                )
                metrics.setdefault("dockq_>0.49", []).append(
                    np.mean(
                        [float(v > 0.49) for v in eval_data["dockq"] if v is not None]
                    )
                )
                metrics.setdefault("len_dockq_", []).append(
                    len([v for v in eval_data["dockq"] if v is not None])
                )

        eval_file = Path(evals) / f"{name}_model_{model_id}_ligand.json"
        with eval_file.open("r") as f:
            eval_data = json.load(f)
            if "lddt_pli" in eval_data:
                lddt_plis = [
                    x["score"] for x in eval_data["lddt_pli"]["assigned_scores"]
                ]
                for _ in eval_data["lddt_pli"][
                    "model_ligand_unassigned_reason"
                ].items():
                    lddt_plis.append(0)
                if not lddt_plis:
                    continue
                lddt_pli = np.mean([x for x in lddt_plis])
                metrics.setdefault("lddt_pli", []).append(lddt_pli)
                metrics.setdefault("len_lddt_pli", []).append(len(lddt_plis))

            if "rmsd" in eval_data:
                rmsds = [x["score"] for x in eval_data["rmsd"]["assigned_scores"]]
                for _ in eval_data["rmsd"]["model_ligand_unassigned_reason"].items():
                    rmsds.append(100)
                if not rmsds:
                    continue
                rmsd2 = np.mean([x < 2.0 for x in rmsds])
                rmsd5 = np.mean([x < 5.0 for x in rmsds])
                metrics.setdefault("rmsd<2", []).append(rmsd2)
                metrics.setdefault("rmsd<5", []).append(rmsd5)
                metrics.setdefault("len_rmsd", []).append(len(rmsds))

    # Get oracle
    oracle = {k: min(v) if k == "rmsd" else max(v) for k, v in metrics.items()}
    avg = {k: sum(v) / len(v) for k, v in metrics.items()}
    top1 = {k: v[top_model] for k, v in metrics.items()}

    results = {}
    for metric_name in metrics:
        if metric_name.startswith("len_"):
            continue
        if metric_name == "lddt_pli":
            l = metrics["len_lddt_pli"][0]
        elif metric_name == "rmsd<2" or metric_name == "rmsd<5":
            l = metrics["len_rmsd"][0]
        elif metric_name == "dockq_>0.23" or metric_name == "dockq_>0.49":
            l = metrics["len_dockq_"][0]
        else:
            l = 1
        results[metric_name] = {
            "oracle": oracle[metric_name],
            "average": avg[metric_name],
            "top1": top1[metric_name],
            "len": l,
        }

    return results


def compute_chai_metrics(preds, evals, name):
    metrics = {}

    top_model = None
    top_confidence = 0
    for model_id in range(5):
        # Load confidence file
        confidence_file = Path(preds) / f"scores.model_idx_{model_id}.npz"
        confidence_data = np.load(confidence_file)
        confidence = confidence_data["aggregate_score"].item()
        if confidence > top_confidence:
            top_model = model_id
            top_confidence = confidence

        # Load eval file
        eval_file = Path(evals) / f"{name}_model_{model_id}.json"
        with eval_file.open("r") as f:
            eval_data = json.load(f)
            for metric_name in METRICS:
                if metric_name in eval_data:
                    metrics.setdefault(metric_name, []).append(eval_data[metric_name])

            if "dockq" in eval_data and eval_data["dockq"] is not None:
                metrics.setdefault("dockq_>0.23", []).append(
                    np.mean(
                        [float(v > 0.23) for v in eval_data["dockq"] if v is not None]
                    )
                )
                metrics.setdefault("dockq_>0.49", []).append(
                    np.mean(
                        [float(v > 0.49) for v in eval_data["dockq"] if v is not None]
                    )
                )
                metrics.setdefault("len_dockq_", []).append(
                    len([v for v in eval_data["dockq"] if v is not None])
                )

        eval_file = Path(evals) / f"{name}_model_{model_id}_ligand.json"
        with eval_file.open("r") as f:
            eval_data = json.load(f)
            if "lddt_pli" in eval_data:
                lddt_plis = [
                    x["score"] for x in eval_data["lddt_pli"]["assigned_scores"]
                ]
                for _ in eval_data["lddt_pli"][
                    "model_ligand_unassigned_reason"
                ].items():
                    lddt_plis.append(0)
                if not lddt_plis:
                    continue
                lddt_pli = np.mean([x for x in lddt_plis])
                metrics.setdefault("lddt_pli", []).append(lddt_pli)
                metrics.setdefault("len_lddt_pli", []).append(len(lddt_plis))

            if "rmsd" in eval_data:
                rmsds = [x["score"] for x in eval_data["rmsd"]["assigned_scores"]]
                for _ in eval_data["rmsd"]["model_ligand_unassigned_reason"].items():
                    rmsds.append(100)
                if not rmsds:
                    continue
                rmsd2 = np.mean([x < 2.0 for x in rmsds])
                rmsd5 = np.mean([x < 5.0 for x in rmsds])
                metrics.setdefault("rmsd<2", []).append(rmsd2)
                metrics.setdefault("rmsd<5", []).append(rmsd5)
                metrics.setdefault("len_rmsd", []).append(len(rmsds))

    # Get oracle
    oracle = {k: min(v) if k == "rmsd" else max(v) for k, v in metrics.items()}
    avg = {k: sum(v) / len(v) for k, v in metrics.items()}
    top1 = {k: v[top_model] for k, v in metrics.items()}

    results = {}
    for metric_name in metrics:
        if metric_name.startswith("len_"):
            continue
        if metric_name == "lddt_pli":
            l = metrics["len_lddt_pli"][0]
        elif metric_name == "rmsd<2" or metric_name == "rmsd<5":
            l = metrics["len_rmsd"][0]
        elif metric_name == "dockq_>0.23" or metric_name == "dockq_>0.49":
            l = metrics["len_dockq_"][0]
        else:
            l = 1
        results[metric_name] = {
            "oracle": oracle[metric_name],
            "average": avg[metric_name],
            "top1": top1[metric_name],
            "len": l,
        }

    return results


def compute_boltz_metrics(preds, evals, name):
    metrics = {}

    top_model = None
    top_confidence = 0
    for model_id in range(5):
        # Load confidence file
        confidence_file = (
            Path(preds) / f"confidence_{Path(preds).name}_model_{model_id}.json"
        )
        with confidence_file.open("r") as f:
            confidence_data = json.load(f)
            confidence = confidence_data["confidence_score"]
            if confidence > top_confidence:
                top_model = model_id
                top_confidence = confidence

        # Load eval file
        eval_file = Path(evals) / f"{name}_model_{model_id}.json"
        with eval_file.open("r") as f:
            eval_data = json.load(f)
            for metric_name in METRICS:
                if metric_name in eval_data:
                    metrics.setdefault(metric_name, []).append(eval_data[metric_name])

            if "dockq" in eval_data and eval_data["dockq"] is not None:
                metrics.setdefault("dockq_>0.23", []).append(
                    np.mean(
                        [float(v > 0.23) for v in eval_data["dockq"] if v is not None]
                    )
                )
                metrics.setdefault("dockq_>0.49", []).append(
                    np.mean(
                        [float(v > 0.49) for v in eval_data["dockq"] if v is not None]
                    )
                )
                metrics.setdefault("len_dockq_", []).append(
                    len([v for v in eval_data["dockq"] if v is not None])
                )

        eval_file = Path(evals) / f"{name}_model_{model_id}_ligand.json"
        with eval_file.open("r") as f:
            eval_data = json.load(f)
            if "lddt_pli" in eval_data:
                lddt_plis = [
                    x["score"] for x in eval_data["lddt_pli"]["assigned_scores"]
                ]
                for _ in eval_data["lddt_pli"][
                    "model_ligand_unassigned_reason"
                ].items():
                    lddt_plis.append(0)
                if not lddt_plis:
                    continue
                lddt_pli = np.mean([x for x in lddt_plis])
                metrics.setdefault("lddt_pli", []).append(lddt_pli)
                metrics.setdefault("len_lddt_pli", []).append(len(lddt_plis))

            if "rmsd" in eval_data:
                rmsds = [x["score"] for x in eval_data["rmsd"]["assigned_scores"]]
                for _ in eval_data["rmsd"]["model_ligand_unassigned_reason"].items():
                    rmsds.append(100)
                if not rmsds:
                    continue
                rmsd2 = np.mean([x < 2.0 for x in rmsds])
                rmsd5 = np.mean([x < 5.0 for x in rmsds])
                metrics.setdefault("rmsd<2", []).append(rmsd2)
                metrics.setdefault("rmsd<5", []).append(rmsd5)
                metrics.setdefault("len_rmsd", []).append(len(rmsds))

    # Get oracle
    oracle = {k: min(v) if k == "rmsd" else max(v) for k, v in metrics.items()}
    avg = {k: sum(v) / len(v) for k, v in metrics.items()}
    top1 = {k: v[top_model] for k, v in metrics.items()}

    results = {}
    for metric_name in metrics:
        if metric_name.startswith("len_"):
            continue
        if metric_name == "lddt_pli":
            l = metrics["len_lddt_pli"][0]
        elif metric_name == "rmsd<2" or metric_name == "rmsd<5":
            l = metrics["len_rmsd"][0]
        elif metric_name == "dockq_>0.23" or metric_name == "dockq_>0.49":
            l = metrics["len_dockq_"][0]
        else:
            l = 1
        results[metric_name] = {
            "oracle": oracle[metric_name],
            "average": avg[metric_name],
            "top1": top1[metric_name],
            "len": l,
        }

    return results


def eval_models(
    chai_preds,
    chai_evals,
    af3_preds,
    af3_evals,
    boltz_preds,
    boltz_evals,
    boltz_preds_x,
    boltz_evals_x,
):
    # Load preds and make sure we have predictions for all models
    chai_preds_names = {
        x.name.lower(): x
        for x in Path(chai_preds).iterdir()
        if not x.name.lower().startswith(".")
    }
    af3_preds_names = {
        x.name.lower(): x
        for x in Path(af3_preds).iterdir()
        if not x.name.lower().startswith(".")
    }
    boltz_preds_names = {
        x.name.lower(): x
        for x in Path(boltz_preds).iterdir()
        if not x.name.lower().startswith(".")
    }
    boltz_preds_names_x = {
        x.name.lower(): x
        for x in Path(boltz_preds_x).iterdir()
        if not x.name.lower().startswith(".")
    }

    print("Chai preds", len(chai_preds_names))
    print("Af3 preds", len(af3_preds_names))
    print("Boltz preds", len(boltz_preds_names))
    print("Boltzx preds", len(boltz_preds_names_x))

    common = (
        set(chai_preds_names.keys())
        & set(af3_preds_names.keys())
        & set(boltz_preds_names.keys())
        & set(boltz_preds_names_x.keys())
    )

    # Remove examples in the validation set
    keys_to_remove = ["t1133", "h1134", "r1134s1", "t1134s2", "t1121", "t1123", "t1159"]
    for key in keys_to_remove:
        if key in common:
            common.remove(key)
    print("Common", len(common))

    # Create a dataframe with the following schema:
    # tool, name, metric, oracle, average, top1
    results = []
    for name in tqdm(common):
        try:
            af3_results = compute_af3_metrics(
                af3_preds_names[name],
                af3_evals,
                name,
            )

        except Exception as e:
            import traceback

            traceback.print_exc()
            print(f"Error evaluating AF3 {name}: {e}")
            continue
        try:
            chai_results = compute_chai_metrics(
                chai_preds_names[name],
                chai_evals,
                name,
            )
        except Exception as e:
            import traceback

            traceback.print_exc()
            print(f"Error evaluating Chai {name}: {e}")
            continue
        try:
            boltz_results = compute_boltz_metrics(
                boltz_preds_names[name],
                boltz_evals,
                name,
            )
        except Exception as e:
            import traceback

            traceback.print_exc()
            print(f"Error evaluating Boltz {name}: {e}")
            continue

        try:
            boltz_results_x = compute_boltz_metrics(
                boltz_preds_names_x[name],
                boltz_evals_x,
                name,
            )
        except Exception as e:
            import traceback

            traceback.print_exc()
            print(f"Error evaluating Boltzx {name}: {e}")
            continue

        for metric_name in af3_results:
            if metric_name in chai_results and metric_name in boltz_results:
                if (
                    (
                        af3_results[metric_name]["len"]
                        == chai_results[metric_name]["len"]
                    )
                    and (
                        af3_results[metric_name]["len"]
                        == boltz_results[metric_name]["len"]
                    )
                    and (
                        af3_results[metric_name]["len"]
                        == boltz_results_x[metric_name]["len"]
                    )
                ):
                    results.append(
                        {
                            "tool": "AF3 oracle",
                            "target": name,
                            "metric": metric_name,
                            "value": af3_results[metric_name]["oracle"],
                        }
                    )
                    results.append(
                        {
                            "tool": "AF3 top-1",
                            "target": name,
                            "metric": metric_name,
                            "value": af3_results[metric_name]["top1"],
                        }
                    )
                    results.append(
                        {
                            "tool": "Chai-1 oracle",
                            "target": name,
                            "metric": metric_name,
                            "value": chai_results[metric_name]["oracle"],
                        }
                    )
                    results.append(
                        {
                            "tool": "Chai-1 top-1",
                            "target": name,
                            "metric": metric_name,
                            "value": chai_results[metric_name]["top1"],
                        }
                    )
                    results.append(
                        {
                            "tool": "Boltz-1 oracle",
                            "target": name,
                            "metric": metric_name,
                            "value": boltz_results[metric_name]["oracle"],
                        }
                    )
                    results.append(
                        {
                            "tool": "Boltz-1 top-1",
                            "target": name,
                            "metric": metric_name,
                            "value": boltz_results[metric_name]["top1"],
                        }
                    )
                    results.append(
                        {
                            "tool": "Boltz-1x oracle",
                            "target": name,
                            "metric": metric_name,
                            "value": boltz_results_x[metric_name]["oracle"],
                        }
                    )
                    results.append(
                        {
                            "tool": "Boltz-1x top-1",
                            "target": name,
                            "metric": metric_name,
                            "value": boltz_results_x[metric_name]["top1"],
                        }
                    )
                else:
                    print(
                        "Different lengths",
                        name,
                        metric_name,
                        af3_results[metric_name]["len"],
                        chai_results[metric_name]["len"],
                        boltz_results[metric_name]["len"],
                        boltz_results_x[metric_name]["len"],
                    )
            else:
                print(
                    "Missing metric",
                    name,
                    metric_name,
                    metric_name in chai_results,
                    metric_name in boltz_results,
                    metric_name in boltz_results_x,
                )

    # Write the results to a file, ensure we only keep the target & metrics where we have all tools
    df = pd.DataFrame(results)
    return df


def eval_validity_checks(df):
    # Filter the dataframe to only include the targets in the validity checks
    name_mapping = {
        "af3": "AF3 top-1",
        "chai": "Chai-1 top-1",
        "boltz1": "Boltz-1 top-1",
        "boltz1x": "Boltz-1x top-1",
    }
    top1 = df[df["model_idx"] == 0]
    top1 = top1[["tool", "pdb_id", "valid"]]
    top1["tool"] = top1["tool"].apply(lambda x: name_mapping[x])
    top1 = top1.rename(columns={"tool": "tool", "pdb_id": "target", "valid": "value"})
    top1["metric"] = "physical validity"
    top1["target"] = top1["target"].apply(lambda x: x.lower())
    top1 = top1[["tool", "target", "metric", "value"]]

    name_mapping = {
        "af3": "AF3 oracle",
        "chai": "Chai-1 oracle",
        "boltz1": "Boltz-1 oracle",
        "boltz1x": "Boltz-1x oracle",
    }
    oracle = df[["tool", "model_idx", "pdb_id", "valid"]]
    oracle = oracle.groupby(["tool", "pdb_id"])["valid"].max().reset_index()
    oracle = oracle.rename(
        columns={"tool": "tool", "pdb_id": "target", "valid": "value"}
    )
    oracle["tool"] = oracle["tool"].apply(lambda x: name_mapping[x])
    oracle["metric"] = "physical validity"
    oracle = oracle[["tool", "target", "metric", "value"]]
    oracle["target"] = oracle["target"].apply(lambda x: x.lower())
    out = pd.concat([top1, oracle])
    return out


def bootstrap_ci(series, n_boot=1000, alpha=0.05):
    """
    Compute 95% bootstrap confidence intervals for the mean of 'series'.
    """
    n = len(series)
    boot_means = []
    # Perform bootstrap resampling
    for _ in range(n_boot):
        sample = series.sample(n, replace=True)
        boot_means.append(sample.mean())

    boot_means = np.array(boot_means)
    mean_val = np.mean(series)
    lower = np.percentile(boot_means, 100 * alpha / 2)
    upper = np.percentile(boot_means, 100 * (1 - alpha / 2))
    return mean_val, lower, upper


def plot_data(desired_tools, desired_metrics, df, dataset, filename):
    filtered_df = df[
        df["tool"].isin(desired_tools) & df["metric"].isin(desired_metrics)
    ]

    # Apply bootstrap to each (tool, metric) group
    boot_stats = filtered_df.groupby(["tool", "metric"])["value"].apply(bootstrap_ci)

    # boot_stats is a Series of tuples (mean, lower, upper). Convert to DataFrame:
    boot_stats = boot_stats.apply(pd.Series)
    boot_stats.columns = ["mean", "lower", "upper"]

    # Unstack to get a DataFrame suitable for plotting
    plot_data = boot_stats["mean"].unstack("tool")
    plot_data = plot_data.reindex(desired_metrics)

    lower_data = boot_stats["lower"].unstack("tool")
    lower_data = lower_data.reindex(desired_metrics)

    upper_data = boot_stats["upper"].unstack("tool")
    upper_data = upper_data.reindex(desired_metrics)

    # If you need a specific order of tools:
    tool_order = [
        "AF3 oracle",
        "AF3 top-1",
        "Chai-1 oracle",
        "Chai-1 top-1",
        "Boltz-1 oracle",
        "Boltz-1 top-1",
        "Boltz-1x oracle",
        "Boltz-1x top-1",
    ]
    plot_data = plot_data[tool_order]
    lower_data = lower_data[tool_order]
    upper_data = upper_data[tool_order]

    # Rename metrics
    renaming = {
        "lddt_pli": "Mean LDDT-PLI",
        "rmsd<2": "L-RMSD < 2A",
        "lddt": "Mean LDDT",
        "dockq_>0.23": "DockQ > 0.23",
        "physical validity": "Physical Validity",
    }
    plot_data = plot_data.rename(index=renaming)
    lower_data = lower_data.rename(index=renaming)
    upper_data = upper_data.rename(index=renaming)
    mean_vals = plot_data.values

    # Colors
    tool_colors = [
        "#994C00",  # AF3 oracle
        "#FFB55A",  # AF3 top-1
        "#931652",  # Chai-1 oracle
        "#FC8AD9",  # Chai-1 top-1
        "#188F52",  # Boltz-1 oracle
        "#86E935",  # Boltz-1 top-1
        "#004D80",  # Boltz-1x oracle
        "#55C2FF",  # Boltz-1x top-1
    ]

    fig, ax = plt.subplots(figsize=(10, 5))

    x = np.arange(len(plot_data.index))
    bar_spacing = 0.015
    total_width = 0.7
    # Adjust width to account for the spacing
    width = (total_width - (len(tool_order) - 1) * bar_spacing) / len(tool_order)

    for i, tool in enumerate(tool_order):
        # Each subsequent bar moves over by width + bar_spacing
        offsets = x - (total_width - width) / 2 + i * (width + bar_spacing)
        # Extract the means and errors for this tool
        tool_means = plot_data[tool].values
        tool_yerr_lower = mean_vals[:, i] - lower_data.values[:, i]
        tool_yerr_upper = upper_data.values[:, i] - mean_vals[:, i]
        # Construct yerr array specifically for this tool
        tool_yerr = np.vstack([tool_yerr_lower, tool_yerr_upper])

        ax.bar(
            offsets,
            tool_means,
            width=width,
            color=tool_colors[i],
            label=tool,
            yerr=tool_yerr,
            capsize=2,
            error_kw={"elinewidth": 0.75},
        )

    ax.set_xticks(x)
    ax.set_xticklabels(plot_data.index, rotation=0)
    ax.set_ylabel("Value")
    ax.set_title(f"Performances on {dataset} with 95% CI (Bootstrap)")

    plt.tight_layout()
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 0.85), ncols=4, frameon=False)

    plt.savefig(filename)
    plt.show()


def main():
    eval_folder = "../../boltz_results_final/"
    output_folder = "../../boltz_results_final/"

    # Eval the test set
    chai_preds = eval_folder + "outputs/test/chai"
    chai_evals = eval_folder + "evals/test/chai"

    af3_preds = eval_folder + "outputs/test/af3"
    af3_evals = eval_folder + "evals/test/af3"

    boltz_preds = eval_folder + "outputs/test/boltz/predictions"
    boltz_evals = eval_folder + "evals/test/boltz"

    boltz_preds_x = eval_folder + "outputs/test/boltzx/predictions"
    boltz_evals_x = eval_folder + "evals/test/boltzx"

    validity_checks = eval_folder + "physical_checks_test.csv"

    df_validity_checks = pd.read_csv(validity_checks)
    df_validity_checks = eval_validity_checks(df_validity_checks)

    df = eval_models(
        chai_preds,
        chai_evals,
        af3_preds,
        af3_evals,
        boltz_preds,
        boltz_evals,
        boltz_preds_x,
        boltz_evals_x,
    )

    df = pd.concat([df, df_validity_checks]).reset_index(drop=True)
    df.to_csv(output_folder + "results_test.csv", index=False)

    desired_tools = [
        "AF3 oracle",
        "AF3 top-1",
        "Chai-1 oracle",
        "Chai-1 top-1",
        "Boltz-1 oracle",
        "Boltz-1 top-1",
        "Boltz-1x oracle",
        "Boltz-1x top-1",
    ]
    desired_metrics = ["lddt", "dockq_>0.23", "lddt_pli", "rmsd<2", "physical validity"]
    plot_data(
        desired_tools, desired_metrics, df, "PDB Test", output_folder + "plot_test.pdf"
    )

    # Eval CASP
    chai_preds = eval_folder + "outputs/casp15/chai"
    chai_evals = eval_folder + "evals/casp15/chai"

    af3_preds = eval_folder + "outputs/casp15/af3"
    af3_evals = eval_folder + "evals/casp15/af3"

    boltz_preds = eval_folder + "outputs/casp15/boltz/predictions"
    boltz_evals = eval_folder + "evals/casp15/boltz"

    boltz_preds_x = eval_folder + "outputs/casp15/boltzx/predictions"
    boltz_evals_x = eval_folder + "evals/casp15/boltzx"

    validity_checks = eval_folder + "physical_checks_casp.csv"

    df_validity_checks = pd.read_csv(validity_checks)
    df_validity_checks = eval_validity_checks(df_validity_checks)

    df = eval_models(
        chai_preds,
        chai_evals,
        af3_preds,
        af3_evals,
        boltz_preds,
        boltz_evals,
        boltz_preds_x,
        boltz_evals_x,
    )

    df = pd.concat([df, df_validity_checks]).reset_index(drop=True)
    df.to_csv(output_folder + "results_casp.csv", index=False)

    plot_data(
        desired_tools, desired_metrics, df, "CASP15", output_folder + "plot_casp.pdf"
    )


if __name__ == "__main__":
    main()
