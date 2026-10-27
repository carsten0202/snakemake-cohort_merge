import csv
import random
from pathlib import Path


def main(snakemake):
    with open(snakemake.input.manifest) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    manifest = {row["merged_sample_id"]: row for row in rows}
    if len(manifest) != len(rows):
        raise ValueError("Sample manifest contains duplicate merged IDs")

    neighbors = {}
    with open(snakemake.input.pairs) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            first, second = row["sample_1"], row["sample_2"]
            if first not in manifest or second not in manifest or first == second:
                raise ValueError(f"Invalid duplicate pair: {first!r}, {second!r}")
            neighbors.setdefault(first, set()).add(second)
            neighbors.setdefault(second, set()).add(first)

    random_generator = random.Random(snakemake.params.random_seed)
    seen = set()
    decisions = []
    removals = []
    group_count = 0
    for sample in sorted(neighbors):
        if sample in seen:
            continue
        members = []
        pending = [sample]
        seen.add(sample)
        while pending:
            current = pending.pop()
            members.append(current)
            for neighbor in neighbors[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    pending.append(neighbor)
        members.sort()
        group_count += 1
        best_priority = min(int(manifest[member]["priority"]) for member in members)
        candidates = [
            member for member in members
            if int(manifest[member]["priority"]) == best_priority
        ]
        representative = random_generator.choice(candidates)
        fully_connected = all(len(neighbors[member]) == len(members) - 1 for member in members)
        for member in members:
            row = manifest[member]
            if member == representative:
                reason = "seeded_random_selection" if len(candidates) > 1 else "best_cohort_priority"
                decision = "keep"
            else:
                reason = "within_cohort_duplicate" if member in candidates else "lower_cohort_priority"
                decision = "remove"
                removals.append(member)
            decisions.append([
                f"group_{group_count:06d}", member, row["original_sample_id"],
                row["cohort"], row["priority"], decision, representative,
                reason, str(fully_connected).lower(),
            ])

    with open(snakemake.output.decisions, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "group_id", "merged_sample_id", "original_sample_id", "cohort",
            "priority", "decision", "retained_representative", "selection_reason",
            "fully_pairwise_connected",
        ])
        writer.writerows(decisions)
    Path(snakemake.output.removals).write_text("".join(f"{sample}\n" for sample in sorted(removals)))
    Path(snakemake.log[0]).write_text(
        f"Selected {group_count} representatives; removed {len(removals)} samples; "
        f"random seed {snakemake.params.random_seed}\n"
    )


if __name__ == "__main__":
    main(snakemake)
