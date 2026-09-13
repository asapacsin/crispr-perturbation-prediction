"""Read-only identity helpers for Norman condition labels and proposed splits."""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Mapping, Sequence

CONTROL_TOKEN = "ctrl"
SOURCE_TO_STORED_TARGET = {"RHOXF2": "RHOXF2BB"}
SPLIT_NAMES = ("train", "validation", "test")
BRIDGE_SPLIT = "excluded_val_test_bridge"
KNOWN_SPLITS = SPLIT_NAMES + (BRIDGE_SPLIT,)
GENE_PARTITION_LABELS = frozenset(SPLIT_NAMES)


class InvalidConditionLabel(ValueError):
    """Raised when a condition label is missing or syntactically invalid."""


def _is_missing(label: object) -> bool:
    if label is None:
        return True
    # pd.NA: `pd.NA != pd.NA` is itself NA and cannot be used as a boolean.
    if type(label).__name__ in {"NAType", "NaTType"}:
        return True
    if isinstance(label, float) and label != label:
        return True
    return False


def parse_condition_tokens(
    label: object,
    *,
    apply_source_alias: bool = False,
) -> list[str]:
    """Return non-control gene tokens from a stored or source condition label."""
    if _is_missing(label):
        raise InvalidConditionLabel("null/missing condition label")
    if not isinstance(label, str):
        raise InvalidConditionLabel(
            f"condition label must be a string, got {type(label).__name__}"
        )
    text = label.strip()
    if not text:
        raise InvalidConditionLabel("empty condition label")
    raw_tokens = [token.strip() for token in text.split("+")]
    if any(token == "" for token in raw_tokens):
        raise InvalidConditionLabel(f"empty token in label {label!r}")
    tokens: list[str] = []
    for token in raw_tokens:
        if token == CONTROL_TOKEN:
            continue
        if apply_source_alias:
            token = SOURCE_TO_STORED_TARGET.get(token, token)
        tokens.append(token)
    return tokens


def canonicalize_condition(
    label: object,
    *,
    apply_source_alias: bool = False,
) -> str:
    """Map a condition label to an order-insensitive target-set key."""
    tokens = parse_condition_tokens(label, apply_source_alias=apply_source_alias)
    if not tokens:
        return CONTROL_TOKEN
    return "+".join(sorted(set(tokens)))


def target_set(
    label: object,
    *,
    apply_source_alias: bool = False,
) -> frozenset[str]:
    return frozenset(parse_condition_tokens(label, apply_source_alias=apply_source_alias))


def n_targets(label: object, *, apply_source_alias: bool = False) -> int:
    return len(target_set(label, apply_source_alias=apply_source_alias))


def gene_partition_errors(gene_split: Mapping[str, str]) -> list[str]:
    """Return errors if gene-partition labels are invalid or not disjoint."""
    errors: list[str] = []
    by_split: dict[str, set[str]] = defaultdict(set)
    for gene, split in gene_split.items():
        if not isinstance(gene, str) or not gene.strip() or gene == CONTROL_TOKEN:
            errors.append(f"invalid gene name in partition: {gene!r}")
            continue
        if split not in GENE_PARTITION_LABELS:
            errors.append(f"unknown gene-partition label {split!r} for {gene}")
            continue
        by_split[split].add(gene)
    names = [name for name in SPLIT_NAMES if name in by_split]
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            overlap = sorted(by_split[left] & by_split[right])
            if overlap:
                errors.append(
                    f"gene partition overlap between {left} and {right}: {overlap}"
                )
    return errors


def combination_split_errors(
    rows: Iterable[tuple[object, str]],
    allowed_splits: Sequence[str] = SPLIT_NAMES,
) -> list[str]:
    """Return errors if non-control target sets appear in more than one split.

    Unknown split labels are rejected even for control rows. Orientation-
    equivalent labels are the same set.
    """
    allowed = set(allowed_splits)
    sets_by_split: dict[str, set[str]] = defaultdict(set)
    errors: list[str] = []
    for label, split in rows:
        if split not in allowed:
            errors.append(f"unknown split {split!r} for condition {label!r}")
            continue
        key = canonicalize_condition(label)
        if key == CONTROL_TOKEN:
            continue
        sets_by_split[split].add(key)
    names = [name for name in SPLIT_NAMES if name in sets_by_split]
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            overlap = sorted(sets_by_split[left] & sets_by_split[right])
            for key in overlap:
                errors.append(
                    f"combination leakage: target set {key} is in both {left} and {right}"
                )
    return errors


def target_gene_split_errors(
    rows: Iterable[tuple[object, str]],
    gene_split: Mapping[str, str],
    *,
    bridge_split: str = BRIDGE_SPLIT,
) -> list[str]:
    """Return errors if held-out genes leak across the gene-holdout protocol.

    Train may not contain validation or test genes. Validation may not contain
    test genes. Test may not contain validation genes. A pair that mixes a
    validation gene and a test gene must be assigned to the bridge exclusion
    and must actually contain both. Validation/test noncontrol conditions must
    contain an assigned held-out gene, not only train genes. Unknown genes and
    unknown split labels (including ctrl rows) are rejected.
    """
    errors = gene_partition_errors(gene_split)
    known_splits = set(SPLIT_NAMES) | {bridge_split}
    forbidden = {
        "train": {gene for gene, split in gene_split.items() if split in {"validation", "test"}},
        "validation": {gene for gene, split in gene_split.items() if split == "test"},
        "test": {gene for gene, split in gene_split.items() if split == "validation"},
    }
    for label, split in rows:
        if split not in known_splits:
            errors.append(f"unknown gene-split assignment {split!r} for {label!r}")
            continue
        genes = target_set(label)
        if not genes:
            if split == bridge_split:
                errors.append("control cannot be a validation/test gene bridge")
            continue
        unknown = sorted(gene for gene in genes if gene not in gene_split)
        if unknown:
            errors.append(
                f"unknown target gene(s) {unknown} not in gene partition "
                f"for {canonicalize_condition(label)}"
            )
        val_genes = {gene for gene in genes if gene_split.get(gene) == "validation"}
        test_genes = {gene for gene in genes if gene_split.get(gene) == "test"}
        if split == bridge_split:
            if not (val_genes and test_genes):
                errors.append(
                    f"excluded_val_test_bridge {canonicalize_condition(label)} "
                    "must contain validation AND test genes"
                )
            continue
        leaked = sorted(genes & forbidden[split])
        if leaked:
            errors.append(
                f"target-gene leakage: {split} condition {canonicalize_condition(label)} "
                f"contains forbidden gene(s) {leaked}"
            )
        if val_genes and test_genes:
            errors.append(
                f"val/test bridge not excluded: {canonicalize_condition(label)} "
                f"assigned to {split}"
            )
        if split == "validation" and not val_genes:
            errors.append(
                f"validation condition {canonicalize_condition(label)} "
                "has no validation-held-out gene"
            )
        if split == "test" and not test_genes:
            errors.append(
                f"test condition {canonicalize_condition(label)} "
                "has no test-held-out gene"
            )
    return errors


def barcode_split_errors(
    assignments: Sequence[tuple[str, str]],
    allowed_splits: Sequence[str] = KNOWN_SPLITS,
) -> list[str]:
    """Return errors if barcodes are duplicated or assigned to unknown splits."""
    seen: dict[str, str] = {}
    errors: list[str] = []
    for barcode, split in assignments:
        if split not in allowed_splits:
            errors.append(f"unknown split {split!r} for barcode {barcode}")
        if barcode in seen:
            errors.append(
                f"barcode {barcode} assigned to both {seen[barcode]} and {split}"
            )
        else:
            seen[barcode] = split
    return errors
